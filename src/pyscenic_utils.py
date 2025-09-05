"""
pySCENIC Utility Functions for Organoid Morphogen Analysis

This module contains utility functions for running pySCENIC analysis on organoid data,
including subsampling, preprocessing, and consensus building functions.
"""

import argparse
import anndata as ad
import pandas as pd
import numpy as np
import os
import random
import pickle
import glob
import scanpy as sc
import warnings
import yaml
import logging
from pathlib import Path
from typing import List, Dict, Optional, Union, Tuple

# pySCENIC imports
from dask.diagnostics import ProgressBar
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from distributed import LocalCluster, Client
import sys
import time
from functools import partial
from multiprocessing import Pool, cpu_count
import loompy as lp
import tqdm
import multiprocessing

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Set system limits for memory optimization
os.environ["MALLOC_ARENA_MAX"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"


def load_config(config_path: str) -> Dict:
    """
    Load configuration from YAML file.
    
    Parameters:
    -----------
    config_path : str
        Path to YAML configuration file
        
    Returns:
    --------
    Dict
        Configuration dictionary
    """
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def subsample_adata_strat(adata: ad.AnnData, 
                         column: str, 
                         ncells: int, 
                         seed: int, 
                         strat: Optional[str] = None, 
                         verbose: bool = False) -> ad.AnnData:
    """
    Subsample cells based on a column in the metadata and a target number of cells,
    optionally ensuring equal representation of strat categories within each group.
    
    Parameters:
    -----------
    adata : AnnData object
        The AnnData object to be subsampled
    column : str
        The column in adata.obs to group cells by
    ncells : int
        The number of cells to sample per group
    seed : int
        Random seed for reproducibility
    strat : str, optional
        The column in adata.obs that has categorical values for balanced sampling
    verbose : bool
        Whether to print warnings
        
    Returns:
    --------
    AnnData object
        Subsampled AnnData object
    """
    np.random.seed(seed)
    unique_groups = adata.obs[column].unique()
    
    subsampled_cells = []
    
    for group in unique_groups:
        group_cells = adata.obs[adata.obs[column] == group]
        
        if strat is not None:
            unique_strat = group_cells[strat].unique()
            cells_per_strat = ncells // len(unique_strat)
            
            for category in unique_strat:
                category_cells = group_cells[group_cells[strat] == category].index
                
                if len(category_cells) > cells_per_strat:
                    sampled_cells = np.random.choice(category_cells, cells_per_strat, replace=False)
                else:
                    sampled_cells = category_cells
                    if verbose:
                        logger.warning(f"Not enough cells in group '{group}', category '{category}'. "
                                     f"Taking all available ({len(category_cells)}).")
                
                subsampled_cells.extend(sampled_cells)
        else:
            if len(group_cells) > ncells:
                sampled_cells = np.random.choice(group_cells.index, ncells, replace=False)
            else:
                sampled_cells = group_cells.index
                if verbose:
                    logger.warning(f"Not enough cells in group '{group}'. "
                                 f"Taking all available ({len(group_cells)}).")
            
            subsampled_cells.extend(sampled_cells)
    
    return adata[subsampled_cells, :].copy()


def hvg_adata(adata: ad.AnnData, 
              n_top_genes: int = 15000, 
              flavor: str = 'seurat_v3') -> ad.AnnData:
    """
    Find highly variable genes in adata and subset adata.
    
    Parameters:
    -----------
    adata : AnnData object
        The AnnData object to be processed
    n_top_genes : int
        Number of top variable genes to return
    flavor : str
        Algorithm to use for finding variable genes
        
    Returns:
    --------
    AnnData object
        AnnData object subset to highly variable genes
    """
    # Make a copy to avoid modifying the original
    adata_hvg = adata.copy()
    
    # Identify highly variable genes
    sc.pp.highly_variable_genes(
        adata_hvg,
        n_top_genes=n_top_genes,
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        flavor=flavor
    )
    
    return adata_hvg[:, adata_hvg.var["highly_variable"]].copy()


def find_free_port() -> int:
    """
    Find a free port for the Dask scheduler.
    
    Returns:
    --------
    int
        Available port number
    """
    import socket
    
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(('', 0))
        return s.getsockname()[1]


def run_infer_partial_network(target_gene_index, gene_names, ex_matrix, 
                            tf_matrix, tf_matrix_gene_names, method_params, seed):
    """
    Run GRNBoost2 inference for a single target gene.
    
    This function is used for multiprocessing-based parallelization.
    """
    from arboreto.core import (
        EARLY_STOP_WINDOW_LENGTH,
        RF_KWARGS,
        SGBM_KWARGS,
        infer_partial_network,
        target_gene_indices,
        to_tf_matrix,
    )

    target_gene_name = gene_names[target_gene_index]
    target_gene_expression = ex_matrix[:, target_gene_index]

    n = infer_partial_network(
        regressor_type=method_params[0],
        regressor_kwargs=method_params[1],
        tf_matrix=tf_matrix,
        tf_matrix_gene_names=tf_matrix_gene_names,
        target_gene_name=target_gene_name,
        target_gene_expression=target_gene_expression,
        include_meta=False,
        early_stop_window_length=EARLY_STOP_WINDOW_LENGTH,
        seed=seed,
    )
    return n


def validate_input_data(adata: ad.AnnData, 
                       required_columns: List[str]) -> bool:
    """
    Validate that the input AnnData object has required metadata columns.
    
    Parameters:
    -----------
    adata : AnnData object
        Input data to validate
    required_columns : List[str]
        List of required column names in adata.obs
        
    Returns:
    --------
    bool
        True if all required columns are present
        
    Raises:
    -------
    ValueError
        If required columns are missing
    """
    missing_columns = [col for col in required_columns if col not in adata.obs.columns]
    
    if missing_columns:
        raise ValueError(f"Missing required columns in adata.obs: {missing_columns}")
    
    logger.info(f"Data validation passed. Shape: {adata.shape}")
    return True


def save_regulons(regulons: List, output_path: str, compress: bool = True):
    """
    Save regulons to pickle file with optional compression.
    
    Parameters:
    -----------
    regulons : List
        List of pySCENIC regulon objects
    output_path : str
        Output file path
    compress : bool
        Whether to compress the output file
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    if compress:
        import gzip
        with gzip.open(output_path + '.gz', 'wb') as f:
            pickle.dump(regulons, f)
    else:
        with open(output_path, 'wb') as f:
            pickle.dump(regulons, f)
    
    logger.info(f"Saved {len(regulons)} regulons to {output_path}")


def load_regulons(input_path: str) -> List:
    """
    Load regulons from pickle file, handling compressed files.
    
    Parameters:
    -----------
    input_path : str
        Input file path
        
    Returns:
    --------
    List
        List of pySCENIC regulon objects
    """
    if input_path.endswith('.gz'):
        import gzip
        with gzip.open(input_path, 'rb') as f:
            regulons = pickle.load(f)
    else:
        with open(input_path, 'rb') as f:
            regulons = pickle.load(f)
    
    logger.info(f"Loaded {len(regulons)} regulons from {input_path}")
    return regulons


def create_output_directories(base_dir: str, 
                            subdirs: List[str] = None) -> Dict[str, str]:
    """
    Create output directory structure.
    
    Parameters:
    -----------
    base_dir : str
        Base output directory
    subdirs : List[str]
        List of subdirectories to create
        
    Returns:
    --------
    Dict[str, str]
        Dictionary mapping subdirectory names to paths
    """
    if subdirs is None:
        subdirs = ['regulons', 'aucell', 'consensus', 'logs', 'plots']
    
    paths = {}
    for subdir in subdirs:
        path = os.path.join(base_dir, subdir)
        os.makedirs(path, exist_ok=True)
        paths[subdir] = path
    
    logger.info(f"Created output directories in {base_dir}")
    return paths


def get_memory_usage() -> str:
    """
    Get current memory usage of the process.
    
    Returns:
    --------
    str
        Memory usage string
    """
    import psutil
    process = psutil.Process(os.getpid())
    memory_mb = process.memory_info().rss / 1024 / 1024
    return f"{memory_mb:.1f} MB"


def log_progress(step: str, start_time: float = None):
    """
    Log progress information with timing.
    
    Parameters:
    -----------
    step : str
        Description of the current step
    start_time : float
        Start time for timing calculation
    """
    current_time = time.time()
    memory_usage = get_memory_usage()
    
    if start_time is not None:
        elapsed = current_time - start_time
        logger.info(f"{step} - Elapsed: {elapsed:.2f}s - Memory: {memory_usage}")
    else:
        logger.info(f"{step} - Memory: {memory_usage}")
    
    return current_time
