#!/usr/bin/env python3
"""
pySCENIC Execution Script for Organoid Morphogen Analysis

This script runs the complete pySCENIC pipeline with consensus building
across multiple random seeds to identify robust gene regulatory networks.
"""

import argparse
import os
import sys
import time
import logging
import yaml
from pathlib import Path
from typing import List, Optional

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from pyscenic_utils import (
    load_config, subsample_adata_strat, hvg_adata, 
    validate_input_data, create_output_directories,
    save_regulons, log_progress
)

import scanpy as sc
import pandas as pd
import numpy as np
import pickle
import glob
import random

# pySCENIC imports
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from dask.diagnostics import ProgressBar

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('pyscenic.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def run_pyscenic_single_seed(adata: sc.AnnData, 
                           config: dict,
                           db_config: dict,
                           output_dir: str, 
                           seed: int) -> bool:
    """
    Run pySCENIC pipeline for a single random seed.
    
    Parameters:
    -----------
    adata : AnnData object
        Expression data
    config : dict
        Configuration parameters
    db_config : dict
        Database paths configuration
    output_dir : str
        Output directory
    seed : int
        Random seed
        
    Returns:
    --------
    bool
        True if successful, False otherwise
    """
    start_time = log_progress(f"Starting pySCENIC with seed {seed}")
    
    try:
        # Set random seeds
        random.seed(seed)
        np.random.seed(seed)
        
        # Create seed-specific output directory
        seed_output_dir = os.path.join(output_dir, f"seed_{seed}")
        os.makedirs(seed_output_dir, exist_ok=True)
        
        # 1. Data preprocessing
        adata_processed = adata.copy()
        
        # Subsampling if specified
        subsample_config = config.get('subsampling', {})
        if subsample_config.get('subsample_column') and subsample_config.get('subsample_ncells'):
            log_progress("Subsampling cells")
            adata_processed = subsample_adata_strat(
                adata_processed,
                subsample_config['subsample_column'],
                subsample_config['subsample_ncells'],
                seed,
                strat=subsample_config.get('stratification_column'),
                verbose=True
            )
        
        # HVG filtering if specified
        if 'hvg_genes' in config.get('pyscenic', {}):
            log_progress("Filtering to highly variable genes")
            adata_processed = hvg_adata(
                adata_processed, 
                n_top_genes=config['pyscenic']['hvg_genes']
            )
        
        log_progress(f"Using data with shape: {adata_processed.shape}")
        
        # Convert to DataFrame
        if hasattr(adata_processed.X, 'toarray'):
            dgem = pd.DataFrame(
                adata_processed.X.toarray(),
                index=adata_processed.obs_names,
                columns=adata_processed.var_names
            )
        else:
            dgem = pd.DataFrame(
                adata_processed.X,
                index=adata_processed.obs_names,
                columns=adata_processed.var_names
            )
        dgem = dgem.astype('float32')
        
        # 2. Load databases
        log_progress("Loading pySCENIC databases")
        
        # Load TF names
        tf_names = load_tf_names(db_config['tf_database'])
        logger.info(f"Loaded {len(tf_names)} transcription factors")
        
        # Load ranking databases
        db_fnames = glob.glob(os.path.join(db_config['ranking_database'], "*.feather"))
        if not db_fnames:
            db_fnames = glob.glob(os.path.join(db_config['ranking_database'], "*.genes_vs_motifs.rankings.feather"))
        
        if not db_fnames:
            raise FileNotFoundError(f"No ranking database files found in {db_config['ranking_database']}")
        
        dbs = [RankingDatabase(fname=fname, name=os.path.splitext(os.path.basename(fname))[0]) 
               for fname in db_fnames]
        logger.info(f"Loaded {len(dbs)} ranking databases")
        
        # 3. GRN inference with GRNBoost2
        log_progress("Running GRNBoost2 for GRN inference")
        
        pyscenic_config = config.get('pyscenic', {})
        method = pyscenic_config.get('method', 'multiproc')
        
        if method == 'dask':
            adjacencies = _run_grnboost2_dask(
                dgem, tf_names, seed, pyscenic_config
            )
        elif method == 'multiproc':
            adjacencies = _run_grnboost2_multiproc(
                adata_processed, tf_names, seed, pyscenic_config
            )
        else:
            raise ValueError(f"Unknown method: {method}")
        
        # Save adjacencies
        adjacencies.to_csv(os.path.join(seed_output_dir, 'adjacencies.csv'))
        
        # 4. Module discovery
        log_progress("Discovering co-expression modules")
        modules = list(modules_from_adjacencies(adjacencies, dgem))
        
        # Save modules
        with open(os.path.join(seed_output_dir, 'modules.p'), 'wb') as f:
            pickle.dump(modules, f)
        
        logger.info(f"Found {len(modules)} co-expression modules")
        
        # 5. Motif enrichment analysis
        log_progress("Running motif enrichment analysis")
        
        with ProgressBar():
            df_motifs = prune2df(dbs, modules, db_config['motif_database'])
        
        # Save motif enrichment results
        df_motifs.to_csv(os.path.join(seed_output_dir, 'motif_enrichment.csv'))
        
        # 6. Create regulons
        log_progress("Creating final regulons")
        
        regulons = df2regulons(df_motifs)
        
        # Filter regulons by minimum size
        min_size = config.get('consensus', {}).get('min_regulon_size', 5)
        regulons_filtered = [r for r in regulons if len(r.gene2weight) >= min_size]
        
        logger.info(f"Created {len(regulons_filtered)} regulons (filtered from {len(regulons)})")
        
        # Save regulons
        save_regulons(regulons_filtered, os.path.join(seed_output_dir, 'regulons.p'))
        
        # 7. AUCell scoring (optional)
        if config.get('run_aucell', True):
            log_progress("Running AUCell scoring")
            
            auc_mtx = aucell(dgem, regulons_filtered, num_workers=pyscenic_config.get('n_workers', 4))
            auc_mtx.to_csv(os.path.join(seed_output_dir, 'aucell_matrix.csv'))
            
            with open(os.path.join(seed_output_dir, 'aucell_matrix.p'), 'wb') as f:
                pickle.dump(auc_mtx, f)
        
        total_time = time.time() - start_time
        logger.info(f"Successfully completed pySCENIC for seed {seed} in {total_time:.2f} seconds")
        
        return True
        
    except Exception as e:
        logger.error(f"Error in pySCENIC with seed {seed}: {str(e)}")
        return False


def _run_grnboost2_dask(dgem: pd.DataFrame, 
                       tf_names: List[str], 
                       seed: int, 
                       config: dict) -> pd.DataFrame:
    """Run GRNBoost2 using Dask backend."""
    from dask.distributed import LocalCluster, Client
    import dask
    
    # Configure Dask
    dask.config.set({
        "distributed.comm.max-message-size": "25GB",
        "distributed.worker.memory.target": 0.8,
        "distributed.worker.memory.spill": 0.9
    })
    
    n_workers = config.get('n_workers', 4)
    memory_limit = config.get('memory_limit', '8GB')
    
    with LocalCluster(
        n_workers=n_workers,
        threads_per_worker=1,
        memory_limit=memory_limit,
        scheduler_port=0,
        silence_logs=logging.ERROR
    ) as cluster:
        with Client(cluster) as client:
            logger.info(f"Dask dashboard: {client.dashboard_link}")
            
            adjacencies = grnboost2(
                expression_data=dgem,
                tf_names=tf_names,
                client_or_address=client,
                seed=seed,
                verbose=True
            )
    
    return adjacencies


def _run_grnboost2_multiproc(adata: sc.AnnData, 
                           tf_names: List[str], 
                           seed: int, 
                           config: dict) -> pd.DataFrame:
    """Run GRNBoost2 using multiprocessing backend."""
    from arboreto.algo import grnboost2
    
    adjacencies = grnboost2(
        expression_data=adata,
        tf_names=tf_names,
        seed=seed,
        verbose=True
    )
    
    return adjacencies


def build_consensus_regulons(output_dir: str, 
                           config: dict) -> List:
    """
    Build consensus regulons from multiple seed runs.
    
    Parameters:
    -----------
    output_dir : str
        Directory containing seed results
    config : dict
        Configuration parameters
        
    Returns:
    --------
    List
        List of consensus regulons
    """
    from collections import defaultdict
    import pyscenic as ps
    
    log_progress("Building consensus regulons")
    
    consensus_config = config.get('consensus', {})
    occurrence_threshold = consensus_config.get('occurrence_threshold', 0)
    
    # Collect all regulons from all seeds
    all_regulons = {}
    seed_dirs = [d for d in os.listdir(output_dir) if d.startswith('seed_')]
    
    for seed_dir in seed_dirs:
        regulon_file = os.path.join(output_dir, seed_dir, 'regulons.p')
        
        if os.path.exists(regulon_file):
            with open(regulon_file, 'rb') as f:
                seed_regulons = pickle.load(f)
            
            seed_num = seed_dir.replace('seed_', '')
            all_regulons[seed_num] = seed_regulons
    
    logger.info(f"Loaded regulons from {len(all_regulons)} seed runs")
    
    # Build consensus
    tf_occurrences = defaultdict(int)
    tf_regulons = defaultdict(list)
    
    # Count TF occurrences across seeds
    for seed, regulons in all_regulons.items():
        for regulon in regulons:
            tf_name = regulon.transcription_factor
            tf_occurrences[tf_name] += 1
            tf_regulons[tf_name].append(regulon)
    
    # Filter TFs by occurrence threshold
    n_seeds = len(all_regulons)
    min_occurrences = max(1, int(occurrence_threshold * n_seeds))
    
    consensus_regulons = []
    
    for tf_name, count in tf_occurrences.items():
        if count >= min_occurrences:
            # Merge regulons for this TF
            all_targets = defaultdict(list)
            
            for regulon in tf_regulons[tf_name]:
                for gene, weight in regulon.gene2weight.items():
                    all_targets[gene].append(weight)
            
            # Calculate consensus weights (mean)
            consensus_targets = {}
            for gene, weights in all_targets.items():
                consensus_targets[gene] = np.mean(weights)
            
            # Create consensus regulon
            consensus_regulon = ps.utils.Regulon(
                name=f"{tf_name}_consensus",
                gene2weight=consensus_targets,
                gene2occurrence={},
                transcription_factor=tf_name
            )
            
            consensus_regulons.append(consensus_regulon)
    
    logger.info(f"Built {len(consensus_regulons)} consensus regulons from {len(tf_occurrences)} total TFs")
    
    return consensus_regulons


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(description="Run pySCENIC pipeline with consensus building")
    
    parser.add_argument("--config", required=True, 
                       help="Path to configuration YAML file")
    parser.add_argument("--db_config", 
                       help="Path to database configuration YAML file")
    parser.add_argument("--input", 
                       help="Input H5AD file (overrides config)")
    parser.add_argument("--output", 
                       help="Output directory (overrides config)")
    parser.add_argument("--seeds", nargs="+", type=int,
                       help="Random seeds to use (overrides config)")
    parser.add_argument("--consensus_only", action="store_true",
                       help="Only build consensus from existing results")
    
    args = parser.parse_args()
    
    # Load configurations
    config = load_config(args.config)
    
    if args.db_config:
        db_config = load_config(args.db_config)
    else:
        # Try to load from config directory
        config_dir = os.path.dirname(args.config)
        db_config_path = os.path.join(config_dir, 'database_paths.yaml')
        if os.path.exists(db_config_path):
            db_config = load_config(db_config_path)
        else:
            raise FileNotFoundError("Database configuration file not found")
    
    # Override config with command line arguments
    if args.input:
        config['data']['input_file'] = args.input
    if args.output:
        config['output']['results_dir'] = args.output
    if args.seeds:
        config['consensus']['seeds'] = args.seeds
    
    # Set up output directory
    output_dir = config['output']['results_dir']
    create_output_directories(output_dir)
    
    if not args.consensus_only:
        # Load and validate input data
        input_file = config['data']['input_file']
        logger.info(f"Loading data from {input_file}")
        
        adata = sc.read_h5ad(input_file)
        
        # Validate input data
        required_columns = []
        if config.get('subsampling', {}).get('subsample_column'):
            required_columns.append(config['subsampling']['subsample_column'])
        if config.get('subsampling', {}).get('stratification_column'):
            required_columns.append(config['subsampling']['stratification_column'])
        
        if required_columns:
            validate_input_data(adata, required_columns)
        
        # Run pySCENIC for multiple seeds
        seeds = config.get('consensus', {}).get('seeds', [42])
        if isinstance(seeds, int):
            seeds = [seeds]
        
        logger.info(f"Running pySCENIC for seeds: {seeds}")
        
        successful_seeds = []
        for seed in seeds:
            success = run_pyscenic_single_seed(adata, config, db_config, output_dir, seed)
            if success:
                successful_seeds.append(seed)
        
        logger.info(f"Successfully completed {len(successful_seeds)}/{len(seeds)} seed runs")
        
        if len(successful_seeds) == 0:
            logger.error("No successful seed runs. Exiting.")
            sys.exit(1)
    
    # Build consensus regulons
    if len(os.listdir(output_dir)) > 0:
        consensus_regulons = build_consensus_regulons(output_dir, config)
        
        # Save consensus regulons
        consensus_dir = os.path.join(output_dir, 'consensus')
        os.makedirs(consensus_dir, exist_ok=True)
        
        save_regulons(consensus_regulons, os.path.join(consensus_dir, 'consensus_regulons.p'))
        
        logger.info("pySCENIC pipeline completed successfully!")
    else:
        logger.warning("No results found for consensus building")


if __name__ == "__main__":
    main()
