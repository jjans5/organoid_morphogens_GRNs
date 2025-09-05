"""
Morphogen Analysis Functions for pySCENIC Results

This module contains functions for analyzing the relationship between
regulon activities and morphogen treatments in organoid data.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy import stats
from scipy.stats import pearsonr, spearmanr
from typing import List, Dict, Optional, Tuple, Union
import logging
import warnings

logger = logging.getLogger(__name__)


def prepare_morphogen_metadata(adata: sc.AnnData, 
                             morphogen_list: List[str],
                             time_variable: str = "Time",
                             medium_variable: str = "medium") -> pd.DataFrame:
    """
    Prepare morphogen metadata for correlation analysis.
    
    Parameters:
    -----------
    adata : AnnData object
        Single-cell data with metadata
    morphogen_list : List[str]
        List of morphogen columns to analyze
    time_variable : str
        Name of the time variable column
    medium_variable : str
        Name of the medium variable column
        
    Returns:
    --------
    pd.DataFrame
        Processed metadata with morphogen and time information
    """
    meta = adata.obs.copy()
    
    # Convert time to quantitative values
    if time_variable in meta.columns:
        meta['Time_q'] = np.nan
        meta.loc[meta[time_variable] == 'late', 'Time_q'] = 1
        meta.loc[meta[time_variable] == 'middle', 'Time_q'] = 0
        meta.loc[meta[time_variable] == 'early', 'Time_q'] = -1
        morphogen_list.append('Time_q')
    
    # One-hot encode medium variable
    if medium_variable in meta.columns:
        for medium in set(meta[medium_variable]):
            meta[medium] = 0
            meta.loc[meta[medium_variable] == medium, medium] = 1
            morphogen_list.append(medium)
    
    # One-hot encode time variable
    if time_variable in meta.columns:
        for time_point in set(meta[time_variable]):
            meta[time_point] = 0
            meta.loc[meta[time_variable] == time_point, time_point] = 1
            morphogen_list.append(time_point)
    
    # Select and normalize morphogen data
    meta_use = meta[morphogen_list].copy()
    
    # Log-transform and normalize
    meta_use = np.log1p(meta_use)
    meta_use = (meta_use - meta_use.min()) / (meta_use.max() - meta_use.min())
    
    return meta_use


def calculate_regulon_morphogen_correlations(aucell_matrix: pd.DataFrame,
                                           morphogen_metadata: pd.DataFrame,
                                           method: str = "pearson",
                                           min_cells: int = 50) -> pd.DataFrame:
    """
    Calculate correlations between regulon activities and morphogen treatments.
    
    Parameters:
    -----------
    aucell_matrix : pd.DataFrame
        AUCell activity matrix (cells x regulons)
    morphogen_metadata : pd.DataFrame
        Morphogen metadata (cells x morphogens)
    method : str
        Correlation method ("pearson" or "spearman")
    min_cells : int
        Minimum number of cells required for correlation calculation
        
    Returns:
    --------
    pd.DataFrame
        Correlation matrix with regulons as rows and morphogens as columns
    """
    # Align indices
    common_cells = aucell_matrix.index.intersection(morphogen_metadata.index)
    
    if len(common_cells) < min_cells:
        raise ValueError(f"Not enough common cells ({len(common_cells)}) for correlation analysis")
    
    aucell_aligned = aucell_matrix.loc[common_cells]
    morphogen_aligned = morphogen_metadata.loc[common_cells]
    
    logger.info(f"Calculating correlations for {len(common_cells)} cells, "
                f"{aucell_aligned.shape[1]} regulons, and {morphogen_aligned.shape[1]} morphogens")
    
    # Calculate correlations
    correlations = []
    p_values = []
    
    for regulon in aucell_aligned.columns:
        regulon_correlations = []
        regulon_p_values = []
        
        for morphogen in morphogen_aligned.columns:
            if method == "pearson":
                corr, p_val = pearsonr(aucell_aligned[regulon], morphogen_aligned[morphogen])
            elif method == "spearman":
                corr, p_val = spearmanr(aucell_aligned[regulon], morphogen_aligned[morphogen])
            else:
                raise ValueError(f"Unknown correlation method: {method}")
            
            regulon_correlations.append(corr)
            regulon_p_values.append(p_val)
        
        correlations.append(regulon_correlations)
        p_values.append(regulon_p_values)
    
    # Create result DataFrames
    corr_df = pd.DataFrame(correlations, 
                          index=aucell_aligned.columns, 
                          columns=morphogen_aligned.columns)
    
    p_val_df = pd.DataFrame(p_values, 
                           index=aucell_aligned.columns, 
                           columns=morphogen_aligned.columns)
    
    return corr_df, p_val_df


def find_morphogen_responsive_regulons(correlation_df: pd.DataFrame,
                                     p_value_df: pd.DataFrame,
                                     morphogen_list: List[str],
                                     correlation_threshold: float = 0.3,
                                     p_value_threshold: float = 0.05) -> pd.DataFrame:
    """
    Identify regulons that are significantly correlated with morphogens.
    
    Parameters:
    -----------
    correlation_df : pd.DataFrame
        Correlation matrix (regulons x morphogens)
    p_value_df : pd.DataFrame
        P-value matrix (regulons x morphogens)
    morphogen_list : List[str]
        List of morphogens to analyze
    correlation_threshold : float
        Minimum absolute correlation threshold
    p_value_threshold : float
        Maximum p-value threshold
        
    Returns:
    --------
    pd.DataFrame
        Summary of morphogen-responsive regulons
    """
    results = []
    
    for morphogen in morphogen_list:
        if morphogen not in correlation_df.columns:
            logger.warning(f"Morphogen {morphogen} not found in correlation matrix")
            continue
        
        # Find significant correlations
        significant_mask = (np.abs(correlation_df[morphogen]) > correlation_threshold) & \
                          (p_value_df[morphogen] < p_value_threshold)
        
        significant_regulons = correlation_df.loc[significant_mask, morphogen]
        significant_p_values = p_value_df.loc[significant_mask, morphogen]
        
        for regulon in significant_regulons.index:
            results.append({
                'regulon': regulon,
                'morphogen': morphogen,
                'correlation': significant_regulons[regulon],
                'p_value': significant_p_values[regulon],
                'abs_correlation': np.abs(significant_regulons[regulon]),
                'direction': 'positive' if significant_regulons[regulon] > 0 else 'negative'
            })
    
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        # Sort by absolute correlation strength
        results_df = results_df.sort_values('abs_correlation', ascending=False)
        logger.info(f"Found {len(results_df)} significant morphogen-regulon associations")
    else:
        logger.warning("No significant morphogen-regulon associations found")
    
    return results_df


def create_morphogen_summary_table(responsive_regulons_df: pd.DataFrame,
                                 regulon_info: Dict = None) -> pd.DataFrame:
    """
    Create a summary table of morphogen-responsive regulons.
    
    Parameters:
    -----------
    responsive_regulons_df : pd.DataFrame
        DataFrame of responsive regulons from find_morphogen_responsive_regulons
    regulon_info : Dict, optional
        Additional information about regulons (e.g., target genes)
        
    Returns:
    --------
    pd.DataFrame
        Summary table with regulon statistics
    """
    summary_stats = []
    
    for regulon in responsive_regulons_df['regulon'].unique():
        regulon_data = responsive_regulons_df[responsive_regulons_df['regulon'] == regulon]
        
        # Count morphogen associations
        n_morphogens = len(regulon_data)
        max_abs_corr = regulon_data['abs_correlation'].max()
        best_morphogen = regulon_data.loc[regulon_data['abs_correlation'].idxmax(), 'morphogen']
        
        # Direction analysis
        n_positive = len(regulon_data[regulon_data['direction'] == 'positive'])
        n_negative = len(regulon_data[regulon_data['direction'] == 'negative'])
        
        summary_stats.append({
            'regulon': regulon,
            'transcription_factor': regulon.split('(')[0] if '(' in regulon else regulon,
            'n_morphogen_associations': n_morphogens,
            'max_abs_correlation': max_abs_corr,
            'best_morphogen': best_morphogen,
            'n_positive_correlations': n_positive,
            'n_negative_correlations': n_negative,
            'morphogens': ', '.join(regulon_data['morphogen'].tolist())
        })
    
    summary_df = pd.DataFrame(summary_stats)
    summary_df = summary_df.sort_values('max_abs_correlation', ascending=False)
    
    return summary_df


def plot_morphogen_correlation_heatmap(correlation_df: pd.DataFrame,
                                     p_value_df: pd.DataFrame,
                                     morphogen_list: List[str],
                                     top_n_regulons: int = 50,
                                     figsize: Tuple[int, int] = (12, 8),
                                     save_path: str = None) -> plt.Figure:
    """
    Create a heatmap of morphogen-regulon correlations.
    
    Parameters:
    -----------
    correlation_df : pd.DataFrame
        Correlation matrix
    p_value_df : pd.DataFrame
        P-value matrix
    morphogen_list : List[str]
        List of morphogens to plot
    top_n_regulons : int
        Number of top regulons to show
    figsize : Tuple[int, int]
        Figure size
    save_path : str, optional
        Path to save the figure
        
    Returns:
    --------
    plt.Figure
        Matplotlib figure object
    """
    # Select morphogens that exist in the data
    available_morphogens = [m for m in morphogen_list if m in correlation_df.columns]
    
    if not available_morphogens:
        raise ValueError("No morphogens found in correlation matrix")
    
    # Select top regulons based on maximum absolute correlation
    max_abs_corr = correlation_df[available_morphogens].abs().max(axis=1)
    top_regulons = max_abs_corr.nlargest(top_n_regulons).index
    
    # Create subsets for plotting
    plot_corr = correlation_df.loc[top_regulons, available_morphogens]
    plot_pval = p_value_df.loc[top_regulons, available_morphogens]
    
    # Create significance mask
    significance_mask = plot_pval < 0.05
    
    # Create the plot
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot heatmap
    sns.heatmap(plot_corr, 
                annot=False,
                cmap='RdBu_r',
                center=0,
                vmin=-1, vmax=1,
                cbar_kws={'label': 'Correlation'},
                ax=ax)
    
    # Add significance markers
    for i, regulon in enumerate(top_regulons):
        for j, morphogen in enumerate(available_morphogens):
            if significance_mask.loc[regulon, morphogen]:
                ax.text(j + 0.5, i + 0.5, '*', 
                       ha='center', va='center', 
                       color='black', fontsize=8, fontweight='bold')
    
    ax.set_title(f'Top {top_n_regulons} Regulon-Morphogen Correlations\\n(* p < 0.05)', 
                fontsize=14, fontweight='bold')
    ax.set_xlabel('Morphogens', fontsize=12)
    ax.set_ylabel('Regulons', fontsize=12)
    
    # Rotate labels for better readability
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved heatmap to {save_path}")
    
    return fig


def compare_cellline_morphogen_responses(cellline_results: Dict[str, pd.DataFrame],
                                       morphogen: str,
                                       min_overlap: int = 5) -> pd.DataFrame:
    """
    Compare morphogen responses across different cell lines.
    
    Parameters:
    -----------
    cellline_results : Dict[str, pd.DataFrame]
        Dictionary mapping cell line names to responsive regulons DataFrames
    morphogen : str
        Morphogen to compare
    min_overlap : int
        Minimum number of cell lines a regulon must appear in
        
    Returns:
    --------
    pd.DataFrame
        Comparison table with consistency metrics
    """
    all_regulons = set()
    cellline_morphogen_data = {}
    
    # Collect data for the specific morphogen
    for cellline, results_df in cellline_results.items():
        morphogen_data = results_df[results_df['morphogen'] == morphogen]
        cellline_morphogen_data[cellline] = morphogen_data.set_index('regulon')
        all_regulons.update(morphogen_data['regulon'].tolist())
    
    # Create comparison table
    comparison_data = []
    
    for regulon in all_regulons:
        regulon_data = {
            'regulon': regulon,
            'n_celllines': 0,
            'consistent_direction': True,
            'mean_correlation': 0,
            'std_correlation': 0,
            'celllines': [],
            'correlations': []
        }
        
        correlations = []
        directions = []
        
        for cellline, data in cellline_morphogen_data.items():
            if regulon in data.index:
                regulon_data['n_celllines'] += 1
                regulon_data['celllines'].append(cellline)
                
                corr = data.loc[regulon, 'correlation']
                correlations.append(corr)
                directions.append('positive' if corr > 0 else 'negative')
        
        if len(correlations) >= min_overlap:
            regulon_data['mean_correlation'] = np.mean(correlations)
            regulon_data['std_correlation'] = np.std(correlations)
            regulon_data['correlations'] = correlations
            regulon_data['consistent_direction'] = len(set(directions)) == 1
            
            comparison_data.append(regulon_data)
    
    comparison_df = pd.DataFrame(comparison_data)
    
    if len(comparison_df) > 0:
        comparison_df = comparison_df.sort_values(['n_celllines', 'mean_correlation'], 
                                                ascending=[False, False])
        logger.info(f"Found {len(comparison_df)} regulons responsive to {morphogen} "
                   f"in multiple cell lines")
    
    return comparison_df


def calculate_morphogen_specificity_scores(correlation_df: pd.DataFrame,
                                         morphogen_list: List[str]) -> pd.DataFrame:
    """
    Calculate specificity scores for regulon-morphogen associations.
    
    Parameters:
    -----------
    correlation_df : pd.DataFrame
        Correlation matrix (regulons x morphogens)
    morphogen_list : List[str]
        List of primary morphogens to analyze
        
    Returns:
    --------
    pd.DataFrame
        DataFrame with specificity scores for each regulon-morphogen pair
    """
    available_morphogens = [m for m in morphogen_list if m in correlation_df.columns]
    specificity_results = []
    
    for regulon in correlation_df.index:
        regulon_correlations = correlation_df.loc[regulon, available_morphogens]
        
        for morphogen in available_morphogens:
            target_corr = regulon_correlations[morphogen]
            other_corrs = regulon_correlations.drop(morphogen)
            
            # Calculate specificity as the difference between target and max other correlation
            max_other_corr = other_corrs.abs().max()
            specificity_score = abs(target_corr) - max_other_corr
            
            # Calculate relative specificity
            if regulon_correlations.abs().sum() > 0:
                relative_specificity = abs(target_corr) / regulon_correlations.abs().sum()
            else:
                relative_specificity = 0
            
            specificity_results.append({
                'regulon': regulon,
                'morphogen': morphogen,
                'correlation': target_corr,
                'specificity_score': specificity_score,
                'relative_specificity': relative_specificity,
                'is_most_specific': abs(target_corr) == regulon_correlations.abs().max()
            })
    
    return pd.DataFrame(specificity_results)
