"""
Visualization Functions for pySCENIC Organoid Analysis

This module contains functions for creating publication-quality plots
for the organoid morphogen pySCENIC analysis.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import pandas as pd
import numpy as np
import networkx as nx
from typing import List, Dict, Optional, Tuple, Union
import logging
from adjustText import adjust_text

logger = logging.getLogger(__name__)

# Set publication-ready style
plt.style.use('default')
sns.set_palette("husl")


def plot_regulon_size_distribution(regulons: List, 
                                 figsize: Tuple[int, int] = (10, 6),
                                 save_path: str = None) -> plt.Figure:
    """
    Plot the distribution of regulon sizes.
    
    Parameters:
    -----------
    regulons : List
        List of pySCENIC regulon objects
    figsize : Tuple[int, int]
        Figure size
    save_path : str, optional
        Path to save the figure
        
    Returns:
    --------
    plt.Figure
        Matplotlib figure object
    """
    # Extract regulon sizes
    sizes = [len(regulon.gene2weight) for regulon in regulons]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Histogram
    ax1.hist(sizes, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    ax1.set_xlabel('Regulon Size (Number of Target Genes)')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Distribution of Regulon Sizes')
    ax1.grid(True, alpha=0.3)
    
    # Add statistics
    mean_size = np.mean(sizes)
    median_size = np.median(sizes)
    ax1.axvline(mean_size, color='red', linestyle='--', label=f'Mean: {mean_size:.1f}')
    ax1.axvline(median_size, color='orange', linestyle='--', label=f'Median: {median_size:.1f}')
    ax1.legend()
    
    # Box plot
    ax2.boxplot(sizes, vert=True)
    ax2.set_ylabel('Regulon Size (Number of Target Genes)')
    ax2.set_title('Regulon Size Distribution')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved regulon size distribution plot to {save_path}")
    
    return fig


def plot_morphogen_network(correlation_df: pd.DataFrame,
                         p_value_df: pd.DataFrame,
                         morphogens: List[str],
                         correlation_threshold: float = 0.3,
                         p_value_threshold: float = 0.05,
                         figsize: Tuple[int, int] = (12, 10),
                         save_path: str = None) -> plt.Figure:
    """
    Create a network plot showing morphogen-regulon relationships.
    
    Parameters:
    -----------
    correlation_df : pd.DataFrame
        Correlation matrix (regulons x morphogens)
    p_value_df : pd.DataFrame
        P-value matrix (regulons x morphogens)
    morphogens : List[str]
        List of morphogens to include
    correlation_threshold : float
        Minimum correlation strength for edges
    p_value_threshold : float
        Maximum p-value for significance
    figsize : Tuple[int, int]
        Figure size
    save_path : str, optional
        Path to save the figure
        
    Returns:
    --------
    plt.Figure
        Matplotlib figure object
    """
    # Filter data
    available_morphogens = [m for m in morphogens if m in correlation_df.columns]
    
    # Create network graph
    G = nx.Graph()
    
    # Add morphogen nodes
    for morphogen in available_morphogens:
        G.add_node(morphogen, node_type='morphogen')
    
    # Add significant regulon-morphogen edges
    significant_connections = []
    
    for morphogen in available_morphogens:
        significant_mask = (np.abs(correlation_df[morphogen]) > correlation_threshold) & \
                          (p_value_df[morphogen] < p_value_threshold)
        
        significant_regulons = correlation_df.loc[significant_mask, morphogen]
        
        for regulon in significant_regulons.index:
            # Add regulon node if not exists
            if regulon not in G.nodes():
                G.add_node(regulon, node_type='regulon')
            
            # Add edge
            correlation = significant_regulons[regulon]
            G.add_edge(morphogen, regulon, 
                      weight=abs(correlation),
                      correlation=correlation,
                      edge_type='positive' if correlation > 0 else 'negative')
            
            significant_connections.append((morphogen, regulon, correlation))
    
    if len(G.nodes()) == 0:
        logger.warning("No significant connections found for network plot")
        return None
    
    # Create layout
    fig, ax = plt.subplots(figsize=figsize)
    
    # Use spring layout with morphogens fixed in center
    morphogen_positions = {}
    angle_step = 2 * np.pi / len(available_morphogens)
    for i, morphogen in enumerate(available_morphogens):
        angle = i * angle_step
        morphogen_positions[morphogen] = (0.3 * np.cos(angle), 0.3 * np.sin(angle))
    
    pos = nx.spring_layout(G, pos=morphogen_positions, fixed=available_morphogens, k=1, iterations=100)
    
    # Draw nodes
    morphogen_nodes = [node for node in G.nodes() if G.nodes[node]['node_type'] == 'morphogen']
    regulon_nodes = [node for node in G.nodes() if G.nodes[node]['node_type'] == 'regulon']
    
    nx.draw_networkx_nodes(G, pos, nodelist=morphogen_nodes, 
                          node_color='red', node_size=500, alpha=0.8, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=regulon_nodes, 
                          node_color='lightblue', node_size=100, alpha=0.6, ax=ax)
    
    # Draw edges with different colors for positive/negative correlations
    positive_edges = [e for e in G.edges(data=True) if e[2]['edge_type'] == 'positive']
    negative_edges = [e for e in G.edges(data=True) if e[2]['edge_type'] == 'negative']
    
    if positive_edges:
        edge_weights = [e[2]['weight'] * 3 for e in positive_edges]
        nx.draw_networkx_edges(G, pos, edgelist=[(e[0], e[1]) for e in positive_edges],
                              width=edge_weights, edge_color='green', alpha=0.6, ax=ax)
    
    if negative_edges:
        edge_weights = [e[2]['weight'] * 3 for e in negative_edges]
        nx.draw_networkx_edges(G, pos, edgelist=[(e[0], e[1]) for e in negative_edges],
                              width=edge_weights, edge_color='red', alpha=0.6, ax=ax)
    
    # Add labels for morphogens only (too many regulons to label clearly)
    morphogen_labels = {node: node for node in morphogen_nodes}
    nx.draw_networkx_labels(G, pos, morphogen_labels, font_size=12, font_weight='bold', ax=ax)
    
    # Add legend
    legend_elements = [
        patches.Patch(color='red', label='Morphogens'),
        patches.Patch(color='lightblue', label='Regulons'),
        plt.Line2D([0], [0], color='green', lw=2, label='Positive correlation'),
        plt.Line2D([0], [0], color='red', lw=2, label='Negative correlation')
    ]
    ax.legend(handles=legend_elements, loc='upper right')
    
    ax.set_title(f'Morphogen-Regulon Network\\n'
                f'(|r| > {correlation_threshold}, p < {p_value_threshold})', 
                fontsize=14, fontweight='bold')
    ax.axis('off')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved network plot to {save_path}")
    
    return fig


def plot_cellline_comparison(comparison_results: Dict[str, pd.DataFrame],
                           morphogen: str,
                           figsize: Tuple[int, int] = (12, 8),
                           save_path: str = None) -> plt.Figure:
    """
    Create a comparison plot across cell lines for a specific morphogen.
    
    Parameters:
    -----------
    comparison_results : Dict[str, pd.DataFrame]
        Dictionary mapping cell line names to results DataFrames
    morphogen : str
        Morphogen to compare
    figsize : Tuple[int, int]
        Figure size
    save_path : str, optional
        Path to save the figure
        
    Returns:
    --------
    plt.Figure
        Matplotlib figure object
    """
    # Extract data for the morphogen
    all_data = []
    
    for cellline, results_df in comparison_results.items():
        morphogen_data = results_df[results_df['morphogen'] == morphogen].copy()
        morphogen_data['cellline'] = cellline
        all_data.append(morphogen_data)
    
    if not all_data:
        logger.warning(f"No data found for morphogen {morphogen}")
        return None
    
    combined_data = pd.concat(all_data, ignore_index=True)
    
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    
    # 1. Correlation distribution by cell line
    ax1 = axes[0, 0]
    celllines = combined_data['cellline'].unique()
    for cellline in celllines:
        cellline_data = combined_data[combined_data['cellline'] == cellline]
        ax1.hist(cellline_data['correlation'], bins=20, alpha=0.6, label=cellline)
    
    ax1.set_xlabel('Correlation with ' + morphogen)
    ax1.set_ylabel('Frequency')
    ax1.set_title('Correlation Distribution by Cell Line')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Number of responsive regulons per cell line
    ax2 = axes[0, 1]
    regulon_counts = combined_data['cellline'].value_counts()
    ax2.bar(regulon_counts.index, regulon_counts.values, color='skyblue')
    ax2.set_xlabel('Cell Line')
    ax2.set_ylabel('Number of Responsive Regulons')
    ax2.set_title(f'Responsive Regulons to {morphogen}')
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45)
    
    # 3. Consistency across cell lines (regulons appearing in multiple lines)
    ax3 = axes[1, 0]
    regulon_counts_per_cellline = combined_data['regulon'].value_counts()
    consistency_counts = regulon_counts_per_cellline.value_counts().sort_index()
    
    ax3.bar(consistency_counts.index, consistency_counts.values, color='lightcoral')
    ax3.set_xlabel('Number of Cell Lines')
    ax3.set_ylabel('Number of Regulons')
    ax3.set_title('Regulon Consistency Across Cell Lines')
    
    # 4. Correlation strength comparison
    ax4 = axes[1, 1]
    combined_data['abs_correlation'] = combined_data['correlation'].abs()
    sns.boxplot(data=combined_data, x='cellline', y='abs_correlation', ax=ax4)
    ax4.set_xlabel('Cell Line')
    ax4.set_ylabel('Absolute Correlation')
    ax4.set_title('Correlation Strength by Cell Line')
    plt.setp(ax4.xaxis.get_majorticklabels(), rotation=45)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved cell line comparison plot to {save_path}")
    
    return fig


def plot_morphogen_timeline(aucell_matrix: pd.DataFrame,
                          metadata: pd.DataFrame,
                          regulons_of_interest: List[str],
                          time_column: str = "Time",
                          morphogen_columns: List[str] = None,
                          figsize: Tuple[int, int] = (15, 10),
                          save_path: str = None) -> plt.Figure:
    """
    Create a timeline plot showing regulon activity during morphogen treatment.
    
    Parameters:
    -----------
    aucell_matrix : pd.DataFrame
        AUCell activity matrix (cells x regulons)
    metadata : pd.DataFrame
        Cell metadata including time and morphogen information
    regulons_of_interest : List[str]
        List of regulon names to plot
    time_column : str
        Name of the time column in metadata
    morphogen_columns : List[str], optional
        List of morphogen columns to include
    figsize : Tuple[int, int]
        Figure size
    save_path : str, optional
        Path to save the figure
        
    Returns:
    --------
    plt.Figure
        Matplotlib figure object
    """
    # Align data
    common_cells = aucell_matrix.index.intersection(metadata.index)
    aucell_aligned = aucell_matrix.loc[common_cells]
    metadata_aligned = metadata.loc[common_cells]
    
    # Filter to regulons of interest
    available_regulons = [r for r in regulons_of_interest if r in aucell_aligned.columns]
    
    if not available_regulons:
        logger.warning("No regulons of interest found in AUCell matrix")
        return None
    
    n_regulons = len(available_regulons)
    n_cols = min(3, n_regulons)
    n_rows = int(np.ceil(n_regulons / n_cols))
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False)
    
    # Get time points
    time_points = sorted(metadata_aligned[time_column].unique())
    time_numeric = {tp: i for i, tp in enumerate(time_points)}
    metadata_aligned['time_numeric'] = metadata_aligned[time_column].map(time_numeric)
    
    for idx, regulon in enumerate(available_regulons):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]
        
        # Prepare data for plotting
        plot_data = pd.DataFrame({
            'regulon_activity': aucell_aligned[regulon],
            'time': metadata_aligned['time_numeric'],
            'time_label': metadata_aligned[time_column]
        })
        
        # Add morphogen information if provided
        if morphogen_columns:
            for morphogen in morphogen_columns:
                if morphogen in metadata_aligned.columns:
                    plot_data[morphogen] = metadata_aligned[morphogen]
        
        # Create timeline plot
        for time_point in time_points:
            time_data = plot_data[plot_data['time_label'] == time_point]
            
            if len(time_data) > 0:
                # Plot mean and confidence interval
                mean_activity = time_data['regulon_activity'].mean()
                std_activity = time_data['regulon_activity'].std()
                
                ax.errorbar(time_numeric[time_point], mean_activity, 
                           yerr=std_activity, fmt='o-', capsize=5)
        
        ax.set_xlabel('Time')
        ax.set_ylabel('Regulon Activity (AUCell)')
        ax.set_title(f'{regulon}')
        ax.set_xticks(list(time_numeric.values()))
        ax.set_xticklabels(time_points, rotation=45)
        ax.grid(True, alpha=0.3)
    
    # Hide empty subplots
    for idx in range(n_regulons, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].set_visible(False)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved timeline plot to {save_path}")
    
    return fig


def create_summary_figure(correlation_df: pd.DataFrame,
                        p_value_df: pd.DataFrame,
                        responsive_regulons_df: pd.DataFrame,
                        morphogen_list: List[str],
                        figsize: Tuple[int, int] = (16, 12),
                        save_path: str = None) -> plt.Figure:
    """
    Create a comprehensive summary figure with multiple panels.
    
    Parameters:
    -----------
    correlation_df : pd.DataFrame
        Correlation matrix
    p_value_df : pd.DataFrame
        P-value matrix
    responsive_regulons_df : pd.DataFrame
        DataFrame of responsive regulons
    morphogen_list : List[str]
        List of morphogens
    figsize : Tuple[int, int]
        Figure size
    save_path : str, optional
        Path to save the figure
        
    Returns:
    --------
    plt.Figure
        Matplotlib figure object
    """
    fig = plt.figure(figsize=figsize)
    
    # Create a grid layout
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # Panel A: Correlation heatmap
    ax1 = fig.add_subplot(gs[0:2, 0:2])
    
    # Select top regulons for heatmap
    available_morphogens = [m for m in morphogen_list if m in correlation_df.columns]
    max_abs_corr = correlation_df[available_morphogens].abs().max(axis=1)
    top_regulons = max_abs_corr.nlargest(30).index
    
    plot_corr = correlation_df.loc[top_regulons, available_morphogens]
    plot_pval = p_value_df.loc[top_regulons, available_morphogens]
    
    # Create significance mask
    significance_mask = plot_pval < 0.05
    
    # Plot heatmap
    sns.heatmap(plot_corr, 
                annot=False,
                cmap='RdBu_r',
                center=0,
                vmin=-1, vmax=1,
                cbar_kws={'label': 'Correlation'},
                ax=ax1)
    
    # Add significance markers
    for i, regulon in enumerate(top_regulons):
        for j, morphogen in enumerate(available_morphogens):
            if significance_mask.loc[regulon, morphogen]:
                ax1.text(j + 0.5, i + 0.5, '*', 
                        ha='center', va='center', 
                        color='white', fontsize=6, fontweight='bold')
    
    ax1.set_title('A. Top Regulon-Morphogen Correlations', fontweight='bold')
    ax1.set_xlabel('')
    ax1.set_ylabel('')
    
    # Panel B: Number of responsive regulons per morphogen
    ax2 = fig.add_subplot(gs[0, 2])
    
    morphogen_counts = responsive_regulons_df['morphogen'].value_counts()
    morphogen_counts = morphogen_counts.reindex([m for m in morphogen_list if m in morphogen_counts.index])
    
    ax2.bar(range(len(morphogen_counts)), morphogen_counts.values, color='skyblue')
    ax2.set_xticks(range(len(morphogen_counts)))
    ax2.set_xticklabels(morphogen_counts.index, rotation=45, ha='right')
    ax2.set_ylabel('Number of\\nResponsive Regulons')
    ax2.set_title('B. Morphogen Responsiveness', fontweight='bold')
    
    # Panel C: Correlation strength distribution
    ax3 = fig.add_subplot(gs[1, 2])
    
    ax3.hist(responsive_regulons_df['abs_correlation'], bins=20, 
             color='lightcoral', alpha=0.7, edgecolor='black')
    ax3.set_xlabel('Absolute Correlation')
    ax3.set_ylabel('Frequency')
    ax3.set_title('C. Correlation Strength\\nDistribution', fontweight='bold')
    
    # Panel D: Top responsive regulons
    ax4 = fig.add_subplot(gs[2, :])
    
    # Get top 15 most responsive regulons
    top_responsive = responsive_regulons_df.nlargest(15, 'abs_correlation')
    
    # Create a horizontal bar plot
    y_pos = np.arange(len(top_responsive))
    correlations = top_responsive['correlation'].values
    colors = ['green' if c > 0 else 'red' for c in correlations]
    
    bars = ax4.barh(y_pos, np.abs(correlations), color=colors, alpha=0.7)
    
    # Add morphogen labels on bars
    for i, (idx, row) in enumerate(top_responsive.iterrows()):
        ax4.text(row['abs_correlation'] + 0.01, i, row['morphogen'], 
                va='center', fontsize=8)
    
    ax4.set_yticks(y_pos)
    ax4.set_yticklabels([r.split('(')[0] if '(' in r else r for r in top_responsive['regulon']], 
                       fontsize=9)
    ax4.set_xlabel('Absolute Correlation')
    ax4.set_title('D. Top 15 Most Responsive Regulons', fontweight='bold')
    ax4.invert_yaxis()
    
    # Add overall title
    fig.suptitle('Organoid Morphogen-Regulon Analysis Summary', 
                fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved summary figure to {save_path}")
    
    return fig
