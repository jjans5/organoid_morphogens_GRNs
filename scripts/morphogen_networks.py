#!/usr/bin/env python3
"""
Morphogen Network Analysis Script for Organoid Study

This script runs the complete morphogen-regulon network analysis pipeline,
combining AUCell regulon activities with morphogen treatment metadata
to identify regulatory relationships using GRNBoost2.

Usage:
    python morphogen_networks.py --cell_line H1 --occur_threshold 0

Author: Research Team
Date: 2025
"""

import sys
import os
import argparse
import logging
from pathlib import Path

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from grnboost_analysis import run_morphogen_network_analysis, MorphogenNetworkAnalyzer
import yaml

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def main():
    """
    Main function for morphogen network analysis.
    """
    parser = argparse.ArgumentParser(description='Morphogen network analysis for organoid study')
    parser.add_argument('--cell_line', type=str, required=True,
                       choices=['H1', 'H9', 'WIBJ2', 'WTC'],
                       help='Cell line to analyze')
    parser.add_argument('--occur_threshold', type=int, default=0,
                       help='Occurrence threshold for consensus regulons')
    parser.add_argument('--config', type=str, default='config/pyscenic_config.yaml',
                       help='Configuration file path')
    parser.add_argument('--n_workers', type=int, default=10,
                       help='Number of workers for GRNBoost2')
    parser.add_argument('--output_dir', type=str, default=None,
                       help='Output directory (default: results/morphogen_networks/{cell_line})')
    
    args = parser.parse_args()
    
    # Load configuration
    config_path = Path(args.config)
    if config_path.exists():
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
    else:
        logger.warning(f"Config file {config_path} not found, using defaults")
        config = {}
    
    # Set up paths
    base_dir = Path.cwd()
    data_dir = base_dir / "data"
    results_dir = base_dir / "results"
    
    # Set output directory
    if args.output_dir is None:
        output_dir = results_dir / "morphogen_networks" / args.cell_line
    else:
        output_dir = Path(args.output_dir)
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Define file paths based on original analysis structure
    aucell_file = f"../pyscenic_post/regulons/consensus_0/aucell_{args.cell_line}_combined.p"
    metadata_file = f"../data/exp1_counts_for_scenic_{args.cell_line}.h5ad"
    
    # Check if files exist
    if not os.path.exists(aucell_file):
        logger.error(f"AUCell file not found: {aucell_file}")
        return 1
    
    if not os.path.exists(metadata_file):
        logger.error(f"Metadata file not found: {metadata_file}")
        return 1
    
    logger.info(f"Starting morphogen network analysis for {args.cell_line}")
    logger.info(f"AUCell file: {aucell_file}")
    logger.info(f"Metadata file: {metadata_file}")
    logger.info(f"Output directory: {output_dir}")
    
    try:
        # Run the complete analysis
        analyzer = run_morphogen_network_analysis(
            cell_line=args.cell_line,
            aucell_file=aucell_file,
            metadata_file=metadata_file,
            output_dir=str(output_dir),
            occur_threshold=args.occur_threshold,
            n_workers=args.n_workers
        )
        
        # Generate additional analysis summary
        create_analysis_summary(analyzer, output_dir)
        
        logger.info(f"Morphogen network analysis completed successfully for {args.cell_line}")
        return 0
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return 1


def create_analysis_summary(analyzer: MorphogenNetworkAnalyzer, output_dir: Path) -> None:
    """
    Create a summary report of the morphogen network analysis.
    
    Args:
        analyzer: MorphogenNetworkAnalyzer instance with results
        output_dir: Output directory for the summary
    """
    logger.info("Creating analysis summary")
    
    summary_lines = [
        f"# Morphogen Network Analysis Summary",
        f"",
        f"## Analysis Parameters",
        f"- Cell line: {analyzer.cell_line}",
        f"- Occurrence threshold: {analyzer.occur_threshold}",
        f"",
        f"## Data Summary",
        f"- AUCell matrix shape: {analyzer.aucell_matrix.shape if analyzer.aucell_matrix is not None else 'N/A'}",
        f"- Morphogen variables: {list(analyzer.metadata.columns) if analyzer.metadata is not None else 'N/A'}",
        f"- Combined matrix shape: {analyzer.morphogen_matrix.shape if analyzer.morphogen_matrix is not None else 'N/A'}",
        f"",
        f"## Network Results",
    ]
    
    if analyzer.network_adjacencies is not None:
        net_adj = analyzer.network_adjacencies
        summary_lines.extend([
            f"- Total network edges: {len(net_adj)}",
            f"- Unique TFs: {net_adj['TF'].nunique()}",
            f"- Unique targets: {net_adj['target'].nunique()}",
            f"- Mean importance: {net_adj['importance'].mean():.4f}",
            f"- Max importance: {net_adj['importance'].max():.4f}",
            f"",
            f"## Top 10 Network Edges",
            f"",
        ])
        
        # Add top edges table
        top_edges = net_adj.nlargest(10, 'importance')[['TF', 'target', 'importance']]
        summary_lines.append("| TF | Target | Importance |")
        summary_lines.append("|---|---|---|")
        
        for _, row in top_edges.iterrows():
            summary_lines.append(f"| {row['TF']} | {row['target']} | {row['importance']:.4f} |")
    
    else:
        summary_lines.append("- Network adjacencies: Not available")
    
    # Save summary
    summary_file = output_dir / "analysis_summary.md"
    with open(summary_file, 'w') as f:
        f.write('\n'.join(summary_lines))
    
    logger.info(f"Analysis summary saved to {summary_file}")


def run_batch_analysis():
    """
    Run morphogen network analysis for all cell lines.
    """
    cell_lines = ['H1', 'H9', 'WIBJ2', 'WTC']
    occur_thresholds = [0]
    
    logger.info("Starting batch morphogen network analysis")
    
    for cell_line in cell_lines:
        for occur_threshold in occur_thresholds:
            logger.info(f"Processing {cell_line} with threshold {occur_threshold}")
            
            # Set up command line arguments
            import sys
            sys.argv = [
                'morphogen_networks.py',
                '--cell_line', cell_line,
                '--occur_threshold', str(occur_threshold),
                '--n_workers', '10'
            ]
            
            try:
                main()
            except Exception as e:
                logger.error(f"Failed to process {cell_line}: {e}")
                continue
    
    logger.info("Batch analysis completed")


if __name__ == "__main__":
    main()
