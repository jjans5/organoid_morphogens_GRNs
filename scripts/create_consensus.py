#!/usr/bin/env python3
"""
Consensus Regulon Generation Script for Organoid Study

This script generates consensus regulons from multiple pySCENIC runs,
following the methodology used in the organoid morphogen response analysis.

Usage:
    python create_consensus.py --sample consensus_0 --occur_threshold 5 --size_threshold 10

Author: Research Team
Date: 2025
"""

import sys
import os
import argparse
import logging
import glob
from pathlib import Path

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from consensus_regulons import ConsensusRegulonGenerator, process_multiple_samples, load_and_convert_motifs_to_regulons
import yaml

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def main():
    """
    Main function for consensus regulon generation.
    """
    parser = argparse.ArgumentParser(description='Generate consensus regulons for organoid study')
    parser.add_argument('--sample', type=str, required=True,
                       help='Sample name (e.g., consensus_0)')
    parser.add_argument('--occur_threshold', type=int, default=5,
                       help='Minimum occurrences for TF-gene links')
    parser.add_argument('--size_threshold', type=int, default=10,
                       help='Minimum regulon size threshold')
    parser.add_argument('--input_dir', type=str, default='../pyscenic/results',
                       help='Input directory with pySCENIC results')
    parser.add_argument('--output_dir', type=str, default='results/regulons',
                       help='Output directory for consensus regulons')
    parser.add_argument('--config', type=str, default='config/pyscenic_config.yaml',
                       help='Configuration file path')
    
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
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info(f"Starting consensus regulon generation for {args.sample}")
    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Occurrence threshold: {args.occur_threshold}")
    logger.info(f"Size threshold: {args.size_threshold}")
    
    try:
        # Step 1: Convert motif files to regulons if needed
        motif_files = list(input_dir.glob("*/enriched_motifs_seed*.csv"))
        if motif_files:
            logger.info(f"Found {len(motif_files)} motif files to convert")
            convert_dir = output_dir / "converted_regulons"
            load_and_convert_motifs_to_regulons(
                [str(f) for f in motif_files], 
                str(convert_dir)
            )
            regulon_input_dir = convert_dir
        else:
            # Look for existing regulon files
            regulon_input_dir = input_dir
        
        # Step 2: Find regulon files
        regulon_files = list(regulon_input_dir.glob("**/regulons_seed*.p"))
        if not regulon_files:
            logger.error("No regulon files found")
            return 1
        
        logger.info(f"Found {len(regulon_files)} regulon files")
        
        # Step 3: Generate consensus regulons
        generator = ConsensusRegulonGenerator(
            sample_name=args.sample,
            occur_threshold=args.occur_threshold,
            size_threshold=args.size_threshold
        )
        
        # Load and process regulon files
        generator.load_regulon_files([str(f) for f in regulon_files])
        generator.plot_regulon_distribution()
        generator.calculate_tf_gene_statistics()
        generator.plot_occurrence_weight_distribution()
        
        # Generate and save consensus regulons
        generator.generate_consensus_regulons()
        generator.save_consensus_regulons(str(output_dir))
        
        # Create analysis summary
        create_consensus_summary(generator, output_dir, args)
        
        logger.info(f"Consensus regulon generation completed for {args.sample}")
        return 0
        
    except Exception as e:
        logger.error(f"Consensus generation failed: {e}")
        import traceback
        traceback.print_exc()
        return 1


def create_consensus_summary(generator: ConsensusRegulonGenerator, 
                           output_dir: Path, args) -> None:
    """
    Create a summary report of the consensus regulon generation.
    
    Args:
        generator: ConsensusRegulonGenerator instance with results
        output_dir: Output directory for the summary
        args: Command line arguments
    """
    logger.info("Creating consensus generation summary")
    
    summary_lines = [
        f"# Consensus Regulon Generation Summary",
        f"",
        f"## Analysis Parameters",
        f"- Sample: {args.sample}",
        f"- Occurrence threshold: {args.occur_threshold}",
        f"- Size threshold: {args.size_threshold}",
        f"",
        f"## Input Data Summary",
    ]
    
    if generator.module_summary_df is not None:
        df = generator.module_summary_df
        summary_lines.extend([
            f"- Total TF-gene relationships: {len(df)}",
            f"- Unique TFs: {df['TF'].nunique()}",
            f"- Unique genes: {df['gene'].nunique()}",
            f"- Unique runs: {df['run'].nunique()}",
            f"",
        ])
    
    if generator.tf_gene_summary is not None:
        tf_gene_df = generator.tf_gene_summary
        summary_lines.extend([
            f"## TF-Gene Link Statistics",
            f"- Unique TF-gene links: {len(tf_gene_df)}",
            f"- Mean occurrences: {tf_gene_df['occurrences'].mean():.2f}",
            f"- Max occurrences: {tf_gene_df['occurrences'].max()}",
            f"- Mean weight: {tf_gene_df['mean_weight'].mean():.4f}",
            f"",
        ])
    
    if generator.consensus_regulons is not None:
        regulons = generator.consensus_regulons
        summary_lines.extend([
            f"## Consensus Regulons",
            f"- Total consensus regulons: {len(regulons)}",
            f"- Mean regulon size: {np.mean([len(r.genes) for r in regulons]):.1f}",
            f"- Largest regulon: {max([len(r.genes) for r in regulons])} genes",
            f"- Smallest regulon: {min([len(r.genes) for r in regulons])} genes",
            f"",
            f"## Top 10 Largest Regulons",
            f"",
        ])
        
        # Sort regulons by size
        regulons_sorted = sorted(regulons, key=lambda x: len(x.genes), reverse=True)[:10]
        
        summary_lines.append("| TF | Genes | Mean Weight |")
        summary_lines.append("|---|---|---|")
        
        for regulon in regulons_sorted:
            mean_weight = np.mean(regulon.weights) if hasattr(regulon, 'weights') and regulon.weights else 0
            summary_lines.append(f"| {regulon.name} | {len(regulon.genes)} | {mean_weight:.4f} |")
    
    # Save summary
    summary_file = output_dir / f"consensus_summary_{args.sample}.md"
    with open(summary_file, 'w') as f:
        f.write('\n'.join(summary_lines))
    
    logger.info(f"Consensus summary saved to {summary_file}")


def run_batch_consensus():
    """
    Run consensus regulon generation for multiple parameter combinations.
    """
    # Parameters based on original analysis
    samples = ['consensus_0']
    occur_thresholds = [0, 5, 10]
    size_thresholds = [5, 10, 20]
    
    logger.info("Starting batch consensus regulon generation")
    
    for sample in samples:
        for occur_thresh in occur_thresholds:
            for size_thresh in size_thresholds:
                logger.info(f"Processing {sample} with thresholds {occur_thresh}/{size_thresh}")
                
                # Set up command line arguments
                import sys
                sys.argv = [
                    'create_consensus.py',
                    '--sample', f"{sample}_occ{occur_thresh}_size{size_thresh}",
                    '--occur_threshold', str(occur_thresh),
                    '--size_threshold', str(size_thresh)
                ]
                
                try:
                    main()
                except Exception as e:
                    logger.error(f"Failed to process {sample}: {e}")
                    continue
    
    logger.info("Batch consensus generation completed")


if __name__ == "__main__":
    # Import numpy for summary statistics
    import numpy as np
    main()
