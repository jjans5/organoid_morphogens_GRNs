#!/usr/bin/env python3
"""
Generate parameter combinations for pySCENIC batch jobs

This script generates parameter combinations for running pySCENIC
across multiple cell regions and random seeds.
"""

import argparse
import itertools
import os
from pathlib import Path


def generate_combinations(cell_regions: list, 
                        seeds: list, 
                        output_file: str = "region_seed_combos.txt"):
    """
    Generate all combinations of cell regions and random seeds.
    
    Parameters:
    -----------
    cell_regions : list
        List of cell region identifiers
    seeds : list
        List of random seeds
    output_file : str
        Output file path
    """
    combinations = list(itertools.product(cell_regions, seeds))
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file) if os.path.dirname(output_file) else '.', exist_ok=True)
    
    with open(output_file, 'w') as f:
        for region, seed in combinations:
            f.write(f"{region} {seed}\\n")
    
    print(f"Generated {len(combinations)} combinations")
    print(f"Saved to: {output_file}")
    
    return combinations


def main():
    parser = argparse.ArgumentParser(description="Generate parameter combinations for pySCENIC batch jobs")
    
    parser.add_argument("--regions", nargs="+", 
                       default=["all", "neurons", "progenitors"],
                       help="List of cell regions to analyze")
    
    parser.add_argument("--seeds", nargs="+", type=int,
                       default=list(range(1, 11)),  # Seeds 1-10
                       help="List of random seeds")
    
    parser.add_argument("--output", default="region_seed_combos.txt",
                       help="Output file path")
    
    parser.add_argument("--custom_regions", 
                       help="Path to file containing custom region list (one per line)")
    
    args = parser.parse_args()
    
    # Load custom regions if provided
    if args.custom_regions:
        with open(args.custom_regions, 'r') as f:
            regions = [line.strip() for line in f if line.strip()]
    else:
        regions = args.regions
    
    print(f"Regions: {regions}")
    print(f"Seeds: {args.seeds}")
    
    combinations = generate_combinations(regions, args.seeds, args.output)
    
    # Print some example combinations
    print("\\nFirst 5 combinations:")
    for i, (region, seed) in enumerate(combinations[:5]):
        print(f"  {i+1}: {region} {seed}")
    
    if len(combinations) > 5:
        print(f"  ... and {len(combinations) - 5} more")


if __name__ == "__main__":
    main()
