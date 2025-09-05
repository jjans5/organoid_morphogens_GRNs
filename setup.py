#!/usr/bin/env python3
"""
Setup script for the Organoid pySCENIC Pipeline

This script helps users set up the pipeline environment and validate their setup.
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path
import shutil
import yaml

def check_conda():
    """Check if conda/mamba is available."""
    for cmd in ['mamba', 'conda']:
        if shutil.which(cmd):
            return cmd
    return None

def check_slurm():
    """Check if SLURM is available."""
    return shutil.which('sbatch') is not None

def create_environment(conda_cmd, env_file):
    """Create the conda environment."""
    print(f"Creating conda environment using {conda_cmd}...")
    try:
        subprocess.run([conda_cmd, 'env', 'create', '-f', env_file], check=True)
        print("✅ Environment created successfully!")
        return True
    except subprocess.CalledProcessError:
        print("❌ Failed to create environment")
        return False

def validate_databases(config_file):
    """Validate that pySCENIC databases are accessible."""
    print("Validating pySCENIC databases...")
    
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        
        databases = {
            'TF database': config.get('tf_database'),
            'Ranking database': config.get('ranking_database'),
            'Motif database': config.get('motif_database')
        }
        
        all_valid = True
        for name, path in databases.items():
            if path and os.path.exists(path):
                print(f"✅ {name}: {path}")
            else:
                print(f"❌ {name}: {path} (not found)")
                all_valid = False
        
        return all_valid
        
    except Exception as e:
        print(f"❌ Error reading database config: {e}")
        return False

def setup_data_directory(data_dir):
    """Set up the data directory structure."""
    data_path = Path(data_dir)
    data_path.mkdir(exist_ok=True)
    
    print(f"Data directory: {data_path.absolute()}")
    
    # Create subdirectories
    subdirs = ['raw', 'processed', 'external']
    for subdir in subdirs:
        (data_path / subdir).mkdir(exist_ok=True)
    
    return True

def main():
    parser = argparse.ArgumentParser(description="Setup Organoid pySCENIC Pipeline")
    
    parser.add_argument("--create-env", action="store_true",
                       help="Create conda environment")
    parser.add_argument("--validate-db", action="store_true",
                       help="Validate database paths")
    parser.add_argument("--setup-data", action="store_true",
                       help="Set up data directory structure")
    parser.add_argument("--all", action="store_true",
                       help="Run all setup steps")
    
    args = parser.parse_args()
    
    if not any(vars(args).values()):
        parser.print_help()
        return
    
    print("🧬 Organoid pySCENIC Pipeline Setup")
    print("=" * 50)
    
    # Check system requirements
    print("Checking system requirements...")
    
    conda_cmd = check_conda()
    if conda_cmd:
        print(f"✅ Conda/Mamba available: {conda_cmd}")
    else:
        print("❌ Conda/Mamba not found. Please install Miniconda or Anaconda.")
        return 1
    
    if check_slurm():
        print("✅ SLURM available")
    else:
        print("⚠️  SLURM not available (local execution only)")
    
    print()
    
    # Create environment
    if args.create_env or args.all:
        if not create_environment(conda_cmd, 'environment.yml'):
            return 1
        print()
    
    # Validate databases
    if args.validate_db or args.all:
        if not validate_databases('config/database_paths.yaml'):
            print("⚠️  Some databases not found. Please download pySCENIC databases.")
            print("See README.md for download instructions.")
        print()
    
    # Setup data directory
    if args.setup_data or args.all:
        setup_data_directory('data')
        print("✅ Data directory structure created")
        print()
    
    print("🎉 Setup complete!")
    print()
    print("Next steps:")
    print("1. Activate the environment: conda activate organoid_pyscenic")
    print("2. Download pySCENIC databases (see README.md)")
    print("3. Update config/database_paths.yaml with correct paths")
    print("4. Place your data in the data/ directory")
    print("5. Run the preprocessing notebook: notebooks/01_data_preprocessing.ipynb")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
