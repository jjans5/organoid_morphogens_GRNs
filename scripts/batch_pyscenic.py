#!/usr/bin/env python3
"""
Batch Job Generator for pySCENIC Pipeline

This script generates SLURM job scripts for running pySCENIC
with different parameter combinations across multiple cell lines.
"""

import argparse
import os
import itertools
from pathlib import Path
import yaml

def load_config(config_file):
    """Load configuration from YAML file."""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def generate_job_script(job_params, template_script="scripts/submit_pyscenic.sh"):
    """Generate a SLURM job script with specific parameters."""
    
    job_script_content = f"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={job_params.get('cpus', 32)}
#SBATCH --mem={job_params.get('memory', '120G')}
#SBATCH --time={job_params.get('time', '48:00:00')}
#SBATCH --job-name=pyscenic_{job_params.get('cell_line', 'unknown')}_{job_params.get('job_id', '0')}
#SBATCH --output=logs/pyscenic_{job_params.get('cell_line', 'unknown')}_{job_params.get('job_id', '0')}_%j.out
#SBATCH --error=logs/pyscenic_{job_params.get('cell_line', 'unknown')}_{job_params.get('job_id', '0')}_%j.err

# Job-specific parameters
CELL_LINE="{job_params.get('cell_line')}"
SEEDS="{' '.join(map(str, job_params.get('seeds', [42])))}"
CONFIG_FILE="{job_params.get('config_file', 'config/pyscenic_config.yaml')}"

# Create logs directory
mkdir -p logs

# Load modules (adjust for your cluster)
module purge
module load gcc/11.2.0
module load python/3.9.7
module load hdf5/1.12.1

# Set environment variables
export MALLOC_ARENA_MAX=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Activate conda environment
source ~/.bashrc
conda activate organoid_pyscenic

echo "Starting pySCENIC for cell line: $CELL_LINE"
echo "Seeds: $SEEDS"
echo "Config: $CONFIG_FILE"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $HOSTNAME"

# Run pySCENIC pipeline
python scripts/run_pyscenic.py \\
    --config "$CONFIG_FILE" \\
    --input "data/exp1_counts_for_scenic_${{CELL_LINE}}.h5ad" \\
    --output "results/cell_line_${{CELL_LINE}}" \\
    --seeds $SEEDS

exit_code=$?

if [[ $exit_code -eq 0 ]]; then
    echo "pySCENIC completed successfully for $CELL_LINE"
else
    echo "pySCENIC failed for $CELL_LINE with exit code $exit_code"
fi

exit $exit_code
"""
    
    return job_script_content

def main():
    parser = argparse.ArgumentParser(description="Generate batch jobs for pySCENIC pipeline")
    
    parser.add_argument("--config", default="config/pyscenic_config.yaml",
                       help="Configuration file")
    parser.add_argument("--cell-lines", nargs="+", 
                       default=["H1", "H9", "WTC", "WIBJ2"],
                       help="Cell lines to process")
    parser.add_argument("--seeds-per-job", type=int, default=5,
                       help="Number of seeds per job")
    parser.add_argument("--total-seeds", type=int, default=10,
                       help="Total number of seeds")
    parser.add_argument("--output-dir", default="batch_jobs",
                       help="Output directory for job scripts")
    parser.add_argument("--dry-run", action="store_true",
                       help="Generate scripts but don't submit")
    parser.add_argument("--submit", action="store_true",
                       help="Submit jobs to SLURM")
    
    args = parser.parse_args()
    
    # Load configuration
    config = load_config(args.config)
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Generate seed combinations
    all_seeds = list(range(1, args.total_seeds + 1))
    seed_chunks = [all_seeds[i:i + args.seeds_per_job] 
                   for i in range(0, len(all_seeds), args.seeds_per_job)]
    
    print(f"Generating jobs for {len(args.cell_lines)} cell lines")
    print(f"Seed chunks: {seed_chunks}")
    print(f"Total jobs: {len(args.cell_lines) * len(seed_chunks)}")
    
    job_scripts = []
    
    for cell_line in args.cell_lines:
        for job_idx, seeds in enumerate(seed_chunks):
            job_params = {
                'cell_line': cell_line,
                'seeds': seeds,
                'job_id': job_idx,
                'config_file': args.config,
                'cpus': config.get('slurm', {}).get('cpus_per_task', 32),
                'memory': config.get('slurm', {}).get('mem_per_cpu', '120G'),
                'time': config.get('slurm', {}).get('time', '48:00:00')
            }
            
            # Generate job script
            job_script_content = generate_job_script(job_params)
            
            # Save job script
            job_script_file = output_dir / f"pyscenic_{cell_line}_{job_idx}.sh"
            with open(job_script_file, 'w') as f:
                f.write(job_script_content)
            
            # Make executable
            os.chmod(job_script_file, 0o755)
            
            job_scripts.append(job_script_file)
            
            print(f"Generated: {job_script_file}")
    
    # Create submission script
    submission_script = output_dir / "submit_all.sh"
    with open(submission_script, 'w') as f:
        f.write("#!/bin/bash\\n")
        f.write("# Submit all pySCENIC jobs\\n\\n")
        
        for script in job_scripts:
            if args.submit:
                f.write(f"sbatch {script}\\n")
            else:
                f.write(f"echo sbatch {script}\\n")
    
    os.chmod(submission_script, 0o755)
    
    print(f"\\nSubmission script: {submission_script}")
    
    if args.submit and not args.dry_run:
        print("Submitting jobs...")
        os.system(f"bash {submission_script}")
    elif args.dry_run:
        print("Dry run - jobs not submitted")
        print(f"To submit jobs, run: bash {submission_script}")
    else:
        print(f"To submit jobs, run: bash {submission_script}")
    
    return 0

if __name__ == "__main__":
    main()
