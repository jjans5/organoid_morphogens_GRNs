#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --job-name=organoid_pyscenic
#SBATCH --output=logs/pyscenic_%j.out
#SBATCH --error=logs/pyscenic_%j.err
#SBATCH --partition=gpu  # Adjust to your cluster's partition

# Script for running the organoid pySCENIC pipeline on SLURM clusters

# Usage:
# sbatch submit_pyscenic.sh [config_file] [additional_args]

# Set default config file
CONFIG_FILE=${1:-"config/pyscenic_config.yaml"}
shift  # Remove first argument so $@ contains only additional arguments

# Create logs directory if it doesn't exist
mkdir -p logs

# Load necessary modules (adjust for your cluster)
module purge
module load gcc/11.2.0
module load python/3.9.7
module load hdf5/1.12.1

# Optional: Load additional modules for your cluster
# module load eth_proxy
# module load eccodes/2.25.0
# module load openblas/0.3.24

# Increase system limits
ulimit -u 32768
ulimit -n 65536

# Set environment variables for memory optimization
export MALLOC_ARENA_MAX=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Dask configuration
export DASK_DISTRIBUTED__COMM__MAX_MESSAGE_SIZE="25GB"
export DASK_DISTRIBUTED__WORKER__MEMORY__TARGET=0.8
export DASK_DISTRIBUTED__WORKER__MEMORY__SPILL=0.9

# Activate conda environment
echo "Activating conda environment..."
source ~/.bashrc

# Check if environment exists, create if not
if ! conda env list | grep -q "organoid_pyscenic"; then
    echo "Creating conda environment..."
    conda env create -f environment.yml
fi

conda activate organoid_pyscenic

# Verify environment is activated
if [[ "$CONDA_DEFAULT_ENV" != "organoid_pyscenic" ]]; then
    echo "ERROR: Failed to activate conda environment"
    exit 1
fi

echo "Running on node: $HOSTNAME"
echo "Job ID: $SLURM_JOB_ID"
echo "Number of CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE MB"
echo "Config file: $CONFIG_FILE"
echo "Additional arguments: $@"

# Change to script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "$SCRIPT_DIR/.."

# Check if config file exists
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Config file $CONFIG_FILE not found"
    exit 1
fi

# Run the pySCENIC pipeline
echo "Starting pySCENIC pipeline at $(date)"
python scripts/run_pyscenic.py --config "$CONFIG_FILE" "$@"

# Check exit status
if [[ $? -eq 0 ]]; then
    echo "pySCENIC pipeline completed successfully at $(date)"
else
    echo "ERROR: pySCENIC pipeline failed at $(date)"
    exit 1
fi

# Optional: Send notification (adjust email if needed)
# echo "pySCENIC job $SLURM_JOB_ID completed on $HOSTNAME" | mail -s "Job Complete" your.email@domain.com

echo "Job finished at $(date)"
