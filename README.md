# Organoid Morphogen pySCENIC Analysis Pipeline

A comprehensive pipeline for analyzing the effects of morphogens on gene regulatory networks in organoid development using pySCENIC.

## Overview

This pipeline performs gene regulatory network (GRN) inference and analysis on single-cell RNA-seq data from organoids treated with different morphogens (FGF8, SHH, CHIR, RA) using the pySCENIC framework. The analysis includes:

1. **Data preprocessing** and quality control
2. **pySCENIC GRN inference** with consensus building across multiple runs
3. **Morphogen correlation analysis** to identify morphogen-responsive regulons
4. **Cross-cell line validation** across multiple human iPSC lines (H1, H9, WTC, WIBJ2)
5. **Downstream visualization** and functional analysis

## Repository Structure

```
organoid_pyscenic_pipeline/
├── README.md                 # This file
├── environment.yml           # Conda environment specification
├── requirements.txt          # Python package requirements
├── config/
│   ├── pyscenic_config.yaml  # pySCENIC configuration parameters
│   └── database_paths.yaml   # Paths to pySCENIC databases
├── src/
│   ├── pyscenic_utils.py     # Utility functions for pySCENIC
│   ├── morphogen_analysis.py # Morphogen correlation analysis functions
│   └── visualization.py     # Plotting and visualization functions
├── scripts/
│   ├── run_pyscenic.py       # Main pySCENIC execution script
│   ├── submit_pyscenic.sh    # SLURM job submission script
│   └── generate_combinations.py # Generate parameter combinations
├── notebooks/
│   ├── 01_data_preprocessing.ipynb    # Data preparation and QC
│   ├── 02_pyscenic_consensus.ipynb    # Consensus regulon building
│   ├── 03_morphogen_analysis.ipynb    # Morphogen correlation analysis
│   ├── 04_cross_cellline_validation.ipynb # Multi-cell line analysis
│   └── 05_visualization.ipynb         # Results visualization
├── data/                     # Input data directory
└── results/                  # Output results directory
```

## Installation

### Prerequisites

- Python 3.8+
- Conda or Mamba package manager
- Access to pySCENIC databases (see Database Setup section)

### Environment Setup

1. Clone this repository:
```bash
git clone <repository-url>
cd organoid_pyscenic_pipeline
```

2. Create the conda environment:
```bash
conda env create -f environment.yml
conda activate organoid_pyscenic
```

Alternatively, use pip:
```bash
pip install -r requirements.txt
```

### Database Setup

Download the required pySCENIC databases:

1. **Transcription factor list:**
   ```bash
   wget https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt
   ```

2. **Ranking databases:**
   ```bash
   wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/
   ```

3. **Motif annotations:**
   ```bash
   wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
   ```

Update the paths in `config/database_paths.yaml` accordingly.

## Usage

### Quick Start

1. **Prepare your data:**
   Place your single-cell count matrix (H5AD format) in the `data/` directory.

2. **Configure the analysis:**
   Edit `config/pyscenic_config.yaml` to set your analysis parameters.

3. **Run the preprocessing:**
   ```bash
   jupyter notebook notebooks/01_data_preprocessing.ipynb
   ```

4. **Execute pySCENIC:**
   ```bash
   python scripts/run_pyscenic.py --config config/pyscenic_config.yaml
   ```

   For SLURM clusters:
   ```bash
   sbatch scripts/submit_pyscenic.sh
   ```

5. **Analyze results:**
   Follow the analysis notebooks in order (02-05).

### Detailed Workflow

#### Step 1: Data Preprocessing
The preprocessing notebook (`01_data_preprocessing.ipynb`) performs:
- Quality control and filtering
- Cell type annotation validation
- Data splitting by cell line
- Preparation for pySCENIC input

#### Step 2: pySCENIC Execution
The pipeline runs pySCENIC with multiple parameter combinations to build consensus regulons:
- GRN inference using GRNBoost2
- Motif enrichment analysis
- AUCell scoring
- Consensus building across multiple runs

#### Step 3: Morphogen Analysis
Correlation analysis between regulon activities and morphogen treatments:
- Temporal correlation analysis
- Morphogen-specific regulon identification
- Statistical significance testing

#### Step 4: Cross-validation
Validation across multiple cell lines:
- Regulon reproducibility assessment
- Cell line-specific effects analysis
- Meta-analysis of morphogen responses

#### Step 5: Visualization
Comprehensive visualization of results:
- Regulatory network plots
- Morphogen response heatmaps
- Cell line comparison plots

## Configuration

### pySCENIC Parameters

Key parameters in `config/pyscenic_config.yaml`:

```yaml
# Data parameters
input_file: "data/exp1_counts_for_scenic_fil.h5ad"
subsample_column: "RNA_snn_res.1"
subsample_ncells: 600
stratification_column: "cell_line"

# pySCENIC parameters
n_workers: 30
memory_limit: "8e9"
method: "multiproc"  # or "dask"

# Consensus parameters
n_seeds: 10  # Number of random seeds for consensus
occurrence_threshold: 0  # Minimum occurrence across runs
```

### Database Paths

Update `config/database_paths.yaml` with your local paths:

```yaml
tf_database: "/path/to/allTFs_hg38.txt"
ranking_database: "/path/to/ranking_databases/"
motif_database: "/path/to/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
```

## Output Files

The pipeline generates several key output files:

- `results/consensus_regulons/`: Consensus regulons from multiple runs
- `results/aucell_scores/`: AUCell activity scores per cell
- `results/morphogen_correlations/`: Morphogen-regulon correlation matrices
- `results/cell_line_comparisons/`: Cross-cell line validation results
- `results/visualizations/`: Generated plots and figures

## Citation

If you use this pipeline in your research, please cite:

```
[Your paper citation here]
```

And the original pySCENIC publication:
```
Aibar, S., González-Blas, C.B., Moerman, T. et al. SCENIC: single-cell regulatory network inference and clustering. Nat Methods 14, 1083–1086 (2017). https://doi.org/10.1038/nmeth.4463
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

For questions and support, please contact:
- [Your name and email]
- [Lab/Institution information]

## Troubleshooting

### Common Issues

1. **Memory errors:** Reduce `n_workers` or `memory_limit` in the configuration
2. **Database access:** Ensure all pySCENIC databases are downloaded and paths are correct
3. **SLURM issues:** Adjust SLURM parameters in `submit_pyscenic.sh` for your cluster

### Performance Optimization

- For large datasets, consider using the Dask backend
- Adjust `subsample_ncells` based on your computational resources
- Use SSD storage for improved I/O performance

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
