# Organoid Morphogen pySCENIC Analysis Pipeline

This repository contains the computational pipeline used for the gene regulatory network analysis of human brain organoid development under different morphogen conditions. This analysis was conducted for our research on morphogen-driven transcriptional regulation in organoid neurogenesis.

## Overview

This pipeline processes single-cell RNA-seq data from human brain organoids treated with different morphogen combinations (FGF8, SHH, CHIR, RA) across multiple time points. The analysis uses pySCENIC to identify gene regulatory networks and their responses to morphogen treatments.

## Key Analysis Components

1. **Consensus Regulon Generation**: Combines multiple pySCENIC runs to create robust consensus regulons
2. **AUCell Activity Calculation**: Computes regulon activity scores for each cell
3. **Morphogen-Regulon Network Inference**: Uses GRNBoost2 to identify regulatory relationships between morphogens and transcription factors
4. **Correlation Analysis**: Quantifies morphogen-regulon associations across experimental conditions

## Pipeline Structure

```
organoid_pyscenic_pipeline/
├── src/                     # Core analysis modules
│   ├── pyscenic_utils.py    # pySCENIC utilities and consensus generation
│   ├── morphogen_analysis.py # Morphogen-regulon correlation analysis
│   ├── grnboost_analysis.py  # GRNBoost2 network inference
│   └── visualization.py     # Plotting and visualization functions
├── scripts/                 # Execution scripts
│   ├── run_pyscenic.py      # Main pySCENIC execution
│   ├── create_consensus.py  # Consensus regulon generation
│   ├── morphogen_networks.py # Morphogen-regulon network analysis
│   └── submit_pyscenic.sh   # SLURM submission script
├── notebooks/               # Analysis notebooks
│   ├── 01_data_preprocessing.ipynb
│   ├── 02_consensus_regulons.ipynb
│   ├── 03_morphogen_analysis.ipynb
│   └── 04_network_visualization.ipynb
├── config/                  # Configuration files
│   ├── pyscenic_config.yaml
│   └── database_paths.yaml
└── results/                 # Output directory
```

## Analysis Workflow

### 1. Data Preprocessing
- Load single-cell RNA-seq data for each cell line (H1, H9, WIBJ2, WTC)
- Filter cells and genes according to quality metrics
- Prepare morphogen treatment metadata

### 2. pySCENIC Analysis
- Run multiple pySCENIC iterations with different parameters
- Generate adjacency matrices, modules, and regulons
- Calculate AUCell activity scores

### 3. Consensus Regulon Generation
- Combine results from multiple pySCENIC runs
- Apply occurrence and size thresholds
- Generate robust consensus regulons

### 4. Morphogen-Regulon Network Analysis
- Combine AUCell matrix with morphogen metadata
- Run GRNBoost2 to infer morphogen-regulon networks
- Identify key regulatory relationships

### 5. Correlation Analysis & Visualization
- Calculate Pearson correlations between morphogens and regulons
- Generate heatmaps and network visualizations
- Statistical significance testing
