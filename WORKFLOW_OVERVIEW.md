# Organoid Morphogen pySCENIC Pipeline - Complete Workflow

This pipeline implements the complete organoid morphogen analysis workflow with 4 main stages:

## Stage 1: pySCENIC Runs
- Multi-array SLURM submission with subsampling
- Individual cell line runs + combined stratified analysis
- Parameter combinations and seed variations

## Stage 2: Consensus Regulon Generation
- Combine multiple pySCENIC runs with occurrence/size thresholds
- Generate consensus regulons for each cell line and combined
- Create AUCell matrices for all cells

## Stage 3: Morphogen Regulon Networks
- Use GRNBoost2 to link morphogens/timing/medium to regulon activities
- Network inference from AUCell + metadata

## Stage 4: Final Analysis & Results
- Combine morphogen regulons across conditions
- Calculate correlations and statistical significance
- Generate final TF-target-morphogen correlation matrix

## Directory Structure

```
organoid_pyscenic_pipeline/
├── 01_pyscenic_runs/           # Stage 1: Multi-run pySCENIC
├── 02_consensus_regulons/      # Stage 2: Consensus generation
├── 03_morphogen_networks/      # Stage 3: GRNBoost morphogen links
├── 04_final_analysis/          # Stage 4: Combined results
├── scripts/                    # Execution scripts
├── config/                     # Configuration files
└── results/                    # Output data
```

Each stage contains clear notebooks and scripts for execution.
