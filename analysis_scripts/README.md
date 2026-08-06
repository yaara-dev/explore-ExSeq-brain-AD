# Analysis Scripts

This directory contains computational analysis pipelines for spatial transcriptomics data from the ExSeq Brain AD study. Each analysis pipeline is self-contained with its own documentation, requirements, and example data.

## Overview

The analysis scripts provide four complementary approaches to analyzing spatial gene expression patterns:

1. **Cell Typing** - Supervised STELLAR encoder label transfer for ExSeq cell-type annotation
2. **Expression, Spatial Distribution & Moran's I** - MATLAB scripts for differential expression, localization/spatial distribution, and Moran’s I
3. **RNA Velocity** - Analysis of RNA velocity and spatial cell-state transitions
4. **SVG Neighborhood Analysis** - Post-CELINA analysis of spatially variable genes and co-clustered genes

## Available Analyses

### 1. Cell Typing (`cell_typing/`)

**Purpose**: Transfer broad cell-type labels from an annotated spatial reference to unlabeled ExSeq tissue using a supervised STELLAR graph encoder.

**Key Features**:
- Supervised training of the STELLAR encoder on reference labels
- KNN/radius spatial graph construction
- Gene alignment between reference and target
- Diagnostic spatial plots and prediction summaries

**Language**: Python  
**Main Script**: `run_label_transfer.py`  
**Documentation**: See [cell_typing/README.md](cell_typing/README.md)

---

### 2. Expression, Spatial Distribution & Moran's I (`expression_spatial_distribution_morans_i/`)

**Purpose**: Compare WT vs 5xFAD spatial transcriptomics (ExSeq and Xenium) for differential expression by region, Moran’s I spatial autocorrelation, and ExSeq spatial distribution/localization on a registered hippocampus map.

**Key Features**:
- Region-wise differential expression between groups
- Moran’s I spatial autocorrelation comparisons
- ExSeq spatial distribution along a hippocampus reference map
- Separate ExSeq and Xenium entry-point scripts

**Language**: MATLAB  
**Main Scripts**: `CompareExpressionBetweenGroups_*.m`, `CompareMoranBetweenGroups_*.m`, `ExSeqSpatialDistributionXY.m`  
**Documentation**: See [expression_spatial_distribution_morans_i/README.md](expression_spatial_distribution_morans_i/README.md)

---

### 3. RNA Velocity (`rna_velocity/`)

**Purpose**: Analyze RNA velocity and spatial cell-state transitions to understand directional changes in gene expression.

**Key Features**:
- Per-cell normalization of spliced/unspliced RNA counts
- KNN pooling in expression/PCA space
- Gene-wise γ parameter estimation
- Distance-stratified analysis of cell-state changes
- Region-stratified analysis within anatomical regions
- Permutation testing for phase analysis

**Language**: Python  
**Main Script**: `src/rna_velocity_all_tissues.py`  
**Documentation**: See [rna_velocity/README.md](rna_velocity/README.md)

---

### 4. SVG Neighborhood Analysis (`svg_neighborhood_analysis/`)

**Purpose**: Post-CELINA analysis to identify and cluster spatially variable genes (SVGs) within spatial neighborhoods.

**Key Features**:
- FDR correction on CELINA p-values (multiple strategies)
- Cell-type-specific expression matrix creation
- Gene clustering to identify co-regulated modules
- K-means and hierarchical clustering methods
- Comprehensive visualizations and statistics

**Language**: Python  
**Main Scripts**: `scripts/apply_fdr_correction.py`, `scripts/create_celltype_matrices.py`, `scripts/cluster_significant_genes.py`  
**Documentation**: See [svg_neighborhood_analysis/README.md](svg_neighborhood_analysis/README.md)

---

## Quick Start

Each analysis has its own detailed README with installation instructions, input requirements, and usage examples. Navigate to the individual analysis directories for:

- Installation and requirements
- Input data format specifications
- Step-by-step usage instructions
- Example workflows
- Output format descriptions

## Common Requirements

Most Python-based analyses require:
- Python 3.8 or higher
- Common scientific Python packages (pandas, numpy, scipy, matplotlib, seaborn)
- See individual `requirements.txt` files for specific dependencies

The cell typing analysis requires:
- Python with the packages listed in `cell_typing/environment.yml`

Expression / spatial distribution / Moran’s I analysis requires:
- MATLAB (R2020b or later recommended; see folder README for toolboxes)

## Directory Structure

```
analysis_scripts/
├── README.md                    # This file
├── cell_typing/                 # Supervised STELLAR label transfer
│   ├── README.md
│   ├── run_label_transfer.py
│   ├── stellar_transfer/
│   ├── examples/
│   └── outputs/
├── expression_spatial_distribution_morans_i/  # Expression, localization & Moran’s I
│   ├── README.md
│   ├── CompareExpressionBetweenGroups_*.m
│   ├── CompareMoranBetweenGroups_*.m
│   ├── ExSeqSpatialDistributionXY.m
│   ├── data/
│   └── results/
├── rna_velocity/                # RNA velocity analysis
│   ├── README.md
│   ├── src/
│   ├── data/
│   └── RESULTS/
└── svg_neighborhood_analysis/   # SVG neighborhood analysis
    ├── README.md
    ├── scripts/
    ├── example_data/
    └── example_output/
```

## Getting Help

For detailed information about each analysis:
1. Navigate to the specific analysis directory
2. Read the individual README.md file
3. Check example data and output directories for format references
4. Review the source code comments for implementation details

