# SpicyFlow: Spatial Proteomics Analysis of cGVHD Spleen Tissue

**Author:** Dr. Michael Zaiken

A comprehensive spatial proteomics analysis pipeline for multiplexed imaging (CellScape/OME-TIFF) data. This pipeline analyzes spleen tissue sections from a chronic graft-versus-host disease (cGVHD) model comparing CAR-Treg versus EGFR-Treg conditions, with a focus on FoxP3+ regulatory T cell spatial organization.

This pipeline is an implementation of the [Spatial Analysis Playbook](https://sydneybiox.github.io/spatialPlaybook/) published by the Sydney Precision Data Science Centre at the University of Sydney, adapted for CellScape multiplex immunofluorescence data.

## Analysis Pipeline

The analysis is implemented in `MAIN/SpicyFlow.Rmd` and flows through six sequential sections, each building on the previous:

| Section | Title | Key Methods |
|---------|-------|-------------|
| 1 | Environment & SpatialExperiment Setup | `renv`, `BiocParallel`, `RBioFormats` |
| 2 | Cell Segmentation & Pre-processing | Cellpose v4 (CP-SAM), GPU-accelerated Python pipeline |
| 3 | Quality Control & Normalization | `simpleSeg` (trim99 + mean + PC1 correction), UMAP |
| 4 | Clustering & Cell Type Annotation | `FuseSOM` unsupervised clustering, manual annotation |
| 5 | Pairwise Cell Co-localization | `spicyR` L-function, pseudotiling for statistical testing |
| 6 | Spatial Domain Detection | `lisaClust` LISA-based domains, Voronoi boundary analysis |

An appendix contains additional exploratory visualizations not used in the final manuscript.

## Methods

**Segmentation**: Nuclear segmentation via Cellpose v4 (CP-SAM model) on GPU with batch size 256, flow threshold 0.4, and cell probability threshold 0.0. Per-cell marker intensities quantified as mean intensity within segmentation masks after adaptive enhancement (percentile normalization, CLAHE, gamma correction).

**Quality Control**: Three-step normalization pipeline — 99th percentile clipping, per-image mean centering, and PC1 regression to remove dominant technical batch effects while preserving biological signal.

**Clustering**: FuseSOM algorithm (Self-Organizing Map + hierarchical clustering) with optimal k estimated via multiple internal validation metrics (Discriminant, Gap, Jump, Slope, WCD, Silhouette). Clusters manually annotated with cell type labels based on dominant marker expression.

**Spatial Analysis**: Pairwise co-localization quantified via L-function (transformation of Ripley's K) at radii of 20, 50, and 100 pixels. Differential co-localization tested using pseudotiling (4x4 grid) with `spicyR::spicy()` linear models. Spatial domains identified by `lisaClust` using Local Indicators of Spatial Association (LISA) functions clustered across cells. Domain boundaries quantified via Voronoi tessellation (`deldir`).

## Prerequisites

- **R** 4.4.3 (pinned in `.Rversion`; install via [rig](https://github.com/r-lib/rig))
- **Python** 3.x (for the segmentation pipeline)
- **Java** JDK (for `RBioFormats` OME-TIFF handling)
- **GPU** with CUDA support (recommended for Cellpose segmentation)
- **QuPath** v0.5.1+ (optional, for image visualization and validation)

## Setup

### 1. Clone the repository

```bash
git clone https://github.com/Anorinoth/P022-SJ-2024-CellScape.git
cd P022-SJ-2024-CellScape
```

### 2. Configure environment

Copy the environment template and set paths for your system:

```bash
cp .Renviron.example .Renviron
```

Edit `.Renviron` to set the required project directory paths:

```
P022_CODE_DIR=/path/to/this/repo
P022_DATA_DIR=/path/to/Data
P022_OUTPUTS_DIR=/path/to/Outputs
P022_REPORTS_DIR=/path/to/Reports
```

### 3. R environment

```r
# renv will bootstrap automatically on first session start.
# Restore all pinned package versions:
renv::restore()
```

### 4. Python dependencies

```bash
pip install -r MAIN/PROG/Scripts/requirements.txt
# For GPU segmentation, install PyTorch with CUDA support:
# https://pytorch.org/get-started/locally/
```

### 5. External tools

- **Bio-Formats CLI** (`bftools`): download from [openmicroscopy.org](https://www.openmicroscopy.org/bio-formats/downloads/) and place in `_vendored/bftools/`
- **ome-tiff-pyramid-tools**: install from [labsyspharm/ome-tiff-pyramid-tools](https://github.com/labsyspharm/ome-tiff-pyramid-tools) into `_vendored/bioio_env/`

## Directory Structure

```
Code/
├── MAIN/
│   ├── SpicyFlow.Rmd              # Main analysis pipeline (6 sections + appendix)
│   ├── CONFIG/                    # Package lists (CRAN, Bioconductor, GitHub)
│   └── PROG/Scripts/              # R and Python helper scripts
│       ├── setup_project.R        # Project environment configuration
│       ├── clustering_analysis.R  # FuseSOM clustering utilities
│       ├── spatial_colocalization.R  # spicyR analysis functions
│       ├── spatial_domain_analysis.R # lisaClust domain detection
│       ├── segment_and_quantify.py   # Cellpose GPU segmentation
│       └── ...                    # Additional utility scripts
├── bioio_correction.py            # OME-TIFF bleedthrough correction pipeline
├── .Renviron.example              # Environment variable template
├── .Rprofile                      # R startup configuration
├── .Rversion                      # Pinned R version (4.4.3)
├── renv.lock                      # Reproducible R package versions
└── README.md
```

## Running the Analysis

1. Open the project in RStudio or start R from the repository root
2. Open `MAIN/SpicyFlow.Rmd`
3. Run `setup_project()` in the first code chunk to initialize paths
4. Execute sections sequentially — each saves an RDS checkpoint to enable resuming from any section
5. Checkpoints are saved to `DATA_DIR/PROG/Checkpoints/` and named by section (e.g., `combined_spe_after_section5.rds`)

## Key R Packages

| Package | Purpose |
|---------|---------|
| `SpatialExperiment` | Central data container for spatial proteomics |
| `simpleSeg` | Cell intensity normalization |
| `FuseSOM` | SOM + hierarchical clustering for cell typing |
| `spicyR` | Pairwise spatial co-localization (L-function) |
| `lisaClust` | Spatial domain detection via LISA functions |
| `Statial` | Spatial statistics framework |
| `scater` | Dimensionality reduction (UMAP) |
| `deldir` | Voronoi tessellation for boundary analysis |
| `BiocParallel` | Parallel processing backend |
| `RBioFormats` | OME-TIFF file I/O |

## Citation

Publication pending.

## License

This project is shared for academic and research purposes. Please cite the associated publication when using this code.
