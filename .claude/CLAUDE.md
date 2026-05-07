# SJ_04.2 — Spatial Proteomics Analysis Pipeline

## Project

Multiplex immunofluorescence (Cellscape platform) analysis of spleen tissue comparing **CAR-Treg** vs **EGFR-Treg** conditions. 26 markers, 4 images (~157K total cells). Primary research question: **Do FoxP3+ regulatory T cells differ spatially between conditions?**

## Role

Claude acts as **assistant bioinformatician** (R + Python expertise). Provide code edits and updates on request. **Do NOT run or test code** unless explicitly instructed — the user handles all execution.

## Hardware

- 72-core Intel Xeon CPU, 1TB RAM, Nvidia A6000 GPU
- WSL2 Linux environment
- Java heap configured at 768GB for RBioFormats

## Coding Standards

**See `MAIN/CONFIG/CODING_STANDARDS.md` for full details.** Key rules:

- **Fail-fast**: NO `tryCatch`, `try`, `suppressWarnings`, `suppressMessages` — let errors halt immediately
- **snake_case** for variables
- **Seed**: `set.seed(51773)` before random operations
- **renv** for reproducible package versions
- Fix root causes, never mask errors

## Pipeline Architecture

**Master file**: `MAIN/SpicyFlow.Rmd` (5,313 lines) — configures parameters, sources scripts from `MAIN/PROG/Scripts/`, generates results, saves checkpoints.

**See `MAIN/PIPELINE_GUIDE.md` for detailed section descriptions.**

### 9 Pipeline Sections

| Section | Purpose | Key Script(s) | Checkpoint |
|---------|---------|---------------|------------|
| 1 | Environment & setup | `setup_project.R`, `configure_parallel_processing.R`, `setup_java_rbioformats.R`, `install_packages.R` | — |
| 2 | Segmentation & pre-processing | `segment_and_quantify.py` (GPU/Cellpose), `import_python_results.R` | section2 |
| 3 | QC & normalization | `multi_image_qc.R` | section3 |
| 4 | Clustering & cell annotation | `clustering_analysis.R` (FuseSOM, k=15) | section4 |
| 5 | Pairwise co-localization | `spatial_colocalization.R` (spicyR L-function) | section5 |
| 6 | Context-aware analysis | Kontextual (Statial package) | section6 |
| 7 | Spatial domain detection | `spatial_domain_analysis.R` (lisaClust, k=6 domains) | section7 |
| 8 | Marker expression changes | SpatioMark (Statial package) | section8 |
| 9 | Feature summary & comparison | Inline in SpicyFlow.Rmd | — |

Checkpoints saved as `PROG/Checkpoints/combined_spe_after_section{N}.rds`.

### Key Parameters

- `n_cores = 48`, parallel mode = fork (socket for Java)
- `optimal_k = 15` clusters, `k_regions = 6` spatial domains
- `radii = c(20, 50, 100)` pixels, `sigma = 50`
- Nuclear marker: `"DNA-Brilliant Violet 421"`
- FoxP3+ threshold: 98th percentile

## Directory Structure

```
MAIN/
├── SpicyFlow.Rmd          # Master orchestration
├── PIPELINE_GUIDE.md       # Detailed section docs
├── CONFIG/                 # Coding standards, package lists
├── PROG/
│   ├── Scripts/            # 32 modular R/Python scripts
│   ├── Checkpoints/        # SPE objects per section (.rds)
│   ├── SEGMENTATION/       # Per-image segmentation outputs
│   │   ├── CAR_TREG_2/
│   │   └── EGFR_TREG_2/
│   ├── SingleChannelImages/ # Extracted channel TIFFs
│   └── tiles/              # Tiled processing intermediates
└── RES/                    # All results
    ├── QC_plots/
    ├── Clustering/
    ├── Colocalization/
    └── SpatialDomains/
```

## Core Packages

**R (Bioconductor):** SpatialExperiment, simpleSeg, FuseSOM, spicyR, Statial, lisaClust, EBImage, cytomapper, scater, ClassifyR, BiocParallel

**R (CRAN):** ggplot2, dplyr, data.table, qs, deldir, igraph

**Python:** cellpose, torch, tifffile, scikit-image, scipy, numpy, pandas

## Central Data Object

`SpatialExperiment` (SPE) — the unified container. All packages read/write to it:
- `assays`: expression matrices (raw + normalized)
- `colData`: cell metadata (cellType, clusters, region, annotation, imageID)
- `spatialCoords`: x, y positions
- `reducedDims`: UMAP, distances, abundances

## Images

| Image | Cells | Condition |
|-------|-------|-----------|
| CAR_TREG_1 | — | CAR-Treg |
| CAR_TREG_2 | 86,199 | CAR-Treg |
| EGFR_TREG_1 | — | EGFR-Treg |
| EGFR_TREG_2 | 71,187 | EGFR-Treg |
