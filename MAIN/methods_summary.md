# Methods Summary — SpicyFlow Spatial Proteomics Pipeline

## Section 1: Environment and SpatialExperiment Setup

The analysis environment was initialized in R using the `renv` package for reproducible dependency management. The `SpatialExperiment` (Bioconductor) object served as the central data container throughout the pipeline, storing expression matrices, cell metadata, and spatial coordinates in a unified structure interoperable with all downstream Bioconductor spatial packages. Parallel processing was configured via `BiocParallel` using a `MulticoreParam` fork-based backend with 48 CPU cores. Java and `RBioFormats` were configured for OME-TIFF file handling using socket-based parallelization to maintain Java compatibility. A global random seed of 51773 was set for reproducibility across all stochastic operations. Visualization defaults were established with `ggplot2` using a classic theme.

**Key packages:** `SpatialExperiment`, `BiocParallel`, `renv`, `rJava`, `RBioFormats`, `ggplot2`

---

## Section 2: Cell Segmentation and Pre-Processing

Cell segmentation was performed using a GPU-accelerated Python pipeline that reads multiplexed OME-TIFF images directly via the `tifffile` library, preserving native float32 precision. Channel identity was determined by parsing OME-XML metadata embedded in the image files, and duplicate channels were removed automatically. A corrected FoxP3 channel (Corrected-FoxP3) was derived by subtracting a blurred and brightness-offset autofluorescence reference channel from the raw FoxP3 signal.

Nuclear segmentation was carried out using Cellpose v4 (CP-SAM model) on an Nvidia A6000 GPU with a batch size of 256, flow threshold of 0.4, cell probability threshold of 0.0, minimum cell size of 15 pixels, and test-time augmentation disabled for performance. Tile normalization with a block size of 100 was applied to the nuclear channel prior to segmentation. A foreground tissue mask was generated from the nuclear channel using Otsu thresholding on a downsampled image, followed by morphological binary dilation and small-hole removal, to exclude background regions from quantification.

Per-cell marker intensities were quantified by computing the mean intensity within each segmentation mask for every channel using `skimage.measure.regionprops`. Prior to quantification, each marker channel underwent adaptive enhancement consisting of percentile normalization to the [0, 1] range, outlier masking at the 98th percentile, a blended CLAHE equalization (70% original, 30% equalized), and gamma correction (gamma = 0.5). The nuclear marker used for segmentation was DNA-Brilliant Violet 421.

Quantified results were exported as CSV files containing per-cell spatial coordinates, cell area, and channel intensities. These were imported into R using `data.table::fread` and assembled into `SpatialExperiment` objects via the `import_python_results` function. Multiple images were batch-merged into a single combined `SpatialExperiment` with harmonized marker names across images.

**Key packages (Python):** `cellpose` (v4, CP-SAM), `torch`, `tifffile`, `numpy`, `pandas`, `scipy.ndimage`, `scikit-image` (exposure, filters, measure, morphology)

**Key packages (R):** `SpatialExperiment`, `data.table`, `qs`

---

## Section 3: Quality Control and Normalization

Quality control was performed to assess and correct technical variation between images. Raw marker intensity distributions were visualized as per-marker density plots across all images to identify batch effects. Dimensionality reduction via Uniform Manifold Approximation and Projection (UMAP) was performed on the raw expression data using `scater::runUMAP`, with DNA and isotype control markers excluded, to determine whether cells clustered by image identity (indicating batch effects) or by biological similarity.

Normalization was applied to the expression matrix using a three-step pipeline. First, the `trim99` step clipped each marker's values at the 99th percentile per image to reduce the influence of extreme outliers. Second, the `mean` step divided each marker's values by its per-image mean to center distributions across samples. Third, the `PC1` step regressed out the first principal component of variation within each image, which typically captures the dominant axis of technical batch effects while preserving biological signal. DNA markers were excluded from normalization. Post-normalization quality was confirmed by regenerating density plots and UMAP visualizations to verify improved cross-image alignment.

FoxP3-positive (FoxP3+) regulatory T cells were identified by thresholding the Corrected-FoxP3 channel at the 98th percentile of normalized expression, classifying the top 2% of cells as FoxP3+.

**Key packages:** `simpleSeg`, `scater`, `tidySingleCellExperiment`, `ggplot2`, `BiocParallel`

---

## Section 4: Unsupervised Clustering and Cell Type Annotation

Unsupervised clustering was used to assign cell type identities based on normalized marker expression profiles. The optimal number of clusters was estimated using `FuseSOM::estimateNumCluster`, which evaluates a range of cluster counts (k = 5 to 25) using multiple internal validation metrics: the Discriminant statistic, Gap statistic, Jump statistic, Slope statistic, Within-Cluster Dissimilarity (WCD), and Silhouette score. These metrics were visualized as elbow plots to identify the inflection point at which additional clusters yield diminishing returns.

Clustering was performed using the FuseSOM algorithm (`FuseSOM::runFuseSOM`), which combines a Self-Organizing Map (SOM) with hierarchical clustering. The final model was run with k = 15 clusters on the normalized expression assay. Marker expression profiles per cluster were visualized as heatmaps to guide biological interpretation.

Each of the 15 clusters was manually annotated with a cell type label based on its dominant marker expression pattern. Annotations included Transgenic Tregs (Corrected-FoxP3+, GFP+), Endogenous Tregs (Foxp3-PE+, GFP−), B Cells (CD19+, MHC-II+), Activated B Cells (CD40+, CD19+), Th2 T Cells (IL-4+), Exhausted T Cells (PD-1+, CD8a+), Activated CD8+ T Cells (Granzyme A+, CD107a+), IL10+ CD4+ T Cells (IL-10+, PD-1+), and Treg-APC Conjugates (CD40+, Corrected-FoxP3+), among others. Two clusters were identified as doublets/multiplets and one as low-quality artifacts based on multi-lineage marker co-expression or isotype control dominance; these were retained but flagged. Cell type composition was compared across images and conditions.

**Key packages:** `FuseSOM`, `MLmetrics`, `scuttle`, `scater`, `SummarizedExperiment`

---

## Section 5: Pairwise Cell Co-Localization Analysis

Spatial co-localization between cell types was quantified using the L-function, a transformation of Ripley's K-function, implemented in the `spicyR` package. For every pairwise combination of cell types, the L-function was calculated at radii of 20, 50, and 100 pixels using `spicyR::getPairwise`, measuring whether cells of type A are found closer to (positive L-value, co-localized) or further from (negative L-value, dispersed) cells of type B than expected under a random Poisson null model.

Co-localization patterns were compared between images using heatmaps of L-function values. To account for confounding due to tissue architecture, the analysis was repeated with an inhomogeneity correction (sigma = 20), which applies kernel density smoothing to adjust for uneven spatial cell density. The effect of this correction was compared against uncorrected results to identify relationships that were artifacts of tissue structure versus genuine biological interactions.

To enable formal statistical testing of differential co-localization between CAR-Treg and EGFR-Treg conditions despite limited biological replicates, a pseudotiling strategy was employed. Each image was divided into a 4 x 4 grid of spatial tiles (minimum 50 cells per tile), generating pseudo-replicates that provide estimates of within-image variance. Differential co-localization between conditions was then tested using `spicyR::spicy()`, which fits linear models comparing L-function values between conditions across pseudo-replicate tiles, yielding coefficients and p-values for each cell type pair.

FoxP3+ Treg-specific spatial relationships were extracted to determine which cell types preferentially co-localize with or are dispersed from regulatory T cells in each condition.

**Key packages:** `spicyR`, `Statial`, `imcRtools`, `dplyr`, `tibble`

---

## Section 7: Spatial Domain Detection and Boundary Analysis

Spatial domains — tissue compartments defined by coordinated cell type arrangements — were identified using the `lisaClust` package. The algorithm computes Local Indicators of Spatial Association (LISA) functions for each cell at multiple radii (20, 50, and 100 pixels) with Gaussian smoothing (sigma = 50), quantifying the degree to which each cell is surrounded by specific cell types relative to random expectation. These LISA vectors were then clustered across all cells to define k = 8 spatial domains, grouping tissue regions with similar local cell type compositions.

Cell type enrichment within each domain was assessed via bubble plots and heatmaps generated by `lisaClust::regionMap`. Based on enrichment patterns, the eight detected regions were manually annotated into six biologically meaningful spatial domains: B Cell Zone (enriched for B cells and activated B cells), T Cell Zone (enriched for CD8+ T cells), Treg Zone (enriched for Transgenic Tregs and Treg-APC conjugates), Th2 T Cell Zone (three regions merged, dominated by Th2 T cells), Resting Lymphoid Zone (predominantly resting/naive lymphocytes), and Low Intensity Zone (indeterminate and low-quality cells). Domain composition was compared between CAR-Treg and EGFR-Treg images, and FoxP3+ cell distribution across domains was assessed.

Domain boundary interfaces were quantified using Voronoi tessellation, computed with the `deldir` package. For each image, a Voronoi diagram was constructed from cell spatial coordinates, partitioning the tissue into polygonal tiles. Shared boundary segments between cells assigned to different spatial domains were identified and their lengths summed to produce total interface lengths for each domain pair. A minimum boundary segment length of 5 units was applied to exclude trivially short Voronoi edges between isolated cells. This analysis quantifies the physical extent of cross-compartment interfaces — for example, measuring the total boundary length between the Treg Zone and B Cell Zone — as a proxy for inter-domain signaling opportunity.

Domain-level co-localization analysis was performed by subsetting cells to individual domains and re-running `spicyR::spicy()` with pseudotiling within each domain to identify condition-specific co-localization patterns that are confined to particular tissue compartments.

**Key packages:** `lisaClust`, `spicyR`, `deldir`, `BiocParallel`, `RColorBrewer`, `ggplot2`
