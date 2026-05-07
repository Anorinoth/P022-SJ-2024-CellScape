# SpicyFlow Pipeline: Step-by-Step Guide

## Overview

SpicyFlow is a spatial proteomics analysis pipeline for multiplexed imaging data (OME-TIFF format). It analyzes spleen tissue sections comparing **CAR-Treg** vs **EGFR-Treg** conditions, with a primary research question: **Do FoxP3+ regulatory T cells differ spatially between conditions?**

The pipeline flows through 9 sections, moving from raw images to comprehensive spatial features.

---

## Section 1: Environment & SpatialExperiment Setup

**What it does:** Initializes the R environment, configures parallel processing, installs packages, and sets up Java/RBioFormats for reading OME-TIFF files.

**Why it's important:** The `SpatialExperiment` (SPE) object is the central data structure used by the entire Bioconductor spatial ecosystem. Setting it up first ensures interoperability across all downstream packages (`simpleSeg`, `FuseSOM`, `spicyR`, `Statial`, `lisaClust`). Parallel processing (via `BiocParallel` with `MulticoreParam` or `SnowParam`) is configured here because many steps are computationally intensive. `renv` manages reproducible package versions.

**Key libraries:** `SpatialExperiment`, `BiocParallel`, `renv`, `rJava`, `RBioFormats`

**What it tells you:** Whether your environment is correctly configured and how many cores are available for computation.

---

## Section 2: Cell Segmentation & Pre-Processing

**What it does:** Runs a Python-based segmentation pipeline (`segment_and_quantify.py`) using **Cellpose v4 (CP-SAM)** on GPU. It reads pyramidal OME-TIFF files directly, discovers channels via OME-XML metadata, segments nuclei from a DNA marker channel, applies optional sharpening/denoising/tile-normalization, and quantifies per-cell marker intensities. Results are imported into R as `SpatialExperiment` objects and batch-merged across images.

**Why it's important:** Segmentation defines the fundamental unit of analysis -- the individual cell. Poor segmentation (over-segmentation splitting cells, under-segmentation merging them) cascades into every downstream result. Cellpose's deep-learning approach (CP-SAM variant) produces more accurate cell boundaries than traditional watershed methods, especially for densely packed tissue. Reading directly from pyramidal OME-TIFF preserves **float32** precision and avoids lossy channel-splitting.

**Key algorithms/libraries:** `Cellpose v4` (deep learning segmentation), `tifffile`/`bioformats` (OME-TIFF I/O), adaptive CLAHE enhancement, tile normalization

**What it tells you:** How many cells were detected per image, their median area, and which markers were quantified. The `cells.csv`, `masks.npz`, and `diagnostics.md` outputs let you assess segmentation quality.

---

## Section 3: Quality Control & Normalization

**What it does:** Assesses and corrects batch effects between images through:

1. **Marker density plots** -- raw intensity distributions per marker per image
2. **UMAP** on raw data -- checks if cells cluster by image (batch effect) or by biology
3. **Normalization** using `simpleSeg::normalizeCells()` with a three-step pipeline: `trim99` (clip top 1% outliers), `mean` (mean-center per image), `PC1` (regress out the first principal component of variation, which typically captures technical batch effects)
4. **Post-normalization QC** -- density plots and UMAP to confirm improvement
5. **FoxP3+ cell identification** using the 98th percentile threshold

**Why it's important:** Raw intensities aren't directly comparable across images due to staining efficiency, autofluorescence, imaging conditions, and sample preparation variability. Without normalization, a marker like CD3 could be bright in one image and dim in another despite identical biology, leading to misclassification. The PC1 correction is particularly powerful because it removes the dominant axis of technical variation while preserving biological signal -- unlike simple scaling, it accounts for correlated batch effects across markers.

**Key libraries:** `simpleSeg`, `scater` (UMAP via `runUMAP`), `tidySingleCellExperiment`

**What it tells you:** Whether batch effects exist (UMAP separation by image), how many markers improved after normalization, and an initial count of FoxP3+ Tregs per image.

---

## Section 4: Unsupervised Clustering & Cell Type Annotation

**What it does:**

1. **Estimates optimal cluster count** using multiple metrics: discriminant analysis, Gap statistic, Jump statistic, Slope statistic, Within-Cluster Dissimilarity, and Silhouette score (via `FuseSOM::estimateNumCluster`)
2. **Runs FuseSOM clustering** -- a Self-Organizing Map combined with hierarchical clustering
3. **Generates marker heatmaps** showing expression profiles per cluster
4. **Manual annotation** of clusters into cell types based on marker profiles (e.g., CD4+FoxP3+ = Tregs, CD19+MHC-II+ = B cells)
5. **Artifact removal** -- removes doublets/multiplets, low-quality clusters, and biologically implausible populations
6. **Composition comparison** between images and FoxP3+ enrichment analysis per cluster

**Why it's important:** Cell type identity is the foundation for all spatial analysis. FuseSOM is preferred over alternatives like FlowSOM for imaging data because it's designed for the lower-dimensional marker panels typical of spatial proteomics (10-40 markers vs thousands of genes in scRNA-seq). The elbow plot approach prevents over- or under-clustering, both of which distort results. Artifact removal prevents technical noise from biasing downstream spatial statistics.

**Key libraries:** `FuseSOM`, `MLmetrics`, `scuttle`, `scater`

**What it tells you:** What cell types are present in your tissue, their relative abundance, how they differ between conditions, and which clusters contain FoxP3+ Tregs.

---

## Section 5: Pairwise Cell Co-localization (spicyR)

**What it does:**

1. **L-function calculation** -- for every pair of cell types, quantifies whether they are found closer together (co-localized) or further apart (dispersed) than expected by random chance, at multiple radii (20, 50, 100 pixels)
2. **Image comparison** -- compares L-function values between CAR and EGFR images
3. **Co-localization heatmaps** -- shows attraction/repulsion for all cell type pairs
4. **Tissue inhomogeneity correction** -- re-runs analysis with a `sigma` smoothing parameter to account for uneven cell density
5. **Pseudotiling + spicyR::spicy()** -- creates pseudo-replicates by tiling images into sub-regions, enabling formal statistical testing (p-values) for differential co-localization between conditions
6. **FoxP3+ Treg spatial relationships** -- specifically examines which cell types co-localize with Tregs

**Why it's important:** Cell-cell proximity is a key indicator of biological interaction. Immune cells that co-localize with tumor cells may indicate active immune surveillance; exhausted T cells dispersed from antigen-presenting cells may indicate immune evasion. The L-function (derived from Ripley's K-function) is the gold standard for point-pattern spatial statistics. The inhomogeneity correction prevents false positives from tissue architecture -- e.g., two cell types might appear co-localized simply because they're both concentrated in the same tissue region, not because they interact.

**Key libraries:** `spicyR` (L-function, `getPairwise`, `spicy`), `imcRtools`

**What it tells you:** Which cell types attract or avoid each other, whether these relationships differ between CAR and EGFR conditions, and specifically how Tregs are positioned relative to other immune populations.

---

## Section 6: Context-Aware Spatial Analysis (Kontextual)

**What it does:**

1. **Defines biological hierarchies** -- groups cell types into lineages (B cells, T cells, Myeloid), functional subsets (Regulatory, Helper, Cytotoxic), and spatial compartments (Germinal Center, T Zone)
2. **Calculates Kontextual metrics** -- evaluates cell-cell spatial relationships *relative to a parent context* (e.g., "Are Tregs near B cells more than expected *within the T cell zone*?")
3. **Compares original vs. Kontextual** -- identifies context-dependent relationships where the conclusion changes when tissue structure is accounted for
4. **Kontextual curves** -- shows how context effects change across spatial scales (50-200 pixel radii)

**Why it's important:** Standard co-localization (Section 5) can be confounded by tissue architecture. For example, if Tregs and B cells are both in follicles, they'll appear co-localized even without direct interaction. Kontextual addresses this by asking: "Given the broader cellular neighborhood, are these two types closer than expected?" A relationship that flips sign between original and Kontextual analysis reveals that tissue structure, not biology, was driving the apparent co-localization. This is analogous to Simpson's paradox in statistics.

**Key libraries:** `Statial` (`Kontextual`, `kontextCurve`, `parentCombinations`)

**What it tells you:** Which spatial relationships are genuinely biological vs. artifacts of tissue structure. Context-dependent relationships (sign changes) are the most informative -- they reveal where tissue architecture creates misleading co-localization signals.

---

## Section 7: Spatial Domain Detection (lisaClust)

**What it does:**

1. **Runs lisaClust** -- identifies spatial domains (tissue compartments) by clustering cells based on their local indicators of spatial association (LISA) at multiple radii (20, 50, 100 pixels) with Gaussian smoothing (sigma=50)
2. **Examines cell type enrichment** per domain via bubble plots and heatmaps
3. **Annotates domains** with biological labels (B Cell Follicle, T Cell Zone, Treg Enriched Zone, etc.)
4. **Compares domain composition** between images
5. **FoxP3+ distribution** across spatial domains
6. **B cell sub-analysis** -- classifies naive vs. activated B cells (CD40 threshold) within domains
7. **Domain boundary analysis** -- quantifies shared boundaries between domains using Voronoi tessellation, specifically measuring Treg Zone <-> B Cell Follicle interface lengths with statistical testing (Wilcoxon)
8. **Domain-level spicyR** -- runs co-localization analysis within specific domains

**Why it's important:** While Sections 5-6 analyze pairwise cell relationships, lisaClust reveals *emergent tissue architecture* -- the higher-order organization of cells into functional compartments. In spleen, this captures germinal centers, follicle mantles, T zones, and interfollicular regions. These domains correspond to known microanatomical structures where specific immune functions occur. The boundary analysis is particularly novel -- it quantifies how much two domains physically interface, which may indicate cross-compartment signaling (e.g., Tregs at the B cell follicle boundary may suppress germinal center reactions).

**Key libraries:** `lisaClust` (LISA-based spatial clustering), `spicyR` (domain co-localization), `deldir` (Voronoi tessellation for boundaries)

**What it tells you:** What tissue compartments exist, which cell types define them, whether domains differ between conditions, where Tregs are concentrated, and how extensively different tissue compartments border each other.

---

## Section 8: Marker Expression Changes (SpatioMark)

**What it does:**

1. **Marker means by cell type and region** -- calculates average marker intensity for each cell type within each spatial domain (via `Statial::getMarkerMeans`)
2. **Distance and abundance calculation** -- computes distance from each cell to nearest cell of every other type (`getDistances`) and local K-function abundance (`getAbundances`)
3. **SpatioMark state changes** -- fits linear models testing whether marker expression in cell type A changes as a function of distance from cell type B (`calcStateChanges`). A negative coefficient means the marker increases when cells are closer together.
4. **Contamination detection** -- uses a random forest classifier (`calcContamination`) to estimate the probability each cell is correctly assigned to its type, detecting lateral marker spillover
5. **Cross-image comparison** -- identifies relationships that flip direction between conditions (sign changes)
6. **FoxP3+ Treg-specific state changes** -- which markers in Tregs are proximity-dependent?

**Why it's important:** Clustering assigns cells a fixed type, but cells exist on a continuum. SpatioMark captures *functional state changes* that are too subtle for clustering -- for instance, T cells that upregulate PD-1 only when near tumor cells, or Tregs that increase IL-10 only within germinal centers. These proximity-dependent expression shifts reveal the microenvironmental signals that drive cell behavior in situ. The contamination correction is critical for imaging data where optical spillover between adjacent cells can create false spatial correlations.

**Key libraries:** `Statial` (`getDistances`, `getAbundances`, `calcStateChanges`, `calcContamination`, `plotStateChanges`)

**What it tells you:** Which markers change based on spatial context, whether cells behave differently depending on their neighbors, which relationships flip between conditions, and whether any detected correlations might be spillover artifacts rather than real biology.

---

## Section 9: Comprehensive Feature Summary & Comparison

**What it does:**

1. **Consolidates all metrics** from Sections 3-8 into structured feature matrices:
   - Cell type proportions (from clustering)
   - Spatial domain proportions (from lisaClust)
   - Co-localization L-function values (from spicyR)
   - Marker means by cell type (from SpatioMark)
   - State change coefficients (from SpatioMark)
   - FoxP3+ Treg-specific features
2. **Ranks features** by discriminatory power between the two images
3. **Compares feature types** to determine which analysis category captures the most inter-image variation
4. **Exports** all feature matrices in RDS/CSV format
5. **Provides a ClassifyR template** for future predictive modeling when more images are available

**Why it's important:** This section transforms the pipeline from a descriptive tool into a quantitative framework. By reducing millions of single-cell measurements into a structured feature matrix (rows = images, columns = features), the data becomes compatible with machine learning approaches. With >=10 images and clinical metadata (treatment response, survival), `ClassifyR` can identify which spatial features predict outcomes. Even with 2 images, the feature comparison reveals what distinguishes the conditions most strongly, guiding hypothesis generation.

**Key libraries:** `ClassifyR` (template for future use), `dplyr`/`tidyr` (feature engineering)

**What it tells you:** Which individual features differ most between conditions, which analysis type (proportions, co-localization, domains, marker means, state changes) has the most discriminatory power, and provides a complete data package ready for expanded analysis with additional samples.

---

## Pipeline Architecture Notes

**Multi-scale spatial analysis:** The workflow examines tissue at progressively larger scales: single-cell markers (Sec 4) -> pairwise relationships (Sec 5) -> context-corrected relationships (Sec 6) -> tissue compartments (Sec 7) -> continuous state gradients (Sec 8) -> integrated features (Sec 9). Each scale reveals different biology.

**The Bioconductor spatial ecosystem:** The SPE object acts as a universal container -- `simpleSeg` writes normalized assays into it, `FuseSOM` adds cluster labels to colData, `spicyR` reads spatial coordinates from spatialCoords, `lisaClust` writes region assignments back to colData, and `Statial` adds distance/abundance matrices to reducedDims. This interoperability is why SPE was chosen as the central framework.

**Checkpoint system:** The pipeline saves RDS checkpoints after every major section. This is critical for spatial workflows because segmentation alone can take 30+ minutes and clustering/spatial analysis can be iterative. You can resume at any section without re-running upstream steps.

---

## Total Pipeline Output

- **30+ CSV/RDS data tables** with quantitative results
- **25+ high-quality plots** for publication/presentation
- **Comprehensive feature matrices** ready for machine learning
- **Complete reproducible workflow** for additional images
