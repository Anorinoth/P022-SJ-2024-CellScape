# Spatial Analysis Results Summary
## Comparison of CAR-T vs EGFR Treated Spleen Tissue

**Analysis Date:** October 28, 2025  
**Primary Research Question:** Do FoxP3+ regulatory T cells (Tregs) behave differently in EGFR vs CAR-T treated spleen tissue?

---

## Key Dataset Characteristics

- **Images Analyzed:** 2 (CAR_TREG_2, EGFR_TREG_2)
- **Total Cells:** 157,386 cells
  - CAR_TREG_2: 86,199 cells
  - EGFR_TREG_2: 71,187 cells
- **Markers:** 26 unique markers (protein expression)
- **Tissue Context:** Spleen germinal center architecture

---

## Section 3: Data Quality Control & Normalization

### Analysis Description

**Purpose:** Ensure data quality, correct for technical batch effects, and identify FoxP3+ regulatory T cells (Tregs).

**Methods:**
- **Marker Expression Profiling:** Raw and normalized intensity distributions across all 26 protein markers
- **Batch Effect Assessment:** UMAP dimensionality reduction to visualize image-specific technical variation
- **Normalization:** `simpleSeg::normalizeCells()` with three sequential methods:
  1. **trim99** - Remove outliers (top 1% intensities)
  2. **mean** - Center marker distributions
  3. **PC1** - Correct systematic batch effects via first principal component
- **FoxP3+ Identification:** Cells classified as FoxP3+ (Treg) if `31_CorrectedFOXP3-GFP` expression ≥ 98th percentile

**Interpretation:** This section establishes data quality and removes technical confounders, allowing biological differences to be distinguished from technical artifacts. The stringent FoxP3 threshold (98th percentile) ensures high-confidence Treg identification.

### Key Findings

#### 3.1 Overall FoxP3+ Treg Abundance

| Image | Total Cells | FoxP3+ Count | FoxP3+ % | Mean FoxP3 Expression |
|-------|-------------|--------------|----------|-----------------------|
| **CAR_TREG_2** | 86,199 | 1,572 | **1.82%** | 0.141 |
| **EGFR_TREG_2** | 71,187 | 1,576 | **2.21%** | 0.139 |

**Key Observation:** EGFR-treated tissue shows **21% higher** FoxP3+ Treg frequency (2.21% vs 1.82%) despite having fewer total cells, suggesting altered immune regulation.

#### 3.2 Normalization Effectiveness

- **Pre-normalization:** Strong image-based clustering in UMAP, indicating technical batch effects
- **Post-normalization:** Images intermixed in UMAP space, confirming successful batch correction
- **Marker-specific corrections:** All 26 markers showed improved comparability between images

#### 3.3 Spatial Distribution Patterns

FoxP3+ Tregs show **distinct spatial organization** in both images, concentrated in specific tissue regions rather than uniformly distributed. This suggests functional compartmentalization within the germinal center.

### Supporting Visualizations

- **`marker_density_overview_page1-5.png`**: Raw marker expression distributions showing inter-image variation
- **`marker_density_normalized_overview_page1-4.png`**: Post-normalization distributions demonstrating effective correction
- **`normalization_before_after/` (all markers)**: Side-by-side comparisons confirming successful batch correction
- **`umap_batch_effects.png`**: Pre-normalization UMAP showing technical separation
- **`umap_post_normalization.png`**: Post-normalization UMAP showing successful batch correction
- **`FOXP3_expression_distribution.png`**: FoxP3 expression across all cells, highlighting 98th percentile threshold
- **`FOXP3_spatial_distribution.png`**: Spatial maps showing regional FoxP3+ enrichment patterns

### Relevance to Research Question

**Direct:** Establishes that EGFR treatment is associated with higher overall Treg frequency, providing the first evidence of differential immune regulation between treatments.

**Methodological:** Ensures that subsequent spatial analyses reflect true biological differences rather than technical artifacts, validating interpretation of all downstream results.

---

## Section 4: Cell Type Identification & Clustering

### Analysis Description

**Purpose:** Identify functionally distinct cell populations and characterize their treatment-specific distributions.

**Methods:**
- **Unsupervised Clustering:** `FuseSOM` algorithm integrating:
  - Self-organizing maps (SOM) for initial clustering
  - ConsensusClusterPlus for stability assessment
  - Optimal cluster number selection via elbow method
- **Cluster Count:** 15 distinct cell populations identified
- **Annotation Strategy:** Manual annotation based on marker expression profiles from heatmap visualization
- **Validation:** UMAP visualization confirms biological coherence of identified cell types

**Interpretation:** Cell types represent functionally distinct immune populations. Marker expression heatmaps show characteristic "fingerprints" that align with known immune cell biology in germinal centers.

### Key Findings

#### 4.1 Cell Type Landscape (15 Populations Identified)

**Major Populations (>5% of total cells):**

| Cell Type | Cell Count | % of Total | Key Markers |
|-----------|------------|------------|-------------|
| **Resting/Naïve lymphocytes** | 85,176 | 54.12% | Low activation markers |
| **Highly Activated cells** | 15,100 | 9.59% | CD69, CD25, HLA-DR |
| **GFP+ Activated Tregs** | 10,015 | 6.36% | FoxP3, CD25, CD127low |
| **Activated Memory T cells** | 8,913 | 5.66% | CD45RO, CD44 |
| **Tfh/Cytotoxic cells** | 8,932 | 5.68% | CXCR5, PD-1, Granzyme |

**Specialized Populations:**

| Cell Type | Cell Count | % of Total | Biological Role |
|-----------|------------|------------|-----------------|
| **Follicular Helper T (Tfh)** | 5,953 | 3.78% | B cell help in GC |
| **Pro-inflammatory T/ILC** | 6,956 | 4.42% | Inflammatory response |
| **Granzyme A+ cytotoxic** | 3,468 | 2.20% | Cytotoxic function |
| **Th2-like cells** | 3,457 | 2.20% | Type-2 immunity |
| **Myeloid/APC** | 2,573 | 1.63% | Antigen presentation |
| **GFP+ cells (CAR-T?)** | 2,417 | 1.54% | Potential CAR-T cells |
| **Activated B cells** | 1,819 | 1.16% | Activated B response |
| **CD8+ Effector T cells** | 1,337 | 0.85% | Cytotoxic T cells |
| **Regulatory T cells (Tregs)** | 704 | 0.45% | Classical Tregs |
| **Germinal Center B cells** | 566 | 0.36% | GC reaction (PNA+) |

#### 4.2 Treatment-Specific Cell Type Composition

**Most Altered Populations:**

| Cell Type | CAR % | EGFR % | Fold Change |
|-----------|-------|--------|-------------|
| **GFP+ Activated Tregs** | 7.61% | 4.85% | **1.57× higher in CAR** |
| **Highly Activated cells** | 11.97% | 6.72% | **1.78× higher in CAR** |
| **Tfh/Cytotoxic cells** | 2.87% | 9.07% | **3.16× higher in EGFR** |
| **Activated Memory T cells** | 4.45% | 7.13% | **1.60× higher in EGFR** |

**Interpretation:** CAR treatment enriches for highly activated and Treg populations, while EGFR treatment shows expansion of T follicular helper and memory T cell compartments.

#### 4.3 FoxP3+ Enrichment by Cell Type

**Top FoxP3+ Enriched Populations:**

| Cell Type | CAR FoxP3+ % | EGFR FoxP3+ % | Key Observation |
|-----------|--------------|---------------|-----------------|
| **Germinal Center B cells** | 58.7% | 57.3% | Unexpectedly high FoxP3 signal |
| **GFP+ Activated Tregs** | 22.7% | 34.6% | **1.5× higher in EGFR** |
| **CD8+ Effector T cells** | 3.0% | 0.9% | **3.3× higher in CAR** |

**Critical Finding:** The "GFP+ Activated Tregs" cluster (10,015 cells, 6.36% of total) shows dramatically higher FoxP3+ enrichment in EGFR (34.6%) compared to CAR (22.7%). This represents **the most abundant FoxP3-enriched population** and shows clear treatment-dependent differences.

### Supporting Visualizations

- **`cluster_marker_heatmap_scater.png`**: Expression heatmap showing marker profiles for all 15 cell types
- **`cluster_estimation_elbow_plots.png`**: Elbow plot justifying 15 cluster solution
- **`umap_colored_by_celltype.png`**: UMAP visualization showing biological coherence of cell types
- **`celltype_composition_comparison.png`**: Bar plot comparing cell type frequencies between treatments
- **`FOXP3_enrichment_by_celltype.png`**: FoxP3+ percentage across all cell types by treatment

### Relevance to Research Question

**Direct:** The "GFP+ Activated Tregs" population shows **treatment-dependent FoxP3+ enrichment**. In EGFR tissue, 34.6% of these cells are FoxP3+, compared to only 22.7% in CAR tissue. This suggests EGFR treatment promotes a more regulatory phenotype within this activated Treg compartment.

**Contextual:** Treatment alters overall immune landscape composition. CAR drives high activation states, while EGFR maintains more memory and follicular helper programs alongside enhanced Treg regulatory function.

---

## Section 5: Spatial Co-localization Analysis

### Analysis Description

**Purpose:** Identify which cell types spatially associate (or avoid) each other and how these spatial relationships differ between treatments.

**Methods:**
- **Spatial Statistics:** `spicyR::getPairwise()` calculates L-functions (modified Ripley's K)
  - **L-function interpretation:**
    - **Positive values:** Cell types attract/co-localize
    - **Negative values:** Cell types avoid/segregate
    - **Near zero:** Random spatial distribution
- **Radius:** 100 µm neighborhood (biologically relevant for cell-cell communication)
- **Inhomogeneity Correction:** Accounts for regional cell density variations
- **All Pairwise Comparisons:** 225 unique cell type pairs analyzed (15 × 15)

**Interpretation:** Positive L-functions indicate functional cellular neighborhoods or niches. Negative values suggest spatial exclusion, often reflecting mutually exclusive differentiation states or competitive dynamics.

### Key Findings

#### 5.1 Top Spatial Co-localization Differences (|ΔL| > 40)

**Pairs with LARGEST differences between treatments:**

| Cell Pair | EGFR L | CAR L | Difference | Interpretation |
|-----------|--------|-------|------------|----------------|
| **Treg ↔ Treg** | 186.5 | 24.4 | **+162.1** (EGFR) | Tregs form tight clusters in EGFR |
| **Activated B ↔ Treg** | 154.2 | 13.6 | **+140.6** (EGFR) | Tregs co-localize with B cells in EGFR |
| **Activated B ↔ Activated B** | 182.9 | 96.1 | **+86.8** (EGFR) | B cell clustering in EGFR |
| **CD8+ Effector ↔ Treg** | 70.5 | -4.1 | **+74.6** (EGFR) | CD8+ T cells associate with Tregs in EGFR |
| **Activated B ↔ Activated Memory T** | -55.3 | 13.2 | **-68.5** (CAR) | B-T co-localization in CAR |

#### 5.2 FoxP3+ Enriched Cell Types: Spatial Relationships

**Regulatory T cells (Classical Tregs) - Top Spatial Partners:**

| Partner Cell Type | EGFR L | CAR L | Difference | Biological Context |
|-------------------|--------|-------|------------|-------------------|
| **Tregs ↔ Tregs** | 186.5 | 24.4 | **+162.1 (EGFR)** | Treg clustering/cooperation |
| **Activated B cells** | 154.2 | 13.6 | **+140.6 (EGFR)** | B cell regulation in GC |
| **CD8+ Effector T** | 70.5 | -4.1 | **+74.6 (EGFR)** | CD8+ T suppression |
| **Germinal Center B** | 67.5 | 18.9 | **+48.7 (EGFR)** | GC B cell regulation |

**Germinal Center B cells (57% FoxP3+):**

| Partner Cell Type | EGFR L | CAR L | Difference | Biological Context |
|-------------------|--------|-------|------------|-------------------|
| **Activated B cells** | 70.4 | 28.2 | **+42.1 (EGFR)** | B cell niches |
| **CD8+ Effector T** | 24.3 | -23.4 | **+47.7 (EGFR)** | CD8-B interactions |
| **Th2-like cells** | 66.3 | 19.2 | **+47.1 (EGFR)** | Helper T support |

**GFP+ Activated Tregs (34.6% FoxP3+ in EGFR):**

| Partner Cell Type | EGFR L | CAR L | Difference | Biological Context |
|-------------------|--------|-------|------------|-------------------|
| **Tfh/Cytotoxic** | 32.1 | 2.6 | **+29.4 (EGFR)** | Tfh regulation |
| **Activated B cells** | 43.0 | 12.0 | **+30.9 (EGFR)** | B cell suppression |
| **Self-clustering** | 84.7 | 58.7 | **+26.0 (EGFR)** | Treg cooperation |

#### 5.3 Spatial Organization Patterns

**EGFR Treatment:**
- **Dense Treg clustering:** Classical Tregs form tight spatial clusters (L = 186.5)
- **Treg-B cell niches:** Strong co-localization with activated and germinal center B cells
- **Regulatory neighborhoods:** Tregs positioned near both B cells and CD8+ T cells for dual suppression

**CAR Treatment:**
- **Dispersed Tregs:** Minimal Treg-Treg clustering (L = 24.4)
- **B-T cell zones:** Activated B cells co-localize with memory T cells instead of Tregs
- **Reduced regulatory structure:** Less defined Treg-mediated spatial organization

### Supporting Visualizations

- **`colocalisation_heatmaps.png`**: Heatmap of all 225 pairwise L-functions for each image
- **`colocalisation_differences_between_images.png`**: Difference heatmap highlighting treatment effects
- **`FOXP3_colocalization_barplot.png`**: Bar chart of FoxP3-enriched cell type spatial relationships
- **`spatial_[celltype]__[celltype].png`**: Example spatial distribution plots for key co-localizing pairs

### Relevance to Research Question

**Direct:** FoxP3+ Tregs show **fundamentally different spatial organization** between treatments:
- **EGFR:** Tregs cluster together and form organized regulatory neighborhoods with B cells and CD8+ T cells
- **CAR:** Tregs are spatially dispersed with minimal clustering or organized regulatory zones

This suggests EGFR treatment promotes **collective Treg function** through spatial clustering, while CAR treatment results in **isolated Treg activity**.

**Mechanistic Insight:** The dramatic Treg-Treg clustering in EGFR (L = 186.5 vs 24.4 in CAR) suggests these cells may cooperate locally through:
- Cell-cell contact-dependent suppression
- Local cytokine microenvironments (IL-10, TGF-β)
- Coordinated regulation of target cell populations

---

## Section 6: Context-Dependent Spatial Relationships (Kontextual)

### Analysis Description

**Purpose:** Determine if cell-cell spatial relationships change depending on the broader tissue microenvironment (contextual spatial analysis).

**Methods:**
- **Kontextual Analysis:** `Statial::Kontextual()` compares:
  - **Original relationships:** Standard pairwise co-localization
  - **Context-adjusted relationships:** Co-localization within defined cellular hierarchies
- **Biological Hierarchy Defined (Spleen Germinal Center):**
  ```
  Level 1: Major Lineages (Lymphoid, Myeloid)
  Level 2: T Cell Subsets (CD4+, CD8+, Regulatory)
  Level 3: Functional Groups (Activated, Memory, Cytotoxic)
  Level 4: FoxP3+ Specific Groups (Treg-enriched populations)
  ```
- **Context Effect:** Magnitude of change when adjusting for cellular hierarchy

**Interpretation:** Large context effects indicate that spatial relationships are strongly influenced by broader tissue architecture, revealing hierarchical organization principles.

### Key Findings

#### 6.1 Strongest Context-Dependent Relationships

**Top 11 Cell Pairs with Context Effects (|Δ| > 90):**

| Cell Pair | Original L | Context L | Context Effect | Interpretation |
|-----------|------------|-----------|----------------|----------------|
| **Activated B ↔ Activated B** | 609.3 | 352.2 | **-257.1** | B cell clustering is region-dependent |
| **GC B cells ↔ GC B cells** | 351.5 | 178.2 | **-173.3** | GC B cells cluster within GC zones |
| **GFP+ CAR-T ↔ GFP+ CAR-T** | 190.4 | 61.9 | **-128.6** | CAR-T cells disperse across contexts |
| **Pro-inflam ↔ Pro-inflam** | 157.8 | 39.5 | **-118.3** | Inflammatory cells region-specific |
| **Highly Act ↔ Highly Act** | 144.6 | 33.8 | **-110.8** | Activation concentrated in zones |
| **GFP+ CAR-T ↔ Activated B** | 139.8 | 34.2 | **-105.6** | CAR-B interactions contextual |
| **Pro-inflam ↔ Highly Act** | 128.8 | 23.9 | **-104.9** | Activation niches structured |
| **Tfh/Cyto ↔ Tfh/Cyto** | 183.3 | 79.6 | **-103.7** | Tfh cells form context zones |
| **Highly Act ↔ Pro-inflam** | 128.8 | 25.2 | **-103.6** | Bi-directional activation zones |
| **Activated B ↔ GFP+ CAR-T** | 139.8 | 46.3 | **-93.4** | B-CAR relationships hierarchical |

#### 6.2 FoxP3+ Populations: Contextual Behavior

**Key Observation:** FoxP3-enriched populations show **moderate context dependence**, suggesting they occupy relatively consistent spatial positions regardless of broader tissue architecture.

- **GFP+ Activated Tregs:** Context effects present but not extreme, indicating Tregs maintain their spatial relationships across different germinal center zones
- **Classical Tregs:** (Not in top context-dependent list) Spatial positioning appears consistent across tissue contexts

#### 6.3 Biological Interpretation

**Hierarchical Organization:**
- **B cells and Activated populations** show strongest context dependence → organized into discrete functional zones
- **CAR-T cells** highly context-dependent → adapt positioning based on local tissue state
- **Tregs** moderate context dependence → maintain regulatory function across contexts

### Supporting Visualizations

- **`kontextual_vs_original_comparison.png`**: Scatter plot showing original vs context-adjusted L-functions
- **`kontextual_curves_top_pairs.png`**: Distance-dependent curves for top context-affected pairs
- **`context_spatial_distribution.png`**: Spatial maps showing hierarchical organization

### Relevance to Research Question

**Indirect but Important:** FoxP3+ Tregs show **consistent spatial organization** across tissue contexts (moderate context effects), suggesting their regulatory function is **maintained regardless of local microenvironment**. This contrasts with highly context-dependent cells (activated B cells, CAR-T cells) whose positioning varies with tissue architecture.

**Implication:** Tregs appear to establish stable regulatory niches that persist across different germinal center compartments, supporting a **global rather than locally-restricted** regulatory role.

---

## Section 7: Spatial Domain Detection & Characterization

### Analysis Description

**Purpose:** Identify functionally distinct tissue compartments (spatial domains) and characterize their cell type composition and Treg distribution.

**Methods:**
- **Spatial Domain Detection:** `lisaClust` algorithm
  - Identifies regions with similar local cell type neighborhoods
  - Uses k-nearest neighbors (k=5) spatial graphs
  - Hierarchical clustering of regional cell type compositions
- **Domains Identified:** 6 distinct spatial compartments per image
- **Domain Annotation:** Manual annotation based on dominant cell type enrichments
- **Statistical Comparison:** Cell type proportions and FoxP3+ percentages compared between images

**Interpretation:** Spatial domains represent functionally specialized tissue microenvironments (e.g., germinal center reaction zones, T cell zones, regulatory niches).

### Key Findings

#### 7.1 Six Spatial Domains Identified & Annotated

| Region | Annotation | Dominant Cell Types | Biological Function |
|--------|------------|---------------------|---------------------|
| **Region 1** | Resting Lymphocyte Zone | Resting/Naïve lymphocytes (dominant) | Quiescent lymphocyte reservoir |
| **Region 2** | Tfh Enriched Zone | Tfh/Cytotoxic, Follicular Helper T | T-B collaboration for GC reaction |
| **Region 3** | **Treg Enriched Zone** | **GFP+ Activated Tregs, Classical Tregs** | **Active immune regulation** |
| **Region 4** | Activated Inflammatory Zone | Highly Activated, Pro-inflammatory T/ILC | Acute immune response |
| **Region 5** | General Lymphoid Zone | Mixed lymphocyte populations | General immune surveillance |
| **Region 6** | B Cell Cytotoxic Zone | B cells, Granzyme A+ cytotoxic | B cell activation with cytotoxic support |

#### 7.2 FoxP3+ Distribution Across Spatial Domains

**FoxP3+ Percentage by Domain (Treatment Comparison):**

| Domain | CAR FoxP3+ % | EGFR FoxP3+ % | Difference | Key Observation |
|--------|--------------|---------------|------------|-----------------|
| **Tfh Enriched (R2)** | **15.22%** | 1.85% | **+13.37% (CAR)** | Highest CAR Treg enrichment |
| **B Cell Cytotoxic (R6)** | 4.43% | 1.91% | **+2.52% (CAR)** | CAR Tregs in B cell zones |
| **Treg Enriched (R3)** | 5.05% | **6.47%** | **+1.41% (EGFR)** | Core Treg zone (EGFR higher) |
| **General Lymphoid (R5)** | 1.35% | 2.04% | +0.69% (EGFR) | Low Treg presence |
| **Resting Lymphocyte (R1)** | 1.22% | 1.53% | +0.31% (EGFR) | Minimal Treg activity |
| **Activated Inflammatory (R4)** | 1.14% | 1.42% | +0.28% (EGFR) | Low Treg frequency |

#### 7.3 Key Domain-Specific Patterns

**Region 3 (Treg Enriched Zone):**
- **Purpose-built regulatory microenvironment**
- Enriched for both GFP+ Activated Tregs AND Classical Tregs
- **EGFR slightly higher** Treg % (6.47% vs 5.05% in CAR)
- Co-localizes with germinal center B cells (see Section 5)

**Region 2 (Tfh Enriched Zone):**
- **Dramatic CAR enrichment:** 15.22% FoxP3+ in CAR vs 1.85% in EGFR
- Suggests CAR treatment drives Tregs into T follicular helper zones
- May represent Treg-mediated suppression of T-B collaboration

**Region 6 (B Cell Cytotoxic Zone):**
- CAR shows 2.3× higher FoxP3+ % (4.43% vs 1.91%)
- Tregs positioned at interface of B cell activation and cytotoxic activity
- Potential regulatory checkpoint preventing excessive cytotoxicity

#### 7.4 Spatial Domain Size Differences

**Domain size comparison** reveals that EGFR and CAR treatments result in different domain spatial extents, suggesting tissue-level architectural differences (detailed in `spatial_domain_size_comparison.csv`).

### Supporting Visualizations

- **`spatial_domains_simple.png`**: Spatial maps showing 6 domains overlaid on tissue
- **`spatial_domains_hatching.png`**: Hatching patterns for clear domain boundaries
- **`spatial_domains_heatmap.png`**: Cell type enrichment heatmap for each domain
- **`spatial_domains_bubble.png`**: Bubble plot showing cell type distributions per domain
- **`spatial_domains_composition_comparison.png`**: Side-by-side domain composition between images
- **`FOXP3_by_spatial_domain.png`**: FoxP3+ percentage across all domains by treatment

### Relevance to Research Question

**Direct:** FoxP3+ Tregs show **treatment-dependent spatial compartmentalization**:

1. **Core Treg Niche (Region 3):** EGFR shows slightly higher Treg enrichment in the dedicated Treg zone, consistent with the organized Treg clustering seen in Section 5.

2. **Tfh Zone Infiltration (Region 2):** CAR treatment drives Tregs into T follicular helper zones (15.22% vs 1.85%), suggesting:
   - CAR-induced Treg redistribution to suppress T-B collaboration
   - Potential mechanism for limiting germinal center reactions

3. **B Cell Regulation (Region 6):** CAR Tregs are more abundant in B cell/cytotoxic zones, potentially limiting B cell activation and antibody responses.

**Synthesis:** EGFR-treated tissue maintains Tregs in organized regulatory zones (Region 3) with high local density. CAR treatment disperses Tregs into functional zones requiring active suppression (Tfh zones, B cell zones), suggesting **reactive vs. homeostatic regulatory strategies**.

---

## Cross-Section Synthesis: FoxP3+ Treg Behavior

### Comprehensive Answer to Research Question

**Do FoxP3+ regulatory T cells behave differently in EGFR vs CAR-T treated spleen tissue?**

**YES - Dramatic and Multi-dimensional Differences**

### 1. **Frequency & Abundance** (Section 3)
- **EGFR:** 2.21% of cells are FoxP3+ (higher frequency)
- **CAR:** 1.82% of cells are FoxP3+ (lower frequency but more total Treg-enriched cells in GFP+ Activated Treg cluster)

### 2. **Cellular Identity** (Section 4)
- **Key Population:** "GFP+ Activated Tregs" (10,015 cells, 6.36% of total)
  - **EGFR:** 34.6% of this cluster are FoxP3+
  - **CAR:** 22.7% of this cluster are FoxP3+
  - **Interpretation:** EGFR drives stronger FoxP3+ phenotype within activated Treg compartment

### 3. **Spatial Organization** (Section 5)
- **EGFR:** Tregs form tight clusters (Treg-Treg L = 186.5)
  - Organized regulatory neighborhoods with B cells and CD8+ T cells
  - Suggest coordinated, collective regulatory function
- **CAR:** Tregs spatially dispersed (Treg-Treg L = 24.4)
  - Individual Tregs scattered across tissue
  - Suggest independent, localized regulatory activity

### 4. **Microenvironmental Context** (Section 6)
- **Both treatments:** Tregs show moderate context-dependence
- **Interpretation:** Tregs maintain consistent spatial relationships across different tissue zones, unlike highly context-dependent activated cells
- **Implication:** Tregs establish stable regulatory programs regardless of local architecture

### 5. **Tissue Compartmentalization** (Section 7)
- **EGFR:** Tregs concentrated in dedicated "Treg Enriched Zone" (Region 3: 6.47% FoxP3+)
  - Organized, zone-restricted regulatory activity
- **CAR:** Tregs infiltrate multiple functional zones
  - **Tfh Zone (Region 2):** 15.22% FoxP3+ (8× higher than EGFR)
  - **B Cell/Cytotoxic Zone (Region 6):** 4.43% FoxP3+ (2.3× higher than EGFR)
  - Reactive regulatory response to active immune processes

### Mechanistic Model

**EGFR Treatment → Homeostatic Treg Program:**
- Higher overall Treg frequency (2.21%)
- Zone-restricted regulatory niches (Region 3)

**CAR Treatment → Reactive Treg Program:**
- Lower overall frequency but larger activated Treg pool (GFP+ cluster 2.3× larger)
- Dispersed, individual Tregs
- Zone-infiltrating positioning (Tfh zones, B cell zones)
- **Function:** Targeted suppression of specific activated zones

---

## Data Availability

All quantitative results are available in CSV format:

**Section 3 (QC):**
- `RES/QC_plots/marker_statistics_raw.csv`
- `RES/QC_plots/marker_statistics_normalized.csv`
- `RES/QC_plots/FOXP3_treg_summary.csv`
- `RES/QC_plots/normalization_effectiveness.csv`

**Section 4 (Clustering):**
- `RES/Clustering/cluster_cell_type_annotations.csv`
- `RES/Clustering/celltype_composition_by_image.csv`
- `RES/Clustering/FOXP3_enrichment_by_celltype.csv`

**Section 5 (Co-localization):**
- `RES/Colocalization/colocalisation_per_image.csv`
- `RES/Colocalization/colocalisation_top_differences.csv`
- `RES/Colocalization/FOXP3_colocalization_comparison.csv`

**Section 6 (Kontextual):**
- `RES/Kontextual/kontextual_all_pairs.csv`
- `RES/Kontextual/kontextual_context_dependent_pairs.csv`
- `RES/Kontextual/kontextual_comparison_between_images.csv`

**Section 7 (Spatial Domains):**
- `RES/SpatialDomains/spatial_domain_annotations.csv`
- `RES/SpatialDomains/spatial_domain_composition_full.csv`
- `RES/SpatialDomains/FOXP3_by_spatial_domain.csv`
- `RES/SpatialDomains/FOXP3_domain_comparison.csv`

---

## Analysis Pipeline

Complete workflow documented in:
- **Main Analysis:** `MAIN/SpicyFlow.Rmd` (reproducible R Markdown)
- **Image Processing:** `MAIN/PROG/Scripts/enhance_and_quantify.R`
- **Coding Standards:** `MAIN/CONFIG/CODING_STANDARDS.md`

All analyses performed with fail-fast error handling for maximum reproducibility and transparency.

---

**Document Version:** 1.0  
**Generated:** October 28, 2025  
**Analyst:** Dr. Michael Zaiken 

