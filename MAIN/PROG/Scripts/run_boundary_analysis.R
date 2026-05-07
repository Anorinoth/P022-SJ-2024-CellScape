# Script to run spatial domain boundary analysis
# Loads data from checkpoint, runs calculation, and generates plots

library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(deldir)

# Define paths
PROG_DIR <- "/mnt/b/Projects/Image_Analysis_Projects/SJ_04.2/MAIN/PROG"
RES_DIR <- "/mnt/b/Projects/Image_Analysis_Projects/SJ_04.2/MAIN/RES"
CHECKPOINT_FILE <- file.path(PROG_DIR, "Checkpoints", "combined_spe_after_section7.rds")

# Load analysis functions
source(file.path(PROG_DIR, "Scripts", "spatial_domain_analysis.R"))
library(BiocParallel)

# Set up parallel backend
# Using SnowParam because MulticoreParam with deldir (Fortran) causes memory instability/hanging
BPPARAM <- SnowParam(workers = 4)

message("Loading checkpoint: ", CHECKPOINT_FILE)
if (!file.exists(CHECKPOINT_FILE)) stop("Checkpoint file not found!")
combined_spe <- readRDS(CHECKPOINT_FILE)

message("Loaded SPE object with ", ncol(combined_spe), " cells")

# Define existing colors for consistency
celltype_colors <- c(
    # Artifacts & low quality (greys)
    "Doublets/Multiplets" = "grey50",
    "Low-quality / Artifacts" = "grey35",
    "Intermediate/Transitional Cells" = "grey70",
    "Resting/Naïve Lymphocytes" = "grey85",
    "IFNg+ Low-intensity Cells" = "grey60",
    # Treg family (purples)
    "Transgenic Tregs (GFP+ FoxP3+)" = "#8B008B", # darkmagenta
    "Endogenous Tregs (FoxP3-high)" = "#9370DB", # mediumpurple
    "CD40+ Treg-APC Conjugates" = "#BA55D3", # mediumorchid
    # B cells / APCs (blues)
    "B cells (CD19+ MHC-II+)" = "#00008B", # darkblue
    "Activated B cells / APCs (CD40+ CD19+)" = "#4169E1", # royalblue
    # CD8+ T cells (oranges/reds)
    "Activated CD8+ T cells" = "#FF4500", # orangered
    "Exhausted T cells (PD1-high)" = "#CD3700", # darkorange3
    # Th2 / cytokine-producing (yellows/greens)
    "Th2/IL4-high T cells" = "#48D1CC", # mediumturquoise
    "Polyfunctional Effectors (IFNg+ IL21+)" = "#00CD00", # green3
    "IL10+ CD4+ T cells (Tr1-like)" = "#FFD700" # gold
)

# Define Domain Colors (Matching SpicyFlow.Rmd / Set2 palette)
domain_colors <- c(
    "B_Cell_Zone"           = "#4A90B8",  # muted steel blue (BCF + GC combined)
    "T_Cell_Zone"           = "#9B59B6",  # soft violet (formerly Treg_Tfh_Zone)
    "Th2_T_Cell_Zone"       = "#2E8B57",  # sea green
    "Marginal_Zone"         = "#E6AB02",  # gold
    "Resting_Lymphoid_Zone" = "#D9D9D9",  # light grey — neutral
    "Low_Intensity_Zone"    = "#BDBDBD"   # medium grey — neutral
)

# Run Calculation
# Run Calculation
domains_dir <- file.path(RES_DIR, "SpatialDomains")

# Force calculation to get plot_data (edges) which are needed for visualization
# logical check for existing file removed to ensure we get the plotting vectors
boundary_results <- calculate_domain_boundaries(
    combined_spe,
    imageID_col = "imageID",
    domain_col = "annotation",
    output_dir = file.path(domains_dir, "Boundaries"),
    min_patch_area = 500,
    BPPARAM = BPPARAM
)

# Generate Plots (Reproducing logic from Rmd)
message("\nGenerating specific visualizations...")
boundary_stats <- boundary_results$stats

if (!is.null(boundary_stats) && nrow(boundary_stats) > 0) {
    # Filter for Treg <-> B Cell Follicle
    target_interaction <- boundary_stats %>%
        filter(
            (source_domain == "T_Cell_Zone" & target_domain == "B_Cell_Zone") |
                (source_domain == "B_Cell_Zone" & target_domain == "T_Cell_Zone")
        )

    if (nrow(target_interaction) > 0) {
        message("Plotting bar chart...")

        # Plot 1: Standard Shared Boundary %
        p_bar <- ggplot(target_interaction, aes(x = imageID, y = pct_shared_boundary, fill = imageID)) +
            geom_col() +
            scale_y_continuous(labels = scales::percent_format(scale = 1)) +
            labs(
                title = "Shared Boundary: T Cell Zone <-> B Cell Zone",
                subtitle = "Pct of Source Domain perimeter shared with Target Domain",
                y = "% Shared Boundary"
            ) +
            facet_wrap(~source_domain, scales = "free_y") +
            theme_classic()

        ggsave(file.path(domains_dir, "Boundaries", "Treg_Bcell_Shared_Boundary.png"),
            p_bar,
            width = 8, height = 6, bg = "white"
        )

        # Plot 2: Enrichment Score
        if ("enrichment_score" %in% colnames(target_interaction)) {
            p_enrich <- ggplot(target_interaction, aes(x = imageID, y = enrichment_score, fill = imageID)) +
                geom_col() +
                geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
                labs(
                    title = "Enrichment: T Cell Zone <-> B Cell Zone",
                    subtitle = "Boundary Enrichment Score (>1 = Preference, <1 = Avoidance)",
                    y = "Enrichment Score (Observed / Expected)"
                ) +
                facet_wrap(~source_domain, scales = "free_y") +
                theme_classic()

            ggsave(file.path(domains_dir, "Boundaries", "Treg_Bcell_Enrichment_Score.png"),
                p_enrich,
                width = 8, height = 6, bg = "white"
            )
        }
    }
}

# Generate Maps
if (!is.null(boundary_results$plot_data)) {
    for (img in names(boundary_results$plot_data)) {
        message("Plotting map for: ", img)
        edges <- boundary_results$plot_data[[img]]
        boundary_edges <- edges %>% filter(type == "Boundary")

        boundary_edges$highlight <- ifelse(
            grepl("T_Cell_Zone", boundary_edges$pair) & grepl("B_Cell_Zone", boundary_edges$pair),
            "Target Interface", "Other Boundary"
        )

        spe_sub <- combined_spe[, combined_spe$imageID == img]
        cell_coords <- as.data.frame(spatialCoords(spe_sub))
        cell_coords$domain <- spe_sub$annotation

        p_map <- ggplot() +
            geom_point(data = cell_coords, aes(x = x, y = y, color = domain), size = 0.1, alpha = 0.3) +
            geom_segment(
                data = boundary_edges,
                aes(x = x1, y = y1, xend = x2, yend = y2, size = highlight, color = highlight),
                alpha = 0.8
            ) +
            scale_color_manual(values = c(
                domain_colors,
                "Target Interface" = "#FF1493",
                "Other Boundary" = "grey40"
            )) +
            scale_size_manual(values = c(
                "Target Interface" = 0.8,
                "Other Boundary" = 0.1
            )) +
            labs(title = paste("Domain Boundaries:", img)) +
            theme_void() +
            guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

        ggsave(file.path(domains_dir, "Boundaries", paste0("boundary_map_", img, ".png")),
            p_map,
            width = 12, height = 10, bg = "white"
        )
    }
}

message("Analysis complete.")
