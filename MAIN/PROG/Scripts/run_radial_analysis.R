# Script to run Radial Enrichment Analysis
# Quantifies Treg density proximal to B Cell Follicle boundaries

library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(FNN)

# Define paths
PROG_DIR <- "/mnt/b/Projects/Image_Analysis_Projects/SJ_04.2/MAIN/PROG"
RES_DIR <- "/mnt/b/Projects/Image_Analysis_Projects/SJ_04.2/MAIN/RES"
CHECKPOINT_FILE <- file.path(PROG_DIR, "Checkpoints", "combined_spe_after_section7.rds")

# Load analysis functions
source(file.path(PROG_DIR, "Scripts", "spatial_radial_analysis.R"))

# Ensure output directory exists
radial_dir <- file.path(RES_DIR, "SpatialDomains", "Radial")
if (!dir.exists(radial_dir)) dir.create(radial_dir, recursive = TRUE)

message("Loading checkpoint: ", CHECKPOINT_FILE)
if (!file.exists(CHECKPOINT_FILE)) stop("Checkpoint file not found!")
combined_spe <- readRDS(CHECKPOINT_FILE)

message("Loaded SPE object with ", ncol(combined_spe), " cells")

# Run Radial Analysis
# Target: B_Cell_Zone
# Interest: Treg (matches "Treg" or "FoxP3")
radial_stats <- calculate_radial_enrichment(
    combined_spe,
    target_domain = "B_Cell_Zone",
    cell_of_interest = "Treg", # Matches "Transgenic Tregs" and "Endogenous Tregs"
    domain_col = "annotation",
    celltype_col = "cellType",
    imageID_col = "imageID",
    bin_width = 20, # 20 micron bins
    max_dist = 200 # Analyze up to 200 microns out
)

# Save Stats
write.csv(radial_stats, file.path(radial_dir, "radial_enrichment_stats.csv"), row.names = FALSE)
message("Saved radial stats to: ", file.path(radial_dir, "radial_enrichment_stats.csv"))


# Visualization
message("Generating Radial Plot...")

# Calculate midpoints for plotting (e.g., bin (0,20] -> 10)
# Extract numeric range from factor bin "(0,20]"
radial_stats$bin_mid <- as.numeric(sub("^\\((.*),.*", "\\1", radial_stats$bin)) + 10

# Plot Treg Proportion vs Distance
p <- ggplot(radial_stats, aes(x = bin_mid, y = treg_prop, color = imageID, group = imageID)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    labs(
        title = "Treg Density vs Distance from B Cell Follicle",
        subtitle = "Radial expansion from follicle boundary",
        x = "Distance from Boundary (µm)",
        y = "Proportion of Tregs (Tregs / Total Cells in Bin)",
        color = "Image / Condition"
    ) +
    theme_classic() +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(breaks = seq(0, 200, 20))

ggsave(file.path(radial_dir, "radial_treg_enrichment.png"), p, width = 8, height = 6, bg = "white")
message("Saved plot to: ", file.path(radial_dir, "radial_treg_enrichment.png"))
