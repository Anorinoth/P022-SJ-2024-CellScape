# Script to generate Spatial Domain Composition Visualizations
# Author: Dr. Michael Zaiken
# Stacked Bar Charts and Pie Charts of domain areas (cell counts)

library(SpatialExperiment)
library(dplyr)
library(ggplot2)

# Load canonical paths from setup_project
source(file.path(dirname(sys.frame(1)$ofile), "setup_project.R"))
dirs <- setup_project(verbose = FALSE)
PROG_DIR        <- dirs$PROG_DIR
RES_DIR         <- dirs$RES_DIR
CHECKPOINT_FILE <- file.path(dirs$CHECKPOINT_DIR, "combined_spe_after_section7_v2.rds")

# Define Domain Colors (Matching SpicyFlow.Rmd / Set2 palette)
domain_colors <- c(
    "B_Cell_Zone"           = "#4A90B8",  # muted steel blue (BCF + GC combined)
    "T_Cell_Zone"           = "#9B59B6",  # soft violet (formerly Treg_Tfh_Zone)
    "Th2_T_Cell_Zone"       = "#2E8B57",  # sea green
    "Marginal_Zone"         = "#E6AB02",  # gold
    "Resting_Lymphoid_Zone" = "#D9D9D9",  # light grey — neutral
    "Low_Intensity_Zone"    = "#BDBDBD"   # medium grey — neutral
)

# Ensure output directory exists
out_dir <- file.path(RES_DIR, "SpatialDomains", "Composition")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message("Loading checkpoint: ", CHECKPOINT_FILE)
if (!file.exists(CHECKPOINT_FILE)) stop("Checkpoint file not found!")
combined_spe <- readRDS(CHECKPOINT_FILE)

message("Loaded SPE object with ", ncol(combined_spe), " cells")

# 1. Calculate Composition
# Count cells per domain per image
# Assumption: Cell density is roughly uniform, so cell count ~ area
composition_data <- as.data.frame(colData(combined_spe)) %>%
    group_by(imageID, annotation) %>%
    summarise(cell_count = n(), .groups = "drop") %>%
    group_by(imageID) %>%
    mutate(
        total_cells = sum(cell_count),
        proportion = cell_count / total_cells,
        percentage = proportion * 100,
        # Create label for pie chart (e.g., "15%")
        label = paste0(round(percentage, 1), "%")
    )

# Save summary table
write.csv(composition_data, file.path(out_dir, "spatial_domain_composition_stats.csv"), row.names = FALSE)
message("Saved composition stats to: ", file.path(out_dir, "spatial_domain_composition_stats.csv"))

# 2. Stacked Bar Chart
message("Generating Stacked Bar Chart...")
p_bar <- ggplot(composition_data, aes(x = imageID, y = proportion, fill = annotation)) +
    geom_col(position = "stack", width = 0.7) +
    scale_fill_manual(values = domain_colors) +
    scale_y_continuous(labels = scales::percent) +
    labs(
        title = "Spatial Domain Composition by Image",
        subtitle = "Proportion of total tissue area (inferred from cell counts)",
        x = "Image / Condition",
        y = "Proportion of Tissue",
        fill = "Spatial Domain"
    ) +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
    )

ggsave(file.path(out_dir, "spatial_domain_composition_bar.png"), p_bar, width = 8, height = 6, bg = "white")

# 3. Pie Charts
message("Generating Pie Charts...")
p_pie <- ggplot(composition_data, aes(x = "", y = proportion, fill = annotation)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = domain_colors) +
    facet_wrap(~imageID) +
    labs(
        title = "Spatial Domain Composition",
        fill = "Spatial Domain"
    ) +
    theme_void() +
    theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold")
    )

ggsave(file.path(out_dir, "spatial_domain_composition_pie.png"), p_pie, width = 10, height = 6, bg = "white")

message("Analysis complete. Visualizations saved to: ", out_dir)
