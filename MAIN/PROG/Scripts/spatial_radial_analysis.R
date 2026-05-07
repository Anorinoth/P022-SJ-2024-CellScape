################################################################################
# spatial_radial_analysis.R
# Author: Dr. Michael Zaiken
#
# Function to calculate cell type enrichment radially from a specific
# spatial domain boundary.
################################################################################

#' Calculate Radial Proximity Enrichment
#'
#' Calculates the density/proportion of specific cell types at increasing
#' distances from a target spatial domain.
#'
#' @param combined_spe SpatialExperiment object
#' @param target_domain Name of the domain to expand from (e.g., "B_Cell_Zone")
#' @param cell_of_interest Name (or regex) of cell type to quantify (e.g., "Treg")
#' @param domain_col Column name for spatial domains (default: "annotation")
#' @param celltype_col Column name for cell types (default: "cellType")
#' @param imageID_col Column name for image IDs
#' @param bin_width Width of distance bins in microns (default: 10)
#' @param max_dist Maximum distance to analyze (default: 200)
#'
#' @return List containing summary stats and plot data
#' @export
calculate_radial_enrichment <- function(combined_spe,
                                        target_domain = "B_Cell_Zone",
                                        cell_of_interest = "Treg",
                                        domain_col = "annotation",
                                        celltype_col = "cellType",
                                        imageID_col = "imageID",
                                        bin_width = 10,
                                        max_dist = 100) {
    require(SpatialExperiment)
    require(dplyr)
    require(FNN) # Fast Nearest Neighbors for distance calc

    message("Starting Radial Enrichment Analysis...")
    message("  Target Domain: ", target_domain)
    message("  Cell of Interest Pattern: ", cell_of_interest)

    # Process each image separately
    image_ids <- unique(combined_spe[[imageID_col]])

    results_list <- list()

    for (img in image_ids) {
        message("  Processing image: ", img)

        # Subset data
        spe_sub <- combined_spe[, combined_spe[[imageID_col]] == img]
        coords <- spatialCoords(spe_sub)
        meta <- colData(spe_sub)

        # 1. Identify Target Cells (The "Anchors")
        # Cells that are INSIDE the target domain
        is_target_domain <- meta[[domain_col]] == target_domain
        target_indices <- which(is_target_domain)

        if (length(target_indices) == 0) {
            message("    ⚠️  No target domain cells found in this image.")
            next
        }

        target_coords <- coords[target_indices, , drop = FALSE]

        # 2. Identify Query Cells (All cells)
        # We want to know distance of EVERY cell to the nearest Target Cell
        query_coords <- coords

        # 3. Calculate Nearest Neighbor Distance
        # knnx.dist returns distance to the k nearest neighbors in 'data' (target) for each point in 'query'
        # k=1 gives standard distance to the set
        nn_result <- FNN::get.knnx(data = target_coords, query = query_coords, k = 1)
        dist_to_domain <- as.vector(nn_result$nn.dist)

        # 4. Create Analysis Dataframe
        df <- data.frame(
            cell_id = rownames(meta),
            imageID = img,
            domain = meta[[domain_col]],
            cellType = meta[[celltype_col]],
            dist_to_domain = dist_to_domain,
            is_target_cell = is_target_domain,
            stringsAsFactors = FALSE
        )

        # Identify "Inside" vs "Outside"
        # If is_target_domain is TRUE, dist should be 0 (distance to self/neighbor in same domain).
        # We focus on the "Outside" or "Boundary" expansion.
        # But user said "From the boundary... expand out".
        # So we look at cells where domain != target_domain, and bin by distance.

        df_outside <- df %>%
            filter(domain != target_domain)

        if (nrow(df_outside) == 0) next

        # 5. Binning
        df_outside$bin <- cut(df_outside$dist_to_domain,
            breaks = seq(0, max_dist + bin_width, by = bin_width),
            include.lowest = TRUE
        )

        # 6. Summarize per Bin
        # Calculate proportion of Tregs in each bin
        bin_stats <- df_outside %>%
            filter(!is.na(bin)) %>%
            group_by(bin) %>%
            summarise(
                n_total = n(),
                n_treg = sum(grepl(cell_of_interest, cellType, ignore.case = TRUE)),
                # Calculate mean distance for plotting x-axis
                mean_dist = mean(dist_to_domain)
            ) %>%
            mutate(
                treg_prop = n_treg / n_total,
                imageID = img
            )

        results_list[[img]] <- bin_stats
    }

    full_stats <- bind_rows(results_list)
    return(full_stats)
}
