################################################################################
# spatial_domain_analysis.R
#
# Functions for analyzing spatial relationships within specific domains
# using spicyR and pseudotiling.
################################################################################

#' Run Spicy Analysis Per Spatial Domain
#'
#' Iterates through each spatial domain, subsets the data, and runs spicy
#' with pseudotiling to identify condition-specific colocalization patterns
#' within that domain.
#'
#' @param combined_spe SpatialExperiment object
#' @param domain_col Column name for spatial domains (default: "annotation")
#' @param condition_col Column name for conditions (default: "condition")
#' @param output_dir Directory to save results
#' @param radii Radii for L-function (default: 100)
#' @param min_cells Minimum cells in a domain to analyze (default: 200)
#'
#' @return List of spicy results per domain
#'
#' @export
run_domain_spicy_analysis <- function(combined_spe,
                                      domain_col = "annotation",
                                      condition_col = "condition",
                                      output_dir = "Results/SpatialDomains/Colocalization",
                                      radii = 100,
                                      min_cells = 200) {
    require(spicyR)
    require(ggplot2)
    require(dplyr)
    require(SpatialExperiment)

    # Ensure output directory exists
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    message("Starting Spatial Domain Analysis...")
    message("  Domain column: ", domain_col)
    message("  Condition column: ", condition_col)
    message("  Output directory: ", output_dir)

    # Get unique domains
    if (!domain_col %in% colnames(colData(combined_spe))) {
        stop("Domain column '", domain_col, "' not found in colData")
    }

    domains <- unique(combined_spe[[domain_col]])
    message("  Domains found: ", paste(domains, collapse = ", "))

    results_list <- list()

    for (dom in domains) {
        clean_dom_name <- gsub("[^A-Za-z0-9_]", "_", dom)
        message("\n=== Analyzing Domain: ", dom, " ===")

        # Subset SPE for this domain
        spe_sub <- combined_spe[, combined_spe[[domain_col]] == dom]

        n_cells <- ncol(spe_sub)
        message("  Cells in domain: ", n_cells)

        if (n_cells < min_cells) {
            message("  ⚠️  Skipping domain (too few cells, < ", min_cells, ")")
            next
        }

        # Run spicy analysis (using existing function which handles tiling)
        # We assume run_spicy_analysis is available in the environment
        # (sourced from spatial_colocalization.R)

        spicy_result <- tryCatch(
            {
                run_spicy_analysis(
                    spe = spe_sub,
                    condition_col = condition_col,
                    subject_col = "imageID", # Use imageID for tiling
                    radii = radii
                )
            },
            error = function(e) {
                message("  ⚠️  Error running spicy for domain ", dom, ": ", e$message)
                return(NULL)
            }
        )

        if (!is.null(spicy_result)) {
            results_list[[dom]] <- spicy_result

            # 1. Save Significance Plot
            p_sig <- spicyR::signifPlot(spicy_result, breaks = c(-1, 0, 1)) +
                labs(
                    title = paste("Significant Interactions in", dom),
                    subtitle = paste("Comparison:", paste(levels(spe_sub[[condition_col]]), collapse = " vs "))
                )

            plot_file <- file.path(output_dir, paste0("spicy_signifPlot_", clean_dom_name, ".png"))
            ggsave(plot_file, p_sig, width = 10, height = 8, bg = "white")
            message("  ✓ Significance plot saved: ", basename(plot_file))

            # 2. Save Statistics Table
            stats <- as.data.frame(spicyR::topPairs(spicy_result, n = Inf))
            stats$domain <- dom
            stats$cell_pair <- paste(stats$from, stats$to, sep = " -> ")

            csv_file <- file.path(output_dir, paste0("spicy_stats_", clean_dom_name, ".csv"))
            write.csv(stats, csv_file, row.names = FALSE)
            message("  ✓ Stats saved: ", basename(csv_file))
        }
    }

    message("\n✓ Domain analysis complete.")
    return(results_list)
}

#' Calculate Shared Boundaries Between Spatial Domains
#'
#' Uses Voronoi tessellation to quantify the length of shared boundaries
#' between different spatial domains.
#'
#' @param combined_spe SpatialExperiment object
#' @param imageID_col Column name for image IDs
#' @param domain_col Column name for spatial domains
#' @param output_dir Directory to save results
#' @param BPPARAM BiocParallel back-end to use (default: SerialParam())
#'
#' @return List containing summary stats and plot data
#' @export
calculate_domain_boundaries <- function(combined_spe,
                                        imageID_col = "imageID",
                                        domain_col = "annotation",
                                        output_dir = "Results/SpatialDomains/Boundaries",
                                        min_patch_area = 0,
                                        min_boundary_length = 0,
                                        BPPARAM = BiocParallel::SerialParam()) {
    require(deldir)
    require(dplyr)
    require(SpatialExperiment)
    require(BiocParallel)

    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    message("Calculating domain boundaries using Voronoi tessellation...")

    # Process each image separately
    image_ids <- unique(combined_spe[[imageID_col]])

    # Prepare lightweight data for each image
    # This avoids exporting the massive combined_spe object to workers
    image_data_list <- list()
    for (img in image_ids) {
        spe_sub <- combined_spe[, combined_spe[[imageID_col]] == img]
        coords <- spatialCoords(spe_sub)

        # Create a lightweight dataframe with only necessary columns
        df <- data.frame(
            x = coords[, 1],
            y = coords[, 2],
            domain = spe_sub[[domain_col]],
            imageID = img,
            stringsAsFactors = FALSE
        )
        image_data_list[[img]] <- df
    }

    # Define worker function that takes a DATAFRAME, not an image ID
    process_single_image_df <- function(df) {
        # Re-load necessary libraries inside the worker
        require(deldir)
        require(dplyr)
        require(igraph)

        img <- unique(df$imageID)
        message("  Processing image: ", img)

        # Remove duplicate coordinates to prevent index mismatch with deldir
        dup_idx <- duplicated(df[, c("x", "y")])
        if (any(dup_idx)) {
            message("    Removing ", sum(dup_idx), " duplicate cells to avoid Voronoi index mismatch")
            df <- df[!dup_idx, ]
        }

        # Fail immediately if too few points for tessellation
        stopifnot(
            "Image has fewer than 10 cells — cannot compute Voronoi tessellation" =
                nrow(df) >= 10
        )

        # Calculate Voronoi Tessellation (fail-fast: no tryCatch)
        vt <- deldir::deldir(df$x, df$y, suppressMsge = TRUE)

        # Always build contiguous domain topological patch identifiers
        # Extract cell areas from Voronoi summary
        cell_areas <- vt$summary$dir.area

        # Extract Delaunay edges to find adjacent cells
        delsgs <- vt$delsgs

        # Create a dataframe of adjacent cells and their domains
        adj_cells <- data.frame(
            ind1 = delsgs$ind1,
            ind2 = delsgs$ind2,
            dom1 = df$domain[delsgs$ind1],
            dom2 = df$domain[delsgs$ind2],
            stringsAsFactors = FALSE
        )

        # Filter to edges connecting the SAME domain
        same_dom_edges <- adj_cells %>% filter(dom1 == dom2)

        # We need all vertices (cells), even isolated ones, for the graph
        nodes <- data.frame(id = 1:nrow(df), domain = df$domain, area = cell_areas, stringsAsFactors = FALSE)
        edges_graph <- data.frame(from = same_dom_edges$ind1, to = same_dom_edges$ind2, stringsAsFactors = FALSE)

        # Create graph and find connected components (patches)
        g <- igraph::graph_from_data_frame(edges_graph, directed = FALSE, vertices = nodes)
        comp <- igraph::components(g)

        # Create unique patch ID per domain and grouping
        nodes$patch_id <- paste0(nodes$domain, "_", comp$membership)

        # Merge individual patch IDs into our primary spatial tracking object!
        df$patch_id <- nodes$patch_id

        # Sum area by patch for potential small-zone exclusion
        patch_areas <- nodes %>%
            group_by(patch_id) %>%
            summarise(total_area = sum(area, na.rm = TRUE), .groups = "drop")

        # Define raw un-filtered patch tracking
        raw_patch_data <- nodes %>%
            group_by(patch_id, domain) %>%
            summarise(total_area = sum(area, na.rm = TRUE), .groups = "drop") %>%
            mutate(imageID = img)

        # Handle min_patch_area filtering
        if (min_patch_area > 0) {
            # Identify small patches
            small_patches <- patch_areas %>%
                filter(total_area < min_patch_area) %>%
                pull(patch_id)

            # Reassign domain for cells in small patches
            cells_to_exclude <- nodes$id[nodes$patch_id %in% small_patches]
            if (length(cells_to_exclude) > 0) {
                message("    Excluding ", length(cells_to_exclude), " cells in patches < ", min_patch_area, " area")
                df$domain[cells_to_exclude] <- "Excluded_Small_Region"
            }
        }

        # Per-cell patch mapping (post min_patch_area filter) for downstream density calculations
        cell_patches_data <- data.frame(
            imageID      = img,
            x            = df$x,
            y            = df$y,
            domain       = df$domain,
            patch_id     = df$patch_id,
            voronoi_area = cell_areas,
            stringsAsFactors = FALSE
        )

        # Get Voronoi edges
        edges <- as.data.frame(vt$dirsgs)

        # Calculate edge lengths
        edges$length <- sqrt((edges$x1 - edges$x2)^2 + (edges$y1 - edges$y2)^2)

        # Track the raw edge properties prior to threshold dropping
        raw_edges_data <- edges %>%
            mutate(
                domain1 = df$domain[ind1],
                domain2 = df$domain[ind2],
                patch1 = df$patch_id[ind1],
                patch2 = df$patch_id[ind2],
                imageID = img
            ) %>%
            select(domain1, domain2, patch1, patch2, length, imageID)

        # Optional: Filter tiny cell-to-cell boundary segments
        if (min_boundary_length > 0) {
            edges <- edges[edges$length >= min_boundary_length, ]
            if (nrow(edges) == 0) {
                return(NULL)
            }
        }

        # Map indices to domains
        dom1 <- df$domain[edges$ind1]
        dom2 <- df$domain[edges$ind2]

        # Classify edges
        edges$domain1 <- dom1
        edges$domain2 <- dom2
        edges$type <- ifelse(dom1 == dom2, "Internal", "Boundary")

        # Create a sorted pair identifier
        edges$pair <- ifelse(edges$type == "Boundary",
            paste(pmin(dom1, dom2), pmax(dom1, dom2), sep = " <-> "),
            paste(dom1, "Internal", sep = " - ")
        )

        # Prepare plotting data
        img_plot_data <- edges %>%
            mutate(imageID = img) %>%
            filter(length < 200)

        # Summarize boundary lengths
        edge_long <- bind_rows(
            edges %>% select(domain = domain1, type, length, pair),
            edges %>% select(domain = domain2, type, length, pair)
        )

        domain_perimeters <- edge_long %>%
            filter(type == "Boundary") %>%
            group_by(domain) %>%
            summarise(total_perimeter = sum(length), .groups = "drop")

        shared_stats <- bind_rows(
            edges %>% filter(type == "Boundary") %>%
                group_by(source_domain = domain1, target_domain = domain2) %>%
                summarise(length = sum(length), .groups = "drop"),
            edges %>% filter(type == "Boundary") %>%
                group_by(source_domain = domain2, target_domain = domain1) %>%
                summarise(length = sum(length), .groups = "drop")
        ) %>%
            group_by(source_domain, target_domain) %>%
            summarise(shared_length = sum(length), .groups = "drop")

        img_stats <- shared_stats %>%
            left_join(domain_perimeters, by = c("source_domain" = "domain")) %>%
            mutate(
                pct_shared_boundary = (shared_length / total_perimeter) * 100,
                imageID = img
            )

        # Calculate Boundary Enrichment Score
        # 1. Total tissue perimeter (sum of all domain perimeters represents "available surface" to touch)
        total_available_perimeter <- sum(domain_perimeters$total_perimeter)

        # 2. Add enrichment score
        img_stats <- img_stats %>%
            left_join(domain_perimeters %>% select(domain, target_perimeter = total_perimeter),
                by = c("target_domain" = "domain")
            ) %>%
            mutate(
                # Available Perimeter = Space available for Source to touch = Total - Source
                denom_perimeter = total_available_perimeter - total_perimeter,

                # Handle edge case where source is the only domain
                denom_perimeter = ifelse(denom_perimeter <= 0, NA, denom_perimeter),

                # Expected Prob = Target's share of the available surface
                expected_prob = target_perimeter / denom_perimeter,
                observed_prob = pct_shared_boundary / 100,
                enrichment_score = observed_prob / expected_prob
            ) %>%
            select(-denom_perimeter, -expected_prob, -observed_prob)

        return(list(
            stats = img_stats,
            plot_data = img_plot_data,
            raw_patches = if (exists("raw_patch_data")) raw_patch_data else NULL,
            raw_edges = raw_edges_data,
            cell_patches = cell_patches_data
        ))
    }

    # Run in parallel using the DATA LIST
    results_list <- BiocParallel::bplapply(image_data_list, process_single_image_df, BPPARAM = BPPARAM)
    names(results_list) <- image_ids

    # Combine results
    boundary_stats_list <- lapply(results_list, function(x) x$stats)
    plot_data_list <- lapply(results_list, function(x) x$plot_data)
    raw_patches_list <- lapply(results_list, function(x) x$raw_patches)
    raw_edges_list <- lapply(results_list, function(x) x$raw_edges)
    cell_patches_list <- lapply(results_list, function(x) x$cell_patches)

    # Filter out NULLs
    boundary_stats_list <- boundary_stats_list[!sapply(boundary_stats_list, is.null)]
    plot_data_list <- plot_data_list[!sapply(plot_data_list, is.null)]
    raw_patches_list <- raw_patches_list[!sapply(raw_patches_list, is.null)]
    raw_edges_list <- raw_edges_list[!sapply(raw_edges_list, is.null)]
    cell_patches_list <- cell_patches_list[!sapply(cell_patches_list, is.null)]

    # Bind stats
    full_stats <- bind_rows(boundary_stats_list)

    # Save statistics
    write.csv(full_stats, file.path(output_dir, "spatial_domain_boundary_stats.csv"), row.names = FALSE)

    return(list(
        stats = full_stats,
        plot_data = plot_data_list,
        raw_patches  = if (length(raw_patches_list)  > 0) bind_rows(raw_patches_list)  else NULL,
        raw_edges    = if (length(raw_edges_list)    > 0) bind_rows(raw_edges_list)    else NULL,
        cell_patches = if (length(cell_patches_list) > 0) bind_rows(cell_patches_list) else NULL
    ))
}
