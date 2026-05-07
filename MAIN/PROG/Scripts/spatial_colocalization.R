################################################################################
# spatial_colocalization.R
#
# Functions for spatial co-localization analysis using spicyR
################################################################################

#' Prepare SpatialExperiment for spatial analysis
#'
#' @param combined_spe SpatialExperiment object
#'
#' @return Updated SpatialExperiment with imageID and cellType columns
#'
#' @export
prepare_spatial_data <- function(combined_spe) {
  require(SpatialExperiment)
  require(SummarizedExperiment)

  message("Checking SpatialExperiment structure...")

  # Verify cell_type column exists (from Section 4 annotation)
  if (!"cell_type" %in% colnames(colData(combined_spe))) {
    stop("No 'cell_type' column found. Please run clustering and annotation first (Section 4.6).")
  }

  # Verify imageID column (should be 'image_name' from our data)
  if ("image_name" %in% colnames(colData(combined_spe))) {
    # Rename to imageID for spicyR compatibility
    colData(combined_spe)$imageID <- colData(combined_spe)$image_name
    message("✓ Renamed 'image_name' to 'imageID'")
  } else if (!"imageID" %in% colnames(colData(combined_spe))) {
    stop("No 'imageID' or 'image_name' column found in colData")
  }

  # Use cell_type (not clusters) for spicyR cellType
  colData(combined_spe)$cellType <- colData(combined_spe)$cell_type
  message("✓ Using cell type annotations (not cluster IDs)")

  # Check spatial coordinates
  coords <- spatialCoords(combined_spe)
  message(
    "✓ Spatial coordinates: ", ncol(coords), " dimensions (",
    paste(colnames(coords), collapse = ", "), ")"
  )

  # Summary
  n_cells <- ncol(combined_spe)
  n_images <- length(unique(combined_spe$imageID))
  n_cell_types <- length(unique(combined_spe$cellType))

  message("\n=== Data Summary ===")
  message("Total cells: ", format(n_cells, big.mark = ","))
  message("Images: ", n_images, " (", paste(unique(combined_spe$imageID), collapse = ", "), ")")
  message("Cell types: ", n_cell_types)

  # Cell type distribution
  cell_type_counts <- table(combined_spe$cellType)
  message("\nCells per cell type:")
  print(cell_type_counts)

  message("\n✓ Data ready for spatial analysis")

  return(combined_spe)
}


#' Calculate L-function for all cell type pairs
#'
#' @param combined_spe SpatialExperiment with cellType and imageID columns
#' @param coloc_dir Character. Colocalization output directory
#' @param radii Numeric vector. Radii to test (default: c(20, 50, 100))
#' @param sigma Numeric or NULL. Tissue inhomogeneity correction (default: NULL)
#' @param seed Integer. Random seed for reproducibility
#'
#' @return List with lfunction_results and output directory
#'
#' @export
calculate_lfunction <- function(combined_spe,
                                coloc_dir,
                                radii = c(20, 50, 100),
                                sigma = NULL,
                                seed = 51773) {
  require(spicyR)
  require(BiocParallel)

  message("Calculating L-function for all cell type pairs...")
  message("  Using radii: ", paste(radii, collapse = ", "), " pixels")
  if (!is.null(sigma)) {
    message("  Using sigma = ", sigma, " for tissue inhomogeneity correction")
  }
  message("  Output directory: ", coloc_dir)

  # Set seed for reproducibility
  set.seed(seed)

  # Use do.call to force evaluation of arguments before calling getPairwise
  # This resolves scoping issues where radii or sigma might not be found by
  # internal parallel workers.
  args <- list(
    cells = combined_spe,
    r = radii,
    sigma = sigma,
    cellType = "cellType",
    imageID = "imageID",
    cores = 1
  )

  message("  Calling spicyR::getPairwise...")
  pairwise_matrix <- do.call(spicyR::getPairwise, args)

  # Convert matrix to list of dataframes (one per image) as expected by downstream
  lfunction_results <- split(as.data.frame(pairwise_matrix), rownames(pairwise_matrix))

  message("  ✓ Calculated co-localization for ", length(lfunction_results), " images")
  message("  ✓ Cell type pairs analyzed: ", ncol(pairwise_matrix))

  message("\n✓ L-function calculation complete")

  return(list(
    lfunction_results = lfunction_results,
    coloc_dir = coloc_dir
  ))
}


#' Compare co-localization between images
#'
#' @param lfunction_results List. Output from calculate_lfunction
#' @param coloc_dir Character. Output directory path
#'
#' @return List with summary data frames and long format data
#'
#' @export
compare_colocalizations <- function(lfunction_results, coloc_dir) {
  require(dplyr)
  require(tidyr)

  # Combine results into a single dataframe
  lfunction_df <- do.call(rbind, lapply(names(lfunction_results), function(img) {
    df <- lfunction_results[[img]]
    df$imageID <- img
    df
  }))

  # Convert to long format for easier comparison
  lfunction_long <- lfunction_df %>%
    pivot_longer(
      cols = -imageID,
      names_to = "cell_pair",
      values_to = "L_function"
    )

  message("Cell type pairs analyzed: ", length(unique(lfunction_long$cell_pair)))

  # Calculate mean L-function for each pair across both images
  lfunction_summary <- lfunction_long %>%
    group_by(cell_pair) %>%
    summarise(
      mean_L = mean(L_function, na.rm = TRUE),
      sd_L = sd(L_function, na.rm = TRUE),
      diff_L = diff(L_function), # Difference between two images
      .groups = "drop"
    ) %>%
    arrange(desc(abs(mean_L)))

  # Top 10 most co-localised pairs (positive L-function)
  message("\nTop 10 most CO-LOCALISED cell type pairs:")
  print(head(lfunction_summary %>% arrange(desc(mean_L)), 10))

  # Top 10 most dispersed pairs (negative L-function)
  message("\nTop 10 most DISPERSED cell type pairs:")
  print(head(lfunction_summary %>% arrange(mean_L), 10))

  # Top 10 pairs with biggest difference between images
  message("\nTop 10 pairs with LARGEST DIFFERENCE between images:")
  top_different_pairs <- lfunction_summary %>% arrange(desc(abs(diff_L)))
  print(head(top_different_pairs, 10))

  # Export per-image L-function values (long format)
  lfunction_long_export <- lfunction_long %>%
    # Parse cell pair into from/to
    mutate(
      from_cellType = sapply(strsplit(cell_pair, "__"), `[`, 1),
      to_cellType = sapply(strsplit(cell_pair, "__"), `[`, 2)
    )

  lfunction_long_file <- file.path(coloc_dir, "colocalisation_per_image.csv")
  write.csv(lfunction_long_export, lfunction_long_file, row.names = FALSE)
  message("\n✓ Per-image L-function values saved to: ", basename(lfunction_long_file))

  return(list(
    lfunction_summary = lfunction_summary,
    lfunction_long = lfunction_long,
    lfunction_long_export = lfunction_long_export,
    top_different_pairs = top_different_pairs
  ))
}


#' Visualize co-localization differences between images
#'
#' @param lfunction_summary Data frame. Summary from compare_colocalizations
#' @param lfunction_long_export Data frame. Long format data with from/to columns
#' @param coloc_dir Character. Output directory
#' @param n_top Integer. Number of top pairs to visualize (default: 20)
#'
#' @return ggplot object
#'
#' @export
visualize_coloc_differences <- function(lfunction_summary,
                                        lfunction_long_export,
                                        coloc_dir,
                                        n_top = 20) {
  require(ggplot2)

  message("\nVisualizing pairs with largest differences between images...")

  # Get top N most different pairs
  top_n_diff <- head(lfunction_summary %>% arrange(desc(abs(diff_L))), n_top)

  # Add direction info
  img_names <- unique(lfunction_long_export$imageID)
  top_n_diff$direction <- ifelse(
    top_n_diff$diff_L > 0,
    paste0("Higher in ", img_names[2]),
    paste0("Higher in ", img_names[1])
  )

  # Create bar plot
  diff_plot <- ggplot(top_n_diff, aes(
    x = reorder(cell_pair, abs(diff_L)),
    y = diff_L,
    fill = direction
  )) +
    geom_col() +
    coord_flip() +
    labs(
      title = paste0("Top ", n_top, " Cell Type Pairs: Largest Co-localization Differences"),
      subtitle = paste0(img_names[1], " vs ", img_names[2]),
      x = "Cell Type Pair",
      y = "Difference in L-function (Image 2 - Image 1)",
      fill = "Enrichment"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 8)
    ) +
    scale_fill_manual(values = c("steelblue", "coral"))

  print(diff_plot)

  # Save
  diff_plot_file <- file.path(coloc_dir, "colocalisation_differences_between_images.png")
  ggsave(diff_plot_file, plot = diff_plot, width = 12, height = 8, dpi = 300, bg = "white")
  message("✓ Difference plot saved to: ", basename(diff_plot_file))

  # Export top differences
  top_diff_file <- file.path(coloc_dir, "colocalisation_top_differences.csv")
  write.csv(top_n_diff, top_diff_file, row.names = FALSE)
  message("✓ Top differences table saved to: ", basename(top_diff_file))

  return(diff_plot)
}


#' Generate co-localization heatmaps
#'
#' @param lfunction_results List. Output from calculate_lfunction
#' @param combined_spe SpatialExperiment object
#' @param coloc_dir Character. Output directory
#'
#' @return List of heatmap grobs
#'
#' @export
generate_coloc_heatmaps <- function(lfunction_results, combined_spe, coloc_dir) {
  require(pheatmap)
  require(gridExtra)

  message("Creating co-localization heatmaps...")

  heatmap_list <- list()
  image_ids <- names(lfunction_results)

  for (img_id in image_ids) {
    # Get L-function matrix for this image
    lf_data <- lfunction_results[[img_id]]

    # Convert to matrix format for heatmap
    # Extract "from" and "to" cell types from column names
    pair_names <- colnames(lf_data)
    cell_types_unique <- unique(combined_spe$cellType)

    # Create a matrix
    n_types <- length(cell_types_unique)
    lf_matrix <- matrix(NA,
      nrow = n_types, ncol = n_types,
      dimnames = list(cell_types_unique, cell_types_unique)
    )

    # Fill matrix
    for (pair in pair_names) {
      parts <- strsplit(pair, "__")[[1]]
      if (length(parts) == 2) {
        from_type <- parts[1]
        to_type <- parts[2]
        if (from_type %in% cell_types_unique && to_type %in% cell_types_unique) {
          lf_matrix[from_type, to_type] <- lf_data[[pair]][1]
        }
      }
    }

    # Create heatmap
    p <- pheatmap::pheatmap(
      lf_matrix,
      main = paste("L-function Heatmap:", img_id),
      color = colorRampPalette(c("blue", "white", "red"))(100),
      breaks = seq(-50, 50, length.out = 101),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      na_col = "grey",
      fontsize = 8,
      silent = TRUE
    )

    heatmap_list[[img_id]] <- p[[4]] # Extract the grob
  }

  # Combine heatmaps
  gridExtra::grid.arrange(grobs = heatmap_list, ncol = 2)

  # Save
  heatmap_file <- file.path(coloc_dir, "colocalisation_heatmaps.png")
  png(heatmap_file, width = 14, height = 12, units = "in", res = 300)
  gridExtra::grid.arrange(grobs = heatmap_list, ncol = 2)
  dev.off()

  message("✓ Co-localisation heatmaps saved to: ", basename(heatmap_file))

  return(heatmap_list)
}


#' Visualize top co-localized pairs spatially
#'
#' @param lfunction_summary Data frame. Summary from compare_colocalizations
#' @param combined_spe SpatialExperiment object
#' @param coloc_dir Character. Output directory
#' @param n_pairs Integer. Number of top pairs to visualize (default: 3)
#'
#' @return NULL (plots are saved)
#'
#' @export
visualize_top_pairs <- function(lfunction_summary, combined_spe, coloc_dir, n_pairs = 3, priority_type = "Transgenic Tregs") {
  require(spicyR)
  require(gridExtra)
  require(dplyr)

  # Split cell_pair into from and to if not already present
  lfunction_summary <- lfunction_summary %>%
    mutate(
      from_type = sapply(strsplit(as.character(cell_pair), "__"), `[`, 1),
      to_type = sapply(strsplit(as.character(cell_pair), "__"), `[`, 2)
    )

  # Filter out self-interactions
  other_pairs <- lfunction_summary %>%
    filter(from_type != to_type) %>%
    arrange(desc(abs(mean_L)))

  # Deduplicate mirror interactions (e.g., A->B and B->A)
  other_pairs <- other_pairs %>%
    mutate(sorted_pair = sapply(1:n(), function(i) {
      paste(sort(c(from_type[i], to_type[i])), collapse = "__")
    })) %>%
    distinct(sorted_pair, .keep_all = TRUE)

  # Select top general non-self pairs
  top_general <- head(other_pairs, n_pairs)

  # Select top interactions with priority cell type (non-self)
  priority_interactions <- other_pairs %>%
    filter(grepl(priority_type, from_type) | grepl(priority_type, to_type)) %>%
    filter(!cell_pair %in% top_general$cell_pair) %>%
    arrange(desc(abs(mean_L))) %>%
    head(3) # Limit to top 3 priority interactions to keep output manageable

  # Combine
  top_pairs <- rbind(top_general, priority_interactions)

  message("Visualizing selected cell type pairs (skipping self-interactions):")
  print(top_pairs[, c("cell_pair", "mean_L", "diff_L")])

  image_ids <- unique(combined_spe$imageID)

  # For each top pair, create spatial plots for both images
  for (i in 1:nrow(top_pairs)) {
    pair_name <- top_pairs$cell_pair[i]
    from_type <- top_pairs$from_type[i]
    to_type <- top_pairs$to_type[i]

    message("\nPlotting: ", from_type, " → ", to_type)

    # Create plots for both images
    plot_list <- list()
    for (img_id in image_ids) {
      p <- spicyR::plotImage(
        combined_spe,
        imageToPlot = img_id,
        from = from_type,
        to = to_type,
        cellType = "cellType",
        imageID = "imageID"
      ) +
        ggtitle(paste0(img_id, ": ", from_type, " → ", to_type))

      plot_list[[img_id]] <- p
    }

    # Display side-by-side
    gridExtra::grid.arrange(grobs = plot_list, ncol = 2)

    # Save
    safe_pair_name <- gsub("[^A-Za-z0-9_]", "_", pair_name)
    pair_file <- file.path(coloc_dir, paste0("spatial_", safe_pair_name, ".png"))
    ggsave(
      filename = pair_file,
      plot = gridExtra::arrangeGrob(grobs = plot_list, ncol = 2),
      width = 14,
      height = 6,
      dpi = 300,
      bg = "white"
    )

    message("  ✓ Saved to: ", basename(pair_file))
  }

  message("\n✓ Top pair visualizations complete")
}


#' Calculate tissue inhomogeneity corrected L-function
#'
#' @param combined_spe SpatialExperiment object
#' @param lfunction_results List. Uncorrected results from calculate_lfunction
#' @param coloc_dir Character. Output directory
#' @param sigma Numeric. Tissue inhomogeneity correction parameter
#' @param radii Numeric vector. Radii to test
#'
#' @return List with corrected results and comparison data
#'
#' @export
calculate_inhomogeneity_corrected <- function(combined_spe,
                                              lfunction_results,
                                              coloc_dir,
                                              sigma = 20,
                                              radii = c(20, 50, 100)) {
  require(spicyR)
  require(BiocParallel)
  require(dplyr)
  require(tidyr)

  message("Recalculating L-function with tissue inhomogeneity correction...")
  message("  Using sigma = ", sigma, " pixels")

  # Use do.call to force evaluation of arguments
  args <- list(
    cells = combined_spe,
    r = radii,
    sigma = sigma,
    cellType = "cellType",
    imageID = "imageID",
    cores = 1
  )

  pairwise_corrected <- do.call(spicyR::getPairwise, args)

  # Convert to list format
  lfunction_corrected <- split(as.data.frame(pairwise_corrected), rownames(pairwise_corrected))

  message("\n✓ Corrected L-function calculation complete")

  # Compare uncorrected vs corrected for top 3 pairs
  # Get top 3 pairs from uncorrected results
  lfunction_df <- do.call(rbind, lapply(names(lfunction_results), function(img) {
    df <- lfunction_results[[img]]
    df$imageID <- img
    df
  }))

  lfunction_long <- lfunction_df %>%
    pivot_longer(cols = -imageID, names_to = "cell_pair", values_to = "L_function")

  lfunction_summary <- lfunction_long %>%
    group_by(cell_pair) %>%
    summarise(mean_L = mean(L_function, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(abs(mean_L)))

  comparison_pairs <- head(lfunction_summary$cell_pair, 3)

  message("\nEffect of tissue inhomogeneity correction (top 3 pairs):")
  correction_comparison <- data.frame()

  for (pair in comparison_pairs) {
    if (pair %in% colnames(lfunction_results[[image_ids[1]]])) {
      message("\n", pair, ":")
      for (img_id in image_ids) {
        uncorr <- lfunction_results[[img_id]][[pair]][1]
        corr <- lfunction_corrected[[img_id]][[pair]][1]
        change <- corr - uncorr
        message(
          "  ", img_id, ": ", round(uncorr, 2), " → ", round(corr, 2),
          " (Δ = ", round(change, 2), ")"
        )

        # Store for export
        correction_comparison <- rbind(correction_comparison, data.frame(
          cell_pair = pair,
          imageID = img_id,
          uncorrected = uncorr,
          corrected = corr,
          change = change,
          percent_change = 100 * change / (abs(uncorr) + 0.01)
        ))
      }
    }
  }

  # Export correction comparison
  correction_file <- file.path(coloc_dir, "colocalisation_inhomogeneity_correction.csv")
  write.csv(correction_comparison, correction_file, row.names = FALSE)
  message("\n✓ Inhomogeneity correction comparison saved to: ", basename(correction_file))

  # Export corrected L-function values
  lfunction_corrected_long <- do.call(rbind, lapply(names(lfunction_corrected), function(img) {
    df <- lfunction_corrected[[img]]
    df_long <- df %>%
      pivot_longer(cols = everything(), names_to = "cell_pair", values_to = "L_function_corrected")
    df_long$imageID <- img
    df_long
  }))

  corrected_file <- file.path(coloc_dir, "colocalisation_per_image_corrected.csv")
  write.csv(lfunction_corrected_long, corrected_file, row.names = FALSE)
  message("✓ Corrected L-function values saved to: ", basename(corrected_file))

  return(list(
    lfunction_corrected = lfunction_corrected,
    correction_comparison = correction_comparison,
    lfunction_corrected_long = lfunction_corrected_long
  ))
}


#' Analyze FoxP3+ Treg spatial relationships
#'
#' @param combined_spe SpatialExperiment object
#' @param lfunction_long_export Data frame. Long format L-function data with from/to columns
#' @param coloc_dir Character. Output directory
#' @param foxp3_marker Character. FoxP3 marker name (default: "31_CorrectedFOXP3-GFP")
#' @param n_top_types Integer. Number of top FoxP3+ cell types to identify (default: 4)
#' @param n_top_pairs Integer. Number of top pairs to visualize (default: 20)
#'
#' @return List with FoxP3 results and plot
#'
#' @export
analyze_foxp3_colocalization <- function(combined_spe,
                                         lfunction_long_export,
                                         coloc_dir,
                                         foxp3_marker = "31_CorrectedFOXP3-GFP",
                                         n_top_types = 4,
                                         n_top_pairs = 20) {
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(SummarizedExperiment)

  message("\n🔬 FOXP3+ TREG CO-LOCALIZATION")

  if (!foxp3_marker %in% rownames(combined_spe)) {
    message("⚠️  FoxP3 marker (", foxp3_marker, ") not found in data.")
    return(NULL)
  }

  # Calculate mean FoxP3 expression per cell type
  foxp3_by_celltype <- colData(combined_spe) %>%
    as.data.frame() %>%
    mutate(foxp3_expression = assay(combined_spe[foxp3_marker, ], "norm")[1, ]) %>%
    group_by(cell_type, imageID) %>%
    summarise(
      mean_foxp3 = mean(foxp3_expression, na.rm = TRUE),
      n_cells = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_foxp3))

  message("\nMean FoxP3 expression by cell type:")
  print(foxp3_by_celltype)

  # Identify FoxP3+ enriched cell types (top N cell types)
  foxp3_celltypes <- foxp3_by_celltype %>%
    group_by(cell_type) %>%
    summarise(overall_mean_foxp3 = mean(mean_foxp3), .groups = "drop") %>%
    arrange(desc(overall_mean_foxp3)) %>%
    head(n_top_types) %>%
    pull(cell_type)

  message(
    "\nFoxP3+ enriched cell types (top ", n_top_types, "): ",
    paste(foxp3_celltypes, collapse = ", ")
  )

  # Extract L-function results for FoxP3+ cell types
  foxp3_coloc_results <- lfunction_long_export %>%
    filter(from_cellType %in% foxp3_celltypes | to_cellType %in% foxp3_celltypes) %>%
    arrange(desc(abs(L_function)))

  message("\nTotal FoxP3+ cell type pairs found: ", nrow(foxp3_coloc_results))
  message("\nTop 10 cell type pairs involving FoxP3+ cell types:")
  print(head(foxp3_coloc_results, 10))

  # Export
  foxp3_coloc_file <- file.path(coloc_dir, "FOXP3_colocalization.csv")
  write.csv(foxp3_coloc_results, foxp3_coloc_file, row.names = FALSE)
  message("✓ FoxP3+ co-localization saved to: ", basename(foxp3_coloc_file))

  # Compare FoxP3+ Treg co-localization between images
  foxp3_coloc_comparison <- foxp3_coloc_results %>%
    mutate(cell_pair_label = paste(from_cellType, to_cellType, sep = " → ")) %>%
    select(imageID, cell_pair_label, L_function) %>%
    tidyr::pivot_wider(names_from = imageID, values_from = L_function, values_fill = 0)

  # Calculate difference if we have 2 images
  if (ncol(foxp3_coloc_comparison) == 3) {
    foxp3_coloc_comparison <- foxp3_coloc_comparison %>%
      mutate(difference = abs(.[[3]] - .[[2]])) %>%
      arrange(desc(difference))
  }

  message("\nTop 5 FoxP3+ cell type pairs with largest differences:")
  print(head(foxp3_coloc_comparison, 5))

  # Export comparison
  foxp3_coloc_comp_file <- file.path(coloc_dir, "FOXP3_colocalization_comparison.csv")
  write.csv(foxp3_coloc_comparison, foxp3_coloc_comp_file, row.names = FALSE)
  message("✓ FoxP3+ co-localization comparison saved to: ", basename(foxp3_coloc_comp_file))

  # Visualize top FoxP3+ Treg interactions (bar plot)
  top_foxp3_pairs <- head(foxp3_coloc_results, n_top_pairs)

  foxp3_coloc_plot <- ggplot(
    top_foxp3_pairs,
    aes(
      x = reorder(paste(from_cellType, to_cellType, sep = " → "), L_function),
      y = L_function, fill = imageID
    )
  ) +
    geom_col(position = "dodge") +
    coord_flip() +
    labs(
      title = "FoxP3+ Cell Type Co-localization",
      subtitle = paste0("Top ", n_top_pairs, " cell type pairs involving FoxP3+ enriched cell types (r = 100 pixels)"),
      x = "Cell Type Pair",
      y = "L-function (Positive = Co-localized, Negative = Dispersed)",
      fill = "Image"
    ) +
    theme_classic() +
    theme(legend.position = "top")

  print(foxp3_coloc_plot)

  ggsave(file.path(coloc_dir, "FOXP3_colocalization_barplot.png"),
    foxp3_coloc_plot,
    width = 12, height = 10, dpi = 300, bg = "white"
  )

  message("\n✓ FoxP3+ Treg co-localization analysis complete")
  message("\n🎯 KEY QUESTION: Which cell types are near FoxP3+ enriched cell types?")
  message("  Positive L = co-localization (cells cluster together)")
  message("  Negative L = dispersion (cells avoid each other)")
  message("  Expected in spleen: Tregs near T zones, GC borders, possibly near B cells")

  return(list(
    foxp3_celltypes = foxp3_celltypes,
    foxp3_coloc_results = foxp3_coloc_results,
    foxp3_coloc_comparison = foxp3_coloc_comparison,
    foxp3_coloc_plot = foxp3_coloc_plot
  ))
}


#' Tile SpatialExperiment object
#'
#' Splits the spatial domain into a grid of tiles, effectively creating
#' pseudo-replicates for analysis. This is crucial for n=1 comparisons.
#'
#' @param spe SpatialExperiment object
#' @param n_rows Number of rows in the grid (default: 4)
#' @param n_cols Number of columns in the grid (default: 4)
#' @param min_cells Minimum cells required to keep a tile (default: 50)
#' @param image_col Column containing original image IDs
#'
#' @return SpatialExperiment with updated imageID column reflecting tiles
#'
#' @export
tile_spatial_experiment <- function(spe, n_rows = 4, n_cols = 4, min_cells = 50, image_col = "imageID") {
  message("Tiling SpatialExperiment: ", n_rows, "x", n_cols, " grid")

  # Get data
  cd <- as.data.frame(SummarizedExperiment::colData(spe))
  coords <- SpatialExperiment::spatialCoords(spe)

  # Check if we have multiple images
  if (!image_col %in% colnames(cd)) {
    stop("Image column '", image_col, "' not found")
  }

  # Process each image separately
  image_ids <- unique(cd[[image_col]])

  new_image_ids <- rep(NA, nrow(cd))
  names(new_image_ids) <- rownames(cd)

  total_tiles_created <- 0
  kept_tiles <- 0

  for (img in image_ids) {
    # Subset cells for this image
    img_mask <- cd[[image_col]] == img
    img_cells <- rownames(cd)[img_mask]

    if (length(img_cells) == 0) next

    # Get bounds for this image
    img_coords <- coords[img_mask, , drop = FALSE]
    x_min <- min(img_coords[, 1])
    x_max <- max(img_coords[, 1])
    y_min <- min(img_coords[, 2])
    y_max <- max(img_coords[, 2])

    x_step <- (x_max - x_min) / n_cols
    y_step <- (y_max - y_min) / n_rows

    # Assign tiles
    for (r in 1:n_rows) {
      for (c in 1:n_cols) {
        total_tiles_created <- total_tiles_created + 1

        # Define tile bounds
        tile_x_min <- x_min + (c - 1) * x_step
        tile_x_max <- x_min + c * x_step
        tile_y_min <- y_min + (r - 1) * y_step
        tile_y_max <- y_min + r * y_step

        # Find cells in this tile
        in_tile <- which(
          img_coords[, 1] >= tile_x_min & img_coords[, 1] <= tile_x_max &
            img_coords[, 2] >= tile_y_min & img_coords[, 2] <= tile_y_max
        )

        if (length(in_tile) >= min_cells) {
          # Create unique tile ID: ImageID_Tile_R_C
          tile_id <- paste0(img, "_Tile_", r, "_", c)

          # Update image IDs for these cells
          new_image_ids[img_cells[in_tile]] <- tile_id
          kept_tiles <- kept_tiles + 1
        }
      }
    }
  }

  # Filter cells that fell into valid tiles
  keep_mask <- !is.na(new_image_ids)

  message("  Tiles created: ", kept_tiles, "/", total_tiles_created)
  message("  Cells retained: ", sum(keep_mask), "/", length(keep_mask))

  # Subset object
  spe_tiled <- spe[, keep_mask]

  # Update imageID in the subsetted object
  SummarizedExperiment::colData(spe_tiled)[[image_col]] <- new_image_ids[keep_mask]

  # Update factor levels
  SummarizedExperiment::colData(spe_tiled)[[image_col]] <- factor(SummarizedExperiment::colData(spe_tiled)[[image_col]])

  return(spe_tiled)
}


#' Run Spicy Colocalization Analysis
#'
#' @param spe SpatialExperiment object
#' @param condition_col Column in colData specifying treatment groups
#' @param subject_col Column in colData specifying pseudo-replicates (tiles)
#' @param radii Radii to test
#' @param sigma Tissue inhomogeneity correction
#' @param min_replicates Minimum replicates per condition to trigger tiling (default: 3)
#'
#' @return Spicy results object
#'
#' @export
run_spicy_analysis <- function(spe,
                               condition_col = "condition",
                               subject_col = "imageID",
                               radii = 100,
                               sigma = NULL,
                               min_replicates = 3) {
  require(spicyR)

  # Check condition column
  if (!condition_col %in% colnames(colData(spe))) {
    stop("Condition column '", condition_col, "' not found in SPE")
  }

  # Check sample size per condition
  sample_counts <- table(spe[[condition_col]][!duplicated(spe[[subject_col]])])
  min_samples <- min(sample_counts)

  message("Samples per condition:")
  print(sample_counts)

  if (min_samples < min_replicates) {
    message("\nWARNING: Low sample size detected (n < ", min_replicates, ").")
    message("  Enabling TILING strategy to generate pseudo-replicates...")

    spe <- tile_spatial_experiment(spe, n_rows = 4, n_cols = 4, min_cells = 50, image_col = subject_col)

    # Re-check counts
    new_counts <- table(spe[[condition_col]][!duplicated(spe[[subject_col]])])
    message("  New '", subject_col, "' counts (tiles):")
    print(new_counts)
  }

  # Run spicy
  message("\nRunning spicy differential colocalization...")
  args <- list(
    cells = spe,
    condition = condition_col,
    subject = subject_col,
    r = radii,
    sigma = sigma,
    cores = 1
  )
  spicy_test <- do.call(spicyR::spicy, args)

  return(spicy_test)
}
