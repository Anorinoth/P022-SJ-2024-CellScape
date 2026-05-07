################################################################################
# clustering_analysis.R
# Author: Dr. Michael Zaiken
#
# Functions for unsupervised clustering and cell type annotation
################################################################################

#' Estimate optimal number of clusters using FuseSOM
#'
#' @param combined_spe SpatialExperiment object with normalized data
#' @param clustering_markers Character vector. Markers to use for clustering
#' @param res_dir Character. Results directory path
#' @param kSeq Numeric vector. Range of k values to test (default: 5:25)
#' @param initial_k Integer. Initial k for FuseSOM model (default: 15)
#' @param seed Integer. Random seed for reproducibility (default: 51773)
#'
#' @return List with cluster_est object and optimal_k_disc value
#'
#' @export
estimate_optimal_clusters <- function(combined_spe,
                                      clustering_markers,
                                      res_dir,
                                      kSeq = 5:25,
                                      initial_k = 15,
                                      seed = 51773) {
  require(FuseSOM)
  require(ggplot2)
  require(gridExtra)
  require(grid)
  require(S4Vectors)

  message("Estimating optimal number of clusters...")
  message("  Testing k = ", min(kSeq), " to ", max(kSeq))
  message("  This may take several minutes...")

  # Run initial FuseSOM model
  set.seed(seed)
  initial_fsom <- FuseSOM::runFuseSOM(
    data = combined_spe,
    markers = clustering_markers,
    assay = "norm",
    numClusters = initial_k,
    verbose = TRUE
  )

  message("✓ Initial FuseSOM model generated")

  # Estimate optimal k
  cluster_est <- FuseSOM::estimateNumCluster(
    data = initial_fsom,
    kSeq = kSeq,
    method = c("Discriminant", "Distance")
  )

  message("✓ Cluster estimation complete")

  # Get discriminant result from metadata
  cluster_estimation_results <- S4Vectors::metadata(cluster_est)$clusterEstimation
  optimal_k_disc <- cluster_estimation_results$Discriminant
  message("\nDiscriminant method suggests: ", optimal_k_disc, " clusters")

  # Generate elbow plots
  message("\nGenerating elbow plots...")

  p_slope <- FuseSOM::optiPlot(cluster_est, method = "slope") +
    labs(title = "Slope Statistic") +
    theme(plot.title = element_text(size = 12, face = "bold"))

  p_jump <- FuseSOM::optiPlot(cluster_est, method = "jump") +
    labs(title = "Jump Statistic") +
    theme(plot.title = element_text(size = 12, face = "bold"))

  p_wcd <- FuseSOM::optiPlot(cluster_est, method = "wcd") +
    labs(title = "Within-Cluster Dissimilarity") +
    theme(plot.title = element_text(size = 12, face = "bold"))

  p_gap <- FuseSOM::optiPlot(cluster_est, method = "gap") +
    labs(title = "Gap Statistic") +
    theme(plot.title = element_text(size = 12, face = "bold"))

  p_sil <- FuseSOM::optiPlot(cluster_est, method = "silhouette") +
    labs(title = "Silhouette Statistic") +
    theme(plot.title = element_text(size = 12, face = "bold"))

  # Combine plots
  combined_elbow <- grid.arrange(
    p_slope, p_jump, p_wcd, p_gap, p_sil,
    ncol = 2,
    top = textGrob("Cluster Number Estimation (Elbow Plots)",
      gp = gpar(fontsize = 14, fontface = "bold")
    )
  )

  # Save
  clustering_dir <- file.path(res_dir, "Clustering")
  dir.create(clustering_dir, showWarnings = FALSE, recursive = TRUE)

  elbow_file <- file.path(clustering_dir, "cluster_estimation_elbow_plots.png")
  ggsave(
    filename = elbow_file,
    plot = combined_elbow,
    width = 12,
    height = 10,
    dpi = 300,
    bg = "white"
  )

  message("✓ Elbow plots saved to: ", elbow_file)

  return(list(
    cluster_est = cluster_est,
    optimal_k_disc = optimal_k_disc,
    cluster_estimation_results = cluster_estimation_results,
    elbow_plot = combined_elbow
  ))
}


#' Run FuseSOM clustering with specified k
#'
#' @param combined_spe SpatialExperiment object
#' @param clustering_markers Character vector. Markers for clustering
#' @param optimal_k Integer. Number of clusters
#' @param seed Integer. Random seed (default: 51773)
#'
#' @return Updated SpatialExperiment with cluster assignments
#'
#' @export
run_fusesom_clustering <- function(combined_spe,
                                   clustering_markers,
                                   optimal_k,
                                   seed = 51773) {
  require(FuseSOM)
  require(SummarizedExperiment)

  message("\nRunning FuseSOM clustering with k = ", optimal_k, "...")

  # Run FuseSOM
  set.seed(seed)
  final_fsom <- FuseSOM::runFuseSOM(
    data = combined_spe,
    markers = clustering_markers,
    assay = "norm",
    numClusters = optimal_k,
    clusterCol = "clusters",
    verbose = TRUE
  )

  message("✓ FuseSOM clustering complete")

  # Check cluster distribution
  cluster_table <- table(colData(final_fsom)$clusters)
  message("\nCluster distribution (proportions):")
  print(round(cluster_table / sum(cluster_table), 3))

  return(final_fsom)
}


#' Generate cluster marker heatmaps
#'
#' @param combined_spe SpatialExperiment with cluster assignments
#' @param clustering_markers Character vector. Markers to visualize
#' @param res_dir Character. Results directory
#'
#' @return List with heatmap_data (long format) and heatmap_data_wide (wide format)
#'
#' @export
generate_cluster_heatmaps <- function(combined_spe,
                                      clustering_markers,
                                      res_dir) {
  require(FuseSOM)
  require(scater)
  require(SummarizedExperiment)
  require(dplyr)
  require(tidyr)

  clustering_dir <- file.path(res_dir, "Clustering")
  dir.create(clustering_dir, showWarnings = FALSE, recursive = TRUE)

  # Extract data for FuseSOM heatmap
  cluster_data <- as.data.frame(t(assay(combined_spe, "norm")))
  colnames(cluster_data) <- clustering_markers
  clusters <- colData(combined_spe)$clusters

  # Generate and save FuseSOM heatmap
  message("Generating FuseSOM marker heatmap...")
  heatmap_file <- file.path(clustering_dir, "cluster_marker_heatmap.png")
  png(heatmap_file, width = 12, height = 8, units = "in", res = 300)
  FuseSOM::markerHeatmap(
    data = cluster_data,
    markers = clustering_markers,
    clusters = clusters,
    clusterMarkers = TRUE
  )
  dev.off()
  message("✓ FuseSOM heatmap saved to: ", heatmap_file)

  # Generate and save scater heatmap
  message("Generating scater grouped heatmap...")
  heatmap_scater_file <- file.path(clustering_dir, "cluster_marker_heatmap_scater.png")
  png(heatmap_scater_file, width = 12, height = 10, units = "in", res = 300)
  heatmap_plot <- scater::plotGroupedHeatmap(
    combined_spe,
    features = clustering_markers,
    group = "clusters",
    exprs_values = "norm",
    center = TRUE,
    scale = TRUE,
    zlim = c(-3, 3),
    cluster_rows = FALSE,
    block = "clusters",
    main = "Marker Expression by Cluster (Z-scored)"
  )
  print(heatmap_plot)
  dev.off()
  message("✓ Scater heatmap saved to: ", heatmap_scater_file)

  # Calculate and export underlying data (mean expression per marker per cluster)
  # Memory-efficient approach: calculate statistics per cluster directly
  message("Calculating heatmap underlying data...")

  # Get cluster assignments (small vector)
  cluster_assignments <- colData(combined_spe)$clusters
  unique_clusters <- sort(unique(cluster_assignments))

  # Calculate stats for each marker-cluster combination
  heatmap_data_list <- list()
  idx <- 1

  for (marker in clustering_markers) {
    # Extract marker values (one marker at a time to save memory)
    marker_values <- assay(combined_spe[marker, ], "norm")

    for (cluster in unique_clusters) {
      # Get indices for this cluster
      cluster_idx <- cluster_assignments == cluster
      cluster_values <- marker_values[cluster_idx]

      # Calculate statistics
      heatmap_data_list[[idx]] <- data.frame(
        clusters = cluster,
        marker = marker,
        mean_expression = mean(cluster_values, na.rm = TRUE),
        median_expression = median(cluster_values, na.rm = TRUE),
        sd_expression = sd(cluster_values, na.rm = TRUE),
        n_cells = length(cluster_values),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }

  # Combine all results
  heatmap_data <- do.call(rbind, heatmap_data_list)

  # Export long format
  heatmap_data_file <- file.path(clustering_dir, "cluster_marker_heatmap_data.csv")
  write.csv(heatmap_data, heatmap_data_file, row.names = FALSE)

  # Also export wide format (clusters as rows, markers as columns)
  heatmap_data_wide <- heatmap_data %>%
    select(clusters, marker, mean_expression) %>%
    pivot_wider(
      names_from = marker,
      values_from = mean_expression
    )

  heatmap_data_wide_file <- file.path(clustering_dir, "cluster_marker_heatmap_data_wide.csv")
  write.csv(heatmap_data_wide, heatmap_data_wide_file, row.names = FALSE)

  message("✓ Heatmap data exported:")
  message("  - ", basename(heatmap_data_file), " (long format)")
  message("  - ", basename(heatmap_data_wide_file), " (wide format)")

  return(list(
    heatmap_data = heatmap_data,
    heatmap_data_wide = heatmap_data_wide
  ))
}


#' Annotate clusters with cell type labels
#'
#' @param combined_spe SpatialExperiment with cluster assignments
#' @param cluster_annotations Named character vector. Cluster -> cell type mapping
#' @param res_dir Character. Results directory
#'
#' @return Updated SpatialExperiment with cell_type column
#'
#' @export
annotate_cell_types <- function(combined_spe,
                                cluster_annotations,
                                res_dir) {
  require(SummarizedExperiment)

  # Check if clusters exist
  if (!"clusters" %in% colnames(colData(combined_spe))) {
    stop("Error: 'clusters' column not found. Please run clustering first.")
  }

  message("\nAnnotating clusters with cell types...")

  # Add annotations (robustly handle both "1" and "cluster_1" keys)
  cluster_ids <- as.character(combined_spe$clusters)

  # Check if annotations use "cluster_1" prefix while ids do not
  if (!any(cluster_ids %in% names(cluster_annotations)) &&
    any(paste0("cluster_", cluster_ids) %in% names(cluster_annotations))) {
    cluster_ids <- paste0("cluster_", cluster_ids)
    message("  Note: Adding 'cluster_' prefix to match annotation keys")
  }

  combined_spe$cell_type <- cluster_annotations[cluster_ids]

  message("✓ Cell type annotations added")
  message("\nCluster to Cell Type mapping:")
  for (i in seq_along(cluster_annotations)) {
    message(sprintf("  %s → %s", names(cluster_annotations)[i], cluster_annotations[i]))
  }

  # Create summary table
  cluster_counts <- table(combined_spe$clusters)
  cluster_annotation_summary <- data.frame(
    cluster = names(cluster_annotations),
    cell_type = unname(cluster_annotations),
    stringsAsFactors = FALSE
  )

  cluster_annotation_summary$n_cells <- sapply(names(cluster_annotations), function(cl) {
    if (cl %in% names(cluster_counts)) {
      return(cluster_counts[cl])
    } else {
      return(0)
    }
  })

  cluster_annotation_summary$percentage <- round(
    100 * cluster_annotation_summary$n_cells / ncol(combined_spe), 2
  )

  # Save annotation table
  clustering_dir <- file.path(res_dir, "Clustering")
  annotation_file <- file.path(clustering_dir, "cluster_cell_type_annotations.csv")
  write.csv(cluster_annotation_summary, annotation_file, row.names = FALSE)
  message("\n✓ Cluster annotations saved to: ", annotation_file)

  # Display summary
  print(cluster_annotation_summary)

  # Check for NA values
  n_na <- sum(is.na(combined_spe$cell_type))
  if (n_na > 0) {
    message("\n⚠️  Warning: ", n_na, " cells have NA cell_type annotation")
  } else {
    message("\n✓ All cells successfully annotated")
  }

  return(combined_spe)
}


#' Visualize cell types on UMAP and in spatial coordinates
#'
#' @param combined_spe SpatialExperiment with cell_type annotations
#' @param res_dir Character. Results directory
#' @param optimal_k Integer. Number of clusters for subtitle
#' @param celltype_colors Named vector. Custom colors for each cell type (optional)
#'
#' @return List with umap_plot and spatial_plots
#'
#' @export
visualize_celltype_umap <- function(combined_spe,
                                    res_dir,
                                    optimal_k,
                                    celltype_colors = NULL) {
  require(scater)
  require(ggplot2)

  message("Creating UMAP colored by cell type...")

  # Validate color mapping covers all cell types present in data
  all_cell_types <- sort(unique(combined_spe$cell_type))
  if (!is.null(celltype_colors)) {
    missing_colors <- setdiff(all_cell_types, names(celltype_colors))
    if (length(missing_colors) > 0) {
      message("  WARNING: No color defined for: ", paste(missing_colors, collapse = ", "))
      fallback_pal <- scales::hue_pal()(length(missing_colors))
      names(fallback_pal) <- missing_colors
      celltype_colors <- c(celltype_colors, fallback_pal)
    }
    # Remove colors for cell types no longer in the data (avoids phantom legend entries)
    celltype_colors <- celltype_colors[names(celltype_colors) %in% all_cell_types]
    message("  Color mapping validated: ", length(celltype_colors), " types with colors")
  } else {
    message("  Using default color palette")
  }

  umap_celltype_plot <- scater::plotReducedDim(
    combined_spe,
    dimred = "UMAP_norm",
    colour_by = "cell_type",
    point_size = 0.5,
    point_alpha = 0.7
  ) +
    labs(
      title = "UMAP: Cells Colored by Cell Type",
      subtitle = paste0(optimal_k, " clusters annotated based on marker expression"),
      color = "Cell Type"
    ) +
    theme_classic() +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray30"),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.4, "cm") # Compact spacing
    ) +
    guides(color = guide_legend(override.aes = list(size = 5))) # Large dots

  # Apply custom colors if provided
  if (!is.null(celltype_colors)) {
    umap_celltype_plot <- umap_celltype_plot + scale_color_manual(values = celltype_colors, drop = FALSE)
  }

  # Save
  clustering_dir <- file.path(res_dir, "Clustering")
  umap_celltype_file <- file.path(clustering_dir, "umap_colored_by_celltype.png")
  ggsave(
    filename = umap_celltype_file,
    plot = umap_celltype_plot,
    width = 12,
    height = 8,
    dpi = 300,
    bg = "white"
  )

  message("✓ UMAP saved to: ", umap_celltype_file)

  # Generate spatial plots colored by cell type (one per image)
  message("\nGenerating spatial cell type maps...")
  require(SpatialExperiment)
  require(gridExtra)

  unique_images <- sort(unique(combined_spe$image_name))
  spatial_plots <- list()
  spatial_plots_combined <- list()

  # Ensure consistent factor levels across all images so no cell types are dropped
  combined_spe$cell_type <- factor(combined_spe$cell_type, levels = all_cell_types)

  for (i in seq_along(unique_images)) {
    img_name <- unique_images[i]

    # Subset cells for this image
    img_idx <- combined_spe$image_name == img_name
    img_coords <- SpatialExperiment::spatialCoords(combined_spe)[img_idx, , drop = FALSE]
    img_celltypes <- combined_spe$cell_type[img_idx]

    # Create data frame
    plot_data <- data.frame(
      x = img_coords[, 1],
      y = img_coords[, 2],
      cell_type = img_celltypes
    )

    # Create spatial plot (for individual files - with legend)
    p_individual <- ggplot(plot_data, aes(x = x, y = y, color = cell_type)) +
      geom_point(size = 0.01, alpha = 0.5) + # Much smaller dots
      labs(
        title = img_name,
        subtitle = paste0(format(nrow(plot_data), big.mark = ","), " cells"),
        x = "X coordinate",
        y = "Y coordinate",
        color = "Cell Type"
      ) +
      theme_classic() +
      theme(
        aspect.ratio = 1,
        legend.position = "right",
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.35, "cm"), # Compact spacing
        plot.title = element_text(size = 11, face = "bold")
      ) +
      guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) + # Large dots
      coord_fixed()

    # Apply custom colors if provided (drop = FALSE keeps absent levels in legend)
    if (!is.null(celltype_colors)) {
      p_individual <- p_individual + scale_color_manual(values = celltype_colors, drop = FALSE)
    }

    spatial_plots[[img_name]] <- p_individual

    # Create version for combined plot (no legend - will add separately)
    p_no_legend <- ggplot(plot_data, aes(x = x, y = y, color = cell_type)) +
      geom_point(size = 0.01, alpha = 0.5) + # Much smaller dots
      labs(
        title = img_name,
        subtitle = paste0(format(nrow(plot_data), big.mark = ","), " cells"),
        x = "X coordinate",
        y = "Y coordinate",
        color = "Cell Type"
      ) +
      theme_classic() +
      theme(
        aspect.ratio = 1,
        legend.position = "none",
        plot.title = element_text(size = 11, face = "bold")
      ) +
      coord_fixed()

    # Apply custom colors if provided
    if (!is.null(celltype_colors)) {
      p_no_legend <- p_no_legend + scale_color_manual(values = celltype_colors, drop = FALSE)
    }

    spatial_plots_combined[[img_name]] <- p_no_legend

    # Save individual image
    spatial_file <- file.path(clustering_dir, paste0(
      "spatial_celltype_",
      gsub("[^A-Za-z0-9_]", "_", img_name), ".png"
    ))
    ggsave(spatial_file, plot = p_individual, width = 10, height = 8, dpi = 300, bg = "white")
    message("  ✓ Saved spatial plot: ", basename(spatial_file))
  }

  # Save combined plot with shared legend
  if (length(spatial_plots_combined) > 1) {
    require(gridExtra)
    require(grid)

    # Extract legend from first individual plot (already has colors applied)
    legend_plot <- p_individual +
      theme(
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.35, "cm") # Compact spacing
      ) +
      guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) # Large dots

    # Function to extract legend
    get_legend <- function(a_plot) {
      tmp <- ggplot_gtable(ggplot_build(a_plot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      if (length(leg) > 0) {
        legend <- tmp$grobs[[leg]]
      } else {
        legend <- NULL
      }
      return(legend)
    }

    legend <- get_legend(legend_plot)

    # Arrange plots: 2 plots in a row + legend on the right
    combined_spatial <- grid.arrange(
      arrangeGrob(grobs = spatial_plots_combined, ncol = 2),
      legend,
      ncol = 2,
      widths = c(10, 2),
      top = textGrob("Spatial Distribution of Cell Types",
        gp = gpar(fontsize = 14, fontface = "bold")
      )
    )

    combined_spatial_file <- file.path(clustering_dir, "spatial_celltype_all_images.png")
    ggsave(combined_spatial_file,
      plot = combined_spatial,
      width = 18, height = 8, dpi = 300, bg = "white"
    )
    message("  ✓ Saved combined spatial plot: ", basename(combined_spatial_file))
  }

  message("✓ Spatial cell type maps complete")

  return(list(
    umap_plot = umap_celltype_plot,
    spatial_plots = spatial_plots
  ))
}


#' Compare cell type composition between images
#'
#' @param combined_spe SpatialExperiment with cell_type and image_name
#' @param res_dir Character. Results directory
#' @param optimal_k Integer. Number of clusters
#'
#' @return List with comparison data frame and plot
#'
#' @export
compare_celltype_composition <- function(combined_spe,
                                         res_dir,
                                         optimal_k) {
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(SummarizedExperiment)

  message("\nComparing cell type composition between images...")

  # Calculate proportions per image
  celltype_by_image <- colData(combined_spe) %>%
    as.data.frame() %>%
    group_by(image_name, cell_type) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(image_name) %>%
    mutate(
      total = sum(n),
      percentage = 100 * n / total
    ) %>%
    ungroup()

  # Pivot for comparison
  celltype_comparison <- celltype_by_image %>%
    select(image_name, cell_type, percentage) %>%
    pivot_wider(
      names_from = image_name,
      values_from = percentage,
      values_fill = 0
    )

  # Calculate difference
  img_names <- colnames(celltype_comparison)[2:3]
  celltype_comparison$difference <- celltype_comparison[[img_names[2]]] -
    celltype_comparison[[img_names[1]]]
  celltype_comparison$abs_difference <- abs(celltype_comparison$difference)
  celltype_comparison <- celltype_comparison %>% arrange(desc(abs_difference))

  message("\nCell types ranked by difference:")
  print(celltype_comparison)

  # Export
  clustering_dir <- file.path(res_dir, "Clustering")
  celltype_comp_file <- file.path(clustering_dir, "celltype_composition_by_image.csv")
  write.csv(celltype_comparison, celltype_comp_file, row.names = FALSE)
  message("\n✓ Saved to: ", basename(celltype_comp_file))

  # Visualize
  celltype_plot <- celltype_by_image %>%
    ggplot(aes(x = reorder(cell_type, -percentage), y = percentage, fill = image_name)) +
    geom_col(position = "dodge", width = 0.7) +
    labs(
      title = "Cell Type Composition: Comparison Between Images",
      subtitle = paste0("Total cell types: ", optimal_k),
      x = "Cell Type",
      y = "Percentage of Cells (%)",
      fill = "Image"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "top"
    ) +
    scale_fill_brewer(palette = "Set2")

  # Save plot
  celltype_comp_plot_file <- file.path(clustering_dir, "celltype_composition_comparison.png")
  ggsave(celltype_comp_plot_file,
    plot = celltype_plot, width = 14, height = 8,
    dpi = 300, bg = "white"
  )
  message("✓ Plot saved to: ", basename(celltype_comp_plot_file))

  # Top differences
  top_diff_celltypes <- head(celltype_comparison, 5)
  message("\nTop 5 cell types with largest differences:")
  for (i in 1:nrow(top_diff_celltypes)) {
    ct <- top_diff_celltypes$cell_type[i]
    diff <- round(top_diff_celltypes$difference[i], 2)
    direction <- ifelse(diff > 0, paste0("higher in ", img_names[2]),
      paste0("higher in ", img_names[1])
    )
    message("  ", ct, ": ", abs(diff), "% ", direction)
  }

  return(list(
    comparison = celltype_comparison,
    plot = celltype_plot
  ))
}


#' Export marker statistics by cell type
#'
#' @param combined_spe SpatialExperiment with cell_type annotations
#' @param clustering_markers Character vector. Markers to analyze
#' @param res_dir Character. Results directory
#'
#' @return Data frame with marker statistics
#'
#' @export
export_marker_statistics <- function(combined_spe,
                                     clustering_markers,
                                     res_dir) {
  require(dplyr)
  require(tidyr)
  require(SummarizedExperiment)

  message("\nCalculating marker statistics per cell type...")

  # Memory-efficient approach: calculate statistics directly
  # Get cell type assignments (small vector)
  celltype_assignments <- colData(combined_spe)$cell_type
  unique_celltypes <- sort(unique(celltype_assignments))

  # Calculate stats for each marker-celltype combination
  marker_celltype_list <- list()
  idx <- 1

  for (marker in clustering_markers) {
    # Extract marker values (one marker at a time to save memory)
    marker_values <- assay(combined_spe[marker, ], "norm")

    for (celltype in unique_celltypes) {
      # Get indices for this cell type
      celltype_idx <- celltype_assignments == celltype
      celltype_values <- marker_values[celltype_idx]

      # Calculate statistics
      marker_celltype_list[[idx]] <- data.frame(
        cell_type = celltype,
        marker = marker,
        mean_expression = mean(celltype_values, na.rm = TRUE),
        median_expression = median(celltype_values, na.rm = TRUE),
        sd_expression = sd(celltype_values, na.rm = TRUE),
        q25 = quantile(celltype_values, 0.25, na.rm = TRUE),
        q75 = quantile(celltype_values, 0.75, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }

  # Combine all results
  marker_by_celltype <- do.call(rbind, marker_celltype_list)

  # Pivot to wide format
  marker_wide <- marker_by_celltype %>%
    select(cell_type, marker, mean_expression) %>%
    pivot_wider(
      names_from = marker,
      values_from = mean_expression
    )

  message("✓ Calculated statistics for ", length(clustering_markers), " markers")

  # Export
  clustering_dir <- file.path(res_dir, "Clustering")
  marker_celltype_file <- file.path(clustering_dir, "marker_expression_by_celltype.csv")
  write.csv(marker_by_celltype, marker_celltype_file, row.names = FALSE)

  marker_wide_file <- file.path(clustering_dir, "marker_expression_by_celltype_wide.csv")
  write.csv(marker_wide, marker_wide_file, row.names = FALSE)

  message("\n✓ Data tables exported:")
  message("  - ", basename(marker_celltype_file), " (long format)")
  message("  - ", basename(marker_wide_file), " (wide format)")

  # Show top markers
  message("\nTop expressed marker per cell type:")
  top_markers <- marker_by_celltype %>%
    group_by(cell_type) %>%
    slice_max(mean_expression, n = 1) %>%
    select(cell_type, marker, mean_expression)
  print(top_markers)

  return(marker_by_celltype)
}


#' Analyze FoxP3+ enrichment by cell type
#'
#' @param combined_spe SpatialExperiment with FoxP3_positive and cell_type
#' @param res_dir Character. Results directory
#'
#' @return List with enrichment data and plot
#'
#' @export
analyze_foxp3_enrichment <- function(combined_spe,
                                     res_dir) {
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(SummarizedExperiment)

  message("\n🔬 FOXP3+ ENRICHMENT BY CELL TYPE")

  # Calculate FoxP3+ percentage per cell type and image
  foxp3_by_celltype <- colData(combined_spe) %>%
    as.data.frame() %>%
    group_by(cell_type, image_name) %>%
    summarise(
      n_cells = n(),
      foxp3_positive = sum(FoxP3_positive),
      foxp3_percentage = 100 * foxp3_positive / n_cells,
      mean_foxp3_expr = mean(FoxP3_expression),
      .groups = "drop"
    ) %>%
    arrange(desc(foxp3_percentage))

  message("\nTop 5 FoxP3+ enriched cell types:")
  print(head(foxp3_by_celltype, 5))

  # Export
  clustering_dir <- file.path(res_dir, "Clustering")
  foxp3_celltype_file <- file.path(clustering_dir, "FOXP3_enrichment_by_celltype.csv")
  write.csv(foxp3_by_celltype, foxp3_celltype_file, row.names = FALSE)
  message("✓ FoxP3 enrichment saved to: ", basename(foxp3_celltype_file))

  # Visualize
  foxp3_enrichment_plot <- ggplot(
    foxp3_by_celltype,
    aes(
      x = reorder(cell_type, foxp3_percentage),
      y = foxp3_percentage, fill = image_name
    )
  ) +
    geom_col(position = "dodge") +
    coord_flip() +
    labs(
      title = "FoxP3+ Treg Enrichment by Cell Type",
      subtitle = "Which cell types contain FoxP3+ Tregs?",
      x = "Cell Type",
      y = "FoxP3+ Percentage",
      fill = "Image"
    ) +
    theme_classic() +
    theme(
      legend.position = "top",
      axis.text.y = element_text(size = 9)
    )

  ggsave(file.path(clustering_dir, "FOXP3_enrichment_by_celltype.png"),
    foxp3_enrichment_plot,
    width = 14, height = 10, dpi = 300, bg = "white"
  )

  # Compare between images
  foxp3_celltype_comparison <- foxp3_by_celltype %>%
    select(cell_type, image_name, foxp3_percentage) %>%
    pivot_wider(
      names_from = image_name, values_from = foxp3_percentage,
      values_fill = 0
    ) %>%
    mutate(difference = abs(.[[3]] - .[[2]])) %>%
    arrange(desc(difference))

  message("\nCell types with largest FoxP3+ % differences:")
  print(head(foxp3_celltype_comparison, 5))

  # Export comparison
  foxp3_celltype_comp_file <- file.path(clustering_dir, "FOXP3_celltype_comparison.csv")
  write.csv(foxp3_celltype_comparison, foxp3_celltype_comp_file, row.names = FALSE)
  message("✓ Comparison saved to: ", basename(foxp3_celltype_comp_file))

  message("\n🎯 KEY FINDING: Identify FoxP3+ Treg cell types")
  message("  Expected: High FoxP3+ % in Regulatory T cells")

  return(list(
    enrichment = foxp3_by_celltype,
    comparison = foxp3_celltype_comparison,
    plot = foxp3_enrichment_plot
  ))
}
