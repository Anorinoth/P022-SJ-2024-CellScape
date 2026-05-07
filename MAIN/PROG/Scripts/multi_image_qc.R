#' Quality Control and Processing Wrapper for Multiple Images
#'
#' This script contains wrapper functions to:
#' 1. Load and Harmonize multiple SPE objects
#' 2. Normalize cells across images
#' 3. Visualize QC metrics
#'
#' @import SpatialExperiment
#' @import scater
#' @import simpleSeg
#' @import ggplot2
#' @import gridExtra
#' @import dplyr
#' @import tidyr
#' @import qs
#' @import BiocParallel

# Load libraries
library(SpatialExperiment)
library(scater)
# simpleSeg might not be loaded yet, we'll require it in function
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(qs)
library(BiocParallel)

# Helper function for marker name harmonization
harmonize_marker_names <- function(spe_list, image_names, debug = FALSE) {
  # Standardize marker name function
  standardize_marker_name <- function(name) {
    name <- gsub(" ", "_", name) # Spaces to underscores
    name <- gsub("\\.", "_", name) # Dots to underscores (including Cy5.5)
    name <- gsub("\\(", "__", name) # ( to __
    name <- gsub("\\)", "_", name) # ) to _
    name <- gsub("[^A-Za-z0-9_-]", "_", name) # Other special chars
    name <- gsub("_{2,}", "_", name) # Clean up multiple underscores
    name <- gsub("_+$", "", name) # Remove trailing underscores
    return(name)
  }

  # Get markers from first SPE (reference)
  markers_ref <- rownames(spe_list[[1]])
  markers_ref_std <- sapply(markers_ref, standardize_marker_name, USE.NAMES = FALSE)

  # Check and harmonize each subsequent SPE
  for (i in 2:length(spe_list)) {
    markers_i <- rownames(spe_list[[i]])

    if (!identical(markers_ref, markers_i)) {
      message(
        "\n⚠️  Marker mismatch detected between ", image_names[1],
        " and ", image_names[i], " - harmonizing..."
      )
      message("  ", image_names[1], " markers: ", length(markers_ref))
      message("  ", image_names[i], " markers: ", length(markers_i))

      # Standardize
      markers_i_std <- sapply(markers_i, standardize_marker_name, USE.NAMES = FALSE)

      # Debug output
      if (debug) {
        message("\n  DEBUG: Standardization examples:")
        mismatched_indices <- which(markers_ref != markers_i)
        for (idx in head(mismatched_indices, 3)) {
          message("    ", image_names[1], ": '", markers_ref[idx], "' -> '", markers_ref_std[idx], "'")
          message("    ", image_names[i], ": '", markers_i[idx], "' -> '", markers_i_std[idx], "'")
          message("    Match: ", markers_ref_std[idx] == markers_i_std[idx])
        }
      }

      # Check if standardized names match
      if (identical(markers_ref_std, markers_i_std)) {
        message("  ✓ Marker names are equivalent after standardization")
        message("  Keeping original names from ", image_names[1])

        # Reorder SPE i and rename to match reference
        # We perform temporary renaming to facilitate reordering
        rownames(spe_list[[i]]) <- markers_i_std
        spe_list[[i]] <- spe_list[[i]][markers_ref_std, ]
        rownames(spe_list[[i]]) <- markers_ref

        message("  ✓ Both SPE objects now use consistent marker names")
      } else {
        # True mismatch - different markers present
        message("  ⚠️  True marker mismatch - subsetting to common markers")

        # Find common markers
        common_markers_std <- intersect(markers_ref_std, markers_i_std)

        # Find markers unique to each
        only_in_ref <- markers_ref[!markers_ref_std %in% common_markers_std]
        only_in_i <- markers_i[!markers_i_std %in% common_markers_std]

        if (length(only_in_ref) > 0) {
          message(
            "  Markers only in ", image_names[1], ": ",
            paste(only_in_ref, collapse = ", ")
          )
        }
        if (length(only_in_i) > 0) {
          message(
            "  Markers only in ", image_names[i], ": ",
            paste(only_in_i, collapse = ", ")
          )
        }

        message("  Common markers: ", length(common_markers_std))

        # Subset both SPE objects to common markers
        # Use logical indexing based on standardized names
        keep_ref <- markers_ref_std %in% common_markers_std
        spe_list[[1]] <- spe_list[[1]][keep_ref, ]

        # For SPE i, we need to subset and reorder
        keep_i <- markers_i_std %in% common_markers_std
        spe_list[[i]] <- spe_list[[i]][keep_i, ]

        # Update current markers after subsetting
        markers_ref_curr <- rownames(spe_list[[1]])
        markers_ref_std_curr <- markers_ref_std[keep_ref]

        markers_i_curr <- rownames(spe_list[[i]])
        markers_i_std_curr <- markers_i_std[keep_i]

        # Rename SPE i temporarily to allow matching
        rownames(spe_list[[i]]) <- markers_i_std_curr

        # Reorder SPE i to match SPE 1 (common markers order)
        spe_list[[i]] <- spe_list[[i]][markers_ref_std_curr, ]

        # Finally, assign the reference names to both to ensure exact match
        rownames(spe_list[[i]]) <- markers_ref_curr

        # Update the reference markers for next iteration
        markers_ref <- markers_ref_curr
        markers_ref_std <- markers_ref_std_curr

        message("  ✓ Both SPE objects subset to ", length(common_markers_std), " common markers")
      }

      message("  ✓ Marker names and order now match")
    } else {
      message("\n✓ Marker names match between ", image_names[1], " and ", image_names[i])
    }
  }

  return(spe_list)
}

#' Load and harmonize multiple SPE objects
#'
#' @param prog_dir Character. Project PROG directory
#' @param n_expected Integer. Number of expected images (default: 2)
#' @param debug Logical. Enable debug messages (default: FALSE)
#'
#' @return Combined SpatialExperiment object
#'
#' @export
load_and_harmonize_spe <- function(prog_dir,
                                   n_expected = 2,
                                   debug = FALSE) {
  # Auto-detect and load SPE files
  # Recursive search to find files in SEGMENTATION/{image_name} subdirs
  cells_files <- list.files(prog_dir, pattern = "_cells\\.qs$", full.names = TRUE, recursive = TRUE)

  # Filter out test files
  test_files <- grep("test", cells_files, ignore.case = TRUE, value = TRUE)
  if (length(test_files) > 0) {
    message("Excluding ", length(test_files), " test files from aggregation:")
    message("  ", paste(basename(test_files), collapse = "\n  "))
    cells_files <- setdiff(cells_files, test_files)
  }

  # Resolution of duplicate files (same basename, different paths)
  # This happens if a file is in PROG and also in PROG/SEGMENTATION/...
  file_basenames <- basename(cells_files)
  if (any(duplicated(file_basenames))) {
    dup_names <- unique(file_basenames[duplicated(file_basenames)])
    message("\n⚠️  Duplicate files detected for ", length(dup_names), " images. Resolving...")

    files_to_keep <- c()
    files_to_remove <- c()

    for (dup in dup_names) {
      candidates <- cells_files[basename(cells_files) == dup]

      # Prioritize SEGMENTATION directory
      seg_candidates <- grep("SEGMENTATION", candidates, value = TRUE)

      if (length(seg_candidates) == 1) {
        # Clear winner
        chosen <- seg_candidates
      } else if (length(seg_candidates) > 1) {
        # Multiple in SEGMENTATION? Pick the one with longer path (deeper)
        chosen <- seg_candidates[which.max(nchar(seg_candidates))]
      } else {
        # None in SEGMENTATION? Pick the one with longer path
        chosen <- candidates[which.max(nchar(candidates))]
      }

      discarded <- setdiff(candidates, chosen)
      files_to_keep <- c(files_to_keep, chosen)
      files_to_remove <- c(files_to_remove, discarded)

      message("  For ", dup, ":")
      message("    Keeping:   ", chosen)
      message("    Discarding: ", paste(discarded, collapse = "\n                "))
    }

    # Reconstruct cells_files: keep non-duplicates + resolved duplicates
    non_dups <- cells_files[!file_basenames %in% dup_names]
    cells_files <- c(non_dups, files_to_keep)
  }

  if (length(cells_files) != n_expected) {
    stop(
      "Expected ", n_expected, " _cells.qs files (excluding tests), found ", length(cells_files),
      "\n  Files: ", paste(cells_files, collapse = "\n         ")
    )
  }

  message("Loading SpatialExperiment objects...")
  for (i in seq_along(cells_files)) {
    message("  File ", i, ": ", basename(cells_files[i]))
  }

  # Load all SPE objects
  spe_list <- lapply(cells_files, qs::qread)

  # Get image names
  image_names <- sapply(spe_list, function(spe) unique(spe$image_name))

  message("\n✓ Loaded:")
  for (i in seq_along(spe_list)) {
    message(
      "  Image ", i, ": ", image_names[i],
      " (", format(ncol(spe_list[[i]]), big.mark = ","), " cells, ",
      nrow(spe_list[[i]]), " markers)"
    )
  }

  # Harmonize marker names across all images
  if (n_expected > 1) {
    spe_list <- harmonize_marker_names(spe_list, image_names, debug = debug)
  }

  # Combine into a single SPE object
  combined_spe <- do.call(cbind, spe_list)

  message("\n✓ Combined SPE:")
  message("  Total cells: ", format(ncol(combined_spe), big.mark = ","))
  message("  Markers: ", nrow(combined_spe))
  message("  Images: ", paste(unique(combined_spe$image_name), collapse = ", "))

  return(combined_spe)
}


#' Generate UMAP and density plots
#'
#' @param combined_spe SpatialExperiment object
#' @param res_dir Character. Results directory
#' @param ct_markers Character vector. Markers to use for cell typing/UMAP
#' @param n_cores Integer. Number of cores
#'
#' @return Updated SpatialExperiment with UMAP
#'
#' @export
generate_qc_plots <- function(combined_spe,
                              res_dir,
                              ct_markers = NULL,
                              n_cores = NULL) {
  require(scater)
  require(ggplot2)

  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 2)
  }

  if (is.null(ct_markers)) {
    message("No cell type markers provided, using all markers for UMAP")
    ct_markers <- rownames(combined_spe)
  }

  message("Running UMAP on ", length(ct_markers), " markers...")

  BPPARAM <- MulticoreParam(workers = n_cores)

  # Run UMAP
  combined_spe <- scater::runUMAP(
    combined_spe,
    subset_row = ct_markers,
    exprs_values = "counts",
    BPPARAM = BPPARAM,
    n_threads = n_cores
  )

  message("✓ UMAP complete")

  # Create and save plot
  qc_dir <- file.path(res_dir, "QC_plots")
  dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

  umap_plot <- scater::plotReducedDim(
    combined_spe,
    dimred = "UMAP",
    colour_by = "image_name",
    point_size = 0.5,
    point_alpha = 0.6
  ) +
    ggplot2::labs(
      title = "UMAP: Cells Colored by Image",
      subtitle = "Clustering by image indicates batch effects"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 11, color = "gray30")
    )

  umap_file <- file.path(qc_dir, "umap_batch_effects.png")
  ggplot2::ggsave(
    filename = umap_file,
    plot = umap_plot,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
  )

  message("✓ UMAP plot saved to: ", umap_file)

  # Store ct_markers as attribute for later use
  attr(combined_spe, "ct_markers") <- ct_markers
  attr(combined_spe, "BPPARAM") <- BPPARAM
  attr(combined_spe, "n_cores") <- n_cores

  return(combined_spe)
}


#' Normalize cell intensities across images
#'
#' Applies normalization steps (trim99, mean, PC1) per-image for batch correction.
#' Implements the three normalization steps directly (no simpleSeg dependency):
#'   - trim99: clip values at the 99th percentile per marker
#'   - mean:   divide each marker by its mean
#'   - PC1:    remove the first principal component
#'
#' @param combined_spe SpatialExperiment object
#' @param exclude_patterns Character vector. Markers to exclude (default: "DNA")
#' @param methods Character vector. Normalization methods (default: c("trim99", "mean", "PC1"))
#' @param n_cores Integer. Number of cores (unused currently, kept for API compat)
#'
#' @return Updated SpatialExperiment with 'norm' assay
#'
#' @export
normalize_cells_wrapper <- function(combined_spe,
                                    exclude_patterns = "DNA",
                                    methods = c("trim99", "mean", "PC1"),
                                    n_cores = NULL) {
    if (is.null(n_cores)) {
        n_cores <- max(1, parallel::detectCores() - 2)
    }

    all_markers <- rownames(combined_spe)

    # Exclude markers
    exclude_regex <- paste(exclude_patterns, collapse = "|")
    exclude_markers <- all_markers[grepl(exclude_regex, all_markers, ignore.case = TRUE)]

    if (length(exclude_markers) > 0) {
        message("Excluding markers from normalization: ", paste(exclude_markers, collapse = ", "))
        use_markers <- all_markers[!all_markers %in% exclude_markers]
    } else {
        message("No markers excluded, normalizing all markers")
        use_markers <- all_markers
    }

    message("Normalizing ", length(use_markers), " markers across images...")

    # Initialize 'norm' assay with original counts
    assay(combined_spe, "norm", withDimnames = FALSE) <- assay(combined_spe, "counts")

    # Check if image_name column exists
    if (!"image_name" %in% names(colData(combined_spe))) {
        stop("Column 'image_name' not found in combined_spe")
    }
    image_ids <- unique(combined_spe$image_name)

    message("  Normalizing per-image (", length(image_ids), " images)...")

    # Work on a matrix copy for used markers
    full_norm_subset <- assay(combined_spe, "counts")[use_markers, ]

    for (img in image_ids) {
        cell_idx <- which(combined_spe$image_name == img)
        if (length(cell_idx) == 0) next

        message("  Processing image: ", img, " (", format(length(cell_idx), big.mark = ","), " cells)")

        # Extract matrix: markers x cells
        mat <- full_norm_subset[, cell_idx, drop = FALSE]

        # Apply each normalization method in order
        for (m in methods) {
            if (m == "trim99") {
                # Clip each marker at its 99th percentile
                for (r in seq_len(nrow(mat))) {
                    q99 <- quantile(mat[r, ], 0.99, na.rm = TRUE)
                    mat[r, mat[r, ] > q99] <- q99
                }
                message("    ✓ trim99 applied")
            } else if (m == "mean") {
                # Divide each marker by its mean (avoid division by zero)
                for (r in seq_len(nrow(mat))) {
                    m_val <- mean(mat[r, ], na.rm = TRUE)
                    if (m_val > 0) {
                        mat[r, ] <- mat[r, ] / m_val
                    }
                }
                message("    ✓ mean applied")
            } else if (m == "PC1") {
                # Remove the first principal component
                # mat is markers x cells, PCA expects samples (cells) in rows
                # Transpose -> cells x markers
                mat_t <- t(mat)

                # Center the data (mean-center each column/marker)
                col_means <- colMeans(mat_t, na.rm = TRUE)
                mat_centered <- sweep(mat_t, 2, col_means, "-")

                # Compute first PC using SVD (more efficient than prcomp for large data)
                # We only need the first component
                svd_result <- svd(mat_centered, nu = 1, nv = 1)

                # Reconstruct the PC1 contribution: u1 * d1 * v1^T
                pc1_contribution <- svd_result$u[, 1, drop = FALSE] %*%
                    (svd_result$d[1] * t(svd_result$v[, 1, drop = FALSE]))

                # Subtract PC1 contribution and add back the means
                mat_cleaned <- mat_centered - pc1_contribution
                mat_cleaned <- sweep(mat_cleaned, 2, col_means, "+")

                # Transpose back to markers x cells
                mat <- t(mat_cleaned)
                message("    ✓ PC1 removed")
            } else {
                warning("Unknown normalization method: ", m, ". Skipping.")
            }
        }

        # Store back
        full_norm_subset[, cell_idx] <- mat
    }

    # Assign the normalized values back to the SPE
    full_norm_assay <- assay(combined_spe, "norm")
    full_norm_assay[use_markers, ] <- full_norm_subset
    assay(combined_spe, "norm", withDimnames = FALSE) <- full_norm_assay

    message("✓ Normalization complete")
    message("  Methods applied: ", paste(methods, collapse = " → "))
    message("  Result stored in 'norm' assay")

    # Store use_markers as attribute
    attr(combined_spe, "use_markers") <- use_markers

    return(combined_spe)
}


#' Generate post-normalization density plots
#'
#' @param combined_spe SpatialExperiment object with 'norm' assay
#' @param res_dir Character. Results directory
#' @param marker_stats_raw Data frame. Raw marker statistics (optional)
#'
#' @return List with stats and plots
#'
#' @export
generate_post_norm_umap <- function(combined_spe,
                                    res_dir,
                                    marker_stats_raw = NULL) {
  require(scater)
  require(dplyr)
  require(SummarizedExperiment)

  # Get stored parameters
  ct_markers <- attr(combined_spe, "ct_markers")
  BPPARAM <- attr(combined_spe, "BPPARAM")
  n_cores <- attr(combined_spe, "n_cores")
  use_markers <- attr(combined_spe, "use_markers")

  message("Running UMAP on normalized data...")

  # Run UMAP on normalized data
  combined_spe <- scater::runUMAP(
    combined_spe,
    subset_row = ct_markers,
    exprs_values = "norm",
    name = "UMAP_norm",
    BPPARAM = BPPARAM,
    n_threads = n_cores
  )

  message("✓ Post-normalization UMAP complete")

  # Create plot
  qc_dir <- file.path(res_dir, "QC_plots")

  umap_norm_plot <- scater::plotReducedDim(
    combined_spe,
    dimred = "UMAP_norm",
    colour_by = "image_name",
    point_size = 0.5,
    point_alpha = 0.6
  ) +
    ggplot2::labs(
      title = "UMAP: After Normalization",
      subtitle = "Better mixing indicates reduced batch effects"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 11, color = "gray30")
    )

  umap_norm_file <- file.path(qc_dir, "umap_post_normalization.png")
  ggplot2::ggsave(
    filename = umap_norm_file,
    plot = umap_norm_plot,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
  )

  message("✓ Post-normalization UMAP saved to: ", umap_norm_file)

  # Calculate statistics
  marker_stats_norm <- lapply(use_markers, function(marker) {
    intensity_data <- as.vector(assay(combined_spe[marker, ], "norm"))
    image_data <- combined_spe$image_name

    data.frame(
      intensity = intensity_data,
      image = image_data
    ) %>%
      group_by(image) %>%
      summarise(
        marker = marker,
        mean_intensity = mean(intensity, na.rm = TRUE),
        median_intensity = median(intensity, na.rm = TRUE),
        sd_intensity = sd(intensity, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      select(marker, image, everything())
  }) %>% bind_rows()

  # Calculate normalization effectiveness
  if (!is.null(marker_stats_raw)) {
    norm_effect <- marker_stats_raw %>%
      inner_join(
        marker_stats_norm,
        by = c("marker", "image"),
        suffix = c("_raw", "_norm")
      ) %>%
      group_by(marker) %>%
      summarise(
        mean_diff_raw = diff(mean_intensity_raw),
        mean_diff_norm = diff(mean_intensity_norm),
        sd_diff_raw = diff(sd_intensity_raw),
        sd_diff_norm = diff(sd_intensity_norm),
        reduction_in_mean_diff = abs(mean_diff_raw) - abs(mean_diff_norm),
        .groups = "drop"
      ) %>%
      arrange(desc(abs(reduction_in_mean_diff)))

    # Save statistics
    norm_stats_file <- file.path(qc_dir, "marker_statistics_normalized.csv")
    write.csv(marker_stats_norm, norm_stats_file, row.names = FALSE)

    norm_effect_file <- file.path(qc_dir, "normalization_effectiveness.csv")
    write.csv(norm_effect, norm_effect_file, row.names = FALSE)

    message("✓ Normalization statistics exported")

    improved <- sum(norm_effect$reduction_in_mean_diff > 0, na.rm = TRUE)
    total_markers <- nrow(norm_effect)
    message(
      "✓ Normalization improved ", improved, "/", total_markers, " markers (",
      round(100 * improved / total_markers, 1), "%)"
    )
  } else {
    norm_effect <- NULL
    improved <- NA
    total_markers <- NA
  }

  return(list(
    spe = combined_spe,
    marker_stats_norm = marker_stats_norm,
    norm_effect = norm_effect,
    improved = improved,
    total_markers = total_markers
  ))
}


#' Analyze FoxP3+ cells
#'
#' @param combined_spe SpatialExperiment object
#' @param res_dir Character. Results directory
#' @param foxp3_channel Character. FoxP3 channel name
#' @param threshold_percentile Numeric. Percentile for FoxP3+ threshold (default: 0.98 for top 2%)
#'
#' @return List with updated SPE and FoxP3 summary
#'
#' @export
analyze_foxp3_cells <- function(combined_spe,
                                res_dir,
                                foxp3_channel = "31_CorrectedFOXP3-GFP",
                                threshold_percentile = 0.98) {
  require(ggplot2)
  require(gridExtra)
  require(dplyr)
  require(SpatialExperiment)
  require(SummarizedExperiment)

  message("\n🔬 FOXP3+ CELL ANALYSIS")
  message(paste(rep("=", 70), collapse = ""))

  # Check channel exists
  if (!foxp3_channel %in% rownames(combined_spe)) {
    stop("FoxP3 channel '", foxp3_channel, "' not found in data")
  }

  foxp3_expression <- assay(combined_spe, "norm")[foxp3_channel, ]
  foxp3_threshold <- quantile(foxp3_expression, threshold_percentile)

  message(
    "FoxP3+ threshold (", threshold_percentile * 100, "th percentile): ",
    round(foxp3_threshold, 4)
  )

  combined_spe$FoxP3_positive <- foxp3_expression > foxp3_threshold
  combined_spe$FoxP3_expression <- foxp3_expression

  # Summary by image
  foxp3_summary <- colData(combined_spe) %>%
    as.data.frame() %>%
    group_by(image_name) %>%
    summarise(
      total_cells = n(),
      foxp3_positive_count = sum(FoxP3_positive),
      foxp3_percentage = 100 * foxp3_positive_count / total_cells,
      mean_foxp3_expression = mean(FoxP3_expression),
      median_foxp3_expression = median(FoxP3_expression),
      .groups = "drop"
    )

  message("\n📊 FoxP3+ TREG SUMMARY:")
  print(foxp3_summary)

  if (nrow(foxp3_summary) == 2) {
    pct_diff <- diff(foxp3_summary$foxp3_percentage)
    message("\n🎯 IMAGE COMPARISON:")
    message("  FoxP3+ % difference: ", round(abs(pct_diff), 2), "%")
    message("  Higher in: ", foxp3_summary$image_name[which.max(foxp3_summary$foxp3_percentage)])
  }

  # Save summary
  qc_dir <- file.path(res_dir, "QC_plots")
  foxp3_file <- file.path(qc_dir, "FOXP3_treg_summary.csv")
  write.csv(foxp3_summary, foxp3_file, row.names = FALSE)

  # Density plot
  foxp3_density_plot <- ggplot(
    colData(combined_spe) %>% as.data.frame(),
    aes(x = FoxP3_expression, color = image_name)
  ) +
    geom_density(linewidth = 1.5) +
    geom_vline(xintercept = foxp3_threshold, linetype = "dashed", color = "red") +
    labs(
      title = "FoxP3 Expression by Image",
      subtitle = paste0("Red line = FoxP3+ threshold: ", round(foxp3_threshold, 3)),
      x = "FoxP3 Expression (Normalized)",
      y = "Density",
      color = "Image"
    ) +
    theme_classic() +
    theme(legend.position = "top", plot.title = element_text(face = "bold"))

  ggsave(file.path(qc_dir, "FOXP3_expression_distribution.png"),
    foxp3_density_plot,
    width = 10, height = 6, dpi = 300, bg = "white"
  )

  # Spatial distribution
  spatial_plots <- list()
  for (img_name in unique(combined_spe$image_name)) {
    img_cells <- combined_spe$image_name == img_name
    coords <- spatialCoords(combined_spe)[img_cells, , drop = FALSE]
    cell_data <- colData(combined_spe)[img_cells, ] %>% as.data.frame()
    img_data <- cbind(as.data.frame(coords), cell_data)

    p <- ggplot(img_data, aes(x = x, y = y, color = FoxP3_positive)) +
      geom_point(size = 0.4, alpha = 0.7) +
      scale_color_manual(
        values = c("FALSE" = "gray85", "TRUE" = "#E31A1C"),
        labels = c("FALSE" = "FoxP3-", "TRUE" = "FoxP3+ (Treg)")
      ) +
      labs(
        title = img_name,
        subtitle = paste0(
          "FoxP3+ Tregs: ", sum(img_data$FoxP3_positive),
          " (", round(100 * mean(img_data$FoxP3_positive), 1), "%)"
        ),
        color = NULL
      ) +
      theme_classic() +
      theme(aspect.ratio = 1, legend.position = "top")

    spatial_plots[[img_name]] <- p
  }

  spatial_combined <- arrangeGrob(
    grobs = spatial_plots, ncol = 2,
    top = "FoxP3+ Treg Spatial Distribution"
  )

  ggsave(
    file.path(qc_dir, "FOXP3_spatial_distribution.png"),
    spatial_combined,
    width = 14, height = 7, dpi = 300, bg = "white"
  )

  message("✓ FoxP3 visualizations saved")

  return(list(
    spe = combined_spe,
    foxp3_summary = foxp3_summary
  ))
}
