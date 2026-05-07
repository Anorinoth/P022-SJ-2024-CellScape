# Process Tiled Segmentation Results
# Author: Dr. Michael Zaiken
#
# This function loads cached tiled segmentation results, applies coordinate
# transformations, displays diagnostics, and creates visualization plots.
#
# @param result_cache_file Path to cached result object (.qs file)
# @param transform_mode Transformation mode: "swap_reflect" or "rotate_reflect"
# @param image_name Name of the image (for output files)
# @param prog_dir PROG directory path (for output)
#
# @return List containing:
#   - cells: Transformed SpatialExperiment object
#   - n_cells: Total cell count
#   - tile_dir: Tile directory path
#   - diagnostics: Segmentation diagnostics data frame

process_tiled_results <- function(result_cache_file,
                                  transform_mode = "swap_reflect",
                                  image_name = NULL,
                                  prog_dir = NULL) {
  # Validate inputs--------------------------------------------------------------
  if (!file.exists(result_cache_file)) {
    stop("Result cache not found: ", result_cache_file)
  }

  if (!transform_mode %in% c("swap_reflect", "rotate_reflect", "none", "unscramble_rot")) {
    stop(
      "Unknown transform_mode: ", transform_mode,
      ". Must be 'swap_reflect', 'rotate_reflect', 'unscramble_rot' or 'none'"
    )
  }

  # Load cached results----------------------------------------------------------
  message("\n=== Loading Cached Results Object ===")
  message("  Loading from: ", result_cache_file)
  result <- qs::qread(result_cache_file)
  message("  ✓ Results object loaded from cache")

  # Extract results--------------------------------------------------------------
  cells <- result$spe
  n_cells <- result$n_cells
  tile_dir <- result$tile_dir
  diagnostics <- result$diagnostics

  # Load tile metadata early (needed for coordinate repair and transforms)
  metadata_file <- file.path(tile_dir, "tile_grid_metadata.rds")
  tile_metadata <- NULL
  tile_size <- 4096 # Default fallback
  n_tiles_x <- 0
  n_tiles_y <- 0
  crop_x_orig <- 0
  crop_y_orig <- 0

  if (file.exists(metadata_file)) {
    tile_metadata <- readRDS(metadata_file)
    tile_size <- tile_metadata$tile_size %||% 4096 # Use fallback if missing
    n_tiles_x <- tile_metadata$n_tiles_x
    n_tiles_y <- tile_metadata$n_tiles_y
    crop_x_orig <- tile_metadata$crop_x
    crop_y_orig <- tile_metadata$crop_y
    message("  ✓ Loaded tile metadata early (tile_size=", tile_size, ")")
  } else {
    message("  ⚠️  Tile metadata not found at: ", metadata_file)
  }

  # AUTO-REPAIR: Check for scrambling bug (Row->X, Col->Y) and fix it
  # The bug in segment_and_quantify_tiled.R caused tiles to be transposed.
  # We can detect this by checking if the tiles match their ID.
  # Or simply force-repair if we are confident.
  # Since the user confirmed "scrambled", we will apply the repair based on imageID.

  if ("imageID" %in% names(SingleCellExperiment::colData(cells))) {
    # Special Handling for 'unscramble_rot' mode (User hypothesis: Global CW, Local CCW)
    if (transform_mode == "unscramble_rot") {
      message("  Applying 'unscramble_rot' logic (Global 90° CW + Local 90° CCW)...")

      image_ids <- SingleCellExperiment::colData(cells)$imageID
      matches <- regmatches(image_ids, regexec("tile_(\\d+)_(\\d+)", image_ids))

      if (length(matches) > 0 && length(matches[[1]]) == 3) {
        matches_mat <- do.call(rbind, matches)
        rows <- as.integer(matches_mat[, 2]) # Original Y index
        cols <- as.integer(matches_mat[, 3]) # Original X index

        c_x <- floor(crop_x_orig)
        c_y <- floor(crop_y_orig)

        current_coords <- SingleCellExperiment::int_colData(cells)$spatialCoords

        # 1. Recover original local coordinates (assuming 1-based indexing passed through)
        # Use original dimensions for recovery
        local_x <- (current_coords[, 1] - c_x - 1) %% tile_size + 1
        local_y <- (current_coords[, 2] - c_y - 1) %% tile_size + 1

        # 2. Apply Local 90° CCW Rotation
        # (x, y) -> (y, size - x)
        # Center of rotation is center of tile. 0-based: (y, size-1-x). 1-based: (y, size-x+1)
        # Let's assume arithmetic:
        # x_new = y
        # y_new = size - x

        # 0-based math is safer:
        lx0 <- local_x - 1
        ly0 <- local_y - 1

        lx0_new <- ly0
        ly0_new <- tile_size - 1 - lx0

        local_x_new <- lx0_new + 1
        local_y_new <- ly0_new + 1

        # 3. Apply Global 90° CW Rotation to Grid
        # Old Grid: Rows (0..Ny-1), Cols (0..Nx-1)
        # New Grid: Row' = Col, Col' = Ny - 1 - Row
        # Wait, (r, c) -> (c, Ny - 1 - r) for CW?
        # Let's re-verify:
        # (0,0) [Top-Left] -> (0, Ny-1) [Top-Right].
        # Is Col index the X-axis? Yes. So (0, Ny-1) means X=0, Y=Ny-1? No.
        # Grid coordinates (row, col) map to (Y, X).
        # Standard Image Rotation 90 CW:
        # Point (x, y) -> (H - y, x).
        # Grid Cell (c, r) -> (Ny - 1 - r, c).
        # Let's check:
        # (0,0) -> (Ny-1, 0) [Far Right X, Top Y]. Correct for 90 CW.

        new_grid_col <- n_tiles_y - 1 - rows # New X index (from old Row)
        new_grid_row <- cols # New Y index (from old Col)

        # 4. Reconstruct Global Position
        # Global X = New_Col * Size + Local_X_New
        # Global Y = New_Row * Size + Local_Y_New

        # Note: We ignore original crops here because the image dimension flipped.
        # Ideally we swap crops too. Crop_X_New = Crop_Y_Orig.

        new_global_x <- new_grid_col * tile_size + c_y + local_x_new
        new_global_y <- new_grid_row * tile_size + c_x + local_y_new

        current_coords[, 1] <- new_global_x
        current_coords[, 2] <- new_global_y
        SingleCellExperiment::int_colData(cells)$spatialCoords <- current_coords
        coords <- SpatialExperiment::spatialCoords(cells)
        message("     ✓ Unscramble rotation applied.")

        # Update global Max limits for plotting (Flip dimensions)
        # We hack the 'tile_metadata' variables locally for the plot loop
        # But those are local variables. We need to update them or create special ones.
        # The plot logic reads 'n_tiles_x' from metadata object directly if it exists.
        # We should probably update the metadata object in memory if we can.

        # Update n_tiles_x/y for subsequent use
        temp <- n_tiles_x
        n_tiles_x <<- n_tiles_y
        n_tiles_y <<- temp
        temp_c <- crop_x_orig
        crop_x_orig <<- crop_y_orig
        crop_y_orig <<- temp_c
      }
    } else {
      # Standard Auto-Repair (Same as before)
      message("  Checking for coordinate scrambling (Tile Transposition Bug)...")

      # Parse tile indices from imageID (format: tile_RR_CC)
      image_ids <- SingleCellExperiment::colData(cells)$imageID
      matches <- regmatches(image_ids, regexec("tile_(\\d+)_(\\d+)", image_ids))

      if (length(matches) > 0 && length(matches[[1]]) == 3) {
        matches_mat <- do.call(rbind, matches)
        rows <- as.integer(matches_mat[, 2])
        cols <- as.integer(matches_mat[, 3])

        # Use floor for crop to match segmentation script
        c_x <- floor(crop_x_orig)
        c_y <- floor(crop_y_orig)

        current_coords <- SingleCellExperiment::int_colData(cells)$spatialCoords

        # Recover local coords (assuming 1-based indexing passed through)
        # We subtract 1 before modulo to handle the 'Size' case mapping to 0
        local_x <- (current_coords[, 1] - c_x - 1) %% tile_size + 1
        local_y <- (current_coords[, 2] - c_y - 1) %% tile_size + 1

        # Reconstruct Global:
        # Correct Logic [v12]: Global_X = Col * Size + CropX + Local_X
        #                Global_Y = Row * Size + CropY + Local_Y

        new_x <- cols * tile_size + c_x + local_x
        new_y <- rows * tile_size + c_y + local_y

        # Detect misalignment: If existing coordinates are far from expected boundaries
        # Or if the repair significantly changes the mean position
        dist_x <- abs(new_x - current_coords[, 1])
        dist_y <- abs(new_y - current_coords[, 2])

        if (mean(dist_x) > tile_size / 2 || mean(dist_y) > tile_size / 2) {
          message("  ⚠️  DETECTED COORDINATE MISALIGNMENT. Reconstructing coordinates from tile IDs...")
          message(sprintf("     Reconstructing %d cells...", length(new_x)))

          current_coords[, 1] <- new_x
          current_coords[, 2] <- new_y
          SingleCellExperiment::int_colData(cells)$spatialCoords <- current_coords

          # Reload accessors
          coords <- SpatialExperiment::spatialCoords(cells)
          message("     ✓ Coordinates reconstructed.")
        }
      }
    }
  }

  # Apply coordinate transformations---------------------------------------------
  message("\n=== Correcting Spatial Coordinates ===")
  message("  Transformation mode: ", transform_mode)
  coords <- SpatialExperiment::spatialCoords(cells)


  message("  Original coordinate range:")
  message("    X: ", round(min(coords[, 1])), " - ", round(max(coords[, 1])))
  message("    Y: ", round(min(coords[, 2])), " - ", round(max(coords[, 2])))

  # Apply "none" transformation
  if (transform_mode == "none") {
    message("  Result: No coordinate transformation applied")
    message("    ✓ Transformation skipped")
  }

  # Apply swap_reflect transformation (CAR_TREG_2)
  if (transform_mode == "swap_reflect") {
    # Step 1: Swap X ↔ Y
    message("  Step 1: Swapping X ↔ Y...")
    temp_coords <- coords
    # Directly modify internal storage to bypass potential S4 replacement issues
    current_coords <- SingleCellExperiment::int_colData(cells)$spatialCoords
    current_coords[, 1] <- temp_coords[, 2]
    current_coords[, 2] <- temp_coords[, 1]
    SingleCellExperiment::int_colData(cells)$spatialCoords <- current_coords

    coords_swapped <- SpatialExperiment::spatialCoords(cells)
    message("    Swapped coordinate range:")
    message("      X: ", round(min(coords_swapped[, 1])), " - ", round(max(coords_swapped[, 1])))
    message("      Y: ", round(min(coords_swapped[, 2])), " - ", round(max(coords_swapped[, 2])))

    # Step 2: Reflect over Y-axis
    message("  Step 2: Reflecting over Y-axis...")
    max_x <- max(coords_swapped[, 1])
    max_x_swapped_for_tiles <- max_x # Capture for tile logic
    current_coords <- SingleCellExperiment::int_colData(cells)$spatialCoords
    current_coords[, 1] <- max_x - current_coords[, 1]
    SingleCellExperiment::int_colData(cells)$spatialCoords <- current_coords

    coords_final <- SpatialExperiment::spatialCoords(cells)
    message("    Final coordinate range:")
    message("      X: ", round(min(coords_final[, 1])), " - ", round(max(coords_final[, 1])))
    message("      Y: ", round(min(coords_final[, 2])), " - ", round(max(coords_final[, 2])))
    message("    ✓ CAR transformation complete")
  }

  # Apply rotate_reflect transformation (EGFR_TREG_2)
  # IMPORTANT: Save transformation max values for tile overlay plotting
  max_x_swapped_for_tiles <- NULL
  max_y_car_for_tiles <- NULL
  max_x_final_for_tiles <- NULL

  if (transform_mode == "rotate_reflect") {
    # Step 1: Apply CAR corrections (swap + reflect)
    message("  Step 1: Applying CAR corrections (swap X↔Y + Y-reflect)...")
    temp_coords <- coords
    current_coords <- SingleCellExperiment::int_colData(cells)$spatialCoords
    current_coords[, 1] <- temp_coords[, 2]
    current_coords[, 2] <- temp_coords[, 1]
    SingleCellExperiment::int_colData(cells)$spatialCoords <- current_coords

    coords_swapped <- SpatialExperiment::spatialCoords(cells)
    max_x_swapped <- max(coords_swapped[, 1])
    max_x_swapped_for_tiles <- max_x_swapped # Save for tile plotting
    current_coords <- SingleCellExperiment::int_colData(cells)$spatialCoords
    current_coords[, 1] <- max_x_swapped - current_coords[, 1]
    SingleCellExperiment::int_colData(cells)$spatialCoords <- current_coords

    coords_after_car <- SpatialExperiment::spatialCoords(cells)
    message("    After CAR corrections:")
    message("      X: ", round(min(coords_after_car[, 1])), " - ", round(max(coords_after_car[, 1])))
    message("      Y: ", round(min(coords_after_car[, 2])), " - ", round(max(coords_after_car[, 2])))

    # Step 2: Rotate 90° clockwise
    message("  Step 2: Rotating 90° clockwise...")
    max_y_car <- max(coords_after_car[, 2])
    max_y_car_for_tiles <- max_y_car # Save for tile plotting
    new_x <- max_y_car - coords_after_car[, 2]
    new_y <- coords_after_car[, 1]

    current_coords <- SingleCellExperiment::int_colData(cells)$spatialCoords
    current_coords[, 1] <- new_x
    current_coords[, 2] <- new_y
    SingleCellExperiment::int_colData(cells)$spatialCoords <- current_coords

    coords_rotated <- SpatialExperiment::spatialCoords(cells)
    message("    After rotation:")
    message("      X: ", round(min(coords_rotated[, 1])), " - ", round(max(coords_rotated[, 1])))
    message("      Y: ", round(min(coords_rotated[, 2])), " - ", round(max(coords_rotated[, 2])))

    # Step 3: Reflect over Y-axis
    message("  Step 3: Final Y-axis reflection...")
    max_x_final <- max(coords_rotated[, 1])
    max_x_final_for_tiles <- max_x_final # Save for tile plotting
    current_coords <- SingleCellExperiment::int_colData(cells)$spatialCoords
    current_coords[, 1] <- max_x_final - current_coords[, 1]
    SingleCellExperiment::int_colData(cells)$spatialCoords <- current_coords

    coords_final <- SpatialExperiment::spatialCoords(cells)
    message("    Final coordinate range:")
    message("      X: ", round(min(coords_final[, 1])), " - ", round(max(coords_final[, 1])))
    message("      Y: ", round(min(coords_final[, 2])), " - ", round(max(coords_final[, 2])))
    message("    ✓ EGFR transformation complete")
  }

  # Display timing summary-------------------------------------------------------
  if (!is.null(result$timing)) {
    message("\n=== Timing Summary ===")
    timing_df <- data.frame(
      Stage = names(result$timing),
      Seconds = round(result$timing, 2)
    )
    print(timing_df, row.names = FALSE)
  }

  # Display segmentation diagnostics---------------------------------------------
  if (!is.null(diagnostics) && nrow(diagnostics) > 0) {
    message("\n=== Segmentation Diagnostics ===")
    message("  Diagnostics file: ", file.path(tile_dir, "segmentation_diagnostics.csv"))
    message("  Total tiles: ", nrow(diagnostics))
    message("  Total cells: ", format(sum(diagnostics$n_cells_final), big.mark = ","))

    # Tile type distribution
    message("\n  Tile Type Distribution:")
    tile_type_summary <- table(diagnostics$tile_type)
    for (tt in names(tile_type_summary)) {
      avg_cells <- mean(diagnostics$n_cells_final[diagnostics$tile_type == tt])
      message(sprintf(
        "    %-12s: %2d tiles (avg %5.0f cells/tile)",
        tt, tile_type_summary[tt], avg_cells
      ))
    }

    # Show top 10 and bottom 10 tiles by cell count
    diagnostics_sorted <- diagnostics[order(-diagnostics$n_cells_final), ]

    message("\n  Top 10 tiles (most cells):")
    top10 <- head(diagnostics_sorted, 10)
    for (i in 1:nrow(top10)) {
      message(sprintf(
        "    %s: %5d cells (%s, gamma=%.2f, bright=%.2f)",
        top10$tile_id[i], top10$n_cells_final[i], top10$tile_type[i],
        top10$gamma_used[i], top10$brightness_used[i]
      ))
    }

    message("\n  Bottom 10 tiles (fewest cells):")
    bottom10 <- tail(diagnostics_sorted, 10)
    for (i in nrow(bottom10):1) {
      message(sprintf(
        "    %s: %5d cells (%s, gamma=%.2f, bright=%.2f)",
        bottom10$tile_id[i], bottom10$n_cells_final[i], bottom10$tile_type[i],
        bottom10$gamma_used[i], bottom10$brightness_used[i]
      ))
    }

    message("\n  ✓ Full diagnostics table saved to CSV file")
  }

  # Display results summary------------------------------------------------------
  message("\n=== Tiled Processing Results ===")
  message("  Total cells detected: ", format(n_cells, big.mark = ","))
  message("  Channels quantified: ", nrow(cells))
  message("  Tile checkpoint directory: ", tile_dir)
  if (!is.null(prog_dir) && !is.null(image_name)) {
    message("  Final SPE file: ", file.path(prog_dir, paste0(image_name, "_cells.qs")))
  }

  # Verify SpatialExperiment structure-------------------------------------------
  message("\n=== Verifying SpatialExperiment Structure ===")
  message("  Spatial coordinate names: ", paste(SpatialExperiment::spatialCoordsNames(cells), collapse = ", "))
  message(
    "  Coordinate range (x): ", round(min(SpatialExperiment::spatialCoords(cells)[, 1])), " - ",
    round(max(SpatialExperiment::spatialCoords(cells)[, 1]))
  )
  message(
    "  Coordinate range (y): ", round(min(SpatialExperiment::spatialCoords(cells)[, 2])), " - ",
    round(max(SpatialExperiment::spatialCoords(cells)[, 2]))
  )
  message("  Assays: ", paste(names(SummarizedExperiment::assays(cells)), collapse = ", "))
  message("  Features (channels): ", nrow(cells))
  message("  Cells: ", ncol(cells))
  message("  colData columns: ", paste(colnames(SummarizedExperiment::colData(cells)), collapse = ", "))

  # Create cell positions plot with tile overlays-------------------------------
  if (!is.null(prog_dir) && !is.null(image_name)) {
    message("\n  Plotting cell positions...")
    coords <- SpatialExperiment::spatialCoords(cells)

    viz_dir <- file.path(prog_dir, "tiles", image_name, "visualizations")
    dir.create(viz_dir, showWarnings = FALSE, recursive = TRUE)

    # Metadata already loaded at start of function
    if (is.null(tile_metadata)) {
      message("  ⚠️  Tile metadata not found, plotting without tile labels")
    }

    # Start the plot with explicit axis limits to ensure all tiles are visible
    plot_file <- file.path(viz_dir, paste0("cell_positions_", transform_mode, ".png"))
    if (file.exists(plot_file)) {
      message("  Replacing existing plot: ", plot_file)
      unlink(plot_file)
    } else {
      message("  Creating new plot: ", plot_file)
    }

    png(plot_file, width = 10, height = 10, units = "in", res = 300)

    # Set x and y limits to include a small margin beyond the data
    # [v12] Invert Y-axis: ylim = c(max_y, 0) maps image 0 to the TOP
    x_range <- range(coords[, 1])
    y_range <- range(coords[, 2])
    x_margin <- (x_range[2] - x_range[1]) * 0.05
    y_margin <- (y_range[2] - y_range[1]) * 0.05

    plot(coords[, 1], coords[, 2],
      pch = ".", col = "black",
      main = sprintf("Cell Positions: %s (n=%s cells)", image_name, format(n_cells, big.mark = ",")),
      xlab = "X coordinate (pixels)",
      ylab = "Y coordinate (pixels)",
      xlim = x_range + c(-x_margin, x_margin),
      ylim = c(y_range[2] + y_margin, y_range[1] - y_margin), # Inverted: Max to Min
      asp = 1 # Keep 1:1 aspect ratio
    )

    # Combined robust tile overlay logic
    if (!is.null(tile_metadata)) {
      message("  Adding tile grid overlays...")

      # Iterate through ORIGINAL grid indices
      # Standard Logic: Col -> X, Row -> Y
      for (row_idx in 0:(n_tiles_y - 1)) {
        for (col_idx in 0:(n_tiles_x - 1)) {
          # 1. Define original tile corners
          # Standard mapping:
          # X comes from col_idx (Width)
          # Y comes from row_idx (Height)
          x1 <- crop_x_orig + col_idx * tile_size
          x2 <- crop_x_orig + (col_idx + 1) * tile_size
          y1 <- crop_y_orig + row_idx * tile_size
          y2 <- crop_y_orig + (row_idx + 1) * tile_size

          corners <- data.frame(
            x = c(x1, x2, x2, x1),
            y = c(y1, y1, y2, y2)
          )

          # 2. Apply transformations to corners based on mode
          trans_corners <- corners

          if (transform_mode == "swap_reflect") {
            # Step 1: Swap
            tx <- trans_corners$y
            ty <- trans_corners$x
            trans_corners$x <- tx
            trans_corners$y <- ty

            # Step 2: Reflect Y (using max_x captured from cell transform)
            if (!is.null(max_x_swapped_for_tiles)) {
              trans_corners$x <- max_x_swapped_for_tiles - trans_corners$x
            }
          } else if (transform_mode == "rotate_reflect") {
            # Step 1: CAR (Swap + Reflect Y)
            tx <- trans_corners$y
            ty <- trans_corners$x
            trans_corners$x <- tx
            trans_corners$y <- ty

            if (!is.null(max_x_swapped_for_tiles)) {
              trans_corners$x <- max_x_swapped_for_tiles - trans_corners$x
            }

            # Step 2: Rotate 90 CW (X' = max_y - y, Y' = x)
            # Use max_y_car_for_tiles
            if (!is.null(max_y_car_for_tiles)) {
              tx <- max_y_car_for_tiles - trans_corners$y
              ty <- trans_corners$x
              trans_corners$x <- tx
              trans_corners$y <- ty
            }

            # Step 3: Reflect Y
            # Use max_x_final_for_tiles
            if (!is.null(max_x_final_for_tiles)) {
              trans_corners$x <- max_x_final_for_tiles - trans_corners$x
            }
          }

          # 3. Draw polygon connecting transformed corners
          polygon(trans_corners$x, trans_corners$y, border = "red", lwd = 0.5)

          # 4. Label with ORIGINAL indices at the visual center of the polygon
          center_x <- mean(trans_corners$x)
          center_y <- mean(trans_corners$y)

          # Default label
          tile_label <- sprintf("%02d_%02d", row_idx, col_idx)

          # If unscramble_rot, calculate original indices for labeling
          if (transform_mode == "unscramble_rot") {
            # Inverse of: new_grid_col = n_tiles_y_orig - 1 - old_row
            #             new_grid_row = old_col
            # Note: n_tiles_x/y were SWAPPED in the repair block.
            # So n_tiles_x (current) = n_tiles_y_orig.

            # c_new (col_idx) = Nx_curr - 1 - r_old  => r_old = Nx_curr - 1 - col_idx
            # r_new (row_idx) = c_old                => c_old = row_idx

            # orig_row = n_tiles_x - 1 - col_idx
            # orig_col = row_idx

            orig_row <- n_tiles_x - 1 - col_idx
            orig_col <- row_idx
            tile_label <- sprintf("%02d_%02d", orig_row, orig_col)
          } else if (length(grep("transpose", transform_mode)) > 0) {
            # If we added other modes, handle here
            # For standard modes, we used row-major or col-major fix.
            # The standard fix used row_idx->X and col_idx->Y logic,
            # so the label row_idx, col_idx (X_Y) matched the visual grid.
          }

          text(center_x, center_y, tile_label,
            col = "red", cex = 0.6, font = 2 # Reduced size
          )
        }
      }
    }

    dev.off()

    if (!is.null(tile_metadata)) {
      message("  ✓ Cell position plot with tile labels saved: ", file.path(viz_dir, "cell_positions.png"))
    } else {
      message("  ✓ Cell position plot saved: ", file.path(viz_dir, "cell_positions.png"))
    }
  }

  message("✓ SpatialExperiment processing complete")

  # Return results---------------------------------------------------------------
  return(list(
    cells = cells,
    n_cells = n_cells,
    tile_dir = tile_dir,
    diagnostics = diagnostics
  ))
}
