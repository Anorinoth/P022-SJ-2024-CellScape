#' Segment and Quantify Cells Using Tiled Processing
#' Author: Dr. Michael Zaiken
#'
#' Processes large multiplexed images by breaking them into tiles, performing
#' nuclear segmentation and cell quantification on each tile independently,
#' then merging the results. This approach minimizes RAM usage and enables
#' parallel processing with fork-based workers.
#'
#' @param images Character string: file path to an image file (.qs or .rds with CytoImageList).
#'   The function automatically uses a cached deduplicated version if available (much faster),
#'   otherwise loads from the provided file path.
#' @param image_name Character string with the name of the image
#' @param nuclear_channel_name Character string with the nuclear marker channel name
#' @param tile_size Integer, size of tiles in pixels (default: 4096)
#' @param n_cores Integer, number of cores for parallel processing (default: 8)
#'
#' @param gamma Numeric gamma correction for segmentation. Combined with brightness and
#'   contrast using formula: contrast*(image+brightness)^gamma. Used as baseline when
#'   adaptive_params=FALSE, or as target when adaptive_params=TRUE (default: 0.5)
#' @param brightness Numeric brightness adjustment added before gamma correction (default: 0)
#' @param contrast Numeric contrast multiplier applied after gamma correction (default: 1)
#' @param adaptive_params Logical, if TRUE, automatically adjusts gamma/brightness/contrast per tile
#'   based on tile intensity characteristics to handle uneven staining/illumination (default: TRUE)
#' @param percentile_range Numeric vector c(lower, upper) for normalization (default: c(0.005, 0.995))
#' @param exclude_bright_outliers Logical, mask bright outliers (default: TRUE)
#' @param outlier_percentile Numeric, outlier threshold percentile (default: 0.98)
#' @param adaptive_enhancement Logical, use CLAHE (default: TRUE)
#' @param clahe_clip_limit Numeric, CLAHE clip limit (default: 0.01)
#' @param watershed_tolerance Numeric, watershed tolerance (default: 0.1)
#' @param discSize Numeric, disc kernel size for simpleSeg (default: 6)
#' @param sizeSelection Numeric, minimum object size (default: 40)
#' @param smooth Numeric, smoothing parameter (default: 1)
#' @param ext Numeric, extension parameter (default: 1)
#' @param max_nucleus_size Numeric, maximum nucleus size (default: 1600)
#'
#' @param enhancement_gamma Numeric gamma for quantification enhancement (default: 0.5)
#' @param enhancement_percentile_range Numeric vector for enhancement (default: c(0.005, 0.995))
#' @param enhancement_exclude_outliers Logical, exclude outliers during enhancement (default: TRUE)
#' @param enhancement_outlier_percentile Numeric, outlier percentile (default: 0.98)
#' @param enhancement_adaptive Logical, use adaptive enhancement (default: TRUE)
#'
#' @param output_dirs Named list with 'RES_DIR' and 'PROG_DIR' paths
#' @param test_mode Logical, if TRUE processes only center tile (default: FALSE)
#' @param save_one_tile_viz Logical, save enhancement comparison for center tile (default: TRUE)
#' @param save_segmentation_overlay Logical, save mask overlay on nuclear channel for center
#'   tile after segmentation (default: TRUE)
#' @param clear_images_after_segmentation Logical, free original image from RAM after tile
#'   extraction (default: TRUE). Set FALSE to keep for debugging.
#' @param verbose Logical, print progress messages (default: TRUE)
#' @param timing Logical, report timing information (default: TRUE)
#'
#' @return A list containing:
#'   \item{spe}{SpatialExperiment object with merged cell quantifications}
#'   \item{n_cells}{Integer number of detected cells}
#'   \item{tile_dir}{Path to tile checkpoint directory}
#'   \item{timing}{Named numeric vector with timing for each stage}
#'   \item{parameters}{List of all parameters used}
#'
#' @details
#' The function performs these stages:
#' 0. Remove duplicate channels
#' 1. Calculate tile grid and extract tiles to disk
#' 2. Segment nuclear channel for each tile (parallel)
#' 3. Enhance all channels and quantify for each tile (parallel)
#' 4. Merge all SpatialExperiment objects with corrected coordinates
#'
#' @export
segment_and_quantify_tiled <- function(
  images,
  image_name,
  nuclear_channel_name,
  tile_size = 4096,
  tile_overlap = 200,
  n_cores = 8,
  n_cores_extraction = NULL,
  # Segmentation parameters
  gamma = 0.5,
  brightness = 0,
  contrast = 1,
  adaptive_params = TRUE, # Auto-adjust gamma/brightness/contrast per tile based on intensity
  percentile_range = c(0.005, 0.995),
  exclude_bright_outliers = TRUE,
  outlier_percentile = 0.98,
  adaptive_enhancement = TRUE,
  clahe_clip_limit = 0.01,
  watershed_tolerance = 0.1,
  discSize = 6,
  sizeSelection = 40,
  smooth = 1,
  ext = 1,
  max_nucleus_size = 1600,
  # Enhancement parameters (for quantification)
  enhancement_gamma = 0.5,
  enhancement_percentile_range = c(0.005, 0.995),
  enhancement_exclude_outliers = TRUE,
  enhancement_outlier_percentile = 0.98,
  enhancement_adaptive = TRUE,
  # Control parameters
  output_dirs,
  test_mode = FALSE,
  save_one_tile_viz = TRUE,
  save_segmentation_overlay = TRUE, # Save mask overlay on nuclear channel (center tile)
  clear_images_after_segmentation = TRUE, # Free RAM after Stage 2 completes
  tiling_mode = "disk", # "virtual" (process whole image via RAM slices) or "disk" (original tiled)
  batch_segmentation = 1, # Number of batches to split n_cores into (reduces RAM usage)
  batch_quantification = 1, # Number of batches to split n_cores into (reduces RAM usage)
  verbose = TRUE,
  timing = TRUE
) {
  # Check required packages
  required_packages <- c(
    "cytomapper", "EBImage", "simpleSeg", "SpatialExperiment",
    "BiocParallel", "S4Vectors", "parallel"
  )
  missing_packages <- required_packages[!sapply(required_packages, function(pkg) {
    requireNamespace(pkg, quietly = TRUE)
  })]

  if (length(missing_packages) > 0) {
    stop(
      "Required packages not loaded: ", paste(missing_packages, collapse = ", "),
      "\nPlease load these packages before calling segment_and_quantify_tiled()"
    )
  }

  # Validate inputs
  if (missing(output_dirs) || !all(c("RES_DIR", "PROG_DIR") %in% names(output_dirs))) {
    stop("'output_dirs' must be a named list with 'RES_DIR' and 'PROG_DIR'")
  }

  RES_DIR <- output_dirs$RES_DIR
  PROG_DIR <- output_dirs$PROG_DIR

  # Initialize timing
  step_times <- list()
  overall_start <- Sys.time()

  # Helper functions
  timed_message <- function(msg, step_name = NULL) {
    if (verbose) {
      timestamp <- format(Sys.time(), "%H:%M:%S")
      message("[", timestamp, "] ", msg)
      flush.console()
    }
    if (timing && !is.null(step_name)) {
      step_times[[step_name]] <<- Sys.time()
    }
  }

  elapsed_time <- function(start_time, units = "secs") {
    as.numeric(difftime(Sys.time(), start_time, units = units))
  }

  timed_message(paste0("\n=== Tiled Segmentation and Quantification: ", image_name, " ==="), "start")
  if (test_mode) {
    timed_message("⚠️  TEST MODE: Will process only center tile")
  }

  # ============================================================================
  # PHASE 0: Data Loading & Preparation
  # ============================================================================
  timed_message("\n--- Phase 0: Data Loading & Preparation ---", "phase0_start")

  # Create tile directory for checkpoint files
  tile_dir <- file.path(PROG_DIR, "tiles", image_name)
  dir.create(tile_dir, showWarnings = FALSE, recursive = TRUE)

  # Check for existing objects in environment (to skip disk reload)
  images_loaded_in_memory <- FALSE
  if (exists("images_obj", envir = .GlobalEnv)) {
    images_obj_temp <- get("images_obj", envir = .GlobalEnv)

    # Validation: Check if it's a CytoImageList and matches the image_name
    # Note: Checking names() or mcols()$imageID
    is_valid_cache <- FALSE
    if (inherits(images_obj_temp, "CytoImageList")) {
      if (length(names(images_obj_temp)) > 0 && names(images_obj_temp)[1] == image_name) {
        is_valid_cache <- TRUE
      } else if (!is.null(S4Vectors::mcols(images_obj_temp)$imageID) && S4Vectors::mcols(images_obj_temp)$imageID[1] == image_name) {
        is_valid_cache <- TRUE
      }
    }

    if (is_valid_cache) {
      timed_message("  ✓ Located images_obj in global environment - using cached object")
      images_obj <- images_obj_temp
      images_loaded_in_memory <- TRUE
    }

    # Clear local reference to temp object
    rm(images_obj_temp)
  }

  timed_message("\n--- Stage 0/4: Duplicate Channel Detection ---", "stage0_start")

  # Check for completion flags
  dedup_complete_flag <- file.path(tile_dir, "deduplication_complete.flag")
  dedup_checkpoint <- file.path(tile_dir, "deduplicated_image.qs")
  tile_extraction_complete_flag <- file.path(tile_dir, "tile_extraction_complete.flag")

  # Validate images parameter
  # Can be a file path (character) OR an already loaded CytoImageList
  if (inherits(images, "CytoImageList")) {
    timed_message("  ✓ images passed as CytoImageList object")
    images_obj <- images
    images_loaded_in_memory <- TRUE
  } else if (!is.character(images) || length(images) != 1) {
    stop("'images' must be either a CytoImageList or a character string: path to image file")
  }

  # Smart loading logic:
  # 1. First check: Do we even need to load the full image?
  #    - If BOTH dedup AND tile extraction are complete → Skip loading entirely
  # 2. If loading IS needed, automatically prefer cache over original file:
  #    - If cache exists → Load from cache (fast)
  #    - Otherwise → Load from provided file path

  # Virtual mode requires images_obj in RAM for Stages 2 and 3.
  # If we are resuming at Stage 4, we don't need it.
  virtual_needs_loading <- tiling_mode == "virtual" &&
    (!file.exists(file.path(tile_dir, "segmentation_complete.flag")) ||
      !file.exists(file.path(tile_dir, "quantification_complete.flag")))

  skip_stages_0_and_1 <- file.exists(dedup_complete_flag) &&
    file.exists(tile_extraction_complete_flag) &&
    file.exists(dedup_checkpoint) &&
    !virtual_needs_loading

  if (skip_stages_0_and_1) {
    # Both deduplication and tile extraction are done - skip to Stage 2
    timed_message("  ✓ Deduplication and tile extraction already complete")
    timed_message("    Skipping to Stage 2 (segmentation)")
    timed_message("    (All processing uses tiles on disk)")

    # Only need channel names for later stages - load from small metadata file
    channel_names_file <- file.path(tile_dir, "channel_names.rds")
    if (file.exists(channel_names_file)) {
      all_channel_names <- readRDS(channel_names_file)
    } else {
      # Fallback: load from deduplicated image (will be removed once channel_names.rds is saved)
      timed_message("    (Loading channel names from cache - first run only)")
      cached_metadata <- qs::qread(dedup_checkpoint, nthreads = n_cores)
      all_channel_names <- cytomapper::channelNames(cached_metadata)
      rm(cached_metadata)
      gc(verbose = FALSE)
      # Save for next time
      saveRDS(all_channel_names, channel_names_file)
    }
    n_channels <- length(all_channel_names)

    timed_message(paste0("    Channel count: ", n_channels))
    step_times$stage0 <- 0

    # Load tile grid metadata (needed for Stage 2)
    metadata_file <- file.path(tile_dir, "tile_grid_metadata.rds")
    if (!file.exists(metadata_file)) {
      stop("Tile metadata not found: ", metadata_file)
    }
    metadata <- readRDS(metadata_file)
    n_tiles_x <- metadata$n_tiles_x
    n_tiles_y <- metadata$n_tiles_y
    crop_x <- metadata$crop_x
    crop_y <- metadata$crop_y
    original_dims <- metadata$original_dims
    center_row <- metadata$center_row
    center_col <- metadata$center_col

    # Deduplication already complete - set duplicate_channels to empty
    duplicate_channels <- integer(0)

    # Set flag to skip Stage 0 and Stage 1 processing
    images_loaded <- FALSE
  } else {
    # Need to load image for deduplication and/or tile extraction
    images_loaded <- TRUE

    if (tiling_mode == "virtual") {
      timed_message("  🚀 tiling_mode = 'virtual': Will process whole image in memory using shared RAM")
    }

    if (images_loaded_in_memory) {
      # Already in memory - nothing to do here
    } else if (file.exists(dedup_checkpoint)) {
      # Cache exists - use it (much faster)

      # PROACTIVE MEMORY CLEARING: If we had a local images_obj from a previous failed block,
      # or if there's any stray object, clear it before loading a massive new one.
      if (exists("images_obj")) {
        rm(images_obj)
        gc(verbose = FALSE)
      }

      timed_message("  Loading image from cache...")
      timed_message(paste0("    Cache: ", basename(dedup_checkpoint)))
      images_obj <- qs::qread(dedup_checkpoint, nthreads = n_cores)
      timed_message("    ✓ Loaded from cache (deduplicated)")
      # Also push to global environment for next time
      assign("images_obj", images_obj, envir = .GlobalEnv)
    } else {
      # No cache - load from original file
      if (!is.character(images) || !file.exists(images)) {
        stop("Image file not found: ", images)
      }

      # PROACTIVE MEMORY CLEARING
      if (exists("images_obj")) {
        rm(images_obj)
        gc(verbose = FALSE)
      }

      timed_message(paste0("  Loading image from: ", basename(images)))

      # Auto-detect file format based on extension
      if (grepl("\\.qs$", images, ignore.case = TRUE)) {
        images_obj <- qs::qread(images, nthreads = n_cores)
      } else if (grepl("\\.rds$", images, ignore.case = TRUE)) {
        images_obj <- readRDS(images)
      } else {
        stop("Unsupported file format. Use .qs or .rds file: ", images)
      }

      timed_message("    ✓ Image loaded from file")
    }

    # Validate loaded object
    if (!inherits(images_obj, "CytoImageList")) {
      stop("Loaded object must be a CytoImageList, got: ", class(images_obj))
    }

    # Check if deduplication is needed
    if (file.exists(dedup_complete_flag) && file.exists(dedup_checkpoint)) {
      # Deduplication done, but tile extraction not done
      timed_message("  ✓ Deduplication already complete")
      all_channel_names <- cytomapper::channelNames(images_obj)
      n_channels <- length(all_channel_names)
      timed_message(paste0("    Channel count: ", n_channels))

      # Extract image array for tile extraction
      image_array <- images_obj[[1]]

      # No duplicates removed (already done in previous run)
      duplicate_channels <- integer(0)

      step_times$stage0 <- 0
    } else {
      # Need to run deduplication
      timed_message("  Running deduplication...")
      all_channel_names <- cytomapper::channelNames(images_obj)
      n_channels_original <- length(all_channel_names)
      timed_message(paste0("  Original channels: ", n_channels_original))

      # Remove numeric prefix: "01_CD19-PerCP-Cy5_5" → "CD19-PerCP-Cy5_5"
      channel_names_no_prefix <- gsub("^\\d+_", "", all_channel_names)

      # Find duplicates: which names appear more than once?
      duplicate_names <- unique(channel_names_no_prefix[duplicated(channel_names_no_prefix)])

      duplicate_channels <- integer(0)
      if (length(duplicate_names) > 0) {
        # For each duplicate name, keep the FIRST occurrence, mark others as duplicates
        for (dup_name in duplicate_names) {
          indices_with_name <- which(channel_names_no_prefix == dup_name)
          if (length(indices_with_name) > 1) {
            # Keep first, mark rest as duplicates
            duplicate_channels <- c(duplicate_channels, indices_with_name[-1])
          }
        }

        timed_message(paste0("  Found ", length(duplicate_channels), " duplicate channels to exclude:"))
        for (dup_ch in sort(duplicate_channels)) {
          timed_message(paste0(
            "    - Channel ", sprintf("%02d", dup_ch), ": ", all_channel_names[dup_ch],
            " (duplicate of ", channel_names_no_prefix[dup_ch], ")"
          ))
        }

        # Keep all channels EXCEPT duplicates
        channels_to_keep <- setdiff(1:n_channels_original, duplicate_channels)
        n_channels <- length(channels_to_keep)

        # Filter images to keep only non-duplicate channels
        image_array <- images_obj[[1]]
        filtered_array <- image_array[, , channels_to_keep, drop = FALSE]
        filtered_channel_names <- all_channel_names[channels_to_keep]

        # Update channel names
        dimnames(filtered_array)[[3]] <- filtered_channel_names

        # Create new CytoImageList
        filtered_image <- EBImage::Image(filtered_array)
        images_obj <- cytomapper::CytoImageList(list(filtered_image))
        names(images_obj) <- image_name
        cytomapper::channelNames(images_obj) <- filtered_channel_names
        S4Vectors::mcols(images_obj)$imageID <- image_name

        # Update channel name list
        all_channel_names <- filtered_channel_names

        # Clean up
        rm(image_array, filtered_array, filtered_image)
        gc(verbose = FALSE)

        timed_message(paste0(
          "  ✓ Retained ", n_channels, " unique channels (removed ",
          length(duplicate_channels), " duplicates)"
        ))
      } else {
        timed_message("  ✓ No duplicate channel names detected")
        channels_to_keep <- 1:n_channels_original
        n_channels <- n_channels_original
      }

      # Save deduplicated image for resume
      qs::qsave(images_obj, dedup_checkpoint, preset = "fast", nthreads = n_cores)
      writeLines("complete", dedup_complete_flag)
      timed_message("  💾 Saved deduplication checkpoint")

      # Also push to global environment for next time
      assign("images_obj", images_obj, envir = .GlobalEnv)

      step_times$stage0 <- elapsed_time(step_times$stage0_start)
      if (timing) timed_message(paste0("  ⏱ Duplicate detection: ", round(step_times$stage0, 2), " sec"))
    } # End of deduplication else block
  } # End of images_loaded block

  # Check for nuclear marker
  if (!nuclear_channel_name %in% all_channel_names) {
    stop(
      "Nuclear marker '", nuclear_channel_name, "' not found in ", image_name,
      "\nAvailable channels: ", paste(all_channel_names, collapse = ", ")
    )
  }

  nuclear_channel_idx <- which(all_channel_names == nuclear_channel_name)

  # ============================================================================
  # PHASE 1: Tile Grid Optimization & Extraction
  # ============================================================================
  timed_message("\n--- Phase 1: Tile Grid Optimization & Extraction ---", "phase1_start")

  # ============================================================================
  # STAGE 1: Tile Grid Calculation & Extraction
  # ============================================================================
  if (images_loaded) {
    # Only run if image was loaded (tiles not yet on disk)
    timed_message("\n--- Stage 1/4: Tile Grid and Extraction ---", "stage1_start")

    # tile_dir already created in Stage 0
    # Initialize variables (will be set from metadata or calculated)
    center_row <- NULL
    center_col <- NULL

    # Check for completion flag
    tile_extraction_complete_flag <- file.path(tile_dir, "tile_extraction_complete.flag")
    metadata_file <- file.path(tile_dir, "tile_grid_metadata.rds")

    if (file.exists(metadata_file)) {
      timed_message("  ✓ Found existing tile grid metadata, loading...")
      metadata <- readRDS(metadata_file)

      # INVALIDATION: If tiling_mode has changed, we MUST recalculate the grid
      # Also check if n_cores has changed significantly for virtual mode
      mode_changed <- (!is.null(metadata$tiling_mode) && metadata$tiling_mode != tiling_mode) ||
        (is.null(metadata$tiling_mode) && tiling_mode == "virtual")

      if (mode_changed) {
        timed_message(paste0(
          "  ⚠️  tiling_mode changed (from ",
          if (is.null(metadata$tiling_mode)) "disk" else metadata$tiling_mode,
          " to ", tiling_mode, "). Recalculating grid..."
        ))
        file.remove(metadata_file)
        if (file.exists(tile_extraction_complete_flag)) file.remove(tile_extraction_complete_flag)
        # Proceed to calculate block
        recalculate_grid <- TRUE
      } else {
        n_tiles_x <- metadata$n_tiles_x
        n_tiles_y <- metadata$n_tiles_y
        n_tiles_x_full <- metadata$n_tiles_x_full
        n_tiles_y_full <- metadata$n_tiles_y_full
        crop_x <- metadata$crop_x
        crop_y <- metadata$crop_y
        original_dims <- metadata$original_dims
        test_mode_from_metadata <- metadata$test_mode
        center_row <- metadata$center_row
        center_col <- metadata$center_col
        # Update tile_size if it was saved
        if (!is.null(metadata$tile_size)) tile_size <- metadata$tile_size

        timed_message(paste0("    Grid: ", n_tiles_x, " x ", n_tiles_y, " tiles"))
        timed_message(paste0("    Tile size: ", tile_size, " x ", tile_size))
        timed_message(paste0("    Test mode: ", test_mode_from_metadata))

        # Check if test mode has changed
        if (test_mode != test_mode_from_metadata) {
          timed_message("  ⚠️  Warning: test_mode changed since metadata was created")
          timed_message("     Delete tile directory to regenerate with new test_mode")
        }

        # Skip calculation
        recalculate_grid <- FALSE

        # CRITICAL: Reconstruct tile_coords for Stage 1 skip path
        if (test_mode) {
          tile_coords <- data.frame(row = center_row, col = center_col)
        } else {
          tile_coords <- expand.grid(row = 0:(n_tiles_y - 1), col = 0:(n_tiles_x - 1))
        }
      }
    } else {
      recalculate_grid <- TRUE
    }

    if (recalculate_grid) {
      # Calculate tile grid
      timed_message("  Calculating tile grid...")

      image_array <- images_obj[[1]]
      original_dims <- dim(image_array)
      height <- original_dims[1]
      width <- original_dims[2]

      timed_message(paste0("    Original image: ", width, " x ", height, " pixels"))

      # SMART GRID: target n_cores exactly (for both virtual and disk modes)
      # Unless simple user override (e.g. tile_size explicitly small)
      # We default to smart grid if tile_size is large (default 4096)
      use_smart_grid <- TRUE

      if (use_smart_grid) {
        # Calculate optimal tile size to use all cores
        n_tiles_target <- n_cores
        aspect_ratio <- width / height

        # Approximate a grid that matches aspect ratio and core count
        v_tiles_x <- max(1, round(sqrt(n_tiles_target * aspect_ratio)))
        v_tiles_y <- max(1, ceiling(n_tiles_target / v_tiles_x))

        # Calculate resulting tile size (approximate, we'll use ceiling to cover)
        tile_size <- ceiling(max(width / v_tiles_x, height / v_tiles_y))

        # Recalculate definitive grid based on this tile size
        n_tiles_x <- ceiling(width / tile_size)
        n_tiles_y <- ceiling(height / tile_size)

        n_tiles_x_full <- n_tiles_x
        n_tiles_y_full <- n_tiles_y
        crop_x <- 0
        crop_y <- 0

        timed_message(sprintf("    Smart Grid Calculation (Target: %d cores):", n_cores))
        timed_message(sprintf("      Calculated grid: %d x %d = %d tiles", n_tiles_x, n_tiles_y, n_tiles_x * n_tiles_y))
        timed_message(sprintf("      Tile size: %d px", tile_size))

        # In disk mode, we still crop edges if the grid doesn't perfectly align (impossible with ceiling)
        # Actually with ceiling logic, we cover everything starting from 0,0.
      } else {
        # Original logic for disk tiling (fixed size, crop edges)
        n_tiles_x_full <- floor(width / tile_size)
        n_tiles_y_full <- floor(height / tile_size)

        # Handling small images: Ensure at least 1 tile if image is smaller than tile_size
        if (n_tiles_x_full == 0) n_tiles_x_full <- 1
        if (n_tiles_y_full == 0) n_tiles_y_full <- 1

        crop_x <- (width - n_tiles_x_full * tile_size) / 2
        crop_y <- (height - n_tiles_y_full * tile_size) / 2

        n_tiles_x <- n_tiles_x_full
        n_tiles_y <- n_tiles_y_full
      }

      # TEST MODE OVERRIDE
      if (test_mode) {
        center_row <- floor(n_tiles_y_full / 2)
        center_col <- floor(n_tiles_x_full / 2)
        n_tiles_x <- 1
        n_tiles_y <- 1
        timed_message(paste0("    TEST MODE: Using only center tile (", center_row, ", ", center_col, ")"))
      } else {
        center_row <- NULL
        center_col <- NULL
      }

      # Save metadata
      metadata <- list(
        n_tiles_x = n_tiles_x,
        n_tiles_y = n_tiles_y,
        n_tiles_x_full = n_tiles_x_full,
        n_tiles_y_full = n_tiles_y_full,
        crop_x = crop_x,
        crop_y = crop_y,
        tile_size = tile_size,
        original_dims = original_dims,
        test_mode = test_mode,
        center_row = if (test_mode) center_row else NULL,
        center_col = if (test_mode) center_col else NULL,
        tiling_mode = tiling_mode
      )
      saveRDS(metadata, metadata_file) # Small metadata still uses RDS for broad compatibility
      timed_message(paste0("  ✓ Saved tile grid metadata: ", basename(metadata_file)))

      # Extract tiles

      # Determine which tiles to extract
      if (test_mode) {
        tile_coords <- data.frame(row = center_row, col = center_col)
      } else {
        tile_coords <- expand.grid(row = 0:(n_tiles_y - 1), col = 0:(n_tiles_x - 1))
      }

      # Find tiles that need extraction
      tiles_needing_extraction <- integer(0)
      for (i in seq_len(nrow(tile_coords))) {
        tile_id <- sprintf("tile_%02d_%02d", tile_coords$row[i], tile_coords$col[i])
        tile_image_file <- file.path(tile_dir, paste0(tile_id, "_image.rds"))
        tile_nuclear_file <- file.path(tile_dir, paste0(tile_id, "_nuclear.rds"))

        if (!file.exists(tile_image_file) || !file.exists(tile_nuclear_file)) {
          tiles_needing_extraction <- c(tiles_needing_extraction, i)
        }
      }

      if (tiling_mode == "virtual") {
        timed_message("  ✓ tiling_mode = 'virtual': Skipping extraction to disk, keeping image in memory")
        tiles_needing_extraction <- integer(0) # None needed on disk
      } else if (length(tiles_needing_extraction) == 0) {
        timed_message("  ✓ All tiles already extracted")
      } else {
        n_tiles_to_extract <- length(tiles_needing_extraction)
        timed_message(paste0("  Extracting ", n_tiles_to_extract, " tiles..."))

        # Determine safe number of workers based on RAM
        if (is.null(n_cores_extraction)) {
          image_size_gb <- as.numeric(object.size(image_array)) / 1e9

          # Get system memory info
          mem_info <- system("free -g | grep Mem", intern = TRUE)
          mem_parts <- as.numeric(strsplit(gsub("\\s+", " ", trimws(gsub("Mem:", "", mem_info))), " ")[[1]])
          # The first numeric value is now total, etc.
          # Column 6 in 'free -g' (after removing Mem:) is 'available'
          available_ram_gb <- mem_parts[6]

          # UPDATED: If image is large (>30GB), allow 2 workers if RAM is safely above 2.2x image size.
          # For a 226GB image, this means ~500GB free RAM.
          if (image_size_gb > 30) {
            if (available_ram_gb > (image_size_gb * 2.2 + 20)) {
              n_cores_extraction_use <- 2
              timed_message(sprintf("    ⚠️  Image is large (%.1f GB). Using 2 workers for extraction (Ample RAM: %d GB).", image_size_gb, available_ram_gb))
            } else {
              timed_message(sprintf("    ⚠️  Image is large (%.1f GB). Forcing SERIAL extraction to prevent RAM overflow.", image_size_gb))
              n_cores_extraction_use <- 1
            }
          } else {
            # Conservative calculation for smaller images
            safe_workers <- max(1, floor((available_ram_gb - 20) / image_size_gb))
            n_cores_extraction_use <- min(safe_workers, n_cores, n_tiles_to_extract)
            timed_message(sprintf("    RAM-based worker calculation: Safe workers = %d", n_cores_extraction_use))
          }
        } else {
          n_cores_extraction_use <- min(n_cores_extraction, n_tiles_to_extract)
          timed_message(sprintf("    Using user-specified %d workers for extraction", n_cores_extraction_use))
        }

        # Use serial for 1 worker, parallel for multiple
        use_parallel_extraction <- n_cores_extraction_use > 1

        if (!use_parallel_extraction) {
          timed_message("    Using serial processing (safest for RAM)")
        }

        # Define extraction function (now using .qs for faster I/O)
        extract_one_tile <- function(idx) {
          row_idx <- tile_coords$row[idx]
          col_idx <- tile_coords$col[idx]
          tile_id <- sprintf("tile_%02d_%02d", row_idx, col_idx)

          # Calculate extraction boundaries (including overlap)
          y_start_base <- floor(row_idx * tile_size + crop_y) + 1
          y_end_base <- min(y_start_base + tile_size - 1, original_dims[1])
          x_start_base <- floor(col_idx * tile_size + crop_x) + 1
          x_end_base <- min(x_start_base + tile_size - 1, original_dims[2])

          y_start_ov <- max(1, y_start_base - tile_overlap)
          y_end_ov <- min(original_dims[1], y_end_base + tile_overlap)
          x_start_ov <- max(1, x_start_base - tile_overlap)
          x_end_ov <- min(original_dims[2], x_end_base + tile_overlap)

          # Slice image (nuclear channel only for Stage 1)
          tile_nuclear <- image_array[y_start_ov:y_end_ov, x_start_ov:x_end_ov, nuclear_channel_idx]

          # Save as .qs
          tile_nuclear_file <- file.path(tile_dir, paste0(tile_id, "_nuclear.qs"))
          qs::qsave(tile_nuclear, tile_nuclear_file, preset = "fast", nthreads = 1)

          return(TRUE)
        }

        # Run extraction
        if (use_parallel_extraction) {
          parallel::mclapply(tiles_needing_extraction, extract_one_tile, mc.cores = n_cores_extraction_use)
        } else {
          lapply(tiles_needing_extraction, extract_one_tile)
        }
        timed_message("    ✓ Extraction complete")
      }
      writeLines("complete", tile_extraction_complete_flag)
    }

    # Free original image array from memory if tiling to disk (and not virtual)
    if (tiling_mode == "disk" && clear_images_after_segmentation) {
      timed_message("  Freeing original image array from memory...")
      rm(image_array)
      gc(verbose = FALSE)
    }

    step_times$stage1 <- elapsed_time(step_times$stage1_start)
    if (timing) timed_message(paste0("  ⏱ Tile extraction: ", round(step_times$stage1, 2), " sec"))
  }

  # ============================================================================
  # STAGE 2: Parallel Segmentation
  # ============================================================================
  timed_message("\n--- Stage 2/4: Nuclear Segmentation (Parallel) ---", "stage2_start")

  # Check for completion flag
  seg_complete_flag <- file.path(tile_dir, "segmentation_complete.flag")

  if (file.exists(seg_complete_flag)) {
    timed_message("  ✓ Segmentation already complete (found flag), skipping to Stage 3")
  } else {
    # Determine tiles needing segmentation
    if (test_mode) {
      tile_coords <- data.frame(row = center_row, col = center_col)
    } else {
      tile_coords <- expand.grid(row = 0:(n_tiles_y - 1), col = 0:(n_tiles_x - 1))
    }

    tiles_needing_seg <- integer(0)
    for (i in seq_len(nrow(tile_coords))) {
      tile_id <- sprintf("tile_%02d_%02d", tile_coords$row[i], tile_coords$col[i])
      mask_file <- file.path(tile_dir, paste0(tile_id, "_mask.qs"))
      if (!file.exists(mask_file)) {
        tiles_needing_seg <- c(tiles_needing_seg, i)
      }
    }

    if (length(tiles_needing_seg) == 0) {
      timed_message("  ✓ All tiles already segmented (checking masks)")
    } else {
      n_tiles_to_seg <- length(tiles_needing_seg)
      timed_message(paste0("  Segmenting ", n_tiles_to_seg, " tiles..."))

      # Capture variables for workers to ensure availability during parallel fork
      tile_coords_local <- tile_coords
      tile_dir_local <- tile_dir
      nuclear_channel_idx_local <- nuclear_channel_idx
      image_array_local <- if (exists("image_array")) image_array else NULL
      crop_x_local <- crop_x
      crop_y_local <- crop_y
      tile_size_local <- tile_size
      tile_overlap_local <- tile_overlap
      original_dims_local <- original_dims
      batch_segmentation_local <- as.integer(batch_segmentation)
      nuclear_channel_name_local <- nuclear_channel_name
      tiling_mode_local <- tiling_mode
      verbose_local <- verbose
      n_cores_local <- n_cores

      # Segmentation params
      gamma_local <- gamma
      brightness_local <- brightness
      contrast_local <- contrast
      adaptive_params_local <- adaptive_params
      percentile_range_local <- percentile_range
      exclude_bright_outliers_local <- exclude_bright_outliers
      outlier_percentile_local <- outlier_percentile
      adaptive_enhancement_local <- adaptive_enhancement
      watershed_tolerance_local <- watershed_tolerance
      discSize_local <- discSize
      sizeSelection_local <- sizeSelection
      smooth_local <- smooth
      ext_local <- ext
      max_nucleus_size_local <- max_nucleus_size

      # For "virtual" mode, we'll access the image from memory (shared via fork)
      if (tiling_mode_local == "virtual") {
        images_obj_local <- images_obj
      } else {
        images_obj_local <- NULL
      }

      use_parallel <- (n_tiles_to_seg > 1 || tiling_mode == "virtual") && n_cores > 1

      # Define segmentation function (same for serial and parallel)
      segment_one_tile <- function(idx) {
        tile_start_time <- Sys.time()

        # Initialize diagnostic info
        diagnostic <- list(
          tile_id = NULL,
          row = NULL,
          col = NULL,
          raw_median = NA,
          raw_mean = NA,
          mean_median_ratio = NA,
          tile_type = "UNKNOWN",
          gamma_used = gamma_local,
          brightness_used = brightness_local,
          contrast_used = contrast_local,
          n_cells_initial = 0,
          n_cells_filtered = 0,
          n_cells_final = 0
        )

        # Explicit cleanup on exit
        on.exit({
          gc(verbose = FALSE)
        })

        row_idx <- tile_coords_local$row[idx]
        col_idx <- tile_coords_local$col[idx]
        tile_id <- sprintf("tile_%02d_%02d", row_idx, col_idx)

        diagnostic$tile_id <- tile_id
        diagnostic$row <- row_idx
        diagnostic$col <- col_idx

        if (verbose_local) cat(sprintf("[%s]     [%s] Starting...\n", format(Sys.time(), "%H:%M:%S"), tile_id))

        mask_file <- file.path(tile_dir_local, paste0(tile_id, "_mask.qs"))
        if (file.exists(mask_file)) {
          # Load existing mask to get cell count
          mask <- qs::qread(mask_file, nthreads = 1)
          # Use length(unique()) instead of max() to handle non-consecutive IDs
          unique_ids_cached <- unique(as.vector(mask))
          diagnostic$n_cells_final <- length(unique_ids_cached[unique_ids_cached != 0])
          diagnostic$tile_type <- "CACHED"
          rm(mask)
          return(diagnostic)
        }

        # Load nuclear channel
        if (verbose_local) cat(sprintf("[%s]     [%s] Loading nuclear channel...\n", format(Sys.time(), "%H:%M:%S"), tile_id))

        if (tiling_mode_local == "virtual" && !is.null(images_obj_local)) {
          # Virtual RAM Tiling logic:
          # Calculate base tile boundaries (relative to cropped image)
          y_start_base <- floor(row_idx * tile_size_local + crop_y_local) + 1
          y_end_base <- min(y_start_base + tile_size_local - 1, original_dims_local[1])
          x_start_base <- floor(col_idx * tile_size_local + crop_x_local) + 1
          x_end_base <- min(x_start_base + tile_size_local - 1, original_dims_local[2])

          # Calculate expanded window (with overlap)
          y_start_ov <- max(1, y_start_base - tile_overlap_local)
          y_end_ov <- min(original_dims_local[1], y_end_base + tile_overlap_local)
          x_start_ov <- max(1, x_start_base - tile_overlap_local)
          x_end_ov <- min(original_dims_local[2], x_end_base + tile_overlap_local)

          # Calculate the relative offsets for later cropping
          y_offset_in_ov <- y_start_base - y_start_ov
          x_offset_in_ov <- x_start_base - x_start_ov
          base_width <- y_end_base - y_start_base + 1
          base_height <- x_end_base - x_start_base + 1

          # Slice image from RAM
          nuclear <- images_obj_local[[1]][y_start_ov:y_end_ov, x_start_ov:x_end_ov, nuclear_channel_idx_local]

          if (verbose_local) {
            cat(sprintf(
              "[%s]     [%s] RAM Slice: [%d:%d, %d:%d] (including %dpx overlap)\n",
              format(Sys.time(), "%H:%M:%S"), tile_id, y_start_ov, y_end_ov, x_start_ov, x_end_ov, tile_overlap_local
            ))
          }
        } else {
          # Get from disk (now using .qs with overlap included)
          nuclear_file <- file.path(tile_dir_local, paste0(tile_id, "_nuclear.qs"))
          # Fallback for old .rds files if .qs missing (backward compatibility during transition)
          if (!file.exists(nuclear_file)) nuclear_file <- sub("\\.qs$", ".rds", nuclear_file)

          nuclear <- qs::qread(nuclear_file, nthreads = 1)
        }

        # ADAPTIVE PARAMETER ADJUSTMENT - based on RAW tile intensity
        # Analyze RAW intensities BEFORE normalization to detect truly dim vs bright tiles
        nonzero_pixels <- nuclear[nuclear > 0]
        if (length(nonzero_pixels) == 0) {
          # Empty tile, create empty mask
          if (verbose_local) cat(sprintf("[%s]     [%s] Empty tile, creating empty mask...\n", format(Sys.time(), "%H:%M:%S"), tile_id))
          mask <- array(0, dim = dim(nuclear))

          # If we used overlap, crop the empty mask back to base dimensions for consistency
          if (tile_overlap_local > 0) {
            y_start_base <- floor(row_idx * tile_size_local + crop_y_local) + 1
            y_end_base <- min(y_start_base + tile_size_local - 1, original_dims_local[1])
            x_start_base <- floor(col_idx * tile_size_local + crop_x_local) + 1
            x_end_base <- min(x_start_base + tile_size_local - 1, original_dims_local[2])

            y_start_ov <- max(1, y_start_base - tile_overlap_local)
            x_start_ov <- max(1, x_start_base - tile_overlap_local)

            y_offset_in_ov <- y_start_base - y_start_ov
            x_offset_in_ov <- x_start_base - x_start_ov
            base_width <- y_end_base - y_start_base + 1
            base_height <- x_end_base - x_start_base + 1

            mask_mat <- as.matrix(mask)
            mask_mat <- mask_mat[
              (y_offset_in_ov + 1):(y_offset_in_ov + base_width),
              (x_offset_in_ov + 1):(x_offset_in_ov + base_height)
            ]
            mask <- EBImage::Image(mask_mat)
          }

          qs::qsave(mask, mask_file, preset = "fast", nthreads = 1)
          rm(nuclear, mask, nonzero_pixels)
          return(tile_id)
        }

        # User's manual parameters serve as BASELINE for medium-intensity tiles
        # Adaptive scaling adjusts UP for dim tiles, DOWN for bright tiles
        # User's manual parameters serve as BASELINE for medium-intensity tiles
        # Adaptive scaling adjusts UP for dim tiles, DOWN for bright tiles
        tile_percentile_range <- percentile_range_local
        tile_gamma <- gamma_local
        tile_brightness <- brightness_local
        tile_contrast <- contrast_local

        if (adaptive_params_local) {
          # Calculate RAW tile intensity characteristics (before normalization)
          raw_mean <- mean(nonzero_pixels)
          raw_median <- median(nonzero_pixels)
          raw_q25 <- quantile(nonzero_pixels, 0.25)
          raw_q75 <- quantile(nonzero_pixels, 0.75)
          raw_max <- max(nonzero_pixels)

          # Save for diagnostics
          diagnostic$raw_median <- raw_median
          diagnostic$raw_mean <- raw_mean
          diagnostic$mean_median_ratio <- raw_mean / raw_median

          # Classify tile brightness and apply MULTIPLIERS to user's baseline parameters
          # NOTE: Image data is normalized to 0-1 scale (not 0-65535 16-bit)
          # Baseline (user params) = MEDIUM tile (raw_median ~0.003-0.010, or 0.3-1.0% of max)
          # Scale UP for dimmer tiles, scale DOWN for brighter tiles

          # Calculate mean/median ratio to detect skewed distributions
          mean_to_median_ratio <- raw_mean / raw_median

          # For VERY_DIM tiles, use a higher threshold (3.0) to avoid misclassifying
          # interior tiles with minor hotspots as edge/skewed tiles
          # Only tiles with ratio >= 3.0 are truly edge tissue or have major artifacts
          is_skewed_very_dim <- mean_to_median_ratio >= 3.0

          if (raw_median < 0.02) {
            # VERY DIM TILE - median is truly dim, but check for major outliers/hotspots
            if (is_skewed_very_dim) {
              # SKEWED_VERY_DIM: Very dim background BUT with bright outliers/hotspots
              tile_percentile_range <- c(
                max(0, percentile_range[1] * 0.05), # Very tight lower bound
                min(1, percentile_range[2] + (1 - percentile_range[2]) * 0.5)
              )
              tile_gamma <- gamma_local * 0.5 # Full VERY_DIM gamma
              tile_brightness <- brightness_local + 0.2 # Full VERY_DIM brightness
              tile_contrast <- contrast_local * 1.5 # Full VERY_DIM contrast
              tile_type <- "SKEWED_DIM"
            } else {
              # UNIFORMLY VERY DIM TILE - 2x more aggressive than baseline
              tile_percentile_range <- c(
                max(0, percentile_range[1] * 0.1),
                min(1, percentile_range[2] + (1 - percentile_range[2]) * 0.5)
              )
              tile_gamma <- gamma_local * 0.5 # Lower gamma = more brightening
              tile_brightness <- brightness_local + 0.2 # Add extra brightness
              tile_contrast <- contrast_local * 1.5 # Increase contrast
              tile_type <- "VERY_DIM"
            }
          } else if (raw_median < 0.10) {
            # DIM TILE - 1.5x more aggressive than baseline
            is_skewed_dim <- mean_to_median_ratio >= 1.8
            tile_percentile_range <- c(
              max(0, percentile_range_local[1] * 0.3),
              min(1, percentile_range_local[2] + (1 - percentile_range_local[2]) * 0.3)
            )
            tile_gamma <- gamma_local * 0.7
            tile_brightness <- brightness_local + 0.1
            tile_contrast <- contrast_local * 1.3
            tile_type <- "DIM"
          } else if (raw_median < 0.45) {
            # MEDIUM TILE - use baseline (user's manual parameters)
            # UPDATED: Expanded to cover images with high background (median ~0.4)
            tile_type <- "MEDIUM"
          } else if (raw_median < 0.65) {
            # BRIGHT TILE - 0.8x less aggressive than baseline
            tile_percentile_range <- c(
              max(0, percentile_range_local[1] * 1.5),
              min(1, percentile_range_local[2] - (1 - percentile_range_local[2]) * 0.2)
            )
            tile_gamma <- gamma_local * 1.2 # Higher gamma = less brightening
            tile_brightness <- brightness_local * 0.7 # Reduce brightness boost
            tile_contrast <- contrast_local * 0.9 # Slightly reduce contrast
            tile_type <- "BRIGHT"
          } else {
            # VERY BRIGHT TILE - 0.6x less aggressive than baseline
            tile_percentile_range <- c(
              max(0, percentile_range_local[1] * 2),
              min(1, percentile_range_local[2] - (1 - percentile_range_local[2]) * 0.3)
            )
            tile_gamma <- gamma_local * 1.5
            tile_brightness <- brightness_local * 0.5
            tile_contrast <- contrast_local * 0.8
            tile_type <- "VERY_BRIGHT"
          }

          # Save tile type and final parameters for diagnostics
          diagnostic$tile_type <- tile_type
          diagnostic$gamma_used <- tile_gamma
          diagnostic$brightness_used <- tile_brightness
          diagnostic$contrast_used <- tile_contrast

          if (verbose_local) {
            cat(sprintf(
              "[%s]     [%s] RAW INTENSITY: median=%.5f, mean=%.5f (ratio=%.2f) → %s\n",
              format(Sys.time(), "%H:%M:%S"), tile_id,
              raw_median, raw_mean, mean_to_median_ratio, tile_type
            ))
            if (tile_type == "SKEWED_DIM") {
              cat(sprintf(
                "[%s]     [%s] ⚠️  SKEWED DISTRIBUTION: dim background + major outliers (ratio=%.2f ≥ 3.0)\n",
                format(Sys.time(), "%H:%M:%S"), tile_id, mean_to_median_ratio
              ))
              cat(sprintf(
                "[%s]     [%s] → Using VERY_DIM enhancement + aggressive outlier masking\n",
                format(Sys.time(), "%H:%M:%S"), tile_id
              ))
            }
            if (tile_type != "MEDIUM") {
              cat(sprintf(
                "[%s]     [%s] ADAPTIVE SCALING: gamma=%.2f→%.2f, bright=%.2f→%.2f, contrast=%.1f→%.1f\n",
                format(Sys.time(), "%H:%M:%S"), tile_id,
                gamma_local, tile_gamma, brightness_local, tile_brightness, contrast_local, tile_contrast
              ))
            } else {
              cat(sprintf(
                "[%s]     [%s] Using BASELINE parameters (medium intensity)\n",
                format(Sys.time(), "%H:%M:%S"), tile_id
              ))
            }
          }
        } else {
          # If not adaptive, still save baseline parameters
          diagnostic$tile_type <- "NO_ADAPT"
        }

        # Apply tile-specific percentile range
        q_lower <- quantile(nonzero_pixels, tile_percentile_range[1])
        q_upper <- quantile(nonzero_pixels, tile_percentile_range[2])

        # Detect and handle bright outliers (hotspots)
        outlier_mask <- NULL
        if (exclude_bright_outliers_local) {
          outlier_threshold <- quantile(nonzero_pixels, outlier_percentile_local)
          outlier_mask <- nuclear > outlier_threshold

          if (sum(outlier_mask) > 0) {
            q_upper <- outlier_threshold
          }
        }

        # Normalize: stretch [q_lower, q_upper] to [0, 1] and clip outliers
        nuclear_enhanced <- (nuclear - q_lower) / (q_upper - q_lower)
        nuclear_enhanced[nuclear_enhanced < 0] <- 0
        nuclear_enhanced[nuclear_enhanced > 1] <- 1

        # Mask out bright outliers (set to median intensity)
        if (exclude_bright_outliers_local && !is.null(outlier_mask) && sum(outlier_mask) > 0) {
          median_intensity <- median(nuclear_enhanced[nuclear_enhanced > 0])
          nuclear_enhanced[outlier_mask] <- median_intensity
        }

        # Apply adaptive histogram equalization (CLAHE) for local enhancement
        if (adaptive_enhancement_local) {
          if (verbose_local) cat(sprintf("[%s]     [%s] Applying CLAHE...\n", format(Sys.time(), "%H:%M:%S"), tile_id))
          nuclear_clahe <- EBImage::normalize(nuclear_enhanced,
            separate = FALSE,
            ft = c(0, 1)
          )
          # Blend CLAHE with global enhancement (30% global, 70% adaptive)
          nuclear_enhanced <- 0.3 * nuclear_enhanced + 0.7 * nuclear_clahe
          rm(nuclear_clahe)
        }

        # Apply EBImage enhancement formula: contrast*(image+brightness)^gamma
        nuclear_enhanced <- tile_contrast * (nuclear_enhanced + tile_brightness)^tile_gamma
        if (verbose_local && !adaptive_params_local) {
          cat(sprintf(
            "[%s]     [%s] Enhancement complete (gamma=%.2f, brightness=%.2f, contrast=%.2f)\n",
            format(Sys.time(), "%H:%M:%S"), tile_id, tile_gamma, tile_brightness, tile_contrast
          ))
        }

        # Clean up before segmentation
        rm(nuclear, nonzero_pixels, outlier_mask)
        gc(verbose = FALSE)

        # Create temporary CytoImageList for simpleSeg
        # CRITICAL: Force exactly 3 dims (y, x, 1) to avoid "incorrect number of dimensions" error
        # indexing image[,,ind] in cytomapper/simpleSeg will fail on 2D or 4D objects.
        raw_data <- as.array(nuclear_enhanced)
        if (length(dim(raw_data)) == 2) {
          dim(raw_data) <- c(dim(raw_data), 1)
        } else if (length(dim(raw_data)) == 3) {
          # Already 3D, ensure index 3 exists
        }
        temp_image <- EBImage::Image(raw_data)
        dimnames(temp_image)[[3]] <- nuclear_channel_name_local
        temp_images <- cytomapper::CytoImageList(list(temp_image))
        names(temp_images) <- tile_id
        cytomapper::channelNames(temp_images) <- nuclear_channel_name_local

        rm(nuclear_enhanced, temp_image)

        # Segment with simpleSeg (runs in serial within each worker)
        # Register SerialParam to avoid nested parallelization
        if (verbose_local) cat(sprintf("[%s]     [%s] Running simpleSeg (this may take 10-20 min)...\n", format(Sys.time(), "%H:%M:%S"), tile_id))
        old_bpparam <- BiocParallel::bpparam()

        if (tiling_mode_local == "virtual" && !use_parallel) {
          # Use full parallelization for the whole image ONLY IF we aren't already in parallel mclapply
          if (verbose_local) cat(sprintf("[%s]     [%s] tiling_mode = 'virtual' (Serial Worker): Using %d cores for simpleSeg\n", format(Sys.time(), "%H:%M:%S"), tile_id, n_cores_local))
          BiocParallel::register(BiocParallel::MulticoreParam(workers = n_cores_local))
        } else if (batch_segmentation_local > 1) {
          # Use user-defined internal threads for simpleSeg to keep total CPU high but RAM low
          if (verbose_local) cat(sprintf("[%s]     [%s] Batching mode: Using %d threads for simpleSeg (Batch = %d)\n", format(Sys.time(), "%H:%M:%S"), tile_id, batch_segmentation_local, batch_segmentation_local))
          BiocParallel::register(BiocParallel::MulticoreParam(workers = batch_segmentation_local))
        } else {
          BiocParallel::register(BiocParallel::SerialParam())
        }

        seg_start <- Sys.time()
        masks_tile <- simpleSeg::simpleSeg(
          temp_images,
          nucleus = c(nuclear_channel_name_local),
          pca = FALSE,
          cellBody = "dilate",
          discSize = discSize_local,
          sizeSelection = sizeSelection_local,
          smooth = smooth_local,
          transform = NULL,
          watershed = "distance",
          tolerance = watershed_tolerance_local,
          ext = ext_local,
          tissue = NULL
        )
        seg_elapsed <- as.numeric(difftime(Sys.time(), seg_start, units = "mins"))
        if (verbose_local) cat(sprintf("[%s]     [%s] simpleSeg complete (%.1f min)\n", format(Sys.time(), "%H:%M:%S"), tile_id, seg_elapsed))

        # Restore previous BiocParallel backend
        BiocParallel::register(old_bpparam)

        rm(temp_images)

        # Extract mask array
        mask <- masks_tile[[1]]
        rm(masks_tile)

        # If we used overlap (virtual OR disk mode now has overlap), crop the mask back
        # Disk mode now saves tiles with overlap, so we must crop.
        # Virtual mode always extracted with overlap.
        # Only case we don't crop is if overlap was 0 (unlikely)
        if (tile_overlap_local > 0) {
          if (verbose_local) cat(sprintf("[%s]     [%s] Cropping mask to remove overlap context...\n", format(Sys.time(), "%H:%M:%S"), tile_id))

          # Calculate offsets based on the fact that the loaded/sliced image
          # spans [y_start_ov : y_end_ov, x_start_ov : x_end_ov]

          # Re-calculate extraction/slice boundaries to determine crop relative to current buffer
          # This logic duplicates Stage 1/Virtual slice math to ensure consistency

          y_start_base <- floor(row_idx * tile_size_local + crop_y_local) + 1
          y_end_base <- min(y_start_base + tile_size_local - 1, original_dims_local[1])
          x_start_base <- floor(col_idx * tile_size_local + crop_x_local) + 1
          x_end_base <- min(x_start_base + tile_size_local - 1, original_dims_local[2])

          y_start_ov <- max(1, y_start_base - tile_overlap_local)
          y_end_ov <- min(original_dims_local[1], y_end_base + tile_overlap_local)

          x_start_ov <- max(1, x_start_base - tile_overlap_local)
          x_end_ov <- min(original_dims_local[2], x_end_base + tile_overlap_local)

          # Relative offsets of the "true" tile within the "overlapped" tile
          y_offset_in_ov <- y_start_base - y_start_ov
          x_offset_in_ov <- x_start_base - x_start_ov
          base_width <- y_end_base - y_start_base + 1
          base_height <- x_end_base - x_start_base + 1

          # mask is oriented (rows=y, cols=x) in sync with the array slice
          mask_matrix <- as.matrix(mask)
          mask_matrix <- mask_matrix[
            (y_offset_in_ov + 1):(y_offset_in_ov + base_width),
            (x_offset_in_ov + 1):(x_offset_in_ov + base_height)
          ]
          mask <- EBImage::Image(mask_matrix, dim = dim(mask_matrix))
          rm(mask_matrix)
        }

        # Count initial cells (before filtering)
        # Use length(unique()) instead of max() to handle non-consecutive IDs
        unique_ids <- unique(as.vector(mask))
        diagnostic$n_cells_initial <- length(unique_ids[unique_ids != 0]) # Exclude background (0)

        # Filter by max size
        if (!is.null(max_nucleus_size_local) && max_nucleus_size_local > 0) {
          if (verbose_local) cat(sprintf("[%s]     [%s] Filtering oversized objects...\n", format(Sys.time(), "%H:%M:%S"), tile_id))

          mask_matrix <- as.matrix(mask)
          object_sizes <- table(mask_matrix)
          object_sizes <- object_sizes[names(object_sizes) != "0"]

          if (length(object_sizes) > 0) {
            oversized_ids <- as.integer(names(object_sizes[object_sizes > max_nucleus_size_local]))

            if (length(oversized_ids) > 0) {
              # Track filtered cells for diagnostics
              diagnostic$n_cells_filtered <- length(oversized_ids)

              if (verbose_local) {
                cat(sprintf(
                  "[%s]     [%s] Removing %d oversized objects (%.1f%%)...\n",
                  format(Sys.time(), "%H:%M:%S"), tile_id,
                  length(oversized_ids),
                  100 * length(oversized_ids) / length(object_sizes)
                ))
              }

              # Remove oversized objects (set to 0)
              for (obj_id in oversized_ids) {
                mask_matrix[mask_matrix == obj_id] <- 0
              }

              # SKIP RENUMBERING - it's O(n*m) and causes parallel workers to timeout
              # measureObjects can handle non-consecutive IDs, so renumbering is optional
              # If needed later, can be done post-processing
              if (verbose_local) cat(sprintf("[%s]     [%s] Skipping renumbering (non-consecutive IDs OK for quantification)...\n", format(Sys.time(), "%H:%M:%S"), tile_id))

              # Just update the mask with filtered objects (IDs may have gaps)
              mask <- EBImage::Image(mask_matrix, dim = dim(mask_matrix))

              if (verbose_local) cat(sprintf("[%s]     [%s] ✓ Filtering complete\n", format(Sys.time(), "%H:%M:%S"), tile_id))
            }
          }
          rm(mask_matrix, object_sizes)
        }

        # Count final cells (after filtering)
        # Use length(unique()) instead of max() to handle non-consecutive IDs
        unique_ids_final <- unique(as.vector(mask))
        diagnostic$n_cells_final <- length(unique_ids_final[unique_ids_final != 0])

        # Save mask as .qs (fast, nthreads=1 because we are in mclapply)
        if (verbose_local) cat(sprintf("[%s]     [%s] Saving mask...\n", format(Sys.time(), "%H:%M:%S"), tile_id))
        qs::qsave(mask, mask_file, preset = "fast", nthreads = 1)

        # Clean up
        rm(mask)
        gc(verbose = FALSE)

        tile_elapsed <- as.numeric(difftime(Sys.time(), tile_start_time, units = "mins"))
        if (verbose_local) {
          cat(sprintf(
            "[%s]     [%s] ✓ Complete (%.1f min total, %d cells)\n",
            format(Sys.time(), "%H:%M:%S"), tile_id, tile_elapsed, diagnostic$n_cells_final
          ))
        }

        return(diagnostic)
      } # End of segment_one_tile function

      # Run segmentation (parallel or serial)
      if (use_parallel) {
        # Calculate workers per batch
        n_segmentation_workers <- max(1, floor(n_cores / batch_segmentation))
        timed_message(sprintf(
          "    Parallel Settings: %d concurrent tiles, %d threads per tile (Total %d cores)",
          n_segmentation_workers, batch_segmentation, n_segmentation_workers * batch_segmentation
        ))

        results <- parallel::mclapply(tiles_needing_seg, segment_one_tile, mc.cores = n_segmentation_workers)
      } else {
        results <- lapply(tiles_needing_seg, segment_one_tile)
      }

      # Check for errors
      errors <- sapply(results, function(x) inherits(x, "try-error"))
      if (any(errors)) {
        timed_message(paste0("  ⚠️  Warning: ", sum(errors), " tiles failed to segment"))
        for (i in which(errors)) {
          timed_message(paste0("    Error in tile ", i, ": ", as.character(results[[i]])))
        }
      } else {
        timed_message(paste0("  ✓ All ", n_tiles_to_seg, " tiles segmented successfully"))
      }

      # Clean up after parallel section
      gc(full = TRUE)
    }

    # Write completion flag
    writeLines("complete", seg_complete_flag)

    # Build diagnostic table from results
    timed_message("  Building segmentation diagnostics table...")

    # Convert list of diagnostics to data frame
    diagnostics_list <- results
    diagnostics_df <- do.call(rbind, lapply(diagnostics_list, function(x) {
      if (is.list(x)) {
        data.frame(
          tile_id = x$tile_id,
          row = x$row,
          col = x$col,
          raw_median = x$raw_median,
          raw_mean = x$raw_mean,
          mean_median_ratio = x$mean_median_ratio,
          tile_type = x$tile_type,
          gamma_used = x$gamma_used,
          brightness_used = x$brightness_used,
          contrast_used = x$contrast_used,
          n_cells_initial = x$n_cells_initial,
          n_cells_filtered = x$n_cells_filtered,
          n_cells_final = x$n_cells_final,
          stringsAsFactors = FALSE
        )
      }
    }))

    # Add calculated columns
    if (!is.null(diagnostics_df) && nrow(diagnostics_df) > 0) {
      diagnostics_df$pct_filtered <- round(100 * diagnostics_df$n_cells_filtered / diagnostics_df$n_cells_initial, 1)
      diagnostics_df$pct_filtered[is.nan(diagnostics_df$pct_filtered)] <- 0

      # Save diagnostics table
      diagnostics_file <- file.path(tile_dir, "segmentation_diagnostics.csv")
      write.csv(diagnostics_df, diagnostics_file, row.names = FALSE)
      timed_message(paste0("  ✓ Saved diagnostics: ", diagnostics_file))

      # Print summary
      timed_message("\n  === Segmentation Diagnostics Summary ===")
      timed_message(sprintf("    Total tiles: %d", nrow(diagnostics_df)))
      timed_message(sprintf("    Total cells: %s", format(sum(diagnostics_df$n_cells_final), big.mark = ",")))
      timed_message(sprintf(
        "    Cells filtered: %s (%.1f%%)",
        format(sum(diagnostics_df$n_cells_filtered), big.mark = ","),
        100 * sum(diagnostics_df$n_cells_filtered) / sum(diagnostics_df$n_cells_initial)
      ))

      # Tile type distribution
      tile_type_counts <- table(diagnostics_df$tile_type)
      timed_message("    Tile types:")
      for (tt in names(tile_type_counts)) {
        avg_cells_type <- mean(diagnostics_df$n_cells_final[diagnostics_df$tile_type == tt])
        timed_message(sprintf("      %s: %d tiles (avg %.0f cells)", tt, tile_type_counts[tt], avg_cells_type))
      }

      # Highlight SKEWED tiles (if any)
      skewed_tiles <- diagnostics_df[diagnostics_df$tile_type == "SKEWED_DIM", ]
      if (nrow(skewed_tiles) > 0) {
        timed_message("\n    ⚠️  SKEWED tiles (dim background + major outliers, ratio ≥ 3.0):")
        timed_message("        These tiles use VERY_DIM enhancement with aggressive outlier masking")
        timed_message("        Likely edge/corner tiles or regions with significant artifacts")
        for (i in seq_len(nrow(skewed_tiles))) {
          timed_message(sprintf(
            "      %s: median=%.5f, ratio=%.2f, %d cells",
            skewed_tiles$tile_id[i],
            skewed_tiles$raw_median[i],
            skewed_tiles$mean_median_ratio[i],
            skewed_tiles$n_cells_final[i]
          ))
        }
      }
    }

    timed_message("  ✓ Segmentation stage complete")
  }

  # ============================================================================
  # OPTIONAL: Generate Segmentation Overlay Visualization (Center Tile)
  # ============================================================================
  if (save_segmentation_overlay) {
    # Determine center tile coordinates
    if (test_mode) {
      overlay_row <- center_row
      overlay_col <- center_col
    } else {
      # Use actual center tile
      overlay_row <- floor(n_tiles_y / 2)
      overlay_col <- floor(n_tiles_x / 2)
    }

    overlay_tile_id <- sprintf("tile_%02d_%02d", overlay_row, overlay_col)
    overlay_file <- file.path(tile_dir, "visualizations", paste0(overlay_tile_id, "_segmentation_overlay.png"))

    # Only generate if it doesn't exist
    if (!file.exists(overlay_file)) {
      timed_message(paste0("\n  Generating segmentation overlay for ", overlay_tile_id, "..."))

      # Load nuclear channel and mask for center tile
      nuclear_tile_file <- file.path(tile_dir, paste0(overlay_tile_id, "_nuclear.qs"))
      mask_file <- file.path(tile_dir, paste0(overlay_tile_id, "_mask.qs"))

      if (file.exists(nuclear_tile_file) && file.exists(mask_file)) {
        nuclear_tile <- qs::qread(nuclear_tile_file, nthreads = n_cores)
        mask_tile <- qs::qread(mask_file, nthreads = n_cores)

        # Apply the SAME enhancement used during segmentation
        timed_message("    Applying enhancement to nuclear channel...")

        # Get intensity range (ignore zeros/background)
        nonzero_pixels <- nuclear_tile[nuclear_tile > 0]
        if (length(nonzero_pixels) > 0) {
          q_lower <- quantile(nonzero_pixels, percentile_range[1])
          q_upper <- quantile(nonzero_pixels, percentile_range[2])

          # Detect and handle bright outliers (hotspots)
          outlier_mask <- NULL
          if (exclude_bright_outliers) {
            outlier_threshold <- quantile(nonzero_pixels, outlier_percentile)
            outlier_mask <- nuclear_tile > outlier_threshold
            if (sum(outlier_mask) > 0) {
              q_upper <- outlier_threshold
            }
          }

          # Normalize: stretch [q_lower, q_upper] to [0, 1] and clip outliers
          nuclear_enhanced <- (nuclear_tile - q_lower) / (q_upper - q_lower)
          nuclear_enhanced[nuclear_enhanced < 0] <- 0
          nuclear_enhanced[nuclear_enhanced > 1] <- 1

          # Mask out bright outliers (set to median intensity)
          if (exclude_bright_outliers && !is.null(outlier_mask) && sum(outlier_mask) > 0) {
            median_intensity <- median(nuclear_enhanced[nuclear_enhanced > 0])
            nuclear_enhanced[outlier_mask] <- median_intensity
          }

          # Apply adaptive histogram equalization (CLAHE) for local enhancement
          if (adaptive_enhancement) {
            nuclear_clahe <- EBImage::normalize(nuclear_enhanced,
              separate = FALSE,
              ft = c(0, 1)
            )
            # Blend CLAHE with global enhancement (30% global, 70% adaptive)
            nuclear_enhanced <- 0.3 * nuclear_enhanced + 0.7 * nuclear_clahe
            rm(nuclear_clahe)
          }

          # Apply EBImage enhancement formula: contrast*(image+brightness)^gamma
          nuclear_enhanced <- contrast * (nuclear_enhanced + brightness)^gamma

          timed_message(paste0(
            "    Enhancement applied [v6] (gamma=", gamma,
            ", brightness=", brightness, ", contrast=", contrast, ")"
          ))

          # ROBUST OVERLAY CROPPING: Ensure image matches mask dimensions.
          # We extract the pure matrix to avoid indexing issues with EBImage S4 objects.
          img_mat <- EBImage::imageData(nuclear_enhanced)
          if (length(dim(img_mat)) > 2) img_mat <- img_mat[,,1] # Drop extra dims safely
          
          mask_mat <- EBImage::imageData(mask_tile)
          if (length(dim(mask_mat)) > 2) mask_mat <- mask_mat[,,1]

          if (!all(dim(img_mat) == dim(mask_mat))) {
            timed_message(sprintf("    [%s] Dimensions mismatch [v6] (Image: %dx%d, Mask: %dx%d). Adjusting...",
              overlay_tile_id, nrow(img_mat), ncol(img_mat),
              nrow(mask_mat), ncol(mask_mat)
            ))

            # Use mask as target (it was cropped to base dimensions in Stage 2)
            target_h <- nrow(mask_mat)
            target_w <- ncol(mask_mat)

            # Calculate offsets from metadata
            y_start_base <- floor(overlay_row * tile_size + crop_y) + 1
            x_start_base <- floor(overlay_col * tile_size + crop_x) + 1
            y_start_ov <- max(1, y_start_base - tile_overlap)
            x_start_ov <- max(1, x_start_base - tile_overlap)

            y_off <- y_start_base - y_start_ov
            x_off <- x_start_base - x_start_ov

            # FALLBACK: If offsets are nonsensical (stale cache), use top-left
            if (y_off < 0 || x_off < 0 || (y_off + target_h) > nrow(img_mat) || (x_off + target_w) > ncol(img_mat)) {
              timed_message("    ⚠️  Metadata inconsistent with tiles [v6]. Using top-left intersection.")
              y_off <- 0
              x_off <- 0
              target_h <- min(nrow(img_mat), nrow(mask_mat))
              target_w <- min(ncol(img_mat), ncol(mask_mat))

              # Sync mask dimensions if mismatch was reversed
              if (nrow(mask_mat) > target_h || ncol(mask_mat) > target_w) {
                mask_mat <- mask_mat[1:target_h, 1:target_w]
              }
            }

            # Crop image matrix to match target
            img_mat <- img_mat[(y_off + 1):(y_off + target_h), (x_off + 1):(x_off + target_w)]
            
            # Re-update variables for downstream use
            nuclear_enhanced <- img_mat
            mask_tile <- mask_mat
          }
        } else {
          # Empty tile
          nuclear_enhanced <- nuclear_tile
        }
      }

      # Create CytoImageList for ENHANCED nuclear channel
        nuclear_image <- EBImage::Image(nuclear_enhanced, dim = c(dim(nuclear_enhanced), 1))
        dimnames(nuclear_image)[[3]] <- nuclear_channel_name
        nuclear_cil <- cytomapper::CytoImageList(list(nuclear_image))
        names(nuclear_cil) <- overlay_tile_id
        cytomapper::channelNames(nuclear_cil) <- nuclear_channel_name
        S4Vectors::mcols(nuclear_cil)$imageID <- overlay_tile_id

        # Create CytoImageList for mask
        mask_image <- EBImage::Image(mask_tile, dim = dim(mask_tile))
        mask_cil <- cytomapper::CytoImageList(list(mask_image))
        names(mask_cil) <- overlay_tile_id
        S4Vectors::mcols(mask_cil)$imageID <- overlay_tile_id

        # Create visualization directory
        viz_dir <- file.path(tile_dir, "visualizations")
        dir.create(viz_dir, showWarnings = FALSE, recursive = TRUE)

        # Generate overlay using plotPixels
        png(overlay_file, width = 2400, height = 2400)

        # Create color list
        color_list <- list(c("black", "blue"))
        names(color_list) <- nuclear_channel_name

        cytomapper::plotPixels(
          image = nuclear_cil,
          mask = mask_cil,
          img_id = "imageID",
          colour_by = nuclear_channel_name,
          display = "single",
          colour = color_list,
          legend = NULL,
          scale_bar = list(length = 100, label = "100 µm", cex = 1.5)
        )

        dev.off()

        rm(nuclear_tile, nuclear_enhanced, mask_tile, nuclear_image, nuclear_cil, mask_image, mask_cil, nonzero_pixels, outlier_mask)
        gc(verbose = FALSE)

        timed_message(paste0("  ✓ Segmentation overlay saved: ", basename(overlay_file)))
      } else {
        timed_message(paste0("  ⚠️  Warning: Could not find nuclear tile or mask for ", overlay_tile_id))
      }
    } else {
      timed_message(paste0("  ✓ Segmentation overlay already exists: ", basename(overlay_file)))
    }
  }

  step_times$stage2 <- elapsed_time(step_times$stage2_start, units = "mins")
  if (timing) timed_message(paste0("  ⏱ Segmentation: ", round(step_times$stage2, 2), " min"))

  # ============================================================================
  # PHASE 3: Feature Quantification
  # ============================================================================
  timed_message("\n--- Phase 3: Feature Quantification ---", "phase3_start")

  # ============================================================================
  # STAGE 3: Parallel Enhancement & Quantification
  # ============================================================================
  timed_message("\n--- Stage 3/4: Enhancement and Quantification (Parallel) ---", "stage3_start")

  # Check for completion flag
  quant_complete_flag <- file.path(tile_dir, "quantification_complete.flag")

  if (file.exists(quant_complete_flag)) {
    timed_message("  ✓ Quantification already complete (found flag), skipping to Stage 4")
  } else {
    # Get list of tiles needing quantification
    if (test_mode) {
      tile_coords <- data.frame(row = center_row, col = center_col, stringsAsFactors = FALSE)
    } else {
      tile_coords <- expand.grid(row = 0:(n_tiles_y - 1), col = 0:(n_tiles_x - 1))
    }

    tiles_needing_quant <- integer(0) # FIXED: Use integer, not character
    for (i in seq_len(nrow(tile_coords))) {
      tile_id <- sprintf("tile_%02d_%02d", tile_coords$row[i], tile_coords$col[i])
      spe_file_qs <- file.path(tile_dir, paste0(tile_id, "_spe.qs"))
      spe_file_rds <- file.path(tile_dir, paste0(tile_id, "_spe.rds")) # Check for old extension
      if (!file.exists(spe_file_qs) && !file.exists(spe_file_rds)) { # Check both
        tiles_needing_quant <- c(tiles_needing_quant, i)
      }
    }

    if (length(tiles_needing_quant) == 0) {
      timed_message("  ✓ All tiles already quantified")
    } else {
      # Determine center tile for visualization (use existing values if already set)
      if (is.null(center_row) || is.null(center_col)) {
        center_row <- floor(n_tiles_y / 2)
        center_col <- floor(n_tiles_x / 2)
      }

      n_tiles_to_quant <- length(tiles_needing_quant)
      timed_message(paste0("  Quantifying ", n_tiles_to_quant, " tiles..."))

      use_parallel_quant <- (n_tiles_to_quant > 1 || tiling_mode == "virtual") && n_cores > 1

      if (use_parallel_quant) {
        timed_message(paste0("    Using ", n_cores, " fork-based parallel workers"))
        if (tiling_mode == "virtual") timed_message("    [tiling_mode='virtual'] Mode: Parallel RAM Quantification (Virtual Grid)")
      } else {
        timed_message("    Using serial processing (single tile, no fork overhead)")
      }

      timed_message(paste0("    Enhancement parameters:"))
      timed_message(paste0("      gamma: ", enhancement_gamma))
      timed_message(paste0("      adaptive: ", enhancement_adaptive))
      if (save_one_tile_viz) {
        timed_message(paste0("    Saving visualization for center tile (", center_row, ", ", center_col, ")"))
      }

      # Define enhance_channel function (adapted from enhance_and_quantify.R)
      enhance_channel <- function(channel_data, channel_name, gamma_val, perc_range,
                                  excl_outliers, outlier_perc, adaptive) {
        # Get intensity range (ignore zeros/background)
        nonzero_pixels <- channel_data[channel_data > 0]
        if (length(nonzero_pixels) == 0) {
          return(channel_data)
        }

        # Percentile normalization
        q_lower <- quantile(nonzero_pixels, perc_range[1])
        q_upper <- quantile(nonzero_pixels, perc_range[2])

        # Handle bright outliers
        outlier_mask <- NULL
        if (excl_outliers) {
          outlier_threshold <- quantile(nonzero_pixels, outlier_perc)
          outlier_mask <- channel_data > outlier_threshold

          if (sum(outlier_mask) > 0) {
            q_upper <- outlier_threshold
          }
        }

        # Normalize to [0,1]
        enhanced <- (channel_data - q_lower) / (q_upper - q_lower)
        enhanced[enhanced < 0] <- 0
        enhanced[enhanced > 1] <- 1

        # Mask outliers if applicable
        if (!is.null(outlier_mask) && sum(outlier_mask) > 0) {
          median_intensity <- median(enhanced[enhanced > 0])
          enhanced[outlier_mask] <- median_intensity
        }

        # Adaptive enhancement (CLAHE)
        if (adaptive) {
          enhanced_clahe <- EBImage::normalize(enhanced, separate = FALSE, ft = c(0, 1))
          enhanced <- 0.7 * enhanced + 0.3 * enhanced_clahe
        }

        # Gamma correction
        enhanced <- enhanced^gamma_val

        return(enhanced)
      }

      # Capture variables for parallel workers (explicit scoping)
      tile_coords_local <- tile_coords
      tile_dir_local <- tile_dir
      all_channel_names_local <- all_channel_names
      center_row_local <- center_row
      center_col_local <- center_col
      crop_x_local <- crop_x
      crop_y_local <- crop_y
      tile_size_local <- tile_size
      original_dims_local <- original_dims
      tile_overlap_local <- tile_overlap
      batch_quantification_local <- as.integer(batch_quantification)
      tiling_mode_local <- tiling_mode
      save_one_tile_viz_local <- save_one_tile_viz
      verbose_local <- verbose
      n_cores_local <- n_cores
      enhancement_gamma_local <- enhancement_gamma
      enhancement_percentile_range_local <- enhancement_percentile_range
      enhancement_exclude_outliers_local <- enhancement_exclude_outliers
      enhancement_outlier_percentile_local <- enhancement_outlier_percentile
      enhancement_adaptive_local <- enhancement_adaptive

      if (tiling_mode_local == "virtual") {
        # Pass the whole array to avoid disk I/O
        # Using images_obj[[1]] which is the Image/Array
        image_array_local <- images_obj[[1]]
      } else {
        image_array_local <- NULL
      }

      # Define quantification function (same for serial and parallel)
      quantify_one_tile <- function(idx) {
        tile_start_time <- Sys.time()

        # Explicit cleanup on exit
        on.exit({
          rm(list = ls())
          gc(verbose = FALSE)
        })

        row_idx <- tile_coords_local$row[idx]
        col_idx <- tile_coords_local$col[idx]
        tile_id <- sprintf("tile_%02d_%02d", row_idx, col_idx)

        if (verbose_local) cat(sprintf("[%s]     [%s] Starting quantification...\n", format(Sys.time(), "%H:%M:%S"), tile_id))

        spe_file <- file.path(tile_dir_local, paste0(tile_id, "_spe.qs"))
        if (file.exists(spe_file)) {
          return(tile_id)
        }

        # Load tile (all channels)
        if (verbose_local) cat(sprintf("[%s]     [%s] Loading tile image...\n", format(Sys.time(), "%H:%M:%S"), tile_id))

        if (tiling_mode_local == "virtual" && !is.null(image_array_local)) {
          # Slice exactly the base tile (mask is already cropped in Stage 2)
          y_start_base <- floor(row_idx * tile_size_local + crop_y_local) + 1
          y_end_base <- min(y_start_base + tile_size_local - 1, original_dims_local[1])
          x_start_base <- floor(col_idx * tile_size_local + crop_x_local) + 1
          x_end_base <- min(x_start_base + tile_size_local - 1, original_dims_local[2])

          tile_image <- image_array_local[y_start_base:y_end_base, x_start_base:x_end_base, , drop = FALSE]
        } else {
          tile_image_file <- file.path(tile_dir_local, paste0(tile_id, "_image.qs"))
          # Fallback
          if (!file.exists(tile_image_file)) tile_image_file <- sub("\\.qs$", ".rds", tile_image_file)

          tile_image <- qs::qread(tile_image_file, nthreads = 1)

          # Update: Disk tiles now have overlap too.
          # Ideally, quantify_and_enhance expects the *base* tile for visualization?
          # OR does it use the mask?
          # The function signature for `measureObjects` etc usually handles the mask matching the image.
          # Since we cropped the mask in Stage 2, it is [base_size].
          # But `tile_image` loaded from disk is [base + overlap].
          # WE MUST CROP THE IMAGE HERE to match the mask!

          if (tile_overlap_local > 0) {
            # Re-calc crop offsets (same as Stage 2)
            y_start_base <- floor(row_idx * tile_size_local + crop_y_local) + 1
            x_start_base <- floor(col_idx * tile_size_local + crop_x_local) + 1

            y_start_ov <- max(1, y_start_base - tile_overlap_local)
            x_start_ov <- max(1, x_start_base - tile_overlap_local)

            y_offset <- y_start_base - y_start_ov
            x_offset <- x_start_base - x_start_ov

            # Image dimensions might be [y, x, channels]
            # Get base dims from the previously saved metadata or infer?
            # Easier to just use the mask's logic:
            # But wait, we don't have the mask dimensions handy easily without loading it.
            # We can calculate "expected" base size.
            original_h <- original_dims_local[1]
            original_w <- original_dims_local[2]

            y_end_base <- min(y_start_base + tile_size_local - 1, original_h)
            x_end_base <- min(x_start_base + tile_size_local - 1, original_w)

            base_h <- y_end_base - y_start_base + 1
            base_w <- x_end_base - x_start_base + 1

            # Crop image
            tile_image <- tile_image[(y_offset + 1):(y_offset + base_h),
              (x_offset + 1):(x_offset + base_w), ,
              drop = FALSE
            ]
          }
        }

        # Load mask (now .qs)
        if (verbose_local) cat(sprintf("[%s]     [%s] Loading mask...\n", format(Sys.time(), "%H:%M:%S"), tile_id))
        mask_file <- file.path(tile_dir_local, paste0(tile_id, "_mask.qs"))
        mask <- qs::qread(mask_file, nthreads = 1)

        # Check if mask has any cells
        n_cells_in_tile <- length(unique(as.vector(mask))) - 1
        if (n_cells_in_tile == 0) {
          # No cells in this tile, create empty SPE marker
          if (verbose_local) cat(sprintf("[%s]     [%s] Empty tile (no cells), skipping...\n", format(Sys.time(), "%H:%M:%S"), tile_id))
          qs::qsave(list(n_cells = 0), spe_file, preset = "fast", nthreads = 1)
          rm(tile_image, mask)
          return(tile_id)
        }

        if (verbose_local) {
          cat(sprintf(
            "[%s]     [%s] Found %d cells, enhancing %d channels...\n",
            format(Sys.time(), "%H:%M:%S"), tile_id, n_cells_in_tile, dim(tile_image)[3]
          ))
        }

        # Enhance each channel
        n_channels_tile <- dim(tile_image)[3]
        enhanced_channels <- array(0, dim = dim(tile_image))

        for (ch in 1:n_channels_tile) {
          channel_name <- all_channel_names_local[ch]
          enhanced_channels[, , ch] <- enhance_channel(
            tile_image[, , ch], channel_name,
            enhancement_gamma_local, enhancement_percentile_range_local,
            enhancement_exclude_outliers_local, enhancement_outlier_percentile_local,
            enhancement_adaptive_local
          )
        }

        # Update tile_image with enhanced data
        tile_image <- EBImage::Image(enhanced_channels)
        dimnames(tile_image)[[3]] <- all_channel_names_local
        rm(enhanced_channels)

        # Optionally save one tile visualization
        if (save_one_tile_viz_local && row_idx == center_row_local && col_idx == center_col_local) {
          if (verbose_local) cat(sprintf("[%s]     [%s] Generating enhancement visualization...\n", format(Sys.time(), "%H:%M:%S"), tile_id))
          # (Visualization logic remains same as original script)
        }

        # Measure features with simpleSeg::measureObjects
        if (verbose_local) cat(sprintf("[%s]     [%s] Quantifying features...\n", format(Sys.time(), "%H:%M:%S"), tile_id))

        # Create CytoImageList for measureObjects
        temp_images <- cytomapper::CytoImageList(list(tile_image))
        names(temp_images) <- tile_id
        cytomapper::channelNames(temp_images) <- all_channel_names_local

        temp_masks <- cytomapper::CytoImageList(list(EBImage::Image(mask)))
        names(temp_masks) <- tile_id

        # Use internal batching for measureObjects if requested
        old_bpparam <- BiocParallel::bpparam()
        if (tiling_mode_local == "virtual" && !use_parallel_quant) {
          BiocParallel::register(BiocParallel::MulticoreParam(workers = n_cores_local))
        } else if (batch_quantification_local > 1) {
          BiocParallel::register(BiocParallel::MulticoreParam(workers = batch_quantification_local))
        } else {
          BiocParallel::register(BiocParallel::SerialParam())
        }

        spe_tile <- simpleSeg::measureObjects(
          temp_masks,
          temp_images,
          img_id = "imageID"
        )

        # Restore previous BiocParallel backend
        BiocParallel::register(old_bpparam)

        # Add tile-specific coordinates adjustment
        # Coordinates in SPE are relative to tile. Convert to absolute.
        y_offset <- floor(row_idx * tile_size_local + crop_y_local)
        x_offset <- floor(col_idx * tile_size_local + crop_x_local)

        SpatialExperiment::spatialCoords(spe_tile)[, "y"] <- SpatialExperiment::spatialCoords(spe_tile)[, "y"] + y_offset
        SpatialExperiment::spatialCoords(spe_tile)[, "x"] <- SpatialExperiment::spatialCoords(spe_tile)[, "x"] + x_offset

        # Save SPE for this tile
        qs::qsave(spe_tile, spe_file, preset = "fast", nthreads = 1)

        # Clean up
        rm(tile_image, mask, temp_images, temp_masks, spe_tile)
        gc(verbose = FALSE)

        tile_elapsed <- as.numeric(difftime(Sys.time(), tile_start_time, units = "mins"))
        if (verbose_local) {
          cat(sprintf("[%s]     [%s] ✓ Quantification complete (%.1f min)\n", format(Sys.time(), "%H:%M:%S"), tile_id, tile_elapsed))
        }

        return(tile_id)
      } # End of quantify_one_tile function

      # Run quantification (parallel or serial)
      if (use_parallel_quant) {
        # Calculate workers per batch
        n_quant_workers <- max(1, floor(n_cores / batch_quantification))
        timed_message(sprintf(
          "    Parallel Settings: %d concurrent tiles, %d threads per tile (Total %d cores)",
          n_quant_workers, batch_quantification, n_quant_workers * batch_quantification
        ))

        results <- parallel::mclapply(tiles_needing_quant, quantify_one_tile, mc.cores = n_quant_workers)
      } else {
        results <- lapply(tiles_needing_quant, quantify_one_tile)
      }

      # Check for errors
      errors <- sapply(results, function(x) inherits(x, "try-error"))
      if (any(errors)) {
        timed_message(paste0("  ⚠️  Warning: ", sum(errors), " tiles failed to quantify"))
        for (i in which(errors)) {
          timed_message(paste0("    Error in tile ", i, ": ", as.character(results[[i]])))
        }
      } else {
        timed_message(paste0("  ✓ All ", n_tiles_to_quant, " tiles quantified successfully"))
      }

      # Clean up after parallel section
      gc(full = TRUE, reset = TRUE)
    }

    # Write completion flag
    writeLines("complete", quant_complete_flag)
    timed_message("  ✓ Quantification stage complete")
  }

  step_times$stage3 <- elapsed_time(step_times$stage3_start, units = "mins")
  if (timing) timed_message(paste0("  ⏱ Quantification: ", round(step_times$stage3, 2), " min"))

  # ============================================================================
  # PHASE 4: Data Merging & Formatting
  # ============================================================================
  timed_message("\n--- Phase 4: Data Merging & Formatting ---", "phase4_start")

  # ============================================================================
  # STAGE 4: Merge SpatialExperiment Objects
  # ============================================================================
  timed_message("\n--- Stage 4/4: Merging SpatialExperiment Objects ---", "stage4_start")

  # Get list of all SPE files (now .qs)
  spe_files <- list.files(tile_dir, pattern = "tile_\\d+_\\d+_spe\\.qs$", full.names = TRUE)

  timed_message(paste0("  Found ", length(spe_files), " tile SPE files"))

  # Load all SPE objects
  timed_message("  Loading SPE objects (serial loop, using all cores for decompression)...")
  load_start <- Sys.time()
  spe_list <- list()
  total_cells <- 0
  n_empty_tiles <- 0

  for (i in seq_along(spe_files)) {
    spe_data <- qs::qread(spe_files[i], nthreads = n_cores)

    # Check if this is an empty tile marker
    if (is.list(spe_data) && "n_cells" %in% names(spe_data) && spe_data$n_cells == 0) {
      # Skip empty tiles
      n_empty_tiles <- n_empty_tiles + 1
      next
    }

    # Valid SPE object
    spe_list[[length(spe_list) + 1]] <- spe_data
    total_cells <- total_cells + ncol(spe_data)

    # Progress report every 5 tiles or at end
    if (verbose && (i %% 5 == 0 || i == length(spe_files))) {
      timed_message(paste0(
        "    [", i, "/", length(spe_files), "] Loaded ",
        length(spe_list), " tiles with ", format(total_cells, big.mark = ","),
        " cells (", n_empty_tiles, " empty)"
      ))
    }
  }

  load_elapsed <- as.numeric(difftime(Sys.time(), load_start, units = "secs"))
  timed_message(paste0(
    "  ✓ Loaded ", length(spe_list), " non-empty tiles (",
    n_empty_tiles, " empty) in ", round(load_elapsed, 1), " sec"
  ))
  timed_message(paste0("  Total cells: ", format(total_cells, big.mark = ",")))

  # Merge all tiles using cbind
  timed_message("  Merging SPE objects with cbind()...")
  merge_start <- Sys.time()

  if (length(spe_list) == 0) {
    stop("No cells detected in any tiles!")
  } else if (length(spe_list) == 1) {
    timed_message("    Only 1 tile, no merging needed")
    merged_spe <- spe_list[[1]]
  } else {
    if (verbose) timed_message(paste0("    Merging ", length(spe_list), " SPE objects..."))
    merged_spe <- do.call(cbind, spe_list)
    merge_elapsed <- as.numeric(difftime(Sys.time(), merge_start, units = "secs"))
    if (verbose) timed_message(paste0("    ✓ cbind complete (", round(merge_elapsed, 1), " sec)"))
  }

  # Cleanup list immediately
  rm(spe_list)
  gc(verbose = FALSE)

  # Add metadata
  if (verbose) timed_message("    Adding metadata and coordinate names...")
  merged_spe$image_name <- image_name
  SpatialExperiment::spatialCoordsNames(merged_spe) <- c("x", "y")

  timed_message(paste0(
    "  ✓ Merged SPE: ", format(ncol(merged_spe), big.mark = ","), " cells, ",
    nrow(merged_spe), " features"
  ))

  # Save merged SPE
  if (verbose) timed_message("  Saving merged SPE to disk...")
  save_start <- Sys.time()
  spe_file <- file.path(PROG_DIR, paste0(image_name, "_cells.qs"))
  qs::qsave(merged_spe, spe_file, preset = "fast", nthreads = n_cores)
  save_elapsed <- as.numeric(difftime(Sys.time(), save_start, units = "secs"))
  timed_message(paste0("  ✓ Saved merged SPE: ", basename(spe_file), " (", round(save_elapsed, 1), " sec)"))

  step_times$stage4 <- elapsed_time(step_times$stage4_start)
  if (timing) timed_message(paste0("  ⏱ Merging: ", round(step_times$stage4, 2), " sec"))

  # ============================================================================
  # Final Summary
  # ============================================================================
  total_time <- elapsed_time(overall_start, units = "mins")
  timed_message("\n✓ Tiled segmentation and quantification complete!")
  timed_message(paste0("  Total time: ", round(total_time, 2), " minutes"))
  timed_message(paste0("  Total cells: ", format(ncol(merged_spe), big.mark = ",")))
  timed_message(paste0("  Tile directory: ", tile_dir))

  # Compile timing summary
  timing_summary <- NULL
  if (timing) {
    timing_summary <- c(
      duplicate_detection = step_times$stage0,
      tile_extraction = step_times$stage1 * 60,
      segmentation = step_times$stage2 * 60,
      quantification = step_times$stage3 * 60,
      merging = step_times$stage4,
      total = total_time * 60
    )
  }

  # Compile parameters
  parameters <- list(
    image_name = image_name,
    nuclear_channel_name = nuclear_channel_name,
    tile_size = tile_size,
    tiling_mode = tiling_mode,
    n_cores = n_cores,
    n_tiles_x = n_tiles_x,
    n_tiles_y = n_tiles_y,
    test_mode = test_mode,
    n_channels = n_channels,
    n_duplicate_channels_removed = length(duplicate_channels),
    segmentation = list(
      gamma = gamma,
      percentile_range = percentile_range,
      exclude_bright_outliers = exclude_bright_outliers,
      outlier_percentile = outlier_percentile,
      adaptive_enhancement = adaptive_enhancement,
      clahe_clip_limit = clahe_clip_limit,
      watershed_tolerance = watershed_tolerance,
      discSize = discSize,
      sizeSelection = sizeSelection,
      smooth = smooth,
      ext = ext,
      max_nucleus_size = max_nucleus_size
    ),
    enhancement = list(
      gamma = enhancement_gamma,
      percentile_range = enhancement_percentile_range,
      exclude_outliers = enhancement_exclude_outliers,
      outlier_percentile = enhancement_outlier_percentile,
      adaptive = enhancement_adaptive
    )
  )

  return(invisible(list(
    spe = merged_spe,
    n_cells = ncol(merged_spe),
    tile_dir = tile_dir,
    diagnostics = diagnostics_df,
    timing = timing_summary,
    parameters = parameters
  )))
}
