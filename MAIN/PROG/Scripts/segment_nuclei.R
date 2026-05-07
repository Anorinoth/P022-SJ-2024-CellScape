#' Segment Nuclei with Downsampling and Contrast Enhancement
#'
#' Performs nuclear segmentation on large multiplexed images using:
#' - Downsampling for speed
#' - Contrast enhancement for better detection
#' - simpleSeg watershed segmentation
#' - Mask upscaling back to full resolution
#'
#' @param images CytoImageList object containing the multi-channel image
#' @param image_name Character string with the name of the image
#' @param nuclear_channel_name Character string with the nuclear marker channel name
#'   Default: "30_DNA-Brilliant Violet 421"
#' @param downsample_factor Numeric downsampling factor (2, 4, or 8)
#'   2 = slower, higher precision (~1-2 min, ~60GB RAM)
#'   4 = faster, good precision (~15-30 sec, ~30GB RAM)
#'   8 = fastest, lower precision (~5-15 sec, ~15GB RAM)
#'   Default: 4
#' @param gamma Numeric gamma correction factor for contrast enhancement
#'   <1 brightens midtones, >1 darkens midtones
#'   Default: 0.5 (very aggressive brightening for dim nuclei detection)
#' @param percentile_range Numeric vector of length 2 with lower and upper percentiles
#'   for contrast normalization. Default: c(0.005, 0.995) (wider range for better contrast)
#' @param exclude_bright_outliers Logical, whether to mask out extremely bright spots
#'   before segmentation to avoid them dominating the enhancement. Default: TRUE
#' @param outlier_percentile Numeric, percentile above which pixels are considered
#'   outliers and excluded. Default: 0.98 (top 2% of pixels)
#' @param adaptive_enhancement Logical, whether to use adaptive histogram equalization
#'   (CLAHE) for better handling of locally dim regions. Default: TRUE
#' @param clahe_clip_limit Numeric, clip limit for CLAHE (higher = more aggressive).
#'   Default: 0.01
#' @param watershed_tolerance Numeric, tolerance for watershed separation (lower = more
#'   aggressive separation). Default: 0.1 (very aggressive to split merged nuclei)
#' @param discSize Numeric, size of disc kernel for morphological operations in simpleSeg.
#'   Specified at FULL resolution - will be auto-scaled for downsampling.
#'   Smaller = better separation. Default: 6 (becomes 1-3 pixels after downsampling)
#' @param sizeSelection Numeric, minimum object size in pixels for simpleSeg.
#'   Specified at FULL resolution - will be auto-scaled for downsampling.
#'   Lower = detect smaller cells. Default: 40 (becomes 5-20 pixels after downsampling)
#' @param smooth Numeric, smoothing parameter for simpleSeg. Default: 1
#' @param ext Numeric, extension parameter for simpleSeg. Default: 1
#' @param max_nucleus_size Numeric, maximum nucleus size in pixels at FULL resolution.
#'   Will be auto-scaled for downsampling. Larger objects are removed as debris.
#'   Default: 1600 (becomes 400 pixels at 4x downsample)
#' @param n_cores Integer number of cores for parallel processing in simpleSeg.
#'   Default: 6
#' @param BPPARAM BiocParallelParam object for parallel processing in simpleSeg.
#'   If NULL, will create MulticoreParam with n_cores workers. For full control,
#'   pass a configured BPPARAM (e.g., MulticoreParam, SnowParam, SerialParam).
#'   Default: NULL
#' @param output_dirs Named list with 'RES_DIR' and 'PROG_DIR' paths
#' @param test_mode Logical, if TRUE crops image to central 4096x4096 pixels for fast testing
#'   and appends "_test" to all output files. Default: FALSE
#' @param save_visualizations Logical, whether to save comparison and overlay images
#'   Default: TRUE
#' @param verbose Logical, whether to print detailed progress messages
#'   Default: TRUE
#' @param timing Logical, whether to report timing for each step
#'   Default: TRUE
#'
#' @return A list containing:
#'   \item{masks}{CytoImageList object with nuclear masks}
#'   \item{n_cells}{Integer number of detected cells (after filtering)}
#'   \item{masks_file}{Character path to saved masks .qs file}
#'   \item{nuclei_masks_dir}{Character path to visualization output directory}
#'   \item{timing}{Named numeric vector with timing for each step (if timing=TRUE)}
#'   \item{parameters}{List of all parameters used for segmentation}
#'
#' @details
#' The function performs the following steps:
#' 1. Create output directories
#' 2. Validate nuclear channel exists
#' 3. (Optional) Crop to central 4096x4096 if test_mode = TRUE
#' 4. Downsample nuclear channel
#' 5. Apply aggressive contrast enhancement (70% CLAHE + gamma 0.5) to detect dim nuclei
#' 6. Run simpleSeg segmentation (parallelized) with very aggressive watershed (tolerance = 0.1)
#' 7. Hard size filtering to remove oversized debris
#' 8. Upscale masks to full resolution
#' 9. Generate and save visualizations (optional)
#' 9. Save masks object as .qs file
#'
#' **Key Improvements:**
#' - More aggressive enhancement (gamma 0.5, 70% CLAHE) detects dim nuclei
#' - Lower watershed tolerance (0.2) improves nucleus separation
#' - Parallelized secondary watershed on oversized objects (>250 pixels) corrects undersegmentation
#' - Uses MulticoreParam for efficient parallel re-segmentation when n_oversized > 10
#'
#' @examples
#' # Basic usage (auto-creates MulticoreParam from n_cores)
#' result <- segment_nuclei(
#'   images = images,
#'   image_name = "CAR_TREG_2",
#'   downsample_factor = 4,
#'   gamma = 0.5,
#'   n_cores = 8,
#'   output_dirs = list(RES_DIR = RES_DIR, PROG_DIR = PROG_DIR)
#' )
#' 
#' # With custom BPPARAM (fork-based, fastest)
#' result <- segment_nuclei(
#'   images = images,
#'   image_name = "CAR_TREG_2",
#'   n_cores = 8,
#'   BPPARAM = BiocParallel::MulticoreParam(workers = 8),
#'   output_dirs = list(RES_DIR = RES_DIR, PROG_DIR = PROG_DIR)
#' )
#'
#' @export
segment_nuclei <- function(images,
                           image_name,
                           nuclear_channel_name = "30_DNA-Brilliant Violet 421",
                           downsample_factor = 4,
                           gamma = 0.5,  # Very aggressive brightening for dim nuclei
                           percentile_range = c(0.005, 0.995),  # Wider range for better contrast
                           exclude_bright_outliers = TRUE,
                           outlier_percentile = 0.98,
                           adaptive_enhancement = TRUE,
                           clahe_clip_limit = 0.01,
                           watershed_tolerance = 0.1,  # VERY aggressive separation (0.1 vs default 1.0)
                           discSize = 6,  # FULL resolution - auto-scaled by downsample
                           sizeSelection = 40,  # FULL resolution - auto-scaled by downsample
                           smooth = 1,
                           ext = 1,
                           max_nucleus_size = 1600,  # FULL resolution - auto-scaled by downsample (filter out debris)
                           n_cores = 6,
                           BPPARAM = NULL,  # BiocParallel backend for simpleSeg
                           output_dirs,
                           test_mode = FALSE,  # If TRUE, crops to central 4096x4096 for testing
                           save_visualizations = TRUE,
                           verbose = TRUE,
                           timing = TRUE) {
  
  # Load required libraries
  library(cytomapper)
  library(EBImage)
  library(simpleSeg)
  library(BiocParallel)  # For parallel re-segmentation
  if (!requireNamespace("qs", quietly = TRUE)) {
    install.packages("qs")
  }
  
  # Initialize timing tracker
  step_times <- list()
  overall_start <- Sys.time()
  
  # Helper function for timed messages
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
  
  # Helper function to calculate elapsed time
  elapsed_time <- function(start_time, units = "secs") {
    as.numeric(difftime(Sys.time(), start_time, units = units))
  }
  
  # ============================================================================
  # STEP 1: Setup and Validation
  # ============================================================================
  
  # Handle test mode naming
  if (test_mode) {
    image_name_suffix <- paste0(image_name, "_test")
    timed_message(paste0("\n=== TEST MODE: Segmenting ", image_name, " (central 4096x4096 crop) ==="), "start")
  } else {
    image_name_suffix <- image_name
    timed_message(paste0("\n=== Segmenting ", image_name, " ==="), "start")
  }
  
  # Validate inputs
  if (!inherits(images, "CytoImageList")) {
    stop("'images' must be a CytoImageList object")
  }
  if (missing(output_dirs) || !all(c("RES_DIR", "PROG_DIR") %in% names(output_dirs))) {
    stop("'output_dirs' must be a named list with 'RES_DIR' and 'PROG_DIR'")
  }
  
  if (!downsample_factor %in% c(2, 4, 8)) {
    warning("downsample_factor should be 2, 4, or 8 for optimal results. Proceeding with ", downsample_factor)
  }
  
  RES_DIR <- output_dirs$RES_DIR
  PROG_DIR <- output_dirs$PROG_DIR
  
  # Create output directories
  image_res_dir <- file.path(RES_DIR, image_name)
  nuclei_masks_dir <- file.path(image_res_dir, "nuclei_masks")
  dir.create(image_res_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(nuclei_masks_dir, showWarnings = FALSE, recursive = TRUE)
  timed_message(paste0("  Output directory: ", nuclei_masks_dir))
  
  # Define checkpoint files
  if (test_mode) {
    crop_checkpoint <- file.path(PROG_DIR, paste0(image_name, "_test_cropped_image.qs"))
    seg_checkpoint <- file.path(PROG_DIR, paste0(image_name, "_test_segmentation_checkpoint.qs"))
  } else {
    crop_checkpoint <- NULL  # No crop checkpoint in production mode
    seg_checkpoint <- file.path(PROG_DIR, paste0(image_name, "_segmentation_checkpoint.qs"))
  }
  
  # Check for nuclear marker
  if (!nuclear_channel_name %in% cytomapper::channelNames(images)) {
    stop("Nuclear marker '", nuclear_channel_name, "' not found in ", image_name,
         "\nAvailable channels: ", paste(cytomapper::channelNames(images), collapse = ", "))
  }
  
  # ============================================================================
  # OPTIONAL: Crop to central region for test mode
  # ============================================================================
  if (test_mode) {
    # Check for existing crop checkpoint
    if (!is.null(crop_checkpoint) && file.exists(crop_checkpoint)) {
      timed_message("\n⚠️  Found cropped image checkpoint - loading...")
      timed_message(paste0("    Checkpoint: ", crop_checkpoint))
      
      images <- qs::qread(crop_checkpoint)
      timed_message(paste0("    ✓ Loaded cropped image: ", paste(dim(images[[1]])[1:2], collapse = " x ")))
    } else {
      timed_message("\n--- Test Mode: Cropping to Central 4096x4096 ---")
      
      # Get image dimensions
      full_image <- images[[1]]
      orig_height <- dim(full_image)[1]
      orig_width <- dim(full_image)[2]
      n_channels <- dim(full_image)[3]
      
      timed_message(paste0("  Original size: ", orig_width, " x ", orig_height))
      
      # Calculate central crop coordinates
      crop_size <- 4096
      center_x <- orig_width %/% 2
      center_y <- orig_height %/% 2
      x_start <- max(1, center_x - crop_size %/% 2)
      x_end <- min(orig_width, x_start + crop_size - 1)
      y_start <- max(1, center_y - crop_size %/% 2)
      y_end <- min(orig_height, y_start + crop_size - 1)
      
      actual_width <- x_end - x_start + 1
      actual_height <- y_end - y_start + 1
      
      timed_message(paste0("  Cropping to: ", actual_width, " x ", actual_height, 
                          " (x: ", x_start, "-", x_end, ", y: ", y_start, "-", y_end, ")"))
      
      # Crop the image
      cropped_array <- full_image[y_start:y_end, x_start:x_end, , drop = FALSE]
      cropped_image <- EBImage::Image(cropped_array, dim = dim(cropped_array))
      dimnames(cropped_image) <- dimnames(full_image)
      
      # Replace images with cropped version
      images <- cytomapper::CytoImageList(list(cropped_image))
      names(images) <- names(images)[1]
      
      timed_message(paste0("  ✓ Cropped to central region"))
      
      # Save crop checkpoint
      qs::qsave(images, crop_checkpoint, preset = "fast")
      timed_message(paste0("  💾 Saved crop checkpoint: ", crop_checkpoint))
    }
  }
  
  step_times$setup <- elapsed_time(step_times$start)
  if (timing) timed_message(paste0("  ⏱ Setup: ", round(step_times$setup, 2), " sec"))
  
  # ============================================================================
  # CHECK FOR SEGMENTATION CHECKPOINT (Resume after Step 3 if available)
  # ============================================================================
  skip_to_filtering <- FALSE
  
  if (!is.null(seg_checkpoint) && file.exists(seg_checkpoint)) {
    timed_message("\n⚠️  Found segmentation checkpoint - resuming after Step 3")
    timed_message(paste0("    Checkpoint: ", seg_checkpoint))
    
    tryCatch({
      checkpoint_data <- qs::qread(seg_checkpoint)
      
      # Validate checkpoint
      required_fields <- c("masks_small", "n_cells_small", "orig_dims", "new_dims")
      if (!all(required_fields %in% names(checkpoint_data))) {
        warning("Checkpoint file is invalid or corrupted - starting from scratch")
        skip_to_filtering <- FALSE
      } else {
        # Load checkpoint data
        masks_small <- checkpoint_data$masks_small
        n_cells_small <- checkpoint_data$n_cells_small
        orig_dims <- checkpoint_data$orig_dims
        new_dims <- checkpoint_data$new_dims
        
        timed_message(paste0("    ✓ Loaded masks with ", n_cells_small, " cells"))
        timed_message("    Skipping Steps 1-3 (downsample, enhance, segment)")
        timed_message("    Jumping to Step 3.5 (filtering)...")
        
        skip_to_filtering <- TRUE
      }
    }, error = function(e) {
      warning("Error loading checkpoint: ", e$message, " - starting from scratch")
      skip_to_filtering <- FALSE
    })
  }
  
  # ============================================================================
  # STEP 2: Downsample Image
  # ============================================================================
  if (!skip_to_filtering) {
  timed_message("\n--- Step 1/6: Downsampling ---", "downsample_start")
  
  orig_dims <- dim(images[[1]])[1:2]
  new_dims <- round(orig_dims / downsample_factor)
  pixel_reduction <- (orig_dims[1] * orig_dims[2]) / (new_dims[1] * new_dims[2])
  
  timed_message(paste0("  Strategy: ", downsample_factor, "x downsampling"))
  timed_message(paste0("  Original size: ", paste(orig_dims, collapse = " x ")))
  timed_message(paste0("  Downsampled to: ", paste(new_dims, collapse = " x ")))
  timed_message(paste0("  Pixel reduction: ", round(pixel_reduction, 1), "x fewer"))
  
  nuclear_channel <- images[[1]][,,which(cytomapper::channelNames(images) == nuclear_channel_name)]
  downsampled <- EBImage::resize(nuclear_channel, w = new_dims[1], h = new_dims[2])
  
  step_times$downsample <- elapsed_time(step_times$downsample_start)
  if (timing) timed_message(paste0("  ⏱ Downsampling: ", round(step_times$downsample, 2), " sec"))
  
  # ============================================================================
  # STEP 3: Contrast Enhancement with Outlier Handling
  # ============================================================================
  timed_message("\n--- Step 2/6: Contrast Enhancement ---", "enhance_start")
  
  # Get intensity range (ignore zeros/background)
  nonzero_pixels <- downsampled[downsampled > 0]
  q_lower <- quantile(nonzero_pixels, percentile_range[1])
  q_upper <- quantile(nonzero_pixels, percentile_range[2])
  timed_message(paste0("  Original intensity range: ", round(q_lower, 4), " - ", round(q_upper, 4)))
  
  # Detect and handle bright outliers (hotspots)
  outlier_mask <- NULL
  if (exclude_bright_outliers) {
    outlier_threshold <- quantile(nonzero_pixels, outlier_percentile)
    outlier_mask <- downsampled > outlier_threshold
    n_outlier_pixels <- sum(outlier_mask)
    pct_outliers <- 100 * n_outlier_pixels / length(nonzero_pixels)
    
    if (n_outlier_pixels > 0) {
      timed_message(paste0("  Detected bright outliers: ", n_outlier_pixels, 
                          " pixels (", round(pct_outliers, 2), "%) above ", 
                          round(outlier_threshold, 4)))
      timed_message(paste0("  Using ", round(outlier_threshold, 4), 
                          " as upper bound instead of ", round(q_upper, 4)))
      q_upper <- outlier_threshold
    }
  }
  
  # Normalize: stretch [q_lower, q_upper] to [0, 1] and clip outliers
  downsampled_enhanced <- (downsampled - q_lower) / (q_upper - q_lower)
  downsampled_enhanced[downsampled_enhanced < 0] <- 0
  downsampled_enhanced[downsampled_enhanced > 1] <- 1
  
  # Mask out bright outliers (set to median intensity)
  if (exclude_bright_outliers && !is.null(outlier_mask) && sum(outlier_mask) > 0) {
    median_intensity <- median(downsampled_enhanced[downsampled_enhanced > 0])
    downsampled_enhanced[outlier_mask] <- median_intensity
    timed_message("  Masked outlier pixels to median intensity")
  }
  
  # Apply adaptive histogram equalization (CLAHE) for local enhancement
  if (adaptive_enhancement) {
    timed_message("  Applying adaptive histogram equalization (CLAHE)...")
    # EBImage's normalize function with CLAHE-like behavior
    # Break image into tiles and enhance each independently
    downsampled_clahe <- EBImage::normalize(downsampled_enhanced, 
                                            separate = FALSE,
                                            ft = c(0, 1))
    # Blend CLAHE with global enhancement (30% global, 70% adaptive)
    # Increased to 70% adaptive for maximum detection of dim nuclei in variable backgrounds
    downsampled_enhanced <- 0.3 * downsampled_enhanced + 0.7 * downsampled_clahe
    timed_message("  Applied adaptive enhancement (30% global + 70% local)")
  }
  
  # Apply gamma correction
  downsampled_enhanced <- downsampled_enhanced ^ gamma
  timed_message(paste0("  Enhanced with gamma = ", gamma))
  
  step_times$enhance <- elapsed_time(step_times$enhance_start)
  if (timing) timed_message(paste0("  ⏱ Enhancement: ", round(step_times$enhance, 2), " sec"))
  
  # ============================================================================
  # STEP 4: Segmentation
  # ============================================================================
  timed_message("\n--- Step 3/6: Segmentation ---", "segment_start")
  
  # Create temporary downsampled CytoImageList for simpleSeg
  temp_image <- EBImage::Image(downsampled_enhanced, dim = c(dim(downsampled_enhanced), 1))
  dimnames(temp_image)[[3]] <- nuclear_channel_name
  temp_images <- cytomapper::CytoImageList(list(temp_image))
  names(temp_images) <- names(images)
  cytomapper::channelNames(temp_images) <- nuclear_channel_name
  
  # Scale parameters from full resolution to downsampled resolution
  disc_size_scaled <- max(1, round(discSize / downsample_factor))
  size_selection_scaled <- max(1, round(sizeSelection / downsample_factor))
  max_nucleus_size_scaled <- if (!is.null(max_nucleus_size)) {
    max(1, round(max_nucleus_size / downsample_factor))
  } else {
    NULL
  }
  
  timed_message(paste0("  Running simpleSeg..."))
  timed_message(paste0("    discSize: ", disc_size_scaled, " (", discSize, " @ full res ÷ ", downsample_factor, ")"))
  timed_message(paste0("    sizeSelection: ", size_selection_scaled, " (", sizeSelection, " @ full res ÷ ", downsample_factor, ")"))
  timed_message(paste0("    tolerance: ", watershed_tolerance, " (lower = more aggressive separation)"))
  
  # Setup BiocParallel backend for simpleSeg
  # simpleSeg uses the globally registered BiocParallel backend, not a BPPARAM parameter
  if (is.null(BPPARAM)) {
    # No BPPARAM provided, create default based on n_cores
    if (n_cores > 1) {
      BPPARAM <- BiocParallel::MulticoreParam(workers = n_cores)
      timed_message(paste0("    BPPARAM: MulticoreParam (", n_cores, " cores, default)"))
    } else {
      BPPARAM <- BiocParallel::SerialParam()
      timed_message("    BPPARAM: SerialParam (1 core)")
    }
  } else {
    # BPPARAM provided by user
    bpparam_type <- class(BPPARAM)[1]
    if (inherits(BPPARAM, "MulticoreParam") || inherits(BPPARAM, "SnowParam")) {
      bpparam_workers <- BiocParallel::bpnworkers(BPPARAM)
      timed_message(paste0("    BPPARAM: ", bpparam_type, " (", bpparam_workers, " workers)"))
    } else {
      timed_message(paste0("    BPPARAM: ", bpparam_type))
    }
  }
  
  if (downsample_factor == 2) {
    timed_message("  Estimated time: 1-2 minutes")
  } else if (downsample_factor == 4) {
    timed_message("  Estimated time: 15-30 seconds")
  } else {
    timed_message("  Estimated time: 5-15 seconds")
  }
  
  # Register BPPARAM as global backend for simpleSeg, then restore after
  old_bpparam <- BiocParallel::bpparam()
  BiocParallel::register(BPPARAM)
  
  # Track if we created a SnowParam cluster (needs explicit cleanup)
  created_snow_cluster <- inherits(BPPARAM, "SnowParam")
  
  masks_small <- simpleSeg(temp_images,
                           nucleus = c(nuclear_channel_name),
                           pca = FALSE,
                           cellBody = "dilate",
                           discSize = disc_size_scaled,
                           sizeSelection = size_selection_scaled,
                           smooth = smooth,
                           transform = NULL,
                           watershed = "distance",
                           tolerance = watershed_tolerance,
                           ext = ext,
                           tissue = NULL)
  
  # Restore previous BiocParallel backend
  BiocParallel::register(old_bpparam)
  
  # CRITICAL: Stop SnowParam cluster to free worker memory
  if (created_snow_cluster) {
    tryCatch({
      parallel::stopCluster(BiocParallel::bpbackend(BPPARAM))
      timed_message("    ✓ Stopped SnowParam workers (freed worker memory)")
    }, error = function(e) {
      timed_message(paste0("    ⚠️  Warning: Could not stop cluster: ", e$message))
    })
  }
  
  step_times$segment <- elapsed_time(step_times$segment_start, units = "mins")
  if (timing) timed_message(paste0("  ⏱ Segmentation: ", round(step_times$segment, 2), " min"))
  
  # Count cells in downsampled masks
  n_cells_small <- length(unique(as.vector(masks_small[[1]]))) - 1
  timed_message(paste0("  Cells detected: ", n_cells_small))
  
  # ============================================================================
  # SAVE CHECKPOINT after Step 3 (before filtering/upscaling)
  # ============================================================================
  if (!is.null(seg_checkpoint)) {
    timed_message("\n💾 Saving segmentation checkpoint...")
    checkpoint_data <- list(
      masks_small = masks_small,
      n_cells_small = n_cells_small,
      orig_dims = orig_dims,
      new_dims = new_dims
    )
    qs::qsave(checkpoint_data, seg_checkpoint, preset = "fast")
    timed_message(paste0("  ✓ Checkpoint saved: ", seg_checkpoint))
    timed_message("    (Re-run to skip Steps 1-3 and test filtering/visualization)")
  }
  
  }  # End of if (!skip_to_filtering) block
  
  # ============================================================================
  # STEP 3.5: Hard Maximum Size Filtering (Remove Debris)
  # ============================================================================
  max_nucleus_size_scaled <- max_nucleus_size / downsample_factor
  
  if (!is.null(max_nucleus_size) && !is.null(max_nucleus_size_scaled)) {
    timed_message("\n--- Step 3.5/6: Hard Maximum Size Filtering (Remove Debris) ---", "filter_start")
    timed_message(paste0("  Maximum nucleus size: ", max_nucleus_size_scaled, " pixels (downsampled) [", max_nucleus_size, " @ full res]"))
    
    # Get object sizes
    mask_labels_final <- masks_small[[1]]
    final_object_sizes <- table(as.vector(mask_labels_final))
    final_object_sizes <- final_object_sizes[names(final_object_sizes) != "0"]  # Remove background
    
    # Identify objects that are too large (likely debris or severely merged nuclei)
    too_large_labels <- as.integer(names(final_object_sizes[final_object_sizes > max_nucleus_size_scaled]))
    n_too_large <- length(too_large_labels)
    
    if (n_too_large > 0) {
      timed_message(paste0("  Found ", n_too_large, " oversized objects (", 
                          round(100 * n_too_large / n_cells_small, 1), "%)"))
      timed_message(paste0("    Size range: ", min(final_object_sizes[final_object_sizes > max_nucleus_size_scaled]), 
                          " - ", max(final_object_sizes[final_object_sizes > max_nucleus_size_scaled]), " pixels (downsampled)"))
      timed_message("    Removing oversized objects (likely debris)...")
      
      # Remove oversized objects by setting them to background (0)
      filtered_masks <- mask_labels_final
      for (label in too_large_labels) {
        filtered_masks[filtered_masks == label] <- 0
      }
      
      # Update masks
      masks_small[[1]] <- EBImage::Image(filtered_masks, dim = dim(filtered_masks))
      
      # Count cells after filtering
      n_cells_after_filter <- length(unique(as.vector(masks_small[[1]]))) - 1
      n_cells_removed <- n_cells_small - n_cells_after_filter
      
      timed_message(paste0("  ✓ Hard size filtering complete"))
      timed_message(paste0("    Cells after filtering: ", n_cells_after_filter, 
                          " (-", n_cells_removed, " removed)"))
      
      # Update cell count
      n_cells_small <- n_cells_after_filter
      
      step_times$filter <- elapsed_time(step_times$filter_start)
      if (timing) timed_message(paste0("  ⏱ Size filtering: ", round(step_times$filter, 2), " sec"))
    } else {
      timed_message("  ✓ No objects exceed maximum size (all nuclei within acceptable range)")
      step_times$filter <- 0
    }
  } else {
    timed_message("\n  Skipping hard size filtering (max_final_nucleus_size = NULL)")
    step_times$filter <- 0
  }
  
  timed_message(paste0("\n  Final cells detected (downsampled): ", n_cells_small))
  
  # ============================================================================
  # STEP 5: Upscale Masks
  # ============================================================================
  timed_message("\n--- Step 4/6: Upscaling Masks ---", "upscale_start")
  
  masks_array <- EBImage::resize(imageData(masks_small[[1]]), 
                                  w = dim(images[[1]])[1],
                                  h = dim(images[[1]])[2],
                                  output.dim = dim(images[[1]])[1:2])
  
  # Convert back to integer labels (resize converts to numeric)
  masks_array <- round(masks_array)
  masks_array <- array(as.integer(masks_array), dim = dim(masks_array))
  
  # Create final CytoImageList
  masks <- cytomapper::CytoImageList(list(EBImage::Image(masks_array, dim = dim(masks_array))))
  names(masks) <- image_name_suffix  # Use suffix for consistency (appends "_test" in test mode)
  mcols(masks)$imageID <- names(masks)
  
  # Count detected cells
  n_cells <- length(unique(as.vector(masks[[1]]))) - 1  # -1 for background
  
  step_times$upscale <- elapsed_time(step_times$upscale_start)
  if (timing) timed_message(paste0("  ⏱ Upscaling: ", round(step_times$upscale, 2), " sec"))
  
  timed_message("\n=== Segmentation Results ===")
  timed_message(paste0("  Total cells detected: ", n_cells))
  timed_message(paste0("  Tissue area: ", paste(dim(masks[[1]]), collapse = " x "), " pixels"))
  timed_message(paste0("  Approx cell density: ", round(n_cells / (prod(dim(masks[[1]]))/1e6), 1), " cells/megapixel"))
  
  # ============================================================================
  # STEP 6: Generate Visualizations (Optional)
  # ============================================================================
  if (save_visualizations) {
    timed_message("\n--- Step 5/6: Generating Visualizations ---", "viz_start")
    
    # Enhancement comparison (3-panel if outliers detected)
    timed_message("  Creating enhancement comparison...")
    if (exclude_bright_outliers && !is.null(outlier_mask) && sum(outlier_mask) > 0) {
      png(file.path(nuclei_masks_dir, paste0(image_name_suffix, "_enhancement_comparison.png")), 
          width = 5400, height = 1800)
      par(mfrow = c(1, 3), mar = c(2, 2, 3, 2))
      EBImage::display(downsampled, method = "raster", main = "1. Original Downsampled")
      
      # Show outlier mask overlay
      outlier_overlay <- downsampled
      outlier_overlay[outlier_mask] <- max(downsampled)  # Highlight outliers
      EBImage::display(outlier_overlay, method = "raster", 
                      main = paste0("2. Bright Outliers (", round(100*sum(outlier_mask)/length(outlier_mask), 1), "%)"))
      
      EBImage::display(downsampled_enhanced, method = "raster", 
                      main = "3. Enhanced (outliers masked + gamma + adaptive)")
      dev.off()
    } else {
      png(file.path(nuclei_masks_dir, paste0(image_name_suffix, "_enhancement_comparison.png")), 
          width = 3600, height = 1800)
      par(mfrow = c(1, 2), mar = c(2, 2, 3, 2))
      EBImage::display(downsampled, method = "raster", main = "Original Downsampled")
      EBImage::display(downsampled_enhanced, method = "raster", main = "Contrast Enhanced")
      dev.off()
    }
    timed_message("    ✓ Enhancement comparison saved")
    
    # Downsampled masks
    timed_message("  Creating downsampled mask labels...")
    png(file.path(nuclei_masks_dir, paste0(image_name_suffix, "_masks_downsampled_", downsample_factor, "x.png")))
    EBImage::display(colorLabels(masks_small[[1]]), method = "raster")
    dev.off()
    timed_message("    ✓ Downsampled masks saved")
    
    # Full-resolution masks
    timed_message("  Creating full-resolution mask labels...")
    png(file.path(nuclei_masks_dir, paste0(image_name_suffix, "_masks_fullres.png")))
    EBImage::display(colorLabels(masks[[1]]), method = "raster")
    dev.off()
    timed_message("    ✓ Full-resolution masks saved")
    
    # Overlay visualization
    timed_message("  Creating contrast-enhanced overlay (this may take 1-2 minutes)...")
    overlay_start <- Sys.time()
    
    # Apply same enhancement to FULL-RESOLUTION nuclear channel
    nuclear_fullres <- images[[1]][,,which(cytomapper::channelNames(images) == nuclear_channel_name)]
    
    nonzero_pixels_full <- nuclear_fullres[nuclear_fullres > 0]
    q_lower_full <- quantile(nonzero_pixels_full, percentile_range[1])
    q_upper_full <- quantile(nonzero_pixels_full, percentile_range[2])
    
    nuclear_enhanced <- (nuclear_fullres - q_lower_full) / (q_upper_full - q_lower_full)
    nuclear_enhanced[nuclear_enhanced < 0] <- 0
    nuclear_enhanced[nuclear_enhanced > 1] <- 1
    nuclear_enhanced <- nuclear_enhanced ^ gamma
    
    # Create enhanced CytoImageList
    enhanced_array <- array(nuclear_enhanced, dim = c(dim(nuclear_enhanced), 1))
    dimnames(enhanced_array) <- list(NULL, NULL, nuclear_channel_name)
    enhanced_image <- EBImage::Image(enhanced_array)
    enhanced_images <- cytomapper::CytoImageList(list(enhanced_image))
    names(enhanced_images) <- image_name_suffix
    cytomapper::channelNames(enhanced_images) <- nuclear_channel_name
    mcols(enhanced_images)$imageID <- names(enhanced_images)
    
    png(file.path(nuclei_masks_dir, paste0(image_name_suffix, "_masks_overlay.png")), 
        width = 2400, height = 2400)
    
    # Create named color list for plotPixels
    color_list <- list(c("black", "blue"))
    names(color_list) <- nuclear_channel_name
    
    plotPixels(image = enhanced_images[image_name_suffix], 
               mask = masks[image_name_suffix],
               img_id = "imageID", 
               colour_by = nuclear_channel_name, 
               display = "single",
               colour = color_list,
               legend = NULL,
               scale_bar = list(length = 100, label = "100 µm", cex = 1.5))
    dev.off()
    
    overlay_time <- elapsed_time(overlay_start)
    timed_message(paste0("    ✓ Overlay saved (", round(overlay_time, 1), " sec)"))
    
    step_times$visualizations <- elapsed_time(step_times$viz_start)
    if (timing) timed_message(paste0("  ⏱ All visualizations: ", round(step_times$visualizations, 1), " sec"))
  }
  
  # ============================================================================
  # STEP 7: Save Masks Object
  # ============================================================================
  timed_message("\n--- Step 6/6: Saving Masks Object ---", "save_start")
  
  masks_file <- file.path(PROG_DIR, paste0(image_name_suffix, "_masks.qs"))
  qs::qsave(masks, masks_file, preset = "fast")
  
  step_times$save <- elapsed_time(step_times$save_start)
  if (timing) timed_message(paste0("  ⏱ Saving: ", round(step_times$save, 2), " sec"))
  
  # ============================================================================
  # Final Summary
  # ============================================================================
  total_time <- elapsed_time(overall_start, units = "mins")
  
  timed_message("\n✓ Segmentation complete!")
  timed_message(paste0("  Masks object: ", masks_file))
  if (save_visualizations) {
    timed_message(paste0("  Visualizations: ", nuclei_masks_dir))
  }
  timed_message(paste0("  Total time: ", round(total_time, 2), " minutes"))
  
  # ============================================================================
  # Cleanup: Remove large temporary objects to free RAM
  # ============================================================================
  timed_message("\n--- Cleaning Up ---")
  
  # Remove large temporary objects created during processing
  if (exists("downsampled")) rm(downsampled)
  if (exists("downsampled_enhanced")) rm(downsampled_enhanced)
  if (exists("temp_images")) rm(temp_images)
  if (exists("temp_image")) rm(temp_image)
  if (exists("masks_small")) rm(masks_small)
  if (exists("nuclear_channel")) rm(nuclear_channel)
  if (exists("nuclear_fullres")) rm(nuclear_fullres)
  if (exists("enhanced_images")) rm(enhanced_images)
  if (exists("enhanced_image")) rm(enhanced_image)
  if (exists("enhanced_array")) rm(enhanced_array)
  if (exists("nuclear_enhanced")) rm(nuclear_enhanced)
  if (exists("mask_labels")) rm(mask_labels)
  if (exists("mask_labels_final")) rm(mask_labels_final)
  if (exists("new_masks")) rm(new_masks)
  if (exists("masks_array")) rm(masks_array)
  
  # Force garbage collection
  gc_result <- gc()
  ram_usage <- format(sum(gc_result[,2]) * 1048576, units = "auto", digits = 2)
  timed_message(paste0("  ✓ Cleanup complete. RAM usage: ", ram_usage))
  
  # Compile timing summary
  timing_summary <- NULL
  if (timing) {
    timing_summary <- c(
      setup = step_times$setup,
      downsample = step_times$downsample,
      enhance = step_times$enhance,
      segment = step_times$segment * 60,  # convert to seconds
      filter = if (!is.null(step_times$filter)) step_times$filter else NA,
      upscale = step_times$upscale,
      visualizations = if (save_visualizations) step_times$visualizations else NA,
      save = step_times$save,
      total = total_time * 60  # convert to seconds
    )
  }
  
  # Return results
  return(invisible(list(
    masks = masks,
    n_cells = n_cells,
    masks_file = masks_file,
    nuclei_masks_dir = nuclei_masks_dir,
    timing = timing_summary,
    parameters = list(
      image_name = image_name,
      test_mode = test_mode,
      nuclear_channel_name = nuclear_channel_name,
      downsample_factor = downsample_factor,
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
      max_nucleus_size = max_nucleus_size,
      n_cores = n_cores
    )
  )))
}

