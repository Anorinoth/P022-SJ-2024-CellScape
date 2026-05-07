#' Enhance Image Channels and Quantify Cell Features
#'
#' Applies contrast enhancement to image channels and performs cell feature
#' quantification. Also generates multi-channel overlay visualizations.
#'
#' @param images_file Character string with path to saved CytoImageList RDS/QS file (images)
#'   Images are only loaded if needed (e.g., for initial channel filtering/downsampling)
#' @param masks_file Character string with path to saved CytoImageList RDS/QS file (masks)
#'   Masks are only loaded when needed (for quantification in Step 4)
#' @param image_name Character string with the name of the image
#' @param enhancement_params List with enhancement parameters:
#'   \item{gamma}{Numeric gamma correction (default 0.65)}
#'   \item{percentile_range}{Numeric vector c(lower, upper) (default c(0.01, 0.99))}
#'   \item{exclude_bright_outliers}{Logical (default TRUE)}
#'   \item{outlier_percentile}{Numeric (default 0.98)}
#'   \item{adaptive_enhancement}{Logical (default TRUE)}
#' @param viz_channels Named character vector of channels to visualize with their colors
#'   Names = channel names, Values = colors
#' @param viz_order Character vector with channel order (bottom to top) for overlay
#' @param output_dirs Named list with 'RES_DIR' and 'PROG_DIR' paths
#' @param BPPARAM BiocParallel parameter object for parallel processing
#'   Use MulticoreParam (fork-based) for memory efficiency with large images
#' @param downsample_factor Numeric, downsample factor for RAM efficiency (e.g., 2 = half resolution)
#'   Enhancement is done on downsampled images, then upsampled back for quantification
#'   Set to 1 for no downsampling (default). Higher values = lower RAM usage.
#' @param cleanup_temp Logical, whether to remove temporary directories after completion (default TRUE)
#'   Set to FALSE to keep temp files for inspection or troubleshooting
#' @param comparison_channels Numeric or character specifying how many channels to compare
#'   Options: numeric (e.g., 5), "all" for all channels, "none" to skip comparison (default 5)
#' @param verbose Logical, whether to print progress messages
#' @param verbose_enhancement Logical, whether to print detailed per-channel enhancement diagnostics (default FALSE)
#'   Set to TRUE to see enhancement mode (minimal/gentle/smooth/normal) for each channel
#' @param generate_visualization Logical, whether to generate multi-channel overlay visualization (default TRUE)
#'   Set to FALSE to skip Step 3 and save time (~1-2 minutes)
#' @param export_tiffs Logical, whether to export single-channel TIFFs after enhancement (default TRUE)
#'   TIFFs are saved to RES_DIR/image_name/enhanced_tiffs/
#' @param timing Logical, whether to report timing
#'
#' @return A list containing:
#'   \item{spe}{SpatialExperiment object with cell features}
#'   \item{enhanced_images}{CytoImageList with enhanced channels (for downstream use)}
#'   \item{viz_file}{Path to saved visualization}
#'   \item{spe_file}{Path to saved SpatialExperiment}
#'   \item{tiff_dir}{Path to directory containing exported single-channel TIFFs (if export_tiffs = TRUE)}
#'   \item{timing}{Timing information}
#'
#' @export
enhance_and_quantify <- function(images_file,
                                  masks_file,
                                  image_name,
                                  enhancement_params = list(
                                    gamma = 0.65,
                                    percentile_range = c(0.01, 0.99),
                                    exclude_bright_outliers = TRUE,
                                    outlier_percentile = 0.98,
                                    adaptive_enhancement = TRUE
                                  ),
                                  viz_channels = c(
                                    "30_DNA-Brilliant Violet 421" = "blue",
                                    "02_Peanut Agglutinin-Fitc" = "red",
                                    "09_CD19-FITC" = "magenta",
                                    "04_CD4-Brilliant Violet_421" = "green",
                                    "31_CorrectedFoxP3-GFP" = "yellow"
                                  ),
                                  viz_order = c("30_DNA-Brilliant Violet 421",
                                               "02_Peanut Agglutinin-Fitc",
                                               "09_CD19-FITC",
                                               "04_CD4-Brilliant Violet_421",
                                               "31_CorrectedFoxP3-GFP"),
                                  output_dirs,
                                  BPPARAM = BiocParallel::SerialParam(),
                                  downsample_factor = 2,
                                  cleanup_temp = TRUE,
                                  comparison_channels = 5,
                                  verbose = TRUE,
                                  verbose_enhancement = FALSE,
                                  generate_visualization = TRUE,
                                 export_tiffs = TRUE,
                                  timing = TRUE) {
  
  # Check required packages are loaded (should be loaded in calling script)
  required_packages <- c("cytomapper", "EBImage", "SpatialExperiment", "BiocParallel", "qs")
  missing_packages <- required_packages[!sapply(required_packages, function(pkg) {
    requireNamespace(pkg, quietly = TRUE)
  })]
  
  if (length(missing_packages) > 0) {
    stop("Required packages not loaded: ", paste(missing_packages, collapse = ", "),
         "\nPlease load these packages before calling enhance_and_quantify()")
  }
  
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
  
  timed_message(paste0("\n=== Enhancing and Quantifying ", image_name, " ==="), "start")
  
  # Check R memory settings for large cell quantification
  expr_limit <- getOption("expressions")
  if (expr_limit < 100000) {
    timed_message("ℹ️  Expression limit is low (", expr_limit, "). Will increase to 500k during quantification.")
  }
  
  # Check if we're likely to hit protection stack overflow
  # Note: --max-ppsize CANNOT be set in .Renviron or within R - it's a startup flag only
  mem_limit <- Sys.getenv("R_MAX_VSIZE")
  if (mem_limit == "") {
    timed_message("⚠️  R_MAX_VSIZE not set in .Renviron (found it set to:", mem_limit, ")")
  } else {
    timed_message("✓ R_MAX_VSIZE set to:", mem_limit)
  }
  
  timed_message("\nℹ️  If 'protection stack overflow' error occurs:")
  timed_message("   1. Exit R completely")
  timed_message("   2. Restart R from terminal with: R --max-ppsize=500000")
  timed_message("      (Valid range: typically 10000-1000000, default ~50000)")
  timed_message("   3. Or in RStudio: Tools → Global Options → Advanced")
  timed_message("   Note: --max-ppsize CANNOT be set in .Renviron\n")
  
  # NOTE: Temporary processing directories (for resume capability)
  # - PROG_DIR/temp_downsampled_{image_name}/  : Downsampled channels
  # - PROG_DIR/temp_enhanced_{image_name}/     : Enhanced channels
  # - PROG_DIR/temp_upsampled_{image_name}/    : Upsampled channels
  # These are cleaned up at the end, but persist if processing is interrupted
  
  # Validate inputs
  if (missing(output_dirs) || !all(c("RES_DIR", "PROG_DIR") %in% names(output_dirs))) {
    stop("'output_dirs' must be a named list with 'RES_DIR' and 'PROG_DIR'")
  }
  
  RES_DIR <- output_dirs$RES_DIR
  PROG_DIR <- output_dirs$PROG_DIR
  
  # Create output directory
  image_res_dir <- file.path(RES_DIR, image_name)
  cell_features_dir <- file.path(image_res_dir, "cell_features")
  dir.create(cell_features_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Extract enhancement parameters
  gamma <- enhancement_params$gamma
  percentile_range <- enhancement_params$percentile_range
  exclude_bright_outliers <- enhancement_params$exclude_bright_outliers
  outlier_percentile <- enhancement_params$outlier_percentile
  adaptive_enhancement <- enhancement_params$adaptive_enhancement
  
  # ============================================================================
  # STEP 1: Define Enhancement Function
  # ============================================================================
  timed_message("\n--- Step 1/5: Setting Up Enhancement Pipeline ---", "setup_start")
  
  #' Enhance a single channel with the optimized pipeline
  enhance_channel <- function(channel_data, channel_name = NULL, verbose_inner = FALSE) {
    # Get intensity range (ignore zeros/background)
    nonzero_pixels <- channel_data[channel_data > 0]
    if (length(nonzero_pixels) == 0) {
      if (verbose_inner) message("    Warning: ", channel_name, " has no non-zero pixels")
      return(channel_data)
    }
    
    # ============================================================================
    # INTELLIGENT SIGNAL PATTERN DETECTION
    # ============================================================================
    # Distinguish between:
    # 1. Rare positive clusters (bright spots = real signal, not outliers)
    # 2. Sparse cellular markers (pinpricks = real signal, not noise)
    # 3. Uniform autofluorescence (low uniform = background, not signal)
    
    # Basic statistics
    signal_max <- max(nonzero_pixels)
    signal_min <- min(nonzero_pixels)
    signal_range <- signal_max - signal_min
    signal_mean <- mean(nonzero_pixels)
    signal_median <- median(nonzero_pixels)
    signal_sd <- sd(nonzero_pixels)
    
    # Quantiles for distribution analysis
    signal_q05 <- quantile(nonzero_pixels, 0.05)
    signal_q25 <- quantile(nonzero_pixels, 0.25)
    signal_q75 <- quantile(nonzero_pixels, 0.75)
    signal_q95 <- quantile(nonzero_pixels, 0.95)
    signal_q99 <- quantile(nonzero_pixels, 0.99)
    signal_q999 <- quantile(nonzero_pixels, 0.999)
    
    # ============================================================================
    # METRIC 1: Coefficient of Variation (CV) - detects uniform vs variable signal
    # ============================================================================
    # High CV = variable signal (good)
    # Low CV = uniform signal (likely autofluorescence)
    coeff_variation <- signal_sd / signal_mean
    is_uniform <- coeff_variation < 0.3  # CV < 0.3 indicates very uniform signal
    
    # ============================================================================
    # METRIC 2: Signal Sparsity - what % of pixels have true signal?
    # ============================================================================
    # Use multiple thresholds to detect different sparsity patterns
    sparse_threshold_conservative <- signal_median * 2  # Conservative (more sensitive)
    sparse_threshold_strict <- signal_mean + 2 * signal_sd  # Strict (less sensitive)
    
    n_pixels_conservative <- sum(channel_data > sparse_threshold_conservative)
    n_pixels_strict <- sum(channel_data > sparse_threshold_strict)
    
    sparsity_conservative <- n_pixels_conservative / length(channel_data)
    sparsity_strict <- n_pixels_strict / length(channel_data)
    
    # Very sparse: rare markers with pinpricks of signal
    is_very_sparse <- sparsity_strict < 0.001  # < 0.1% of pixels above mean+2SD
    
    # ============================================================================
    # METRIC 3: Dynamic Range Ratio - how much of the possible range is used?
    # ============================================================================
    # Low dynamic range = compressed signal (likely autofluorescence)
    # High dynamic range = good separation between background and signal
    dynamic_range_ratio <- signal_range / signal_max
    is_low_dynamic_range <- dynamic_range_ratio < 0.5  # Less than 50% of max
    
    # ============================================================================
    # METRIC 4: Tail Heaviness - are there rare bright outliers?
    # ============================================================================
    # Compare top 0.1% to top 5% to detect rare bright clusters
    tail_ratio <- signal_q999 / signal_q95
    has_bright_clusters <- tail_ratio > 2  # Top 0.1% is >2x brighter than top 5%
    
    # ============================================================================
    # METRIC 5: Signal-to-Noise Ratio
    # ============================================================================
    # Use lower quartile as background, upper quartile as signal
    background_level <- signal_q25
    signal_level <- signal_q95
    background_noise <- sd(nonzero_pixels[nonzero_pixels < signal_median])
    if (is.na(background_noise) || background_noise == 0) background_noise <- 1
    
    signal_to_noise <- (signal_level - background_level) / background_noise
    is_low_snr <- signal_to_noise < 2
    
    # ============================================================================
    # PATTERN CLASSIFICATION
    # ============================================================================
    # CRITICAL: Check patterns in order from most specific to least specific
    # Priority: Autofluorescence > Sparse Cellular > Rare Cluster > Normal
    
    # NOTE: Data is normalized to [0,1] range by CytoImageList
    
    # DEBUG: Add diagnostic info if verbose
    if (verbose_inner) {
      message("    ", channel_name, ": DIAGNOSTICS")
      message("      • max=", round(signal_max, 4), ", range=", round(signal_range, 4), 
              ", median=", round(signal_median, 4))
      message("      • CV=", round(coeff_variation, 3), ", dynamic_ratio=", round(dynamic_range_ratio, 3))
      message("      • sparsity(strict)=", round(sparsity_strict * 100, 3), 
              "%, tail_ratio=", round(tail_ratio, 2))
      message("      • SNR=", round(signal_to_noise, 2), ", uniform=", is_uniform, 
              ", low_dyn_range=", is_low_dynamic_range)
    }
    
    # PATTERN 1: Uniform Autofluorescence (PDL1, CTLA-4)
    # Characteristics:
    # - Low CV (< 1.5) - very little variation
    # - High sparsity (> 15%) - "signal" is widespread everywhere
    # - Low tail ratio (< 4) - no real bright spots
    # Key: Widespread + uniform + no bright spots = background autofluorescence
    is_autofluorescence <- coeff_variation < 1.5 &&  # Very uniform
                          sparsity_conservative > 0.15 &&  # At least 15% of pixels
                          tail_ratio < 4  # No bright clusters
    
    # PATTERN 2: Sparse Cellular Signal - PINPRICKS (CXCR5, EGFR, Granzyme B, IL4, IL-21, FoxP3)
    # Characteristics:
    # - VERY sparse (< 2% strict threshold) - only a few positive cells
    # - Has variation (CV > 1.5 OR tail_ratio > 2) - signal exists when present
    # - Good dynamic range
    # Key: VERY LOW sparsity + some signal variation = sparse cellular markers
    # NOTE: Check this BEFORE rare cluster to catch moderate tail ratios (5-15)
    is_sparse_cellular <- sparsity_strict < 0.02 &&  # < 2% of pixels above mean+2SD
                         (coeff_variation > 1.5 || tail_ratio > 2) &&  # Some signal variation
                         !is_low_dynamic_range  # Good range
    
    # PATTERN 3: Rare Positive Clusters - BRIGHT REGIONS (Peanut Agglutinin)
    # Characteristics:
    # - VERY high tail ratio (>= 15x) - a few pixels are MASSIVELY brighter
    # - Moderate sparsity (> 1%) - larger areas, not just pinpricks
    # - High absolute max (> 0.8) - strong bright signal
    # Key: EXTREME tail + moderate coverage = large bright clusters
    is_rare_cluster <- tail_ratio >= 15 &&  # Extreme tail (>= 15x)
                      sparsity_strict > 0.01 &&  # More than 1% (bigger than pinpricks)
                      signal_max > 0.8 &&  # Very strong signal
                      !is_low_dynamic_range  # Good separation
    
    # ============================================================================
    # DETERMINE ENHANCEMENT STRATEGY
    # ============================================================================
    
    if (is_autofluorescence) {
      enhancement_mode <- "autofluorescence"
      disable_outlier_removal <- TRUE  # Keep all signal
      use_strong_clahe <- TRUE  # Gentle CLAHE only
      gamma_adjustment <- 0.35  # VERY conservative gamma (0.65 + 0.35 = 1.0 = no gamma)
      smoothing_sigma <- 0  # No smoothing
      
      if (verbose_inner) {
        message("    ", channel_name, ": AUTOFLUORESCENCE mode")
        message("      → Uniform signal (CV=", round(coeff_variation, 2), 
                "), sparsity=", round(sparsity_conservative * 100, 1), "%")
        message("      → Strategy: NO gamma correction (mostly background)")
      }
      
    } else if (is_rare_cluster) {
      enhancement_mode <- "rare_cluster"
      disable_outlier_removal <- TRUE  # DON'T remove bright spots - they're real signal!
      use_strong_clahe <- FALSE  # Gentle enhancement
      gamma_adjustment <- 0  # Normal gamma
      smoothing_sigma <- 0  # NO smoothing - preserve bright clusters!
      
      if (verbose_inner) {
        message("    ", channel_name, ": RARE CLUSTER mode")
        message("      → Bright clusters detected (tail ratio=", round(tail_ratio, 2), 
                "), sparse=", round(sparsity_strict * 100, 3), "%")
        message("      → Strategy: Preserve bright spots, NO outlier removal, NO smoothing")
      }
      
    } else if (is_sparse_cellular) {
      enhancement_mode <- "sparse_cellular"
      disable_outlier_removal <- TRUE  # Keep the sparse signal!
      use_strong_clahe <- FALSE  # Minimal CLAHE
      gamma_adjustment <- 0.35  # NO gamma correction (0.65 + 0.35 = 1.0)
      smoothing_sigma <- 2.0  # Strong smoothing to suppress background noise
      
      if (verbose_inner) {
        message("    ", channel_name, ": SPARSE CELLULAR mode")
        message("      → Very sparse signal (", round(sparsity_strict * 100, 3), 
                "%), tail=", round(tail_ratio, 2))
        message("      → Strategy: NO gamma, strong smoothing to suppress background")
      }
      
    } else {
      # Normal channel with good signal distribution
      enhancement_mode <- "normal"
      disable_outlier_removal <- FALSE
      use_strong_clahe <- FALSE
      gamma_adjustment <- 0
      smoothing_sigma <- 0
      
      if (verbose_inner) {
        message("    ", channel_name, ": NORMAL mode")
        message("      → Max=", round(signal_max, 4), ", range=", round(signal_range, 4))
        message("      → Sparsity=", round(sparsity_conservative * 100, 2), 
                "%, tail_ratio=", round(tail_ratio, 2))
        message("      → CV=", round(coeff_variation, 3), ", dynamic_ratio=", round(dynamic_range_ratio, 3))
      }
    }
    
    # ============================================================================
    # APPLY PRE-ENHANCEMENT SMOOTHING
    # ============================================================================
    
    if (smoothing_sigma > 0) {
      channel_data <- EBImage::gblur(channel_data, sigma = smoothing_sigma)
      # Recalculate nonzero pixels after smoothing
      nonzero_pixels <- channel_data[channel_data > 0]
      if (length(nonzero_pixels) == 0) {
        return(channel_data)  # All smoothed to zero somehow
      }
    }
    
    # ============================================================================
    # NORMALIZATION WITH ADAPTIVE OUTLIER HANDLING
    # ============================================================================
    
    q_lower <- quantile(nonzero_pixels, percentile_range[1])
    q_upper <- quantile(nonzero_pixels, percentile_range[2])
    
    # Handle bright outliers ONLY if not disabled by pattern detection
    outlier_mask <- NULL
    if (exclude_bright_outliers && !disable_outlier_removal) {
      outlier_threshold <- quantile(nonzero_pixels, outlier_percentile)
      outlier_mask <- channel_data > outlier_threshold
      
      if (sum(outlier_mask) > 0) {
        if (verbose_inner) {
          pct <- 100 * sum(outlier_mask) / length(nonzero_pixels)
          message("      → Masking ", round(pct, 2), "% outliers (threshold=", 
                  round(outlier_threshold, 1), ")")
        }
        q_upper <- outlier_threshold
      }
    } else if (disable_outlier_removal && verbose_inner) {
      message("      → Outlier removal DISABLED (preserving bright signal)")
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
    
    # ============================================================================
    # ADAPTIVE ENHANCEMENT (CLAHE)
    # ============================================================================
    
    if (adaptive_enhancement) {
      if (use_strong_clahe) {
        # Autofluorescence: MINIMAL CLAHE (mostly background, avoid false signal)
        enhanced_clahe <- EBImage::normalize(enhanced, separate = FALSE, ft = c(0, 1))
        enhanced <- 0.95 * enhanced + 0.05 * enhanced_clahe  # 5% CLAHE (minimal)
        if (verbose_inner) {
          message("      → Applying MINIMAL CLAHE (5% blend) - avoid false signal")
        }
      } else if (enhancement_mode == "sparse_cellular") {
        # Sparse cellular: MINIMAL CLAHE to avoid amplifying background
        enhanced_clahe <- EBImage::normalize(enhanced, separate = FALSE, ft = c(0, 1))
        enhanced <- 0.95 * enhanced + 0.05 * enhanced_clahe  # 5% CLAHE (minimal)
        if (verbose_inner) {
          message("      → Applying MINIMAL CLAHE (5% blend) - suppress background")
        }
      } else if (enhancement_mode == "rare_cluster") {
        # Rare cluster: Gentle CLAHE
        enhanced_clahe <- EBImage::normalize(enhanced, separate = FALSE, ft = c(0, 1))
        enhanced <- 0.85 * enhanced + 0.15 * enhanced_clahe  # 15% CLAHE
        if (verbose_inner) {
          message("      → Applying GENTLE CLAHE (15% blend) to preserve bright clusters")
        }
      } else {
        # Normal: Standard CLAHE
        enhanced_clahe <- EBImage::normalize(enhanced, separate = FALSE, ft = c(0, 1))
        enhanced <- 0.7 * enhanced + 0.3 * enhanced_clahe  # 30% CLAHE
      }
    }
    
    # ============================================================================
    # ADAPTIVE GAMMA CORRECTION
    # ============================================================================
    
    gamma_adjusted <- gamma + gamma_adjustment
    
    if (verbose_inner && gamma_adjustment != 0) {
      message("      → Gamma: ", round(gamma_adjusted, 2), " (adjusted from ", gamma, ")")
    }
    
    enhanced <- enhanced ^ gamma_adjusted
    
    return(enhanced)
  }
  
  timed_message("  ✓ Enhancement function ready")
  timed_message(paste0("    gamma = ", gamma))
  timed_message(paste0("    outlier exclusion = ", exclude_bright_outliers))
  timed_message(paste0("    adaptive enhancement = ", adaptive_enhancement))
  
  step_times$setup <- elapsed_time(step_times$setup_start)
  if (timing) timed_message(paste0("  ⏱ Setup: ", round(step_times$setup, 2), " sec"))
  
  # ============================================================================
  # STEP 0.5: Determine if we need to load images/masks (lazy loading)
  # ============================================================================
  timed_message("\n--- Step 0.5/5: Checking if images need to be loaded ---", "lazy_load_start")
  
  # Validate file paths
  if (!file.exists(images_file)) {
    # Provide helpful error message suggesting both formats
    base_name <- sub("\\.(qs|rds)$", "", images_file, ignore.case = TRUE)
    qs_path <- paste0(base_name, ".qs")
    rds_path <- paste0(base_name, ".rds")
    stop("Images file not found: ", images_file, "\n",
         "  Tried: ", qs_path, " (", ifelse(file.exists(qs_path), "EXISTS", "missing"), ")\n",
         "  Tried: ", rds_path, " (", ifelse(file.exists(rds_path), "EXISTS", "missing"), ")\n",
         "  Make sure you've saved the images using qs::qsave() or saveRDS()")
  }
  if (!file.exists(masks_file)) {
    # Provide helpful error message suggesting both formats
    base_name <- sub("\\.(qs|rds)$", "", masks_file, ignore.case = TRUE)
    qs_path <- paste0(base_name, ".qs")
    rds_path <- paste0(base_name, ".rds")
    stop("Masks file not found: ", masks_file, "\n",
         "  Tried: ", qs_path, " (", ifelse(file.exists(qs_path), "EXISTS", "missing"), ")\n",
         "  Tried: ", rds_path, " (", ifelse(file.exists(rds_path), "EXISTS", "missing"), ")\n",
         "  Make sure you've run segment_nuclei() first")
  }
  
  timed_message(paste0("  Images file: ", images_file, " (", 
                      ifelse(grepl("\\.qs$", images_file, ignore.case = TRUE), "QS", "RDS"), ")"))
  timed_message(paste0("  Masks file: ", masks_file, " (", 
                      ifelse(grepl("\\.qs$", masks_file, ignore.case = TRUE), "QS", "RDS"), ")"))
  
  # Check if we have existing temp files that would let us skip loading images
  temp_upsampled_check <- file.path(PROG_DIR, paste0("temp_upsampled_", image_name))
  temp_downsampled_check <- file.path(PROG_DIR, paste0("temp_downsampled_", image_name))
  temp_fullres_check <- file.path(PROG_DIR, paste0("temp_fullres_", image_name))
  
  existing_upsampled <- list.files(temp_upsampled_check, pattern = "^channel_\\d{2}\\.rds$", full.names = FALSE)
  existing_downsampled <- list.files(temp_downsampled_check, pattern = "^channel_\\d{2}\\.rds$", full.names = FALSE)
  existing_fullres <- list.files(temp_fullres_check, pattern = "^channel_\\d{2}\\.rds$", full.names = FALSE)
  
  has_upsampled_temps <- length(existing_upsampled) > 0
  has_downsampled_temps <- length(existing_downsampled) > 0
  has_fullres_temps <- length(existing_fullres) > 0
  
  # Determine which temp directory to use for metadata extraction (prefer earliest stage)
  metadata_temp_dir <- NULL
  metadata_temp_type <- NULL
  metadata_files <- NULL
  
  if (has_fullres_temps) {
    metadata_temp_dir <- temp_fullres_check
    metadata_temp_type <- "full-res"
    metadata_files <- existing_fullres
  } else if (has_downsampled_temps) {
    metadata_temp_dir <- temp_downsampled_check
    metadata_temp_type <- "downsampled"
    metadata_files <- existing_downsampled
  } else if (has_upsampled_temps) {
    metadata_temp_dir <- temp_upsampled_check
    metadata_temp_type <- "upsampled"
    metadata_files <- existing_upsampled
  }
  
  # We need images ONLY if no temp files exist at all
  need_images <- is.null(metadata_temp_dir)
  
  if (!need_images) {
    timed_message(paste0("  ✓ Found existing ", metadata_temp_type, " temp files (", 
                        length(metadata_files), " channels)"))
    timed_message("  → Images NOT needed (will extract metadata from temp files)")
  } else {
    timed_message("  → Images NEEDED (no temp files found, will process from scratch)")
  }
  
  # Load images or extract metadata from temp files
  images <- NULL
  all_channel_names <- NULL
  n_channels_original <- NULL
  original_dims <- NULL
  original_image_name <- NULL
  
  if (need_images) {
    # CASE 1: No temp files exist - must load full images
    timed_message("\n  Loading images from disk...")
    
    # Detect file format and load
    if (grepl("\\.qs$", images_file, ignore.case = TRUE)) {
      images <- qs::qread(images_file)
      timed_message("    ✓ Loaded from QS format")
    } else {
      images <- readRDS(images_file)
      timed_message("    ✓ Loaded from RDS format")
    }
    
    all_channel_names <- cytomapper::channelNames(images)
    n_channels_original <- length(all_channel_names)
    original_dims <- dim(images[[1]])
    original_image_name <- names(images)[1]
    
    image_ram <- prod(original_dims) * 8 / 1e9
    timed_message(paste0("    Channels: ", n_channels_original))
    timed_message(paste0("    Dimensions: ", paste(original_dims[1:2], collapse = " x ")))
    timed_message(paste0("    RAM: ~", round(image_ram, 1), " GB"))
    
    # Save channel metadata for future runs (tiny file, saves loading 8.3 GB later)
    metadata_file <- file.path(PROG_DIR, paste0(image_name, "_channel_metadata.rds"))
    channel_metadata <- list(
      channel_names = all_channel_names,
      n_channels = n_channels_original,
      dimensions = original_dims,
      image_name = original_image_name,
      created = Sys.time()
    )
    saveRDS(channel_metadata, metadata_file, compress = FALSE)
    timed_message(paste0("    ✓ Saved channel metadata: ", basename(metadata_file)))
    
  } else {
    # CASE 2: Temp files exist - extract metadata WITHOUT loading 8.3 GB images!
    timed_message("\n  Extracting metadata from temp files (avoiding 8.3 GB load)...")
    
    # Extract channel numbers from filenames (e.g., "channel_01.rds" → 1)
    channel_numbers <- as.integer(gsub("channel_(\\d{2})\\.rds", "\\1", metadata_files))
    n_channels_original <- length(channel_numbers)
    
    timed_message(paste0("    Found ", n_channels_original, " channels in ", metadata_temp_type, " temps"))
    
    # Load ONE channel to get dimensions
    first_channel_file <- file.path(metadata_temp_dir, metadata_files[1])
    first_channel <- readRDS(first_channel_file)
    temp_dims <- dim(first_channel)
    rm(first_channel)
    
    # Calculate original dimensions based on temp type
    if (metadata_temp_type == "downsampled") {
      # Upscale dimensions by downsample_factor
      original_dims <- c(
        round(temp_dims[1] * downsample_factor),
        round(temp_dims[2] * downsample_factor),
        n_channels_original
      )
      timed_message(paste0("    Temp dimensions: ", paste(temp_dims[1:2], collapse = " x ")))
      timed_message(paste0("    Original dimensions (×", downsample_factor, "): ", 
                          paste(original_dims[1:2], collapse = " x ")))
    } else {
      # Full-res or upsampled temps already have original dimensions
      original_dims <- c(temp_dims[1], temp_dims[2], n_channels_original)
      timed_message(paste0("    Dimensions: ", paste(original_dims[1:2], collapse = " x ")))
    }
    
    # Try to load channel names from metadata file (if it exists)
    metadata_file <- file.path(PROG_DIR, paste0(image_name, "_channel_metadata.rds"))
    
    if (file.exists(metadata_file)) {
      # Load saved metadata
      channel_metadata <- readRDS(metadata_file)
      all_channel_names <- channel_metadata$channel_names
      original_image_name <- channel_metadata$image_name
      timed_message(paste0("    Channels: ", n_channels_original, " (names loaded from metadata file)"))
    } else {
      # No metadata file - create placeholder names
      # We'll save metadata for next time
      all_channel_names <- paste0("CH_", sprintf("%02d", channel_numbers))
      original_image_name <- image_name  # Use the passed parameter as fallback
      timed_message(paste0("    Channels: ", n_channels_original, " (using placeholder names)"))
      timed_message("    ⚠️  No channel metadata file found - will save for next run")
    }
    
    timed_message("    ✓ Metadata extracted WITHOUT loading full images (saved ~8.3 GB load!)")
  }
  
  step_times$lazy_load <- elapsed_time(step_times$lazy_load_start)
  if (timing) timed_message(paste0("  ⏱ Lazy load check: ", round(step_times$lazy_load, 2), " sec"))
  
  # ============================================================================
  # STEP 1: Exclude Duplicate Channels (ONLY when starting fresh)
  # ============================================================================
  # Identify duplicate channels by comparing names after removing numeric prefix
  # Expected: 31 original channels → 26 unique channels (5 duplicates removed)
  #
  # IMPORTANT: This ONLY runs when loading fresh images. When resuming from temp
  # files, we skip this entirely and process all existing temp files.
  
  timed_message("\n--- Step 1/5: Channel Selection ---", "channel_selection_start")
  timed_message(paste0("  Original channels: ", n_channels_original))
  
  # Use temp file info from lazy load check
  has_existing_temps <- !is.null(metadata_temp_dir)
  
  if (has_existing_temps) {
    # -------------------------------------------------------------------------
    # RESUMING FROM TEMP FILES: Skip duplicate detection entirely
    # -------------------------------------------------------------------------
    timed_message(paste0("  ℹ️  Found ", length(metadata_files), " existing temp files from previous run"))
    timed_message("  → Skipping duplicate detection (already filtered when temp files were created)")
    
    # Process all temp files as-is
    n_channels <- n_channels_original  # This is the count from temp files
    channels_to_keep <- 1:n_channels
    duplicate_channels <- integer(0)
    
    timed_message(paste0("  → Will process all ", n_channels, " channels from temp files"))
    
  } else {
    # -------------------------------------------------------------------------
    # STARTING FRESH: Detect and remove duplicates
    # -------------------------------------------------------------------------
    timed_message("  → No existing temp files found, starting fresh")
    timed_message("  Detecting duplicate channels by comparing names...")
    
    # Remove numeric prefix: "01_CD19-PerCP-Cy5_5" → "CD19-PerCP-Cy5_5"
    channel_names_no_prefix <- gsub("^\\d+_", "", all_channel_names)
    
    # Find duplicates: which names appear more than once?
    duplicate_names <- unique(channel_names_no_prefix[duplicated(channel_names_no_prefix)])
    
    if (length(duplicate_names) > 0) {
      # For each duplicate name, keep the FIRST occurrence, mark others as duplicates
      duplicate_channels <- integer(0)
      for (dup_name in duplicate_names) {
        indices_with_name <- which(channel_names_no_prefix == dup_name)
        if (length(indices_with_name) > 1) {
          # Keep first, mark rest as duplicates
          duplicate_channels <- c(duplicate_channels, indices_with_name[-1])
        }
      }
      
      timed_message(paste0("  Found ", length(duplicate_channels), " duplicate channels to exclude:"))
      for (dup_ch in sort(duplicate_channels)) {
        timed_message(paste0("    - Channel ", sprintf("%02d", dup_ch), ": ", all_channel_names[dup_ch], 
                            " (duplicate of ", channel_names_no_prefix[dup_ch], ")"))
      }
      
      # Keep all channels EXCEPT duplicates
      channels_to_keep <- setdiff(1:n_channels_original, duplicate_channels)
      n_channels <- length(channels_to_keep)
      
      timed_message(paste0("  Retained ", n_channels, " unique channels (removed ", 
                          length(duplicate_channels), " duplicates)"))
      
    } else {
      timed_message("  → No duplicate channel names detected")
      duplicate_channels <- integer(0)
      channels_to_keep <- 1:n_channels_original
      n_channels <- n_channels_original
    }
  }
  
  # Calculate RAM (original_dims already set in lazy load section)
  ram_per_channel <- prod(original_dims[1:2]) * 8 / 1e9
  original_ram <- n_channels_original * ram_per_channel
  reduced_ram <- n_channels * ram_per_channel
  
  if (length(duplicate_channels) > 0) {
    timed_message(paste0("  Retained channels: ", n_channels, " (saved ~", round(original_ram - reduced_ram, 1), " GB RAM)"))
    timed_message(paste0("  RAM: ", round(original_ram, 1), " GB → ", round(reduced_ram, 1), " GB"))
  } else {
    timed_message(paste0("  Using all ", n_channels, " channels (no filtering)"))
    timed_message(paste0("  RAM: ", round(reduced_ram, 1), " GB"))
  }
  
  # Store original image name (if images were loaded, otherwise already set in lazy load section)
  if (!is.null(images)) {
    original_image_name <- names(images)[1]
  }
  
  # Only perform channel extraction if we're actually filtering
  if (length(duplicate_channels) > 0 && !has_existing_temps) {
    # Extract channels WITHOUT creating intermediate copies (in-place operation)
    # NOTE: We keep channels_to_keep for mapping original→filtered indices later
    timed_message("  Extracting non-duplicate channels...")
    
    # Get the original array
    original_full_array <- images[[1]]
    
    # Create filtered array directly
    filtered_array <- original_full_array[,,channels_to_keep, drop = FALSE]
    
    # Update channel names to reflect the filtered set
    filtered_channel_names <- all_channel_names[channels_to_keep]
    dimnames(filtered_array)[[3]] <- filtered_channel_names
    
    # Free original array IMMEDIATELY before creating new image
    rm(original_full_array, images)
    gc()
    timed_message("    ✓ Freed original array")
    
    # Now create the filtered image
    filtered_image <- EBImage::Image(filtered_array)
    images <- cytomapper::CytoImageList(list(filtered_image))
    names(images) <- original_image_name
    cytomapper::channelNames(images) <- filtered_channel_names
    S4Vectors::mcols(images)$imageID <- names(images)
    
    # Update all_channel_names to be the filtered set
    all_channel_names <- filtered_channel_names
    
    # Clean up
    rm(filtered_array, filtered_image)
    gc()
  } else {
    if (has_existing_temps && n_channels == n_channels_original) {
      timed_message("  ✓ Skipping channel extraction (using existing 31-channel temp files)")
    }
    # No filtering needed - keep images as-is
  }
  
  timed_message("  ✓ Channel filtering complete")
  step_times$channel_selection <- elapsed_time(step_times$channel_selection_start)
  if (timing) timed_message(paste0("  ⏱ Channel selection: ", round(step_times$channel_selection, 2), " sec"))
  
  # ============================================================================
  # STEP 1.5: Optional Downsampling for RAM Efficiency
  # ============================================================================
  needs_upsampling <- FALSE
  
  if (downsample_factor > 1) {
    timed_message("\n--- Step 1.5/5: Downsampling for RAM Efficiency ---", "downsample_start")
    timed_message(paste0("  Original dimensions: ", paste(original_dims[1:2], collapse = " x ")))
    timed_message(paste0("  Downsample factor: ", downsample_factor, "x"))
    
    new_height <- round(original_dims[1] / downsample_factor)
    new_width <- round(original_dims[2] / downsample_factor)
    timed_message(paste0("  New dimensions: ", new_height, " x ", new_width))
    
    original_ram <- prod(original_dims) * 8 / 1e9
    downsampled_ram <- new_height * new_width * n_channels * 8 / 1e9
    timed_message(paste0("  RAM reduction: ", round(original_ram, 1), " GB → ", 
                        round(downsampled_ram, 1), " GB (", 
                        round(100 * downsampled_ram / original_ram, 1), "%)"))
    
    # STRATEGY: Two-phase downsampling to avoid socket-based memory duplication
    # Phase 1: Save full-res channels to disk (serial, RAM-safe)
    # Phase 2: Load+downsample in parallel from disk (workers only see 1 channel each)
    # RESUMABLE: Skip channels that are already processed
    
    timed_message("  Downsampling strategy: Disk-based with resume capability...")
    
    # Create temporary directories
    temp_fullres_dir <- file.path(PROG_DIR, paste0("temp_fullres_", image_name))
    temp_downsampled_dir <- file.path(PROG_DIR, paste0("temp_downsampled_", image_name))
    dir.create(temp_fullres_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(temp_downsampled_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Check which channels are already saved at full-res (for Phase 1)
    # NOTE: Temp files are named with FILTERED indices (1-26), matching the filtered array
    existing_fullres <- list.files(temp_fullres_dir, pattern = "^channel_\\d{2}\\.rds$")
    existing_fullres_nums <- as.integer(gsub("channel_(\\d{2})\\.rds", "\\1", existing_fullres))
    channels_needing_phase1 <- setdiff(1:n_channels, existing_fullres_nums)
    
    # Check which channels are already downsampled (for Phase 2)
    existing_downsampled <- list.files(temp_downsampled_dir, pattern = "^channel_\\d{2}\\.rds$")
    existing_downsampled_nums <- as.integer(gsub("channel_(\\d{2})\\.rds", "\\1", existing_downsampled))
    channels_needing_phase2 <- setdiff(1:n_channels, existing_downsampled_nums)
    
    if (length(existing_fullres_nums) > 0) {
      timed_message(paste0("  ✓ Found ", length(existing_fullres_nums), " full-res channels on disk"))
    }
    if (length(existing_downsampled_nums) > 0) {
      timed_message(paste0("  ✓ Found ", length(existing_downsampled_nums), " downsampled channels on disk"))
    }
    
    # CRITICAL SANITY CHECK: If downsampled count doesn't match full-res count, re-process ALL
    if (length(existing_downsampled_nums) > 0 && 
        length(existing_downsampled_nums) != length(existing_fullres_nums)) {
      timed_message(paste0("  ⚠️  MISMATCH DETECTED: ", length(existing_fullres_nums), 
                          " full-res files but only ", length(existing_downsampled_nums), 
                          " downsampled files"))
      timed_message("     This indicates incomplete processing from a previous run")
      timed_message("     → Will DELETE inconsistent downsampled files and re-process ALL from full-res")
      
      # Delete the inconsistent downsampled directory
      if (dir.exists(temp_downsampled_dir)) {
        unlink(temp_downsampled_dir, recursive = TRUE)
        dir.create(temp_downsampled_dir, showWarnings = FALSE, recursive = TRUE)
        timed_message("     ✓ Cleared inconsistent downsampled directory")
      }
      
      # Re-process ALL full-res files (use full-res count, not filtered count)
      channels_needing_phase2 <- existing_fullres_nums
      existing_downsampled_nums <- integer(0)  # Reset for logging
      
      # CRITICAL: Update n_channels to match full-res count
      # This ensures subsequent steps (enhancement, upsampling) process ALL channels
      n_channels <- length(existing_fullres_nums)
      timed_message(paste0("     ✓ Reset n_channels to ", n_channels, " (matching full-res count)"))
    }
    
    # PHASE 1: Save full-res channels to disk (only those not already saved)
    if (length(channels_needing_phase1) > 0) {
      # Get the filtered array (26 channels after duplicate removal)
      filtered_array_local <- images[[1]]
      
      timed_message(paste0("  Phase 1: Saving ", length(channels_needing_phase1), " full-res channels to disk..."))
      for (i in seq_along(channels_needing_phase1)) {
        # ch_idx is the FILTERED index (1-26)
        ch_idx <- channels_needing_phase1[i]
        timed_message(paste0("    [", i, "/", length(channels_needing_phase1), "] Saving full-res channel ", 
                            ch_idx, ": ", all_channel_names[ch_idx]))
        # Access filtered array with filtered index, save with filtered index
        saveRDS(filtered_array_local[,,ch_idx], 
                file.path(temp_fullres_dir, paste0("channel_", sprintf("%02d", ch_idx), ".rds")),
                compress = FALSE)
      }
      timed_message("  ✓ Phase 1 complete: All full-res channels saved to disk")
      
      # FREE original images NOW - this is critical!
      timed_message("  Freeing filtered full-resolution images...")
      rm(images, filtered_array_local)
      gc()
      timed_message("    ✓ Freed ~57 GB (filtered)")
    } else {
      timed_message("  ✓ Phase 1: All full-res channels already on disk, skipping")
      # Still need to free images
      rm(images)
      gc()
      timed_message("    ✓ Freed ~57 GB (filtered)")
    }
    
    # PHASE 2: Load+downsample in parallel from disk (only those not already downsampled)
    if (length(channels_needing_phase2) > 0) {
      # CRITICAL: Use SERIAL processing for downsampling
      # BiocParallel has persistent issues on WSL2
      BPPARAM_downsample <- BiocParallel::SerialParam()
      n_workers <- BiocParallel::bpnworkers(BPPARAM_downsample)
      
      timed_message(paste0("  Phase 2: Downsampling ", length(channels_needing_phase2), " channels in parallel..."))
      timed_message(paste0("    Using ", n_workers, " workers"))
      timed_message(paste0("    Batch processing ", n_workers, " channels at a time..."))
      timed_message("    ⚠️  Using SerialParam (BiocParallel parallel issues on WSL2)")
      
      # CRITICAL: temp_fullres_dir ALWAYS contains filtered channels (1-26)
      # Direct mapping: ch_idx → file index (no channels_to_keep needed)
      # Phase 1 saved filtered channels with filtered indices, so Phase 2 loads with filtered indices
      
      # Each worker loads ONE channel from disk, downsamples it, saves result
      results <- BiocParallel::bplapply(
        channels_needing_phase2,
        function(ch_idx) {
          # Direct mapping: ch_idx is the filtered index (1-26)
          # temp_fullres already contains filtered channels, so load directly
          source_idx <- ch_idx
          
          # Load this ONE full-res channel from disk
          channel <- readRDS(
            file.path(temp_fullres_dir, paste0("channel_", sprintf("%02d", source_idx), ".rds"))
          )
          
          # Downsample it (suppress package startup messages on first EBImage call)
          downsampled <- suppressPackageStartupMessages({
            EBImage::resize(channel, w = new_width, h = new_height)
          })
          
          # Save downsampled channel using FILTERED index
          saveRDS(downsampled,
                  file.path(temp_downsampled_dir, paste0("channel_", sprintf("%02d", ch_idx), ".rds")),
                  compress = FALSE)
          
          # Return channel index for confirmation
          return(ch_idx)
        },
        BPPARAM = BPPARAM_downsample
      )
      
      # Report completion
      completed_channels <- unlist(results)
      timed_message(paste0("  ✓ Phase 2 complete: ", length(completed_channels), " channels downsampled"))
      timed_message(paste0("    Total downsampled: ", length(existing_downsampled_nums) + length(completed_channels), "/", n_channels))
    } else {
      timed_message("  ✓ Phase 2: All channels already downsampled, skipping")
    }
    
    # Clean up full-res temp files (after both phases)
    if (cleanup_temp) {
      timed_message("  Cleaning up temporary full-res files...")
      unlink(temp_fullres_dir, recursive = TRUE)
    } else {
      timed_message(paste0("  ℹ Keeping temp files: ", temp_fullres_dir))
    }
    
    # Note: Downsampled channels remain on disk in temp_downsampled_dir
    # They will be loaded directly by enhancement workers in Phase 2
    # No need to load them into RAM here!
    
    needs_upsampling <- TRUE
    timed_message("  ✓ Downsampling complete - channels saved to disk")
    step_times$downsample <- elapsed_time(step_times$downsample_start)
    if (timing) timed_message(paste0("  ⏱ Downsample: ", round(step_times$downsample, 2), " sec"))
  }
  
  # ============================================================================
  # STEP 2: Enhance ALL Channels for Quantification
  # ============================================================================
  timed_message("\n--- Step 2/5: Enhancing All Channels ---", "enhance_all_start")
  
  # Get current image dimensions (downsampled if applicable)
  # If downsampled, images object won't exist - use calculated dimensions
  if (needs_upsampling) {
    img_dims <- c(new_height, new_width, n_channels)
  } else {
    img_dims <- original_dims
  }
  
  timed_message(paste0("  Processing ", n_channels, " channels..."))
  timed_message(paste0("  Image dimensions: ", paste(img_dims[1:2], collapse = " x ")))
  timed_message(paste0("  RAM per channel: ~", round(prod(img_dims[1:2]) * 8 / 1e9, 2), " GB"))
  
  # STRATEGY: Disk-based enhancement to avoid socket-based memory duplication
  # Workers load downsampled channels directly from temp_downsampled_dir (already on disk),
  # enhance them, and save results. This avoids passing large objects to workers.
  
  timed_message("  Strategy: Disk-based parallel enhancement (RAM-safe)")
  
  # Point to the downsampled channels directory (created during downsampling phase)
  temp_downsampled_channels_dir <- file.path(PROG_DIR, paste0("temp_downsampled_", image_name))
  
  # Create temporary directory for enhanced channels
  temp_enhanced_dir <- file.path(PROG_DIR, paste0("temp_enhanced_", image_name))
  dir.create(temp_enhanced_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Check which channels are already enhanced
  # SIMPLE: Both downsampled and enhanced directories use the SAME indexing (1-26 after filtering)
  # Just check which indices 1:n_channels are missing from enhanced directory
  existing_enhanced <- list.files(temp_enhanced_dir, pattern = "^channel_\\d{2}\\.rds$")
  existing_enhanced_nums <- as.integer(gsub("channel_(\\d{2})\\.rds", "\\1", existing_enhanced))
  
  # Direct comparison: which indices from 1:n_channels don't have enhanced files?
  channels_to_enhance <- setdiff(1:n_channels, existing_enhanced_nums)
  
  if (length(existing_enhanced_nums) > 0) {
    timed_message(paste0("  ✓ Found ", length(existing_enhanced_nums), " already-enhanced channels"))
    if (length(channels_to_enhance) > 0) {
      timed_message(paste0("  → Resuming enhancement from channel ", min(channels_to_enhance)))
      timed_message(paste0("     Channels to enhance: ", paste(channels_to_enhance, collapse=", ")))
    }
  }
  
  # CRITICAL SANITY CHECK: Enhanced count should match downsampled count
  # Get downsampled count for comparison
  temp_downsampled_channels_dir_check <- file.path(PROG_DIR, paste0("temp_downsampled_", image_name))
  if (dir.exists(temp_downsampled_channels_dir_check)) {
    downsampled_files <- list.files(temp_downsampled_channels_dir_check, pattern = "^channel_\\d{2}\\.rds$")
    downsampled_nums <- as.integer(gsub("channel_(\\d{2})\\.rds", "\\1", downsampled_files))
    
    if (length(existing_enhanced_nums) > 0 && length(existing_enhanced_nums) != length(downsampled_nums)) {
      timed_message(paste0("  ⚠️  MISMATCH DETECTED: ", length(downsampled_nums), 
                          " downsampled files but only ", length(existing_enhanced_nums), 
                          " enhanced files"))
      timed_message("     This indicates incomplete processing from a previous run")
      timed_message("     → Will DELETE inconsistent enhanced files and re-process ALL from downsampled")
      
      # Delete the inconsistent enhanced directory
      if (dir.exists(temp_enhanced_dir)) {
        unlink(temp_enhanced_dir, recursive = TRUE)
        dir.create(temp_enhanced_dir, showWarnings = FALSE, recursive = TRUE)
        timed_message("     ✓ Cleared inconsistent enhanced directory")
      }
      
      # Re-process ALL downsampled files
      channels_to_enhance <- downsampled_nums
      existing_enhanced_nums <- integer(0)  # Reset for logging
      
      # CRITICAL: Update n_channels to match downsampled count
      n_channels <- length(downsampled_nums)
      timed_message(paste0("     ✓ Reset n_channels to ", n_channels, " (matching downsampled count)"))
    }
  }
  
  if (length(channels_to_enhance) > 0) {
    # FREE downsampled images if they exist in RAM
    # (Workers will load channels from disk instead)
    if (exists("images")) {
      timed_message("  Freeing downsampled images from RAM...")
      rm(images)
      gc()
      timed_message("    ✓ Freed ~17 GB")
    }
    
    # CRITICAL: Use SERIAL processing for enhancement
    # BiocParallel SnowParam has persistent "wrong args for environment subassignment" 
    # errors on WSL2, even with fresh BPPARAM objects. Serial is slower but reliable.
    BPPARAM_enhance <- BiocParallel::SerialParam()
    n_workers <- BiocParallel::bpnworkers(BPPARAM_enhance)
    
    # Load+enhance in parallel from disk
    timed_message(paste0("  Enhancing ", length(channels_to_enhance), " channels in parallel..."))
    timed_message(paste0("    Using ", n_workers, " workers"))
    timed_message(paste0("    Processing ", n_workers, " channels at a time..."))
    timed_message("    ⚠️  Using SerialParam (BiocParallel parallel issues on WSL2)")
    
    # Pass channel names and mapping info to workers
    ch_names <- all_channel_names
    
    # CRITICAL: Downsampled directory ALWAYS uses filtered indices (1-23)
    # Both downsampled and enhanced use the SAME indexing after duplicate removal
    # So ch_idx maps directly to the file index (no need for channels_to_keep mapping)
    
    # Each worker loads ONE downsampled channel from disk, enhances it, saves result
    results <- BiocParallel::bplapply(
      channels_to_enhance,
      function(ch_idx) {
        # Direct mapping: ch_idx → file index
        # Both downsampled and enhanced directories use filtered indices (1-23)
        source_idx <- ch_idx
        
        # Load this ONE downsampled channel from disk
        channel <- readRDS(
          file.path(temp_downsampled_channels_dir, paste0("channel_", sprintf("%02d", source_idx), ".rds"))
        )
        
        # Enhance it (use verbose_enhancement parameter for diagnostics)
        enhanced <- enhance_channel(channel, ch_names[ch_idx], verbose_inner = verbose_enhancement)
        
        # Save enhanced channel using FILTERED index
        saveRDS(enhanced,
                file.path(temp_enhanced_dir, paste0("channel_", sprintf("%02d", ch_idx), ".rds")),
                compress = FALSE)
        
        # Return channel index for confirmation
        return(ch_idx)
      },
      BPPARAM = BPPARAM_enhance
    )
    
    # Report completion
    completed_channels <- unlist(results)
    timed_message(paste0("  ✓ Enhancement complete: ", length(completed_channels), " channels enhanced"))
    timed_message(paste0("    Total enhanced: ", length(existing_enhanced_nums) + length(completed_channels), "/", n_channels))
  } else {
    timed_message("  ✓ All channels already enhanced, skipping")
    # Still need to free images if we're skipping enhancement
    if (exists("images")) {
      rm(images)
      gc()
    }
  }
  
  # Note: Enhanced channels remain on disk in temp_enhanced_dir
  # If upsampling is needed, they'll be loaded directly by upsampling workers
  # If no upsampling, they'll be loaded for final assembly after Step 2.5
  timed_message("  ✓ Enhanced channels saved to disk")
  
  # ============================================================================
  # STEP 2.5: Upsample Back to Full Resolution (if downsampled)
  # ============================================================================
  if (needs_upsampling) {
    timed_message("\n--- Step 2.5/5: Upsampling to Full Resolution ---", "upsample_start")
    timed_message(paste0("  Current dimensions: ", paste(img_dims[1:2], collapse = " x ")))
    timed_message(paste0("  Target dimensions: ", paste(original_dims[1:2], collapse = " x ")))
    
    # STRATEGY: Parallel disk-based upsampling
    # Workers load enhanced channels from disk, upsample, and save results
    # This avoids loading all enhanced channels into RAM
    
    # Point to the enhanced channels directory
    temp_enhanced_channels_dir <- file.path(PROG_DIR, paste0("temp_enhanced_", image_name))
    
    # Create temporary directory for upsampled channels
    temp_upsampled_dir <- file.path(PROG_DIR, paste0("temp_upsampled_", image_name))
    dir.create(temp_upsampled_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Check which channels are already upsampled
    existing_upsampled <- list.files(temp_upsampled_dir, pattern = "^channel_\\d{2}\\.rds$")
    existing_upsampled_nums <- as.integer(gsub("channel_(\\d{2})\\.rds", "\\1", existing_upsampled))
    channels_to_upsample <- setdiff(1:n_channels, existing_upsampled_nums)
    
    if (length(existing_upsampled_nums) > 0) {
      timed_message(paste0("  ✓ Found ", length(existing_upsampled_nums), " already-upsampled channels"))
      if (length(channels_to_upsample) > 0) {
        timed_message(paste0("  → Resuming upsampling from channel ", min(channels_to_upsample)))
      }
    }
    
    # CRITICAL SANITY CHECK: Upsampled count should match enhanced count
    temp_enhanced_channels_dir_check <- file.path(PROG_DIR, paste0("temp_enhanced_", image_name))
    if (dir.exists(temp_enhanced_channels_dir_check)) {
      enhanced_files <- list.files(temp_enhanced_channels_dir_check, pattern = "^channel_\\d{2}\\.rds$")
      enhanced_nums <- as.integer(gsub("channel_(\\d{2})\\.rds", "\\1", enhanced_files))
      
    if (length(existing_upsampled_nums) > 0 && length(existing_upsampled_nums) != length(enhanced_nums)) {
      timed_message(paste0("  ⚠️  MISMATCH DETECTED: ", length(enhanced_nums), 
                          " enhanced files but only ", length(existing_upsampled_nums), 
                          " upsampled files"))
      timed_message("     This indicates incomplete processing from a previous run")
      timed_message("     → Will DELETE inconsistent upsampled files and re-process ALL from enhanced")
      
      # Delete the inconsistent upsampled directory
      if (dir.exists(temp_upsampled_dir)) {
        unlink(temp_upsampled_dir, recursive = TRUE)
        dir.create(temp_upsampled_dir, showWarnings = FALSE, recursive = TRUE)
        timed_message("     ✓ Cleared inconsistent upsampled directory")
      }
      
      # Re-process ALL enhanced files
      channels_to_upsample <- enhanced_nums
      existing_upsampled_nums <- integer(0)  # Reset for logging
      
      # CRITICAL: Update n_channels to match enhanced count
      n_channels <- length(enhanced_nums)
      timed_message(paste0("     ✓ Reset n_channels to ", n_channels, " (matching enhanced count)"))
    }
    }
    
    if (length(channels_to_upsample) > 0) {
      # CRITICAL: Use SERIAL processing for upsampling
      # BiocParallel has persistent issues on WSL2
      BPPARAM_upsample <- BiocParallel::SerialParam()
      n_workers <- BiocParallel::bpnworkers(BPPARAM_upsample)
      
      # Parallel upsampling from disk
      timed_message(paste0("  Upsampling ", length(channels_to_upsample), " channels in parallel..."))
      timed_message(paste0("    Using ", n_workers, " workers"))
      timed_message(paste0("    Processing ", n_workers, " channels at a time..."))
      timed_message("    ⚠️  Using SerialParam (BiocParallel parallel issues on WSL2)")
      
      # Pass dimensions to workers
      orig_dims <- original_dims
      
      # CRITICAL: Enhanced directory uses filtered indices (1-23)
      # Direct mapping: ch_idx → file index (no channels_to_keep needed)
      
      # Each worker loads ONE enhanced channel from disk, upsamples it, saves result
      results <- BiocParallel::bplapply(
        channels_to_upsample,
        function(ch_idx) {
          # Direct mapping: enhanced files use filtered indices (1-23)
          source_idx <- ch_idx
          
          # Load this ONE enhanced channel from disk
          channel <- readRDS(
            file.path(temp_enhanced_channels_dir, paste0("channel_", sprintf("%02d", source_idx), ".rds"))
          )
          
          # Upsample it (suppress package startup messages)
          upsampled <- suppressPackageStartupMessages({
            EBImage::resize(channel,
                           w = orig_dims[2],
                           h = orig_dims[1])
          })
          
          # Save upsampled channel using FILTERED index
          saveRDS(upsampled,
                  file.path(temp_upsampled_dir, paste0("channel_", sprintf("%02d", ch_idx), ".rds")),
                  compress = FALSE)
          
          # Return channel index for confirmation
          return(ch_idx)
        },
        BPPARAM = BPPARAM_upsample
      )
      
      # Report completion
      completed_channels <- unlist(results)
      timed_message(paste0("  ✓ Upsampling complete: ", length(completed_channels), " channels upsampled"))
      timed_message(paste0("    Total upsampled: ", length(existing_upsampled_nums) + length(completed_channels), "/", n_channels))
    } else {
      timed_message("  ✓ All channels already upsampled, skipping")
    }
    
    timed_message("  ✓ Upsampling complete")
    step_times$upsample <- elapsed_time(step_times$upsample_start)
    if (timing) timed_message(paste0("  ⏱ Upsample: ", round(step_times$upsample, 2), " sec"))
    
    # SANITY CHECK: Create before/after comparison images
    # Determine how many channels to compare based on comparison_channels parameter
    if (is.character(comparison_channels) && tolower(comparison_channels) == "none") {
      timed_message("\n  Skipping enhancement comparison (comparison_channels = 'none')")
    } else {
      timed_message("\n  Creating before/after enhancement comparisons...")
      comparison_dir <- file.path(image_res_dir, "enhancement_comparison")
      dir.create(comparison_dir, showWarnings = FALSE, recursive = TRUE)
      
      # Compare full-res pre-enhanced vs post-enhanced (upsampled)
      temp_fullres_dir <- file.path(PROG_DIR, paste0("temp_fullres_", image_name))
      temp_upsampled_dir <- file.path(PROG_DIR, paste0("temp_upsampled_", image_name))
      
      # Determine number of channels to compare
      if (is.character(comparison_channels) && tolower(comparison_channels) == "all") {
        max_channels_to_compare <- n_channels
      } else {
        max_channels_to_compare <- min(as.numeric(comparison_channels), n_channels)
      }
      
      # Check which comparison images already exist (resume capability)
      existing_comparisons <- list.files(comparison_dir, pattern = "^comparison_ch\\d{2}_.*\\.png$")
      existing_comparison_nums <- as.integer(gsub("comparison_ch(\\d{2})_.*\\.png", "\\1", existing_comparisons))
      channels_to_visualize <- setdiff(1:max_channels_to_compare, existing_comparison_nums)
      
      if (length(existing_comparison_nums) > 0) {
        timed_message(paste0("  ✓ Found ", length(existing_comparison_nums), " existing comparison images"))
        if (length(channels_to_visualize) > 0) {
          timed_message(paste0("  → Creating ", length(channels_to_visualize), " remaining comparisons"))
        }
      }
      
      if (length(channels_to_visualize) > 0) {
        # CRITICAL: Use SERIAL processing for comparison generation
        # BiocParallel has persistent issues on WSL2
        BPPARAM_compare <- BiocParallel::SerialParam()
        n_workers <- BiocParallel::bpnworkers(BPPARAM_compare)
        
        timed_message(paste0("  Generating comparison PNGs in parallel..."))
        timed_message(paste0("    Using ", n_workers, " workers"))
        timed_message(paste0("    Processing ", n_workers, " channels at a time..."))
        timed_message("    ⚠️  Using SerialParam (BiocParallel parallel issues on WSL2)")
        
        # Pass channel names to workers
        ch_names <- all_channel_names
        
        # CRITICAL: BOTH temp_fullres AND temp_upsampled use FILTERED indices (1-26)
        # Phase 1 saves filtered channels with filtered indices
        # Direct mapping: ch_idx → file index (no mapping needed)
        
        # Parallel comparison generation
        results <- BiocParallel::bplapply(
          channels_to_visualize,
          function(ch_idx) {
            # Direct mapping: both directories use filtered indices (1-26)
            # ch_idx is the filtered index, use it directly for both
            
            # Load pre-enhanced (full-res, filtered) with FILTERED index
            pre_enhanced <- readRDS(
              file.path(temp_fullres_dir, paste0("channel_", sprintf("%02d", ch_idx), ".rds"))
            )
            
            # Load post-enhanced (upsampled, filtered) with FILTERED index
            post_enhanced <- readRDS(
              file.path(temp_upsampled_dir, paste0("channel_", sprintf("%02d", ch_idx), ".rds"))
            )
            
            # Create comparison image
            channel_name <- ch_names[ch_idx]
            safe_channel_name <- gsub("[^A-Za-z0-9_-]", "_", channel_name)
            comparison_file <- file.path(comparison_dir, 
                                         paste0("comparison_ch", sprintf("%02d", ch_idx), 
                                               "_", safe_channel_name, ".png"))
            
            # Save side-by-side comparison
            png(comparison_file, width = 2400, height = 1200)
            par(mfrow = c(1, 2), mar = c(2, 2, 3, 1))
            
            # Normalize for display
            pre_norm <- (pre_enhanced - min(pre_enhanced)) / (max(pre_enhanced) - min(pre_enhanced))
            post_norm <- (post_enhanced - min(post_enhanced)) / (max(post_enhanced) - min(post_enhanced))
            
            # Display pre-enhanced
            EBImage::display(pre_norm, method = "raster")
            title(main = paste0("BEFORE Enhancement\n", channel_name), 
                  cex.main = 1.5, col.main = "red")
            
            # Display post-enhanced
            EBImage::display(post_norm, method = "raster")
            title(main = paste0("AFTER Enhancement\n", channel_name), 
                  cex.main = 1.5, col.main = "green")
            
            dev.off()
            
            # Return channel index for confirmation
            return(ch_idx)
          },
          BPPARAM = BPPARAM_compare
        )
        
        # Report completion
        completed_comparisons <- unlist(results)
        timed_message(paste0("  ✓ Comparison generation complete: ", length(completed_comparisons), " images created"))
        timed_message(paste0("    Total comparisons: ", length(existing_comparison_nums) + length(completed_comparisons), "/", max_channels_to_compare))
      } else {
        timed_message(paste0("  ✓ All ", max_channels_to_compare, " comparison images already exist"))
      }
      
      timed_message(paste0("  ✓ Comparison images saved to: ", comparison_dir))
    }
    
    # SKIP full-resolution CytoImageList assembly
    # Channels will be loaded on-demand during quantification (Step 4)
    timed_message("\n  Skipping full-resolution assembly (channels loaded on-demand in Step 4)")
    
    # Keep upsampled temp directory for quantification
    if (cleanup_temp) {
      timed_message(paste0("  ℹ Keeping upsampled temp files for quantification: ", temp_upsampled_dir))
      timed_message("     Will be cleaned after Step 4")
    } else {
      timed_message(paste0("  ℹ Keeping temp files: ", temp_upsampled_dir))
    }
    
  } else {
    # No upsampling needed - channels will be loaded on-demand during quantification
    timed_message("\n  Skipping enhanced assembly (channels loaded on-demand in Step 4)")
    
    # Point to enhanced channels for quantification
    temp_upsampled_dir <- file.path(PROG_DIR, paste0("temp_enhanced_", image_name))
    
    # Keep enhanced temp directory for quantification
    if (cleanup_temp) {
      timed_message(paste0("  ℹ Keeping enhanced temp files for quantification: ", temp_enhanced_dir))
      timed_message("     Will be cleaned after Step 4")
    } else {
      timed_message(paste0("  ℹ Keeping temp files: ", temp_enhanced_dir))
    }
  }
  
  step_times$enhance_all <- elapsed_time(step_times$enhance_all_start, units = "mins")
  if (timing) timed_message(paste0("  ⏱ Total enhancement + resampling: ", round(step_times$enhance_all, 2), " min"))
  
  # ============================================================================
  # STEP 2.75: Export Single-Channel TIFFs (Optional)
  # ============================================================================
  if (export_tiffs) {
    timed_message("\n--- Step 2.75/5: Exporting Single-Channel TIFFs ---", "export_tiffs_start")
    
    # Determine which directory contains the final enhanced channels
    # If upsampling was done, use upsampled directory (full resolution)
    # Otherwise, use enhanced directory (already at full resolution)
    if (needs_upsampling) {
      source_channel_dir <- file.path(PROG_DIR, paste0("temp_upsampled_", image_name))
      timed_message("  Source: Upsampled enhanced channels (full resolution)")
    } else {
      source_channel_dir <- file.path(PROG_DIR, paste0("temp_enhanced_", image_name))
      timed_message("  Source: Enhanced channels (full resolution)")
    }
    
    # CRITICAL: Get channel names from source directory to ensure correct mapping
    # The temp files are the source of truth for which channels actually exist
    source_files <- list.files(source_channel_dir, pattern = "^channel_\\d{2}\\.rds$", full.names = FALSE)
    source_channel_nums <- as.integer(gsub("channel_(\\d{2})\\.rds", "\\1", source_files))
    source_channel_nums <- sort(source_channel_nums)
    
    timed_message(paste0("  Found ", length(source_channel_nums), " channels in source directory"))
    
    # Load channel metadata to get correct names
    metadata_file <- file.path(PROG_DIR, paste0(image_name, "_channel_metadata.rds"))
    if (file.exists(metadata_file)) {
      channel_metadata <- readRDS(metadata_file)
      # Use metadata channel names, filtered to match source files
      all_channel_names_export <- channel_metadata$channel_names[source_channel_nums]
      timed_message("  ✓ Using channel names from metadata file")
    } else {
      # Fallback: use current all_channel_names (may be incorrect if resuming)
      all_channel_names_export <- all_channel_names
      timed_message("  ⚠️  No metadata file found, using current channel names (may be incorrect)")
    }
    
    # Create output directory for TIFFs
    tiff_output_dir <- file.path(image_res_dir, "enhanced_tiffs")
    dir.create(tiff_output_dir, showWarnings = FALSE, recursive = TRUE)
    timed_message(paste0("  Output directory: ", tiff_output_dir))
    
    # Check which TIFFs already exist (resume capability)
    existing_tiffs <- list.files(tiff_output_dir, pattern = "^.*\\.tif$")
    existing_tiff_nums <- as.integer(gsub(".*_ch(\\d{2})_.*\\.tif", "\\1", existing_tiffs))
    channels_to_export <- setdiff(source_channel_nums, existing_tiff_nums)
    
    if (length(existing_tiff_nums) > 0) {
      timed_message(paste0("  ✓ Found ", length(existing_tiff_nums), " existing TIFF files"))
      if (length(channels_to_export) > 0) {
        timed_message(paste0("  → Exporting ", length(channels_to_export), " remaining channels"))
      }
    }
    
    if (length(channels_to_export) > 0) {
      timed_message(paste0("  Exporting ", length(channels_to_export), " channels as TIFFs..."))
      
      # Export each channel
      for (i in seq_along(channels_to_export)) {
        ch_idx <- channels_to_export[i]
        
        # Map channel index to position in export names array
        # source_channel_nums contains the sorted indices (e.g., 1,2,3,...,26)
        # Find position of ch_idx in source_channel_nums
        name_position <- which(source_channel_nums == ch_idx)
        
        if (length(name_position) == 0) {
          timed_message(paste0("    ⚠️  WARNING: Channel ", ch_idx, " not found in source channel list"))
          next
        }
        
        channel_name <- all_channel_names_export[name_position]
        
        if (i %% 5 == 1 || i == length(channels_to_export)) {
          timed_message(paste0("    [", i, "/", length(channels_to_export), "] Exporting channel ", 
                              ch_idx, ": ", channel_name))
        }
        
        # Load channel from temp directory
        channel_file <- file.path(source_channel_dir, paste0("channel_", sprintf("%02d", ch_idx), ".rds"))
        
        if (!file.exists(channel_file)) {
          timed_message(paste0("    ⚠️  WARNING: Channel file not found: ", channel_file))
          next
        }
        
        channel <- readRDS(channel_file)
        
        # If channel is an EBImage object, extract the raw array
        if (inherits(channel, "Image")) {
          channel_data <- EBImage::imageData(channel)
        } else {
          channel_data <- channel
        }
        
        # Create safe filename from channel name
        safe_channel_name <- gsub("[^A-Za-z0-9_-]", "_", channel_name)
        
        # Construct TIFF filename with channel index and name
        tiff_filename <- file.path(tiff_output_dir, 
                                   paste0(image_name, "_ch", sprintf("%02d", ch_idx), 
                                         "_", safe_channel_name, ".tif"))
        
        # Convert to EBImage Image if needed and write TIFF
        if (!inherits(channel_data, "Image")) {
          channel_image <- EBImage::Image(channel_data)
        } else {
          channel_image <- channel_data
        }
        
        # Write TIFF (16-bit format to preserve dynamic range)
        EBImage::writeImage(channel_image, tiff_filename, type = "tiff", bits.per.sample = 16)
        
        # Clean up
        rm(channel, channel_data, channel_image)
      }
      
      timed_message(paste0("  ✓ Export complete: ", length(channels_to_export), " TIFFs created"))
      timed_message(paste0("    Total TIFFs: ", length(existing_tiff_nums) + length(channels_to_export), "/", length(source_channel_nums)))
    } else {
      timed_message(paste0("  ✓ All ", length(source_channel_nums), " TIFF files already exist, skipping export"))
    }
    
    timed_message(paste0("  ✓ Single-channel TIFFs saved to: ", tiff_output_dir))
    
    step_times$export_tiffs <- elapsed_time(step_times$export_tiffs_start)
    if (timing) timed_message(paste0("  ⏱ TIFF export: ", round(step_times$export_tiffs, 2), " sec"))
  } else {
    timed_message("\n--- Step 2.75/5: Skipping TIFF Export (export_tiffs = FALSE) ---")
    step_times$export_tiffs <- 0
  }
  
  # ============================================================================
  # STEP 3: Generate Multi-Channel Visualization (Optional)
  # ============================================================================
  viz_file <- NULL  # Initialize in case visualization is skipped
  
  if (generate_visualization) {
    timed_message("\n--- Step 3/5: Creating Multi-Channel Overlay ---", "viz_start")
    timed_message("  NOTE: Visualization step skipped - enhanced_images loaded later for quantification")
    timed_message("  To enable: Load enhanced channels from disk, extract viz channels, then generate overlay")
    
    # Visualization would require loading enhanced_images from disk here
    # For now, skip this step to save time and RAM
    # Future implementation: Load only viz channels from temp_upsampled_dir
    
    viz_file <- NULL
    step_times$viz <- 0
    
    timed_message("  ✓ Skipped (implement if needed)")
  } else {
    timed_message("\n--- Step 3/5: Skipping Multi-Channel Overlay (generate_visualization = FALSE) ---")
    step_times$viz <- 0
  }
  
  # ============================================================================
  # STEP 4: Quantify Cell Features (Using Downsampled Images for RAM Efficiency)
  # ============================================================================
  timed_message("\n--- Step 4/5: Quantifying Cell Features ---", "quantify_start")
  timed_message("  NOTE: Using DOWNSAMPLED images and masks to stay within RAM limits")
  timed_message(paste0("  Downsampled by factor of ", downsample_factor, "x"))
  
  # Load masks now (only when needed for quantification)
  # DOWNSAMPLE masks to match downsampled images if needed
  masks <- NULL
  if (!exists("masks") || is.null(masks)) {
    timed_message("\n  Loading masks...")
    
    # Detect file format and load
    if (grepl("\\.qs$", masks_file, ignore.case = TRUE)) {
      masks_loaded <- qs::qread(masks_file)
      timed_message("    ✓ Loaded masks from QS format")
    } else {
      masks_loaded <- readRDS(masks_file)
      timed_message("    ✓ Loaded masks from RDS format")
    }
    
    mask_array_loaded <- masks_loaded[[1]]
    loaded_height <- dim(mask_array_loaded)[1]
    loaded_width <- dim(mask_array_loaded)[2]
    
    # Calculate expected downsampled dimensions
    expected_height <- round(original_dims[1] / downsample_factor)
    expected_width <- round(original_dims[2] / downsample_factor)
    
    timed_message(paste0("    Loaded dimensions: ", loaded_height, " x ", loaded_width))
    timed_message(paste0("    Expected dimensions: ", expected_height, " x ", expected_width, 
                        " (", downsample_factor, "x downsampled)"))
    
    # Check if we need to downsample
    if (loaded_height == expected_height && loaded_width == expected_width) {
      # Masks are already at the right resolution
      timed_message("    ✓ Masks are already at downsampled resolution - no resizing needed")
      masks <- masks_loaded
    } else if (loaded_height == original_dims[1] && loaded_width == original_dims[2]) {
      # Masks are at full resolution - need to downsample
      timed_message(paste0("    Downsampling masks by ", downsample_factor, "x..."))
      
      # Downsample using nearest neighbor (preserves cell IDs)
      mask_downsampled <- EBImage::resize(mask_array_loaded, 
                                         w = expected_width, 
                                         h = expected_height, 
                                         filter = "none")  # "none" = nearest neighbor
      
      # Wrap in CytoImageList
      masks <- cytomapper::CytoImageList(list(mask_downsampled))
      names(masks) <- names(masks_loaded)
      S4Vectors::mcols(masks) <- S4Vectors::mcols(masks_loaded)
      
      timed_message(paste0("    ✓ Downsampled to: ", paste(dim(masks[[1]])[1:2], collapse = " x ")))
      
      # Clean up
      rm(mask_downsampled)
    } else {
      # Unexpected dimensions - warn and try to resize anyway
      timed_message("    ⚠️  Unexpected mask dimensions - attempting to resize...")
      mask_resized <- EBImage::resize(mask_array_loaded, 
                                      w = expected_width, 
                                      h = expected_height, 
                                      filter = "none")
      masks <- cytomapper::CytoImageList(list(mask_resized))
      names(masks) <- names(masks_loaded)
      S4Vectors::mcols(masks) <- S4Vectors::mcols(masks_loaded)
      rm(mask_resized)
    }
    
    # Clean up loaded masks
    rm(masks_loaded, mask_array_loaded)
    gc()
    
    masks_ram <- prod(dim(masks[[1]])) * 8 / 1e9
    timed_message(paste0("    Final mask dimensions: ", paste(dim(masks[[1]])[1:2], collapse = " x ")))
    timed_message(paste0("    RAM: ~", round(masks_ram, 1), " GB"))
  }
  
  # Get number of cells for estimation
  mask_array <- masks[[1]]
  n_cells <- max(mask_array)
  
  timed_message(paste0("  Cells to quantify: ", format(n_cells, big.mark = ",")))
  
  # NOTE: Duplicate filtering already happened in Step 1
  # n_channels already reflects the correct filtered count (26 unique channels)
  # No additional duplicate detection needed here
  duplicate_channels <- integer(0)
  
  timed_message(paste0("  Channels for quantification: ", n_channels, " (duplicates already removed in Step 1)"))
  
  # Load all n_channels (duplicates already removed in Step 1)
  channels_for_quant <- 1:n_channels
  n_channels_for_quant <- n_channels
  
  # RAM CALCULATION: Based on DOWNSAMPLED dimensions (not full-res!)
  # Downsampled dimensions
  downsampled_height <- round(original_dims[1] / downsample_factor)
  downsampled_width <- round(original_dims[2] / downsample_factor)
  ram_per_channel_downsampled <- downsampled_height * downsampled_width * 8 / 1e9
  
  reduced_ram <- n_channels_for_quant * ram_per_channel_downsampled
  masks_ram <- prod(dim(mask_array)) * 8 / 1e9  # Actual downsampled mask RAM
  per_worker_ram <- reduced_ram + masks_ram
  
  # CRITICAL FIX: Force serial processing to avoid protection stack overflow
  # The large CytoImageList (56+ GB) causes BiocParallel workers to hit R's
  # internal protection stack limit, even when RAM is sufficient.
  # Serial processing is slower but reliable for objects this large.
  
  timed_message("\n  Creating BiocParallel backend for quantification...")
  timed_message("  ⚠️  Using SERIAL processing (BiocParallel::SerialParam())")
  timed_message("     Large objects (56+ GB) cause 'protection stack overflow' in parallel workers")
  timed_message("     Serial is slower but avoids R stack limitations")
  
  BPPARAM_quant <- BiocParallel::SerialParam()
  
  # Calculate RAM estimate for serial processing
  ram_limit <- 250  # GB
  serial_ram <- reduced_ram + masks_ram
  
  # Increase R's protection stack and max nested expressions
  # This helps prevent stack overflows during complex operations
  old_expressions <- getOption("expressions")
  options(expressions = 500000)  # Default is 5000
  
  timed_message(paste0("  ✓ Set options(expressions = 500000) [was ", old_expressions, "]"))
  
  timed_message("\n  Strategy: Serial processing (protection stack overflow workaround)")
  timed_message(paste0("     Image stack (downsampled): ", round(reduced_ram, 1), " GB (", 
                      n_channels_for_quant, " channels at ", downsampled_height, "x", downsampled_width, ")"))
  timed_message(paste0("     Masks (downsampled): ", round(masks_ram, 1), " GB"))
  timed_message(paste0("     Total RAM estimate: ", round(serial_ram, 1), " GB (< ", ram_limit, " GB ✓)"))
  timed_message(paste0("     Estimated time: ~", round(n_cells / 1500, 1), 
                      " minutes (serial processing)"))
  
  # Point to DOWNSAMPLED enhanced channels directory (not upsampled)
  # This uses the enhanced but still downsampled images for RAM efficiency
  temp_enhanced_dir <- file.path(PROG_DIR, paste0("temp_enhanced_", image_name))
  
  # FREE the original images object if it still exists in RAM
  if (exists("images")) {
    timed_message("\n  Freeing original images object from RAM...")
    original_image_size_gb <- prod(dim(images[[1]])) * 8 / 1e9
    rm(images)
    gc()
    timed_message(paste0("    ✓ Freed ~", round(original_image_size_gb, 1), " GB"))
  }
  
  # Load DOWNSAMPLED enhanced channels from disk (ONLY non-duplicate channels)
  timed_message("\n  Loading DOWNSAMPLED enhanced channels from disk...")
  timed_message(paste0("    Directory: ", temp_enhanced_dir))
  timed_message(paste0("    Using enhanced channels at ", downsample_factor, "x downsampled resolution"))
  
  # CRITICAL: Always exclude duplicates for quantification, even if temp files exist for them
  # This ensures we only quantify unique channels
  # Use channels_for_quant which was calculated earlier with dynamic duplicate detection
  channels_to_load <- channels_for_quant
  n_channels_quant <- n_channels_for_quant
  all_channel_names_quant <- all_channel_names[channels_to_load]
  
  if (has_existing_temps && n_channels == n_channels_original && length(duplicate_channels) > 0) {
    timed_message("  ⚠️  FILTERING duplicates for quantification (even though temp files have all 31)")
    timed_message("     Will load ONLY non-duplicate channels from temp directory")
    timed_message(paste0("     Loading ", n_channels_quant, " non-duplicate channels (skipping ", 
                        length(duplicate_channels), " duplicates: ", paste(duplicate_channels, collapse=", "), ")"))}
  else {
    timed_message(paste0("     Loading ", n_channels_quant, " channels"))
  }
  
  # Get mask dimensions (downsampled)
  target_height <- dim(mask_array)[1]
  target_width <- dim(mask_array)[2]
  
  timed_message(paste0("  Target dimensions: ", target_height, " x ", target_width, 
                      " (downsampled by ", downsample_factor, "x)"))
  
  # Create array for NON-DUPLICATE channels only (26 channels) at downsampled resolution
  enhanced_array <- array(0, dim = c(target_height, target_width, n_channels_quant))
  
  # Load each non-duplicate channel
  for (i in 1:n_channels_quant) {
    if (i %% 5 == 1 || i == n_channels_quant) {
      timed_message(paste0("    Loading ", i, "/", n_channels_quant, ": ", all_channel_names_quant[i]))
    }
    
    # CRITICAL: temp_enhanced uses FILTERED indices (1-23), not original
    # Direct mapping: position i → file i
    channel_idx <- i
    channel_file <- file.path(temp_enhanced_dir, paste0("channel_", sprintf("%02d", channel_idx), ".rds"))
    
    if (!file.exists(channel_file)) {
      stop("Missing channel file: ", channel_file, 
           "\n  Position ", i, " → Channel file ", channel_idx, 
           " (", all_channel_names_quant[i], ")")
    }
    
    channel <- readRDS(channel_file)
    
    # If channel is an EBImage object, extract the raw array
    if (inherits(channel, "Image")) {
      channel <- EBImage::imageData(channel)
    }
    
    # Check if dimensions match mask (allowing for transpose)
    if (dim(channel)[1] != target_height || dim(channel)[2] != target_width) {
      # Try transpose if dimensions are swapped
      if (dim(channel)[1] == target_width && dim(channel)[2] == target_height) {
        channel <- t(channel)
      } else {
        stop("Channel ", i, " (", all_channel_names_quant[i], ") dimensions ", 
             paste(dim(channel), collapse="x"), 
             " don't match mask dimensions ", target_height, "x", target_width)
      }
    }
    
    # Store in array at filtered position
    enhanced_array[,,i] <- channel
    rm(channel)
  }
  
  # Update variables for quantification
  n_channels <- n_channels_quant
  all_channel_names <- all_channel_names_quant
  
  timed_message("  ✓ All channels loaded")
  
  # Set dimension names (use filtered channel names)
  dimnames(enhanced_array)[[3]] <- all_channel_names
  
  # Create EBImage and CytoImageList
  timed_message("  Creating CytoImageList...")
  
  # Verify original_image_name is set
  if (is.null(original_image_name) || is.na(original_image_name)) {
    timed_message("    ⚠️  original_image_name not set, using image_name parameter")
    original_image_name <- image_name
  }
  timed_message(paste0("    Image name: ", original_image_name))
  
  # Create Image and wrap in list with name
  enhanced_image <- EBImage::Image(enhanced_array)
  enhanced_list <- list(enhanced_image)
  names(enhanced_list) <- original_image_name
  
  # Create CytoImageList
  enhanced_images <- cytomapper::CytoImageList(enhanced_list)
  
  # Set channel names and metadata
  cytomapper::channelNames(enhanced_images) <- all_channel_names
  S4Vectors::mcols(enhanced_images)$imageID <- names(enhanced_images)
  
  # Verify metadata is set
  timed_message(paste0("    Names: ", paste(names(enhanced_images), collapse=", ")))
  timed_message(paste0("    imageID in mcols: ", "imageID" %in% colnames(S4Vectors::mcols(enhanced_images))))
  
  # Clean up
  rm(enhanced_array, enhanced_image)
  gc()
  
  timed_message("  ✓ CytoImageList ready for quantification")
  timed_message(paste0("     Dimensions: ", paste(dim(enhanced_images[[1]]), collapse = " x ")))
  # Calculate actual RAM usage based on loaded channels
  actual_ram_usage <- n_channels_quant * ram_per_channel_downsampled
  timed_message(paste0("     RAM usage: ~", round(actual_ram_usage, 1), " GB (", n_channels_quant, " channels)"))
  
  # CRITICAL: Aggressive garbage collection before measureObjects
  # This reduces R's protection stack usage by clearing unused objects
  timed_message("\n  Preparing for quantification...")
  timed_message("     Running aggressive garbage collection...")
  gc(full = TRUE, reset = TRUE, verbose = FALSE)
  timed_message("     ✓ Memory cleared")
  
  # Quantify using SerialParam (protection stack overflow workaround)
  timed_message(paste0("\n  ⏳ Quantifying ", n_cells, " cells across ", n_channels_quant, " channels..."))
  timed_message(paste0("     Using ", class(BPPARAM_quant)[1], " (serial processing)"))
  
  quant_start <- Sys.time()
  
  cells <- cytomapper::measureObjects(
    masks,
    enhanced_images,
    img_id = "imageID",
    return_as = "spe",
    BPPARAM = BPPARAM_quant  # Use fresh BPPARAM (avoids worker corruption)
  )
  
  quant_elapsed <- as.numeric(difftime(Sys.time(), quant_start, units = "mins"))
  timed_message(paste0("  ✓ Quantification complete in ", round(quant_elapsed, 1), " minutes"))
  
  # Restore original expression limit
  options(expressions = old_expressions)
  
  # Clean up
  if (cleanup_temp) {
    timed_message("  Cleaning up temporary channel files...")
    unlink(temp_upsampled_dir, recursive = TRUE)
  } else {
    timed_message(paste0("  ℹ Keeping upsampled files: ", temp_upsampled_dir))
  }
  
  # Note about enhanced images
  enhanced_images_file <- file.path(PROG_DIR, paste0(image_name, "_images_enhanced.qs"))
  
  # FREE masks and enhanced_images
  timed_message("  Freeing masks and enhanced_images from RAM...")
  rm(masks, enhanced_images)
  gc()
  timed_message("    ✓ Freed ~60 GB")
  
  # Set spatial coordinate names
  SpatialExperiment::spatialCoordsNames(cells) <- c("x", "y")
  
  # Add image name to colData
  cells$image_name <- image_name
  
  n_cells_measured <- ncol(cells)
  n_features <- nrow(cells)
  
  timed_message("\n=== Quantification Results ===")
  timed_message(paste0("  Cells measured: ", n_cells_measured))
  timed_message(paste0("  Features per cell: ", n_features))
  timed_message(paste0("  Spatial coordinates: ", paste(SpatialExperiment::spatialCoordsNames(cells), collapse = ", ")))
  
  step_times$quantify <- elapsed_time(step_times$quantify_start, units = "mins")
  if (timing) timed_message(paste0("  ⏱ Quantification: ", round(step_times$quantify, 2), " min"))
  
  # ============================================================================
  # STEP 5: Save Results
  # ============================================================================
  timed_message("\n--- Step 5/5: Saving Results ---", "save_start")
  
  # Save SpatialExperiment as .qs (faster than RDS)
  spe_file <- file.path(PROG_DIR, paste0(image_name, "_cells.qs"))
  qs::qsave(cells, spe_file, preset = "fast")
  timed_message(paste0("  ✓ SpatialExperiment saved: ", spe_file))
  timed_message(paste0("  ℹ Enhanced images already saved: ", enhanced_images_file))
  
  step_times$save <- elapsed_time(step_times$save_start)
  if (timing) timed_message(paste0("  ⏱ Saving: ", round(step_times$save, 2), " sec"))
  
  # Final summary
  total_time <- elapsed_time(overall_start, units = "mins")
  timed_message("\n✓ Enhancement and quantification complete!")
  timed_message(paste0("  Total time: ", round(total_time, 2), " minutes"))
  
  gc_result <- gc()
  ram_usage <- format(sum(gc_result[,2]) * 1048576, units = "auto", digits = 2)
  timed_message(paste0("  RAM usage: ", ram_usage))
  
  # Temp files summary
  if (!cleanup_temp) {
    timed_message("\nℹ️  Temporary files were kept (cleanup_temp = FALSE):")
    temp_dirs <- c(
      paste0("temp_fullres_", image_name),
      paste0("temp_downsampled_", image_name),
      paste0("temp_enhanced_", image_name),
      paste0("temp_upsampled_", image_name)
    )
    for (temp_dir in temp_dirs) {
      full_path <- file.path(PROG_DIR, temp_dir)
      if (dir.exists(full_path)) {
        n_files <- length(list.files(full_path))
        if (n_files > 0) {
          size_mb <- sum(file.info(list.files(full_path, full.names = TRUE))$size, na.rm = TRUE) / 1e6
          timed_message(paste0("  • ", temp_dir, ": ", n_files, " files (", round(size_mb, 1), " MB)"))
        }
      }
    }
    timed_message("  To remove: Set cleanup_temp = TRUE or manually delete directories")
  }
  
  # Compile timing summary
  timing_summary <- NULL
  if (timing) {
    timing_summary <- c(
      setup = step_times$setup,
      enhance_all = step_times$enhance_all * 60,
      visualization = step_times$viz,
      quantification = step_times$quantify * 60,
      save = step_times$save,
      total = total_time * 60
    )
  }
  
  # Prepare TIFF directory path for return value
  tiff_dir_path <- if (export_tiffs) {
    file.path(image_res_dir, "enhanced_tiffs")
  } else {
    NULL
  }
  
  return(invisible(list(
    spe = cells,
    enhanced_images = NULL,  # Freed after quantification; load from enhanced_images_file if needed
    viz_file = viz_file,
    spe_file = spe_file,
    enhanced_images_file = enhanced_images_file,
    tiff_dir = tiff_dir_path,
    timing = timing_summary,
    parameters = enhancement_params
  )))
}

