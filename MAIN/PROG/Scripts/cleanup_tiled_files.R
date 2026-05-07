################################################################################
# cleanup_tiled_files.R
# Author: Dr. Michael Zaiken
#
# Utility functions for cleaning up intermediate tiled segmentation files
# Allows selective deletion based on what needs to be re-run
################################################################################

#' Clean up tiled segmentation intermediate files
#'
#' @param tile_dir Character. Path to the tile directory (e.g., "MAIN/PROG/tiles/CAR_TREG_2")
#' @param image_name Character. Name of the image (e.g., "CAR_TREG_2")
#' @param prog_dir Character. Path to the PROG directory (e.g., "MAIN/PROG")
#' @param cleanup_level Character. One of:
#'   - "quantification": Remove only quantification outputs (SPE files, quant flag, merged SPE)
#'   - "segmentation": Remove segmentation AND quantification (masks, SPE files, visualizations, flags, merged SPE)
#'   - "all": Remove EVERYTHING including tile images (full reset)
#'
#' @return Invisible NULL. Prints messages about deleted files.
#'
#' @details
#' Cleanup levels:
#' 
#' **"quantification"** - Use to re-test enhancement/quantification parameters
#'   Removes:
#'     - SpatialExperiment objects (*_spe.rds)
#'     - Enhancement visualizations (pre/post comparisons)
#'     - Quantification completion flag
#'     - Final merged SPE file
#'   Keeps:
#'     - Tile images (*_image.rds, *_nuclear.rds)
#'     - Segmentation masks (*_mask.rds)
#'     - Segmentation visualizations (overlays)
#'     - All earlier flags
#'
#' **"segmentation"** - Use to re-test segmentation parameters (gamma, watershed, etc.)
#'   Removes:
#'     - Everything from "quantification" level PLUS
#'     - Segmentation masks (*_mask.rds)
#'     - Segmentation visualizations (overlays)
#'     - Segmentation completion flag
#'   Keeps:
#'     - Tile images (*_image.rds, *_nuclear.rds)
#'     - Deduplication checkpoint
#'     - Tile grid metadata
#'     - Tile extraction flag
#'
#' **"all"** - Complete cleanup (CAUTION!)
#'   Removes:
#'     - ALL tile files and directories
#'     - All checkpoints and flags
#'     - All metadata
#'   Use only after confirming results are final
#'
#' @examples
#' \dontrun{
#' # Re-test quantification only
#' cleanup_tiled_files(
#'   tile_dir = "MAIN/PROG/tiles/CAR_TREG_2",
#'   image_name = "CAR_TREG_2",
#'   prog_dir = "MAIN/PROG",
#'   cleanup_level = "quantification"
#' )
#'
#' # Re-test segmentation (will also remove quantification)
#' cleanup_tiled_files(
#'   tile_dir = "MAIN/PROG/tiles/CAR_TREG_2",
#'   image_name = "CAR_TREG_2",
#'   prog_dir = "MAIN/PROG",
#'   cleanup_level = "segmentation"
#' )
#' }
cleanup_tiled_files <- function(tile_dir, 
                                image_name, 
                                prog_dir, 
                                cleanup_level = c("quantification", "segmentation", "all")) {
  
  # Validate cleanup level
  cleanup_level <- match.arg(cleanup_level)
  
  # Check tile directory exists
  if (!dir.exists(tile_dir)) {
    stop("Tile directory not found: ", tile_dir)
  }
  
  message("\n=== Cleaning up tiled files ===")
  message("  Tile directory: ", tile_dir)
  message("  Cleanup level: ", cleanup_level)
  message("")
  
  # Track what gets deleted
  deleted_count <- list(
    spe = 0,
    masks = 0,
    viz_files = 0,
    flags = 0,
    tiles = 0,
    metadata = 0,
    dirs = 0
  )
  
  # LEVEL 1: Remove quantification outputs (SPE files, visualizations, flags)
  if (cleanup_level %in% c("quantification", "segmentation", "all")) {
    message("  Removing quantification outputs...")
    
    # Delete SPE files
    spe_files <- list.files(tile_dir, pattern = "*_spe\\.rds$", full.names = TRUE)
    if (length(spe_files) > 0) {
      unlink(spe_files)
      deleted_count$spe <- length(spe_files)
      message("    ✓ Deleted ", length(spe_files), " SPE files")
    }
    
    # Delete enhancement visualizations
    viz_dir <- file.path(tile_dir, "visualizations")
    if (dir.exists(viz_dir)) {
      enhancement_viz <- list.files(viz_dir, pattern = "enhancement_comparison", full.names = TRUE)
      if (length(enhancement_viz) > 0) {
        unlink(enhancement_viz)
        deleted_count$viz_files <- deleted_count$viz_files + length(enhancement_viz)
        message("    ✓ Deleted ", length(enhancement_viz), " enhancement visualizations")
      }
    }
    
    # Delete quantification completion flag
    quant_flag <- file.path(tile_dir, "quantification_complete.flag")
    if (file.exists(quant_flag)) {
      unlink(quant_flag)
      deleted_count$flags <- deleted_count$flags + 1
      message("    ✓ Deleted quantification_complete.flag")
    }
    
    # Delete final merged SPE file
    spe_file <- file.path(prog_dir, paste0(image_name, "_cells.qs"))
    if (file.exists(spe_file)) {
      unlink(spe_file)
      message("    ✓ Deleted final merged SPE file")
    }
    
    # Delete cached result file
    result_cache <- file.path(prog_dir, paste0(image_name, "_tiled_result.rds"))
    if (file.exists(result_cache)) {
      unlink(result_cache)
      message("    ✓ Deleted cached result file")
    }
  }
  
  # LEVEL 2: Remove segmentation outputs (masks, overlays, flags)
  if (cleanup_level %in% c("segmentation", "all")) {
    message("  Removing segmentation outputs...")
    
    # Delete masks
    mask_files <- list.files(tile_dir, pattern = "*_mask\\.rds$", full.names = TRUE)
    if (length(mask_files) > 0) {
      unlink(mask_files)
      deleted_count$masks <- length(mask_files)
      message("    ✓ Deleted ", length(mask_files), " mask files")
    }
    
    # Delete segmentation visualizations (overlays, diagnostics CSV)
    viz_dir <- file.path(tile_dir, "visualizations")
    if (dir.exists(viz_dir)) {
      seg_viz <- list.files(viz_dir, pattern = "segmentation_overlay|cell_positions\\.png", full.names = TRUE)
      if (length(seg_viz) > 0) {
        unlink(seg_viz)
        deleted_count$viz_files <- deleted_count$viz_files + length(seg_viz)
        message("    ✓ Deleted ", length(seg_viz), " segmentation visualizations")
      }
    }
    
    # Delete segmentation diagnostics CSV
    diag_csv <- file.path(tile_dir, "segmentation_diagnostics.csv")
    if (file.exists(diag_csv)) {
      unlink(diag_csv)
      message("    ✓ Deleted segmentation_diagnostics.csv")
    }
    
    # Delete segmentation completion flag
    seg_flag <- file.path(tile_dir, "segmentation_complete.flag")
    if (file.exists(seg_flag)) {
      unlink(seg_flag)
      deleted_count$flags <- deleted_count$flags + 1
      message("    ✓ Deleted segmentation_complete.flag")
    }
  }
  
  # LEVEL 3: Remove EVERYTHING (tile images, metadata, all flags)
  if (cleanup_level == "all") {
    message("  Removing ALL tile files and metadata...")
    
    # Calculate total size before deletion
    all_files <- list.files(tile_dir, recursive = TRUE, full.names = TRUE)
    total_size_gb <- sum(file.info(all_files)$size, na.rm = TRUE) / 1e9
    
    message("    Total files: ", length(all_files))
    message("    Disk space to free: ", round(total_size_gb, 2), " GB")
    
    # Delete entire tile directory
    unlink(tile_dir, recursive = TRUE)
    deleted_count$dirs <- 1
    
    message("    ✓ Deleted entire tile directory")
  }
  
  # Summary
  message("\n=== Cleanup Summary ===")
  if (cleanup_level == "all") {
    message("  ✓ Complete cleanup - all tile files deleted")
    message("  Next run will start from Stage 0 (deduplication)")
  } else if (cleanup_level == "segmentation") {
    message("  ✓ Segmentation reset complete")
    message("  Deleted: ", deleted_count$masks, " masks, ", 
            deleted_count$spe, " SPE files, ",
            deleted_count$viz_files, " visualizations")
    
    # Count remaining tile files
    tile_image_files <- list.files(tile_dir, pattern = "*_image\\.rds$", full.names = TRUE)
    tile_nuclear_files <- list.files(tile_dir, pattern = "*_nuclear\\.rds$", full.names = TRUE)
    message("  Kept: ", length(tile_image_files), " image tiles, ", 
            length(tile_nuclear_files), " nuclear tiles")
    message("  Next run will start from Stage 2 (segmentation)")
  } else if (cleanup_level == "quantification") {
    message("  ✓ Quantification reset complete")
    message("  Deleted: ", deleted_count$spe, " SPE files, ",
            deleted_count$viz_files, " enhancement visualizations")
    
    # Count remaining masks and tiles
    mask_files <- list.files(tile_dir, pattern = "*_mask\\.rds$", full.names = TRUE)
    tile_image_files <- list.files(tile_dir, pattern = "*_image\\.rds$", full.names = TRUE)
    message("  Kept: ", length(mask_files), " masks, ",
            length(tile_image_files), " image tiles")
    message("  Next run will start from Stage 3 (quantification)")
  }
  
  message("")
  
  invisible(NULL)
}

