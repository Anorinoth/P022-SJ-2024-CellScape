#' Delete Intermediate Files for Specific Channels
#' Author: Dr. Michael Zaiken
#'
#' Selectively deletes temporary processing files for specific channels to force
#' re-processing with updated settings. Allows fine-grained control over which
#' processing steps to redo.
#'
#' @param channels_to_reprocess Numeric vector of channel indices (e.g., c(1, 2, 3))
#' @param image_name Character string with the name of the image
#' @param steps_to_redo Named list specifying which steps to redo:
#'   \item{fullres}{Logical, delete full-resolution channel files (default FALSE)}
#'   \item{downsampled}{Logical, delete downsampled channel files (default FALSE)}
#'   \item{enhanced}{Logical, delete enhanced channel files (default TRUE)}
#'   \item{upsampled}{Logical, delete upsampled channel files (default TRUE)}
#'   \item{comparison}{Logical, delete comparison PNG files (default TRUE)}
#' @param prog_dir Character string with path to PROG directory containing temp files
#' @param res_dir Character string with path to RES directory containing comparison images
#' @param verbose Logical, whether to print progress messages (default TRUE)
#'
#' @return Invisible list with:
#'   \item{channels_deleted}{Vector of channel indices that had files deleted}
#'   \item{files_deleted}{Total count of files deleted}
#'   \item{steps_processed}{Named vector showing counts per step}
#'
#' @examples
#' # Delete enhancement, upsampling, and comparison files for channels 1-3
#' delete_channel_intermediates(
#'   channels_to_reprocess = c(1, 2, 3),
#'   image_name = "CAR_TREG_2",
#'   steps_to_redo = list(
#'     fullres = FALSE,
#'     downsampled = FALSE,
#'     enhanced = TRUE,
#'     upsampled = TRUE,
#'     comparison = TRUE
#'   ),
#'   prog_dir = PROG_DIR,
#'   res_dir = RES_DIR
#' )
#'
#' @export
delete_channel_intermediates <- function(channels_to_reprocess,
                                          image_name,
                                          steps_to_redo = list(
                                            fullres = FALSE,
                                            downsampled = FALSE,
                                            enhanced = TRUE,
                                            upsampled = TRUE,
                                            comparison = TRUE
                                          ),
                                          prog_dir,
                                          res_dir,
                                          verbose = TRUE) {
  
  # Validate inputs
  if (missing(channels_to_reprocess) || length(channels_to_reprocess) == 0) {
    stop("'channels_to_reprocess' must be a non-empty numeric vector")
  }
  if (missing(image_name) || !is.character(image_name)) {
    stop("'image_name' must be a character string")
  }
  if (missing(prog_dir) || !dir.exists(prog_dir)) {
    stop("'prog_dir' must be a valid directory path")
  }
  if (missing(res_dir) || !dir.exists(res_dir)) {
    stop("'res_dir' must be a valid directory path")
  }
  
  # Ensure steps_to_redo has all required fields
  default_steps <- list(
    fullres = FALSE,
    downsampled = FALSE,
    enhanced = TRUE,
    upsampled = TRUE,
    comparison = TRUE
  )
  
  for (step_name in names(default_steps)) {
    if (is.null(steps_to_redo[[step_name]])) {
      steps_to_redo[[step_name]] <- default_steps[[step_name]]
    }
  }
  
  # Print header
  if (verbose) {
    message("\n=== Deleting Temp Files for Specific Channels ===")
    message(paste0("  Channels: ", paste(sprintf("%02d", channels_to_reprocess), collapse = ", ")))
    message("  Steps to redo:")
  }
  
  # Build list of temp directories to clean based on config
  temp_dirs_to_clean <- character(0)
  
  if (steps_to_redo$fullres) {
    temp_dirs_to_clean <- c(temp_dirs_to_clean, paste0("temp_fullres_", image_name))
    if (verbose) message("    ✓ Full-resolution")
  }
  if (steps_to_redo$downsampled) {
    temp_dirs_to_clean <- c(temp_dirs_to_clean, paste0("temp_downsampled_", image_name))
    if (verbose) message("    ✓ Downsampled")
  }
  if (steps_to_redo$enhanced) {
    temp_dirs_to_clean <- c(temp_dirs_to_clean, paste0("temp_enhanced_", image_name))
    if (verbose) message("    ✓ Enhanced")
  }
  if (steps_to_redo$upsampled) {
    temp_dirs_to_clean <- c(temp_dirs_to_clean, paste0("temp_upsampled_", image_name))
    if (verbose) message("    ✓ Upsampled")
  }
  if (steps_to_redo$comparison) {
    if (verbose) message("    ✓ Comparison images")
  }
  
  if (verbose) message("")
  
  # Track deletion statistics
  total_files_deleted <- 0
  step_counts <- list(
    fullres = 0,
    downsampled = 0,
    enhanced = 0,
    upsampled = 0,
    comparison = 0
  )
  channels_with_deletions <- integer(0)
  
  # Delete channel files from selected temp directories
  for (ch_idx in channels_to_reprocess) {
    ch_file <- paste0("channel_", sprintf("%02d", ch_idx), ".rds")
    if (verbose) message(paste0("  Channel ", sprintf("%02d", ch_idx), ":"))
    
    deleted_count <- 0
    
    # Delete from temp directories
    for (temp_dir in temp_dirs_to_clean) {
      full_path <- file.path(prog_dir, temp_dir, ch_file)
      if (file.exists(full_path)) {
        file.remove(full_path)
        if (verbose) message(paste0("    ✓ Removed from ", basename(temp_dir)))
        deleted_count <- deleted_count + 1
        total_files_deleted <- total_files_deleted + 1
        
        # Track which step
        if (grepl("fullres", temp_dir)) step_counts$fullres <- step_counts$fullres + 1
        if (grepl("downsampled", temp_dir)) step_counts$downsampled <- step_counts$downsampled + 1
        if (grepl("enhanced", temp_dir)) step_counts$enhanced <- step_counts$enhanced + 1
        if (grepl("upsampled", temp_dir)) step_counts$upsampled <- step_counts$upsampled + 1
      }
    }
    
    # Delete comparison images if requested
    if (steps_to_redo$comparison) {
      comparison_dir <- file.path(res_dir, image_name, "enhancement_comparison")
      if (dir.exists(comparison_dir)) {
        comparison_files <- list.files(comparison_dir, 
                                       pattern = paste0("^comparison_ch", sprintf("%02d", ch_idx), "_.*\\.png$"),
                                       full.names = TRUE)
        for (comp_file in comparison_files) {
          file.remove(comp_file)
          if (verbose) message(paste0("    ✓ Removed ", basename(comp_file)))
          deleted_count <- deleted_count + 1
          total_files_deleted <- total_files_deleted + 1
          step_counts$comparison <- step_counts$comparison + 1
        }
      }
    }
    
    if (deleted_count == 0) {
      if (verbose) message("    (no files found to delete)")
    } else {
      channels_with_deletions <- c(channels_with_deletions, ch_idx)
    }
  }
  
  # Print summary
  if (verbose) {
    message("\n✓ Temp files deleted - selected steps will be re-processed")
    message(paste0("  Total files deleted: ", total_files_deleted))
    if (total_files_deleted > 0) {
      message("  Breakdown:")
      for (step_name in names(step_counts)) {
        if (step_counts[[step_name]] > 0) {
          message(paste0("    - ", step_name, ": ", step_counts[[step_name]]))
        }
      }
    }
  }
  
  # Return invisibly
  return(invisible(list(
    channels_deleted = channels_with_deletions,
    files_deleted = total_files_deleted,
    steps_processed = unlist(step_counts)
  )))
}

