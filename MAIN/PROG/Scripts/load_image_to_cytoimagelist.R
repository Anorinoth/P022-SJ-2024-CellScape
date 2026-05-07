#' Load Multi-Channel Image as CytoImageList (On-Disk)
#'
#' This function loads single-channel TIFF files, combines them into a multi-channel
#' HDF5-backed image, and creates a CytoImageList object for use with cytomapper.
#' The entire operation is memory-efficient, processing one channel at a time.
#'
#' @param tiff_dir Character, path to directory containing single-channel TIFF files
#' @param image_name Character, name for the image (used for HDF5 dataset name)
#' @param h5_output_dir Character, directory where H5 file will be saved
#' @param channel_names Character vector, names for each channel (optional, inferred from filenames if NULL)
#' @param force_rebuild Logical, if TRUE, rebuild H5 file even if it exists. Default: FALSE
#' @param BPPARAM BiocParallel parameter object for parallel channel writing. Default: SerialParam()
#' @param verbose Logical, if TRUE, display detailed progress messages. Default: TRUE
#'
#' @return A CytoImageList object containing the multi-channel image backed by HDF5 on disk
#'
#' @details
#' The function performs these steps:
#' 1. Checks if H5 file already exists (skips creation if found and force_rebuild=FALSE)
#' 2. Reads TIFF files one at a time from tiff_dir
#' 3. Creates HDF5 file with optimized chunking (stays under 4GB chunk limit)
#' 4. Writes each channel directly to H5 file (memory-efficient)
#' 5. Loads the H5 file as a DelayedArray
#' 6. Creates CytoImageList with proper metadata
#'
#' @note
#' - Requires packages: rhdf5, HDF5Array, cytomapper, tiff, S4Vectors
#' - Peak memory usage: ~2-3 GB (one channel at a time)
#' - Chunk size automatically optimized to stay under HDF5's 4GB limit
#' - Full spatial width preserved for efficient access
#'
#' @examples
#' \dontrun{
#' images <- load_image_to_cytoimagelist(
#'   tiff_dir = "MAIN/PROG/SingleChannelImages/CAR_TREG_2",
#'   image_name = "CAR_TREG_2",
#'   h5_output_dir = "MAIN/PROG/cytomapper_h5"
#' )
#' 
#' # Verify
#' length(channelNames(images))  # Number of channels
#' dim(images[[1]])               # Image dimensions
#' }
#'
#' @export
load_image_to_cytoimagelist <- function(tiff_dir, 
                                         image_name, 
                                         h5_output_dir,
                                         channel_names = NULL,
                                         force_rebuild = FALSE,
                                         BPPARAM = BiocParallel::SerialParam(),
                                         verbose = TRUE) {
  
  # Validate inputs
  if (!dir.exists(tiff_dir)) {
    stop("TIFF directory does not exist: ", tiff_dir)
  }
  
  if (!is.character(image_name) || length(image_name) != 1) {
    stop("image_name must be a single character string")
  }
  
  # Create output directory
  dir.create(h5_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Define H5 file path
  final_h5_path <- file.path(h5_output_dir, paste0(image_name, ".h5"))
  
  message("\n=== Loading Image: ", image_name, " ===")
  message("TIFF directory: ", tiff_dir)
  message("H5 output: ", h5_output_dir)
  
  # Check if already processed
  if (file.exists(final_h5_path) && !force_rebuild) {
    message("✓ Found existing H5 file: ", basename(final_h5_path))
    message("  Loading from disk (use force_rebuild=TRUE to regenerate)")
    multi_channel_image <- HDF5Array::HDF5Array(final_h5_path, name = image_name)
    
  } else {
    # If force_rebuild is TRUE, delete any existing file
    if (force_rebuild && file.exists(final_h5_path)) {
      message("Removing existing H5 file for rebuild...")
      H5close()
      Sys.sleep(0.2)
      unlink(final_h5_path, force = TRUE)
      Sys.sleep(0.2)
      
      if (file.exists(final_h5_path)) {
        stop("Failed to delete existing H5 file: ", final_h5_path)
      }
    }
    
    message("Creating on-disk multi-channel image...")
    
    # Get list of TIFF files
    tiff_files <- list.files(tiff_dir, pattern = "\\.tif$", full.names = TRUE)
    
    if (length(tiff_files) == 0) {
      stop("No TIFF files found in directory: ", tiff_dir)
    }
    
    tiff_files <- sort(tiff_files)  # Ensure correct order
    n_channels <- length(tiff_files)
    
    message("  Found ", n_channels, " TIFF files")
    message("  Processing channels ONE AT A TIME (memory-efficient)")
    
    # Read first channel to get dimensions
    message("  Reading first channel to determine dimensions...")
    first_img <- tiff::readTIFF(tiff_files[1], native = FALSE, as.is = TRUE)
    img_height <- nrow(first_img)
    img_width <- ncol(first_img)
    
    if (verbose) {
      message("  Image dimensions: ", img_height, " x ", img_width, " x ", n_channels)
      message("  Peak RAM per channel: ~", round((img_height * img_width * 8) / 1e9, 2), " GB")
    }
    
    # Load rhdf5
    if (!requireNamespace("rhdf5", quietly = TRUE)) {
      stop("rhdf5 package required. Install with: BiocManager::install('rhdf5')")
    }
    library(rhdf5)
    
    # Create HDF5 file
    h5createFile(final_h5_path)
    
    # Calculate chunk dimensions to stay under HDF5's 4GB limit
    # Each chunk: rows × width × 1 channel × 8 bytes < 4GB
    # 4GB / (width × 8) gives max rows
    max_chunk_rows <- floor(4e9 / (img_width * 8))
    chunk_rows <- min(10000, max_chunk_rows, img_height)
    chunk_dims <- c(chunk_rows, img_width, 1)
    
    # Calculate total chunks needed
    chunks_per_channel_rows <- ceiling(img_height / chunk_rows)
    total_chunks <- chunks_per_channel_rows * n_channels
    
    chunk_size_gb <- (chunk_rows * img_width * 1 * 8) / 1e9
    
    # Calculate estimated file size (avoiding integer overflow by using numeric)
    est_file_size_gb <- (as.numeric(img_height) * as.numeric(img_width) * as.numeric(n_channels) * 8) / 1e9 * 0.7
    
    if (verbose) {
      message("  Chunk dimensions: ", paste(chunk_dims, collapse = " x "))
      message("  Chunk size: ", round(chunk_size_gb, 2), " GB (under 4GB limit ✓)")
      message("  Chunks per channel: ", chunks_per_channel_rows, " row chunks")
      message("  Total chunks: ", total_chunks, " (", chunks_per_channel_rows, " × ", n_channels, " channels)")
      message("  Estimated file size: ", round(est_file_size_gb, 2), " GB (compressed)")
    }
    
    # Create HDF5 file and dataset
    # Note: h5createFile returns FALSE if file exists, so we check file existence instead
    h5createFile(final_h5_path)
    
    # Verify file was actually created
    if (!file.exists(final_h5_path)) {
      stop("Failed to create H5 file (file does not exist after creation): ", final_h5_path)
    }
    
    if (verbose) {
      message("  ✓ H5 file created: ", basename(final_h5_path))
    }
    
    # Create dataset
    dataset_created <- h5createDataset(
      file = final_h5_path,
      dataset = image_name,
      dims = c(img_height, img_width, n_channels),
      chunk = chunk_dims,
      level = 6,
      storage.mode = "double"
    )
    
    if (verbose) {
      message("  Dataset creation returned: ", dataset_created)
    }
    
    if (!dataset_created) {
      # Check if dataset already exists
      existing_datasets <- h5ls(final_h5_path)$name
      stop("Failed to create dataset '", image_name, "' in H5 file: ", final_h5_path, 
           "\nDataset creation returned FALSE.",
           "\nExisting datasets in file: ", paste(existing_datasets, collapse = ", "))
    }
    
    if (verbose) {
      message("  ✓ Dataset created: ", image_name)
    }
    
    # Close any open HDF5 handles before forking (important for parallel processing)
    H5close()
    
    # Determine if parallel processing should be used
    n_workers <- BiocParallel::bpnworkers(BPPARAM)
    use_parallel <- n_workers > 1 && !inherits(BPPARAM, "SerialParam")
    
    # DEBUG: Show what BPPARAM we actually received
    if (verbose) {
      message("  DEBUG: BPPARAM class = ", class(BPPARAM)[1])
      message("  DEBUG: BPPARAM workers = ", n_workers)
    }
    
    # Check parallelization type
    is_fork <- inherits(BPPARAM, "MulticoreParam")
    is_socket <- inherits(BPPARAM, "SnowParam")
    
    if (verbose) {
      message("  DEBUG: is_fork = ", is_fork, ", is_socket = ", is_socket)
    }
    
    # ENFORCE: Socket-based parallelization is NOT allowed for HDF5 writes
    # It's too slow and causes SIGPIPE errors with large data
    # Check for SnowParam specifically, but NOT if it's MulticoreParam (which may inherit from it)
    if (use_parallel && is_socket && !is_fork) {
      stop("\n",
           "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n",
           "ERROR: Socket-based parallelization is NOT supported for HDF5 writes\n",
           "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n",
           "Reason: Socket parallelization causes SIGPIPE errors with large images.\n\n",
           "SOLUTION: Switch to FORK-based parallelization (fast and reliable)\n\n",
           "Run these commands in your R session:\n\n",
           "  # Configure fork-based parallelization with ", n_workers, " workers\n",
           "  BPPARAM <- BiocParallel::MulticoreParam(workers = ", n_workers, ")\n\n",
           "  # Re-run the image loading\n",
           "  images <- load_image_to_cytoimagelist(\n",
           "    tiff_dir = result$tiff_dir,\n",
           "    image_name = result$image_name,\n",
           "    h5_output_dir = dirname(result$h5_dir),\n",
           "    channel_names = result$channel_names,\n",
           "    force_rebuild = TRUE,\n",
           "    BPPARAM = BPPARAM,\n",
           "    verbose = TRUE\n",
           "  )\n\n",
           "Alternative: Use serial processing (slower):\n",
           "  BPPARAM <- BiocParallel::SerialParam()\n\n",
           "Current BPPARAM type: ", class(BPPARAM)[1], "\n",
           "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
    }
    
    if (verbose) {
      if (use_parallel) {
        if (is_fork) {
          message("  Writing channels in parallel using ", n_workers, " workers (FORK-based)...")
          message("  ✓ Fork-based parallelization is FAST and reliable for HDF5 writes")
        } else {
          message("  Writing channels in parallel using ", n_workers, " workers...")
        }
        est_time_min <- ceiling((n_channels / n_workers) * 0.3)  # ~18 seconds per channel
        message("  Estimated time: ", est_time_min, "-", est_time_min * 2, " minutes")
      } else {
        message("  Writing channels serially...")
        est_time_min <- ceiling(n_channels * 0.3)  # ~18 seconds per channel serially
        message("  Estimated time: ", est_time_min, "-", est_time_min * 2, " minutes")
      }
    }
    
    start_time <- Sys.time()
    
    if (use_parallel) {
      # Parallel processing: split channels across workers
      # Fork-based parallelization is fast and reliable for HDF5 writes
      message("  Starting parallel channel writing...")
      
      write_channel <- function(i, tiff_files, final_h5_path, image_name) {
        # Load rhdf5 (fork inherits parent's loaded packages, but explicit load is safer)
        library(rhdf5)
        
        # Read one channel from TIFF
        ch_img <- tiff::readTIFF(tiff_files[i], native = FALSE, as.is = TRUE)
        
        # HDF5 file handles don't survive fork() - each worker must open independently
        # Use low-level H5 operations to ensure proper file handle management
        # Open file in read-write mode
        h5_fid <- H5Fopen(final_h5_path, flags = "H5F_ACC_RDWR")
        
        # Open dataset
        h5_did <- H5Dopen(h5_fid, image_name)
        
        # Create dataspace for the slice we're writing to
        # Write to slice [,, i] (all rows, all cols, channel i)
        h5_sid <- H5Dget_space(h5_did)
        dims <- H5Sget_simple_extent_dims(h5_sid)$size
        
        # Select hyperslab for this channel
        H5Sselect_hyperslab(h5_sid, 
                           start = c(0, 0, i-1),  # 0-indexed
                           count = c(dims[1], dims[2], 1))
        
        # Create memory dataspace
        h5_mid <- H5Screate_simple(c(dims[1], dims[2], 1))
        
        # Write data
        H5Dwrite(h5_did, ch_img, h5spaceMem = h5_mid, h5spaceFile = h5_sid)
        
        # Close all handles
        H5Sclose(h5_mid)
        H5Sclose(h5_sid)
        H5Dclose(h5_did)
        H5Fclose(h5_fid)
        
        # Clean up
        rm(ch_img)
        gc(verbose = FALSE)
        
        return(i)
      }
      
      # Process channels in parallel with fork-based parallelization
      results <- BiocParallel::bplapply(
        1:n_channels,
        write_channel,
        tiff_files = tiff_files,
        final_h5_path = final_h5_path,
        image_name = image_name,
        BPPARAM = BPPARAM
      )
      
      if (verbose) {
        message("  ✓ All ", n_channels, " channels written in parallel")
      }
      
    } else {
      # Serial processing
      for (i in 1:n_channels) {
        if (verbose && (i %% 5 == 0 || i == n_channels)) {
          elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
          pct_complete <- (i / n_channels) * 100
          est_total <- elapsed / i * n_channels
          est_remaining <- est_total - elapsed
          
          message(sprintf("    Channel %d/%d (%.1f%%) - Elapsed: %.1f min, Est. remaining: %.1f min", 
                          i, n_channels, pct_complete, elapsed, est_remaining))
        }
        
        # Read one channel
        ch_img <- tiff::readTIFF(tiff_files[i], native = FALSE, as.is = TRUE)
        
        # Write directly to HDF5 file at the correct position
        h5write(
          obj = ch_img,
          file = final_h5_path,
          name = image_name,
          index = list(NULL, NULL, i)
        )
        
        rm(ch_img)
        if (i %% 10 == 0) gc(verbose = FALSE)
      }
    }
    
    # Close HDF5 file
    H5close()
    
    total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    
    if (verbose) {
      message("  ✓ Written to: ", basename(final_h5_path))
      message("  File size: ", format(file.info(final_h5_path)$size / 1e9, digits = 2), " GB")
      message("  Total time: ", round(total_time, 1), " minutes")
    }
    
    # Load back as HDF5Array
    multi_channel_image <- HDF5Array::HDF5Array(final_h5_path, name = image_name)
  }
  
  if (verbose) {
    message("✓ Multi-channel image loaded (on-disk)")
    message("  Dimensions: ", paste(dim(multi_channel_image), collapse = " x "))
  }
  
  # Extract or use provided channel names
  if (is.null(channel_names)) {
    # Extract from filenames (remove extension and path)
    tiff_files <- list.files(tiff_dir, pattern = "\\.tif$")
    tiff_files <- sort(tiff_files)
    channel_names <- gsub("\\.tif$", "", tiff_files)
    message("  Channel names extracted from filenames")
  }
  
  if (length(channel_names) != dim(multi_channel_image)[3]) {
    warning("Number of channel names (", length(channel_names), 
            ") does not match number of channels (", dim(multi_channel_image)[3], ")")
  }
  
  # Set dimension names
  dimnames(multi_channel_image) <- list(
    NULL,  # row names
    NULL,  # column names
    channel_names  # channel names (3rd dimension)
  )
  
  if (verbose) {
    message("  Channel names set: ", length(channel_names), " channels")
    message("  First 3: ", paste(head(channel_names, 3), collapse = ", "))
    message("  Last 3: ", paste(tail(channel_names, 3), collapse = ", "))
  }
  
  # Create CytoImageList
  if (verbose) message("\n=== Creating CytoImageList ===")
  
  # Create named list for CytoImageList
  image_list <- list(multi_channel_image)
  names(image_list) <- image_name
  
  images <- cytomapper::CytoImageList(image_list)
  
  # Assign metadata
  mcols(images) <- S4Vectors::DataFrame(imageID = image_name)
  
  if (verbose) {
    message("✓ CytoImageList created successfully")
    message("  Class: ", class(images)[1])
    message("  Number of images: ", length(images))
    message("  Number of channels: ", length(cytomapper::channelNames(images)))
    message("  Image dimensions: ", paste(dim(images[[1]]), collapse = " x "))
  }
  
  return(images)
}

# Auto-run example if sourced directly
if (sys.nframe() == 0) {
  message("load_image_to_cytoimagelist.R sourced")
  message("Usage: images <- load_image_to_cytoimagelist(tiff_dir, image_name, h5_output_dir)")
}

