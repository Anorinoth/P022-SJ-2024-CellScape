#' R Wrapper for Python merge_channels_to_ome.py Script
#'
#' This function calls the Python script to merge single-channel TIFFs
#' into a multi-channel OME-TIFF for QuPath v6.
#'
#' @param tiff_dir Character string, path to directory containing single-channel TIFFs
#' @param output_file Character string, path for output OME-TIFF file
#' @param python_path Character string, path to Python executable (default: "python3")
#'   Use "python" on Windows or specify full path if needed
#' @param compression Character string, compression type: "lzw", "zlib", or "none" (default: "lzw")
#' @param check_dependencies Logical, whether to check for required Python packages (default: TRUE)
#'
#' @return Invisible list with:
#'   \item{success}{Logical, whether merge was successful}
#'   \item{output_file}{Path to created OME-TIFF}
#'   \item{exit_code}{Exit code from Python script}
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' merge_tiffs_to_ome(
#'   tiff_dir = "MAIN/RES/sample_01/enhanced_tiffs",
#'   output_file = "MAIN/RES/sample_01/sample_01_multichannel.ome.tiff"
#' )
#'
#' # Specify Python path
#' merge_tiffs_to_ome(
#'   tiff_dir = "MAIN/RES/sample_01/enhanced_tiffs",
#'   output_file = "MAIN/RES/sample_01/sample_01_multichannel.ome.tiff",
#'   python_path = "/usr/bin/python3"
#' )
#' }
#'
#' @export
merge_tiffs_to_ome <- function(tiff_dir,
                                output_file,
                                python_path = "python3",
                                compression = "lzw",
                                check_dependencies = TRUE) {
  
  # Validate inputs
  if (!dir.exists(tiff_dir)) {
    stop("TIFF directory does not exist: ", tiff_dir)
  }
  
  # Get script path (relative to this file)
  script_dir <- dirname(sys.frame(1)$ofile)
  if (is.null(script_dir) || script_dir == "") {
    # Fallback: use current working directory
    script_path <- file.path(getwd(), "merge_channels_to_ome.py")
  } else {
    script_path <- file.path(script_dir, "merge_channels_to_ome.py")
  }
  
  if (!file.exists(script_path)) {
    stop("Python script not found: ", script_path,
         "\nExpected location: MAIN/PROG/Scripts/merge_channels_to_ome.py")
  }
  
  # Check Python installation
  python_check <- system2(python_path, args = "--version", 
                          stdout = TRUE, stderr = TRUE)
  
  if (!is.null(attr(python_check, "status")) && attr(python_check, "status") != 0) {
    stop("Python not found at: ", python_path,
         "\nPlease install Python 3 or specify correct path with python_path argument")
  }
  
  message("Python version: ", python_check[1])
  
  # Check for required Python packages
  if (check_dependencies) {
    message("Checking Python dependencies...")
    
    # Check for tifffile
    tifffile_check <- system2(python_path, 
                              args = c("-c", shQuote("import tifffile; print(tifffile.__version__)")),
                              stdout = TRUE, stderr = TRUE)
    
    if (!is.null(attr(tifffile_check, "status")) && attr(tifffile_check, "status") != 0) {
      stop("Required Python package 'tifffile' not found.\n",
           "Install with: ", python_path, " -m pip install tifffile numpy")
    }
    
    message("  ✓ tifffile version: ", tifffile_check[1])
    
    # Check for numpy
    numpy_check <- system2(python_path,
                          args = c("-c", shQuote("import numpy; print(numpy.__version__)")),
                          stdout = TRUE, stderr = TRUE)
    
    if (!is.null(attr(numpy_check, "status")) && attr(numpy_check, "status") != 0) {
      stop("Required Python package 'numpy' not found.\n",
           "Install with: ", python_path, " -m pip install numpy")
    }
    
    message("  ✓ numpy version: ", numpy_check[1])
  }
  
  # Prepare command
  message("\n=== Merging Single-Channel TIFFs to OME-TIFF ===")
  message("Input directory: ", tiff_dir)
  message("Output file: ", output_file)
  message("Compression: ", compression)
  
  # Build command arguments
  args <- c(
    script_path,
    shQuote(tiff_dir),
    shQuote(output_file),
    compression
  )
  
  # Run Python script
  message("\nRunning Python script...")
  exit_code <- system2(python_path, args = args, stdout = "", stderr = "")
  
  # Check result
  success <- (exit_code == 0)
  
  if (success) {
    if (file.exists(output_file)) {
      file_size <- file.info(output_file)$size / 1e9
      message("\n✓ SUCCESS: Multi-channel OME-TIFF created")
      message("  File: ", output_file)
      message("  Size: ", round(file_size, 2), " GB")
    } else {
      success <- FALSE
      warning("Python script reported success but output file not found: ", output_file)
    }
  } else {
    warning("Python script failed with exit code: ", exit_code)
  }
  
  return(invisible(list(
    success = success,
    output_file = if(success) output_file else NULL,
    exit_code = exit_code
  )))
}

