#' Configure Parallel Processing for BiocParallel
#' Author: Dr. Michael Zaiken
#'
#' Sets up BPPARAM for parallel processing with flexible options for
#' serial, fork-based (MulticoreParam), or socket-based (SnowParam) parallelization.
#'
#' @param n_cores Number of cores to use. Can be:
#'   - Integer: Specific number of cores (e.g., 8)
#'   - "serial": Serial processing (no parallelization)
#'   - "max": Use all available cores
#'   Default: 9 (conservative for I/O operations)
#'
#' @param mode Character, parallelization type:
#'   - "socket": Socket-based parallelization (SnowParam) - Java-compatible, more memory
#'   - "fork": Fork-based parallelization (MulticoreParam) - Unix only, incompatible with Java, faster
#'   Default: "socket"
#'
#' @param for_java Logical, if TRUE and mode="fork", will warn about Java incompatibility.
#'   Default: TRUE.
#'
#' @param verbose Logical, if TRUE, displays configuration messages. Default: TRUE
#'
#' @return A list containing:
#'   - BPPARAM: BiocParallel parameter object
#'   - n_cores: Number of cores configured (integer)
#'   - mode: Parallelization mode used ("socket" or "fork")
#'   - type: Type of parallelization ("serial", "socket", "fork")
#'
#' @details
#' Parallelization Types:
#' 
#' - **Serial** (n_cores="serial"): No parallelization, most stable, slowest
#' - **Socket** (mode="socket"): Creates separate R processes, Java-compatible, uses more memory
#' - **Fork** (mode="fork"): Faster, shared memory, Unix only, INCOMPATIBLE with Java/rJava
#'
#' Guidelines:
#' - For Java/RBioFormats: Use mode="socket" (default)
#' - For non-Java operations: mode="fork" is faster
#' - For debugging: n_cores="serial"
#' - For max speed: n_cores="max" with appropriate mode
#'
#' @examples
#' \dontrun{
#' # Serial processing (no parallelization)
#' config <- configure_parallel_processing(n_cores = "serial")
#'
#' # Maximum parallelization (all cores, socket-based)
#' config <- configure_parallel_processing(n_cores = "max", mode = "socket")
#'
#' # Manual core specification (socket-based for Java)
#' config <- configure_parallel_processing(n_cores = 8, mode = "socket")
#'
#' # Fork-based for non-Java operations (faster)
#' config <- configure_parallel_processing(n_cores = 16, mode = "fork", for_java = FALSE)
#'
#' # Extract BPPARAM for use
#' BPPARAM <- config$BPPARAM
#' }
#'
#' @export
configure_parallel_processing <- function(n_cores = 9,
                                           mode = "socket",
                                           for_java = TRUE,
                                           verbose = TRUE) {
  
  # Validate mode
  valid_modes <- c("socket", "fork")
  if (!mode %in% valid_modes) {
    stop("Invalid mode. Must be one of: ", paste(valid_modes, collapse = ", "))
  }
  
  # Detect system
  is_windows <- .Platform$OS.type == "windows"
  detected_cores <- parallel::detectCores()
  
  if (verbose) {
    message("\n=== Configuring Parallel Processing ===")
    message("  System: ", ifelse(is_windows, "Windows", "Unix"))
    message("  Detected cores: ", detected_cores)
  }
  
  # Parse n_cores parameter
  cores_value <- NULL
  
  if (is.character(n_cores)) {
    # Handle special string values
    if (n_cores == "serial") {
      if (verbose) {
        message("  ✓ Configured: Serial processing (no parallelization)")
      }
      return(list(
        BPPARAM = BiocParallel::SerialParam(),
        n_cores = 1,
        mode = "serial",
        type = "serial"
      ))
    } else if (n_cores == "max") {
      cores_value <- detected_cores
      if (verbose) {
        message("  Cores requested: ALL (max)")
      }
    } else {
      stop("n_cores must be a number, 'serial', or 'max', got: ", n_cores)
    }
  } else if (is.numeric(n_cores)) {
    # Validate numeric input
    if (n_cores < 1) {
      stop("n_cores must be a positive integer, got: ", n_cores)
    }
    cores_value <- as.integer(n_cores)
    if (cores_value > detected_cores) {
      warning("Requested ", cores_value, " cores but only ", detected_cores, 
              " detected. Using ", detected_cores, " cores.")
      cores_value <- detected_cores
    }
    if (verbose) {
      message("  Cores requested: ", cores_value, " (manual)")
    }
  } else {
    stop("n_cores must be a number, 'serial', or 'max'")
  }
  
  # If only 1 core, use serial
  if (cores_value == 1) {
    if (verbose) {
      message("  ✓ Configured: Serial processing (only 1 core)")
    }
    return(list(
      BPPARAM = BiocParallel::SerialParam(),
      n_cores = 1,
      mode = "serial",
      type = "serial"
    ))
  }
  
  # Determine parallelization type based on mode and system constraints
  param_type <- mode
  
  if (mode == "fork") {
    if (is_windows) {
      warning("Fork-based parallelization not available on Windows. Using socket-based.")
      param_type <- "socket"
    } else if (for_java) {
      warning("Fork-based parallelization incompatible with Java/rJava. Using socket-based.")
      param_type <- "socket"
    }
  }
  
  # Create BPPARAM
  if (param_type == "socket") {
    BPPARAM <- BiocParallel::SnowParam(workers = cores_value, type = "SOCK")
    if (verbose) {
      message("  Mode: socket")
      message("  ✓ Configured: Socket-based parallelization (SnowParam)")
      message("    Workers: ", cores_value)
      message("    Memory: ~2-3 GB per worker")
      message("    Java-compatible: YES")
      if (is.character(n_cores) && n_cores == "max") {
        message("    ⚠ Using all cores may impact system responsiveness")
      }
    }
  } else if (param_type == "fork") {
    BPPARAM <- BiocParallel::MulticoreParam(workers = cores_value)
    if (verbose) {
      message("  Mode: fork")
      message("  ✓ Configured: Fork-based parallelization (MulticoreParam)")
      message("    Workers: ", cores_value)
      message("    Memory: Shared memory (more efficient)")
      message("    Java-compatible: NO")
      if (for_java) {
        message("    ⚠ WARNING: This will NOT work with Java/rJava operations!")
      }
    }
  }
  
  # Return configuration
  config <- list(
    BPPARAM = BPPARAM,
    n_cores = cores_value,
    mode = param_type,
    type = param_type
  )
  
  if (verbose) {
    message("  Total cores available: ", detected_cores)
    message("  Cores configured: ", cores_value, " (", 
            round(cores_value/detected_cores * 100, 1), "% utilization)")
  }
  
  return(invisible(config))
}

# Auto-run example if sourced directly
if (sys.nframe() == 0) {
  message("configure_parallel_processing.R sourced")
  message("\nUsage examples:")
  message("  # Serial processing:")
  message("  config <- configure_parallel_processing(n_cores = 'serial')")
  message("")
  message("  # Maximum cores (socket-based):")
  message("  config <- configure_parallel_processing(n_cores = 'max', mode = 'socket')")
  message("")
  message("  # Manual cores (socket-based for Java):")
  message("  config <- configure_parallel_processing(n_cores = 8, mode = 'socket')")
  message("")
  message("  # Fork-based for non-Java operations:")
  message("  config <- configure_parallel_processing(n_cores = 16, mode = 'fork', for_java = FALSE)")
  message("")
  message("  # Extract BPPARAM:")
  message("  BPPARAM <- config$BPPARAM")
}

