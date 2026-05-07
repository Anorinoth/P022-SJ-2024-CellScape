#' Setup Project Environment
#'
#' Configures renv, creates required directories, and sets up temp directory paths.
#' This should be run at the start of SpicyFlow.Rmd or any analysis scripts.
#'
#' @param root_dir Character, path to the project root directory.
#'   Default: "/mnt/b/Projects/Image_Analysis_Projects/SJ_04.2/MAIN"
#' @param auto_snapshot Logical, if TRUE, automatically updates renv.lock when
#'   packages differ from lockfile. Default: FALSE (manual control)
#' @param verbose Logical, if TRUE, displays setup messages. Default: TRUE
#'
#' @return A list containing paths:
#'   - ROOT_DIR: Project root directory
#'   - RES_DIR: Results directory
#'   - PROG_DIR: Programs/scripts directory
#'   - CONFIG_DIR: Configuration directory
#'   - TEMP_DIR: Temporary directory on B drive
#'
#' @details
#' This function handles three main setup tasks:
#' 
#' 1. **renv Configuration**: Sets options for dependency management and checks
#'    if the project is synced with renv.lock
#' 
#' 2. **Directory Structure**: Creates standard project directories if they
#'    don't exist
#' 
#' 3. **Temp Directory**: Configures temporary directory on B drive for large
#'    intermediate files
#'
#' @examples
#' \dontrun{
#' # Standard setup
#' dirs <- setup_project()
#' 
#' # Custom root directory
#' dirs <- setup_project(root_dir = "/path/to/project")
#' 
#' # Auto-update renv lockfile
#' dirs <- setup_project(auto_snapshot = TRUE)
#' 
#' # Extract directories for use
#' ROOT_DIR <- dirs$ROOT_DIR
#' RES_DIR <- dirs$RES_DIR
#' }
#'
#' @export
setup_project <- function(root_dir = "/mnt/b/Projects/Image_Analysis_Projects/SJ_04.2/MAIN",
                          auto_snapshot = FALSE,
                          verbose = TRUE) {
  
  if (verbose) {
    message("\n=== Project Setup ===")
  }
  
  # ============================================================================
  # 1. Configure renv
  # ============================================================================
  
  if (verbose) {
    message("\n1. Configuring renv...")
  }
  
  # Set renv options for better performance and reliability
  options(renv.config.dependencies.limit = 10000L)  # Use large integer instead of Inf
  options(renv.verbose = FALSE)  # Reduce verbosity
  
  # Use explicit dependency discovery - ONLY look at library() and require() calls
  # This prevents renv from trying to parse ALL code and hitting example output
  options(renv.config.dependency.errors = "ignored")  # IGNORE parsing errors entirely
  options(renv.config.dependencies.errors = "ignored")  # Alternative setting
  
  # Additional safety: don't infer dependencies from code, only explicit calls
  Sys.setenv(RENV_CONFIG_DEPENDENCY_ERRORS = "ignored")
  
  # Initialize renv if not already initialized
  if (!file.exists("renv.lock")) {
    if (verbose) {
      message("  Initializing renv for the first time...")
      message("  Note: This may take a while on first run. A .renvignore file will help with future runs.")
    }
    
    # Initialize with settings to avoid common issues
    renv::init(bare = FALSE, restart = FALSE)
    
  } else {
    if (verbose) {
      message("  ✓ renv already initialized")
    }
    
    # Check status - will fail if there are issues (fail-fast approach)
    status <- renv::status()
    
    # If not synced, inform user
    if (!is.null(status) && length(status) > 0) {
      if (verbose) {
        message("  ⚠ Project libraries differ from renv.lock")
        message("    Run renv::restore() to revert to lockfile versions")
        message("    Or run renv::snapshot() to UPDATE lockfile with current versions")
      }
      
      if (auto_snapshot) {
        if (verbose) {
          message("  → Auto-updating renv.lock...")
        }
        renv::snapshot(prompt = FALSE)
      }
      
    } else {
      if (verbose) {
        message("  ✓ Project is synced with renv.lock")
      }
    }
  }
  
  # ============================================================================
  # 2. Setup Directory Structure
  # ============================================================================
  
  if (verbose) {
    message("\n2. Setting up directories...")
  }
  
  # Define standard directories
  ROOT_DIR <- root_dir
  RES_DIR <- file.path(ROOT_DIR, "RES")
  PROG_DIR <- file.path(ROOT_DIR, "PROG")
  CONFIG_DIR <- file.path(ROOT_DIR, "CONFIG")
  
  # Create directories if they don't exist
  dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)
  dir.create(PROG_DIR, showWarnings = FALSE, recursive = TRUE)
  dir.create(CONFIG_DIR, showWarnings = FALSE, recursive = TRUE)
  
  # Print directory structure for confirmation
  if (verbose) {
    message("  Directory structure:")
    message("    ROOT_DIR:   ", ROOT_DIR)
    message("    RES_DIR:    ", RES_DIR)
    message("    PROG_DIR:   ", PROG_DIR)
    message("    CONFIG_DIR: ", CONFIG_DIR)
  }
  
  # ============================================================================
  # 3. Configure Temp Directory
  # ============================================================================
  
  if (verbose) {
    message("\n3. Configuring temp directory...")
  }
  
  # Create temp directory on B drive
  TEMP_DIR <- file.path(ROOT_DIR, "temp")
  dir.create(TEMP_DIR, showWarnings = FALSE, recursive = TRUE)
  
  # Note: Setting environment variables here won't change tempdir() for the current session
  # R's tempdir() is set at startup and cannot be changed mid-session
  # The TMPDIR must be set BEFORE R starts (see .Renviron file in project root)
  Sys.setenv(TMPDIR = TEMP_DIR)
  Sys.setenv(TEMP = TEMP_DIR)
  Sys.setenv(TMP = TEMP_DIR)
  
  # Check if temp directory is on B drive
  current_tempdir <- tempdir()
  if (verbose) {
    message("  R temp directory: ", current_tempdir)
  }
  
  if (!grepl("^/mnt/b/", current_tempdir)) {
    if (verbose) {
      message("  ⚠ WARNING: R temp directory is NOT on B drive")
      message("    To fix: Restart R session after creating .Renviron file")
      message("    For now, using explicit TEMP_DIR path where needed: ", TEMP_DIR)
    }
  } else {
    if (verbose) {
      message("  ✓ SUCCESS: R temp directory is on B drive")
    }
  }
  
  # ============================================================================
  # 4. Check R Memory Settings
  # ============================================================================
  
  if (verbose) {
    message("\n4. Checking R memory settings...")
  }
  
  # Check R_MAX_VSIZE
  max_vsize <- Sys.getenv("R_MAX_VSIZE")
  if (max_vsize != "") {
    if (verbose) {
      message("  ✓ R_MAX_VSIZE:     ", max_vsize)
    }
  } else {
    if (verbose) {
      message("  ⚠ R_MAX_VSIZE:     Not set (add to ~/.Renviron)")
    }
  }
  
  # Check expression limit
  expr_limit <- getOption("expressions")
  if (verbose) {
    if (expr_limit >= 100000) {
      message("  ✓ Expression limit: ", format(expr_limit, big.mark = ","))
    } else {
      message("  ⚠ Expression limit: ", format(expr_limit, big.mark = ","), " (recommend 500,000)")
      message("    Add to .Rprofile: options(expressions = 500000)")
    }
  }
  
  # Check Java heap (if Java configured)
  java_params <- getOption("java.parameters")
  if (!is.null(java_params)) {
    if (verbose) {
      message("  ✓ Java heap:       ", java_params)
    }
  }
  
  # Note about --max-ppsize (cannot be checked programmatically)
  if (verbose) {
    message("\n  Note: --max-ppsize cannot be verified within R")
    message("  If quantification fails with 'protection stack overflow':")
    message("    • Start R with: R --max-ppsize=500000")
    message("    • Or configure in .vscode/settings.json")
  }
  
  if (verbose) {
    message("\n=== Project Setup Complete ===\n")
  }
  
  # Return directory paths and settings
  return(invisible(list(
    ROOT_DIR = ROOT_DIR,
    RES_DIR = RES_DIR,
    PROG_DIR = PROG_DIR,
    CONFIG_DIR = CONFIG_DIR,
    TEMP_DIR = TEMP_DIR,
    settings = list(
      R_MAX_VSIZE = max_vsize,
      expressions = expr_limit,
      java_heap = java_params
    )
  )))
}

# Auto-run example if sourced directly
if (sys.nframe() == 0) {
  message("setup_project.R sourced")
  message("\nUsage:")
  message("  dirs <- setup_project()")
  message("  ROOT_DIR <- dirs$ROOT_DIR")
  message("  RES_DIR <- dirs$RES_DIR")
  message("  PROG_DIR <- dirs$PROG_DIR")
}

