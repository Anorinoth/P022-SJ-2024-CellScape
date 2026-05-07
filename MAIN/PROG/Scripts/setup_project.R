#' Setup Project Environment
#' Author: Dr. Michael Zaiken
#'
#' Configures renv, creates required directories, and sets up temp directory paths.
#' This should be run at the start of SpicyFlow.Rmd or any analysis scripts.
#'
#' Paths default to environment variables (P022_CODE_DIR, P022_DATA_DIR, etc.).
#' Set these in your .Renviron file — see .Renviron.example for a template.
#'
#' @param code_dir Character, path to the Code root (repo root).
#'   Default: P022_CODE_DIR env var, or getwd() if unset.
#' @param data_dir Character, path to the Data root (raw + derived inputs).
#'   Default: P022_DATA_DIR env var. Required.
#' @param outputs_dir Character, path to the Outputs root (figures, CSVs, RDS results).
#'   Default: P022_OUTPUTS_DIR env var. Required.
#' @param reports_dir Character, path to the Reports root (knit HTML, run logs).
#'   Default: P022_REPORTS_DIR env var. Required.
#' @param paper_dir Character, path to the Paper root (manuscript files).
#'   Default: P022_PAPER_DIR env var. Optional.
#' @param auto_snapshot Logical, if TRUE, automatically updates renv.lock when
#'   packages differ from lockfile. Default: FALSE (manual control)
#' @param verbose Logical, if TRUE, displays setup messages. Default: TRUE
#'
#' @return A list containing paths:
#'   - CODE_DIR:         Code repository root
#'   - ROOT_DIR:         MAIN subdirectory (Code/MAIN)
#'   - PROG_DIR:         MAIN/PROG
#'   - CONFIG_DIR:       MAIN/CONFIG
#'   - DATA_DIR:         Data root (raw inputs + derived intermediates)
#'   - INPUT_DIR:        Data/data/input (OME-TIFFs)
#'   - SEGMENTATION_DIR: Data/PROG/SEGMENTATION (per-image cells.csv, masks.npz)
#'   - CHECKPOINT_DIR:   Data/PROG/Checkpoints (RDS checkpoints)
#'   - OUTPUTS_DIR:      Outputs root
#'   - RES_DIR:          Outputs/RES_v2 (current-run results)
#'   - REPORTS_DIR:      Reports root (knit HTML, pipeline logs)
#'   - PAPER_DIR:        Paper root (manuscript)
#'   - TEMP_DIR:         Temp directory for large intermediate files
#'
#' @examples
#' \dontrun{
#' # Standard setup (reads paths from .Renviron)
#' dirs <- setup_project()
#'
#' # Override specific paths
#' dirs <- setup_project(data_dir = "/path/to/Data")
#'
#' # Extract directories for use
#' ROOT_DIR       <- dirs$ROOT_DIR
#' RES_DIR        <- dirs$RES_DIR
#' PROG_DIR       <- dirs$PROG_DIR
#' CHECKPOINT_DIR <- dirs$CHECKPOINT_DIR
#' }
#'
#' @export
setup_project <- function(
    code_dir     = Sys.getenv("P022_CODE_DIR", unset = getwd()),
    data_dir     = Sys.getenv("P022_DATA_DIR", unset = ""),
    outputs_dir  = Sys.getenv("P022_OUTPUTS_DIR", unset = ""),
    reports_dir  = Sys.getenv("P022_REPORTS_DIR", unset = ""),
    paper_dir    = Sys.getenv("P022_PAPER_DIR", unset = ""),
    auto_snapshot = FALSE,
    verbose = TRUE) {

  # Validate required directories
  missing <- character(0)
  if (data_dir == "")    missing <- c(missing, "P022_DATA_DIR (data_dir)")
  if (outputs_dir == "") missing <- c(missing, "P022_OUTPUTS_DIR (outputs_dir)")
  if (reports_dir == "") missing <- c(missing, "P022_REPORTS_DIR (reports_dir)")
  if (length(missing) > 0) {
    stop(
      "Required directory paths not set. Missing:\n  ",
      paste(missing, collapse = "\n  "),
      "\n\nSet these in your .Renviron file (see .Renviron.example for a template),",
      "\nor pass them directly: setup_project(data_dir = '/path/to/Data', ...)",
      call. = FALSE
    )
  }

  if (verbose) {
    message("\n=== Project Setup ===")
  }

  # ============================================================================
  # 1. Configure renv
  # ============================================================================

  if (verbose) {
    message("\n1. Configuring renv...")
  }

  # Activate this project's renv regardless of which directory R started in.
  # Without this, running from a parent workspace uses the wrong renv library.
  if (!identical(normalizePath(getwd()), normalizePath(code_dir))) {
    if (verbose) {
      message("  → Switching to project root: ", code_dir)
    }
    setwd(code_dir)
  }
  renv::load(code_dir)

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
  # 2. Define canonical directory structure
  # ============================================================================

  if (verbose) {
    message("\n2. Setting up directories...")
  }

  # Code tree (read-only from pipeline perspective)
  CODE_DIR   <- code_dir
  ROOT_DIR   <- file.path(CODE_DIR, "MAIN")
  PROG_DIR   <- file.path(ROOT_DIR, "PROG")
  CONFIG_DIR <- file.path(ROOT_DIR, "CONFIG")

  # Data tree (raw inputs + derived intermediates on A drive)
  DATA_DIR         <- data_dir
  INPUT_DIR        <- file.path(DATA_DIR, "data", "input")
  SEGMENTATION_DIR <- file.path(DATA_DIR, "PROG", "SEGMENTATION")
  CHECKPOINT_DIR   <- file.path(DATA_DIR, "PROG", "Checkpoints")

  # Outputs tree (figures, CSVs, RDS results on B drive)
  OUTPUTS_DIR <- outputs_dir
  RES_DIR     <- file.path(OUTPUTS_DIR, "RES_v2")

  # Reports tree (knit HTML, run logs on B drive)
  REPORTS_DIR <- reports_dir

  # Paper root (manuscript on B drive)
  PAPER_DIR <- paper_dir

  # Create writable output directories if they do not exist
  dir.create(RES_DIR,          showWarnings = FALSE, recursive = TRUE)
  dir.create(CHECKPOINT_DIR,   showWarnings = FALSE, recursive = TRUE)
  dir.create(SEGMENTATION_DIR, showWarnings = FALSE, recursive = TRUE)
  dir.create(REPORTS_DIR,      showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(REPORTS_DIR, "pipeline_logs"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(REPORTS_DIR, "run_logs"),      showWarnings = FALSE, recursive = TRUE)

  if (verbose) {
    message("  Directory structure:")
    message("    CODE_DIR:         ", CODE_DIR)
    message("    ROOT_DIR:         ", ROOT_DIR)
    message("    PROG_DIR:         ", PROG_DIR)
    message("    CONFIG_DIR:       ", CONFIG_DIR)
    message("    DATA_DIR:         ", DATA_DIR)
    message("    INPUT_DIR:        ", INPUT_DIR)
    message("    SEGMENTATION_DIR: ", SEGMENTATION_DIR)
    message("    CHECKPOINT_DIR:   ", CHECKPOINT_DIR)
    message("    OUTPUTS_DIR:      ", OUTPUTS_DIR)
    message("    RES_DIR:          ", RES_DIR)
    message("    REPORTS_DIR:      ", REPORTS_DIR)
    message("    PAPER_DIR:        ", PAPER_DIR)
  }

  # ============================================================================
  # 3. Configure Temp Directory
  # ============================================================================

  if (verbose) {
    message("\n3. Configuring temp directory...")
  }

  # Temp directory for large intermediate image files
  # Set TMPDIR in .Renviron to use a drive with sufficient space
  TEMP_DIR <- Sys.getenv("P022_TEMP_DIR", unset = tempdir())
  dir.create(TEMP_DIR, showWarnings = FALSE, recursive = TRUE)

  if (verbose) {
    message("  Temp directory: ", TEMP_DIR)
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

  # Return all canonical directory paths and settings
  return(invisible(list(
    CODE_DIR         = CODE_DIR,
    ROOT_DIR         = ROOT_DIR,
    PROG_DIR         = PROG_DIR,
    CONFIG_DIR       = CONFIG_DIR,
    DATA_DIR         = DATA_DIR,
    INPUT_DIR        = INPUT_DIR,
    SEGMENTATION_DIR = SEGMENTATION_DIR,
    CHECKPOINT_DIR   = CHECKPOINT_DIR,
    OUTPUTS_DIR      = OUTPUTS_DIR,
    RES_DIR          = RES_DIR,
    REPORTS_DIR      = REPORTS_DIR,
    PAPER_DIR        = PAPER_DIR,
    TEMP_DIR         = TEMP_DIR,
    settings = list(
      R_MAX_VSIZE = max_vsize,
      expressions = expr_limit,
      java_heap   = java_params
    )
  )))
}

# Auto-run example if sourced directly
if (sys.nframe() == 0) {
  message("setup_project.R sourced")
  message("\nUsage:")
  message("  dirs <- setup_project()")
  message("  ROOT_DIR       <- dirs$ROOT_DIR")
  message("  RES_DIR        <- dirs$RES_DIR")
  message("  PROG_DIR       <- dirs$PROG_DIR")
  message("  CHECKPOINT_DIR <- dirs$CHECKPOINT_DIR")
}
