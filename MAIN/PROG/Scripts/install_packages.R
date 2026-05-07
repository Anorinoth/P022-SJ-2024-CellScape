#' Install Required Packages
#'
#' Installs packages from Bioconductor, CRAN, and GitHub based on configuration files.
#' Only installs packages that are not already present in the library.
#'
#' @param config_dir Character, path to the directory containing package list files.
#'   Default: "MAIN/CONFIG"
#' @param install_bioc Logical, if TRUE, installs Bioconductor packages. Default: TRUE
#' @param install_cran Logical, if TRUE, installs CRAN packages. Default: TRUE
#' @param install_github Logical, if TRUE, installs GitHub packages. Default: TRUE
#' @param force_reinstall Logical, if TRUE, reinstalls packages even if present. Default: FALSE
#' @param verbose Logical, if TRUE, displays installation messages. Default: TRUE
#'
#' @return A list containing:
#'   - bioc_installed: Names of Bioconductor packages installed
#'   - cran_installed: Names of CRAN packages installed
#'   - github_installed: Names of GitHub packages installed
#'   - all_present: TRUE if all packages were already installed
#'
#' @details
#' Package lists are read from configuration files in CONFIG/:
#' - `packages_bioconductor.txt`: Bioconductor packages (one per line)
#' - `packages_cran.txt`: CRAN packages (one per line)
#' - `packages_github.txt`: GitHub packages in format "user/repo" (one per line)
#' 
#' Lines starting with # are treated as comments and ignored.
#' Blank lines are also ignored.
#'
#' @examples
#' \dontrun{
#' # Install all packages
#' result <- install_packages()
#' 
#' # Install only CRAN packages
#' result <- install_packages(install_bioc = FALSE, install_github = FALSE)
#' 
#' # Force reinstall all packages
#' result <- install_packages(force_reinstall = TRUE)
#' }
#'
#' @export
install_packages <- function(config_dir = "MAIN/CONFIG",
                             install_bioc = TRUE,
                             install_cran = TRUE,
                             install_github = TRUE,
                             force_reinstall = FALSE,
                             verbose = TRUE) {
  
  if (verbose) {
    message("\n=== Package Installation ===")
  }
  
  # Helper function to read package lists from config files
  read_package_list <- function(file_path) {
    if (!file.exists(file_path)) {
      if (verbose) {
        message("  ⚠ Config file not found: ", file_path)
      }
      return(character(0))
    }
    
    # Read file
    lines <- readLines(file_path, warn = FALSE)
    
    # Remove comments (lines starting with #) and blank lines
    lines <- lines[!grepl("^\\s*#", lines)]
    lines <- lines[!grepl("^\\s*$", lines)]
    
    # Trim whitespace
    lines <- trimws(lines)
    
    return(lines)
  }
  
  # Initialize results
  results <- list(
    bioc_installed = character(0),
    cran_installed = character(0),
    github_installed = character(0),
    all_present = TRUE
  )
  
  # ============================================================================
  # 1. Install Bioconductor Packages
  # ============================================================================
  
  if (install_bioc) {
    if (verbose) {
      message("\n1. Bioconductor Packages")
    }
    
    # Install BiocManager if not already installed
    if (!require("BiocManager", quietly = TRUE)) {
      if (verbose) {
        message("  → Installing BiocManager...")
      }
      install.packages("BiocManager")
    }
    
    # Read Bioconductor package list
    bioc_file <- file.path(config_dir, "packages_bioconductor.txt")
    bioc_packages <- read_package_list(bioc_file)
    
    if (length(bioc_packages) > 0) {
      if (verbose) {
        message("  Found ", length(bioc_packages), " Bioconductor packages in config")
      }
      
      # Check which packages are missing or force reinstall
      if (force_reinstall) {
        bioc_missing <- bioc_packages
      } else {
        bioc_missing <- bioc_packages[!bioc_packages %in% rownames(installed.packages())]
      }
      
      # Install missing packages
      if (length(bioc_missing) > 0) {
        if (verbose) {
          message("  → Installing ", length(bioc_missing), " packages:")
          message("    ", paste(bioc_missing, collapse = ", "))
        }
        BiocManager::install(bioc_missing, ask = FALSE, update = TRUE)
        results$bioc_installed <- bioc_missing
        results$all_present <- FALSE
      } else {
        if (verbose) {
          message("  ✓ All Bioconductor packages already installed")
        }
      }
    } else {
      if (verbose) {
        message("  ⚠ No Bioconductor packages listed in config")
      }
    }
  }
  
  # ============================================================================
  # 2. Install CRAN Packages
  # ============================================================================
  
  if (install_cran) {
    if (verbose) {
      message("\n2. CRAN Packages")
    }
    
    # Read CRAN package list
    cran_file <- file.path(config_dir, "packages_cran.txt")
    cran_packages <- read_package_list(cran_file)
    
    if (length(cran_packages) > 0) {
      if (verbose) {
        message("  Found ", length(cran_packages), " CRAN packages in config")
      }
      
      # Check which packages are missing or force reinstall
      if (force_reinstall) {
        cran_missing <- cran_packages
      } else {
        cran_missing <- cran_packages[!cran_packages %in% rownames(installed.packages())]
      }
      
      # Install missing packages
      if (length(cran_missing) > 0) {
        if (verbose) {
          message("  → Installing ", length(cran_missing), " packages:")
          message("    ", paste(cran_missing, collapse = ", "))
        }
        install.packages(cran_missing)
        results$cran_installed <- cran_missing
        results$all_present <- FALSE
      } else {
        if (verbose) {
          message("  ✓ All CRAN packages already installed")
        }
      }
    } else {
      if (verbose) {
        message("  ⚠ No CRAN packages listed in config")
      }
    }
  }
  
  # ============================================================================
  # 3. Install GitHub Packages
  # ============================================================================
  
  if (install_github) {
    if (verbose) {
      message("\n3. GitHub Packages")
    }
    
    # Ensure remotes is installed
    if (!require("remotes", quietly = TRUE)) {
      if (verbose) {
        message("  → Installing remotes package...")
      }
      install.packages("remotes")
    }
    
    # Read GitHub package list
    github_file <- file.path(config_dir, "packages_github.txt")
    github_packages <- read_package_list(github_file)
    
    if (length(github_packages) > 0) {
      if (verbose) {
        message("  Found ", length(github_packages), " GitHub packages in config")
      }
      
      # For GitHub packages, extract the package name from "user/repo"
      # and check if already installed
      for (gh_pkg in github_packages) {
        # Extract package name (last part after /)
        pkg_name <- sub(".*/", "", gh_pkg)
        pkg_name <- sub("@.*", "", pkg_name)  # Remove branch if present
        
        # Check if package is installed
        needs_install <- force_reinstall || !pkg_name %in% rownames(installed.packages())
        
        if (needs_install) {
          if (verbose) {
            message("  → Installing ", gh_pkg, "...")
          }
          remotes::install_github(gh_pkg)
          results$github_installed <- c(results$github_installed, gh_pkg)
          results$all_present <- FALSE
        } else {
          if (verbose) {
            message("  ✓ ", pkg_name, " already installed")
          }
        }
      }
    } else {
      if (verbose) {
        message("  ⚠ No GitHub packages listed in config")
      }
    }
  }
  
  # ============================================================================
  # Summary
  # ============================================================================
  
  if (verbose) {
    message("\n=== Package Installation Complete ===")
    
    total_installed <- length(results$bioc_installed) + 
                      length(results$cran_installed) + 
                      length(results$github_installed)
    
    if (total_installed > 0) {
      message("  Total packages installed: ", total_installed)
      if (length(results$bioc_installed) > 0) {
        message("    Bioconductor: ", length(results$bioc_installed))
      }
      if (length(results$cran_installed) > 0) {
        message("    CRAN: ", length(results$cran_installed))
      }
      if (length(results$github_installed) > 0) {
        message("    GitHub: ", length(results$github_installed))
      }
    } else {
      message("  ✓ All required packages already installed")
    }
    message("")
  }
  
  return(invisible(results))
}

# Auto-run example if sourced directly
if (sys.nframe() == 0) {
  message("install_packages.R sourced")
  message("\nUsage:")
  message("  # Install all packages")
  message("  result <- install_packages()")
  message("")
  message("  # Install only CRAN packages")
  message("  result <- install_packages(install_bioc = FALSE, install_github = FALSE)")
  message("")
  message("  # Force reinstall")
  message("  result <- install_packages(force_reinstall = TRUE)")
}

