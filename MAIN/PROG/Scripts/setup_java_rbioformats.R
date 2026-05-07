#' Setup Java and RBioFormats for OME-TIFF Processing
#' Author: Dr. Michael Zaiken
#'
#' This script configures Java, installs rJava and RBioFormats if needed,
#' downloads required JAR files, and verifies everything is working correctly.
#'
#' @return A list with status information:
#'   - java_available: Logical, TRUE if Java is configured
#'   - rbioformats_available: Logical, TRUE if RBioFormats is ready
#'   - java_version: Character, Java version string
#'   - java_heap_gb: Numeric, Java max heap size in GB
#'   - rbioformats_version: Character, RBioFormats package version
#'   - jar_count: Integer, number of JAR files found
#'
#' @examples
#' \dontrun{
#' source("MAIN/PROG/Scripts/setup_java_rbioformats.R")
#' status <- setup_java_rbioformats()
#' if (status$rbioformats_available) {
#'   message("✓ Ready to process OME-TIFF files!")
#' }
#' }
setup_java_rbioformats <- function() {
  # Initialize return status
  status <- list(
    java_available = FALSE,
    rbioformats_available = FALSE,
    java_version = NA_character_,
    java_heap_gb = NA_real_,
    rbioformats_version = NA_character_,
    jar_count = 0L
  )

  # Configure Java for RBioFormats (required for OME-TIFF processing)
  message("\n=== Configuring Java environment ===")

  # CRITICAL: Java heap space options are set in .Rprofile (before renv loads)
  # This is essential for renv + RBioFormats compatibility
  # See: https://github.com/aoles/RBioFormats#the-javalangoutofmemoryerror-error
  java_params <- getOption("java.parameters")
  if (!is.null(java_params) && grepl("-Xmx", java_params)) {
    message("  Java heap space configured: ", java_params)
    if (!grepl("768g", java_params)) {
      message("  ATTENTION: Java parameters in current session (", java_params, ") do not match 768g.")
      message("  Because Java is already initialized, you MUST restart R for changes to take effect.")
    }
  } else {
    message("  WARNING: Java heap space not configured in .Rprofile!")
    message("  This may cause RBioFormats installation to fail.")
    if (!"rJava" %in% loadedNamespaces()) {
      message("  Attempting to set now (may not work with renv)...")
      options(java.parameters = "-Xmx768g")
    }
  }

  # == CHECK: Verify Java configuration ==
  java_home <- Sys.getenv("JAVA_HOME")

  if (java_home == "") {
    message("WARNING: JAVA_HOME not set!")
    message("  RBioFormats (for OME-TIFF files) requires Java.")
    message("  To fix:")
    message("    1. Verify .Renviron contains: JAVA_HOME=/usr/lib/jvm/default-java")
    message("    2. Run from terminal: sudo R CMD javareconf")
    message("    3. Restart R completely (important!)")
    return(status)
  }

  message("  JAVA_HOME: ", java_home)

  # == STEP 1: Install and initialize rJava ==
  message("\n=== Step 1: Initialize rJava ===")

  # Install rJava if needed
  if (!requireNamespace("rJava", quietly = TRUE)) {
    message("  Installing rJava from source...")
    install.packages("rJava", type = "source", repos = "https://cloud.r-project.org")
  }

  # Try to initialize Java
  status$java_available <- tryCatch(
    {
      suppressPackageStartupMessages(library(rJava))
      .jinit() # Initialize Java Virtual Machine
      java_version <- .jcall("java/lang/System", "S", "getProperty", "java.runtime.version")
      message("  ✓ Java initialized successfully (version ", java_version, ")")
      status$java_version <- java_version

      # Display Java heap space settings
      runtime <- .jcall("java/lang/Runtime", "Ljava/lang/Runtime;", "getRuntime")
      max_memory <- .jcall(runtime, "J", "maxMemory")
      heap_gb <- max_memory / 1024^3
      message("  ✓ Java max heap size: ", format(heap_gb, digits = 2), " GB")
      status$java_heap_gb <- heap_gb

      TRUE
    },
    error = function(e) {
      message("  ✗ Failed to initialize Java: ", conditionMessage(e))
      message("  Solution:")
      message("    1. Close R completely")
      message("    2. Ensure you ran: sudo R CMD javareconf")
      message("    3. Restart R and try again")
      FALSE
    }
  )

  if (!status$java_available) {
    return(status)
  }

  # == STEP 2: Install RBioFormats ==
  message("\n=== Step 2: Install RBioFormats ===")

  if (!requireNamespace("RBioFormats", quietly = TRUE)) {
    message("  Installing RBioFormats from Bioconductor...")
    message("  Note: This may take several minutes...")
    message("  Using --no-test-load to avoid renv + Java interaction issues...")

    # Save and modify environment
    old_check <- Sys.getenv("_R_CHECK_INSTALL_DEPENDS_", unset = NA)
    Sys.setenv("_R_CHECK_INSTALL_DEPENDS_" = "false")

    # Install RBioFormats
    BiocManager::install("RBioFormats",
      ask = FALSE,
      update = FALSE,
      INSTALL_opts = "--no-test-load",
      force = TRUE
    )

    # Restore environment
    if (is.na(old_check)) {
      Sys.unsetenv("_R_CHECK_INSTALL_DEPENDS_")
    } else {
      Sys.setenv("_R_CHECK_INSTALL_DEPENDS_" = old_check)
    }
  } else {
    message("  RBioFormats already installed")
  }

  # == STEP 3: Find RBioFormats JAR files ==
  message("\n=== Step 3: Locate Bio-Formats JAR files ===")

  jar_files <- character(0)
  rbioformats_path <- ""

  if (requireNamespace("RBioFormats", quietly = TRUE)) {
    rbioformats_path <- system.file(package = "RBioFormats")

    if (rbioformats_path == "") {
      message("  ✗ RBioFormats package not found in library")
      return(status)
    }

    message("  RBioFormats installed at: ", rbioformats_path)

    # Check multiple possible JAR locations
    jar_locations <- c(
      file.path(rbioformats_path, "java"),
      file.path(rbioformats_path, "inst", "java"),
      file.path(rbioformats_path, "jars")
    )

    for (loc in jar_locations) {
      if (dir.exists(loc)) {
        jars <- list.files(loc, pattern = "\\.jar$", full.names = TRUE, recursive = TRUE)
        if (length(jars) > 0) {
          jar_files <- c(jar_files, jars)
        }
      }
    }

    message("  Found ", length(jar_files), " JAR file(s)")
    status$jar_count <- length(jar_files)
  }

  # == STEP 4: Download Bio-Formats package JAR if missing ==
  message("\n=== Step 4: Verify Bio-Formats library ===")

  if (length(jar_files) == 1 && basename(jar_files[1]) == "RBioFormats.jar") {
    message("  WARNING: Only RBioFormats.jar found - Bio-Formats library JAR is missing!")
    message("  Downloading bioformats_package.jar from GitHub releases...")

    jar_dir <- file.path(rbioformats_path, "java")
    message("  Download directory: ", jar_dir)

    # Download the complete Bio-Formats package JAR from GitHub releases
    bioformats_version <- "v8.3.0"
    jar_url <- sprintf(
      "https://github.com/ome/bioformats/releases/download/%s/bioformats_package.jar",
      bioformats_version
    )
    jar_file <- file.path(jar_dir, "bioformats_package.jar")

    if (file.exists(jar_file)) {
      message("    bioformats_package.jar already exists")
    } else {
      message("    Downloading bioformats_package.jar from ", bioformats_version, "...")
      download_success <- tryCatch(
        {
          download.file(jar_url, jar_file, mode = "wb", quiet = FALSE)
          message("    ✓ Download successful!")
          TRUE
        },
        error = function(e) {
          message("    ✗ Download failed: ", conditionMessage(e))
          message("    Alternative: remove.packages('RBioFormats'); BiocManager::install('RBioFormats')")
          FALSE
        }
      )

      if (!download_success) {
        return(status)
      }
    }

    # Refresh JAR file list after successful download
    jar_files <- list.files(jar_dir, pattern = "\\.jar$", full.names = TRUE)
    message("  ✓ Bio-Formats JAR ready (total: ", length(jar_files), " files)")
    status$jar_count <- length(jar_files)
  }

  # == STEP 5: Load RBioFormats with all JARs ==
  message("\n=== Step 5: Load RBioFormats ===")

  if (length(jar_files) > 0) {
    message("  Loading RBioFormats with ", length(jar_files), " JAR file(s)...")

    status$rbioformats_available <- tryCatch(
      {
        # Add all JARs to classpath
        for (jar in jar_files) {
          .jaddClassPath(jar)
        }
        message("    ✓ Added all JARs to Java classpath")

        # Load RBioFormats
        suppressPackageStartupMessages(library(RBioFormats))
        pkg_version <- as.character(packageVersion("RBioFormats"))
        message("    ✓ RBioFormats loaded (version ", pkg_version, ")")
        status$rbioformats_version <- pkg_version

        # Verify Bio-Formats is accessible
        tryCatch(
          {
            mem_info <- checkJavaMemory()
            message("    ✓ Bio-Formats library accessible")
            message("    ✓ Available Java memory: ", mem_info["available"], " MB")
          },
          error = function(e) {
            message("    ⚠ Bio-Formats check failed: ", conditionMessage(e))
          }
        )

        TRUE
      },
      error = function(e) {
        message("  ✗ Failed to load RBioFormats: ", conditionMessage(e))
        message("  Try: remove.packages('RBioFormats'); restart R; reinstall")
        FALSE
      }
    )
  }

  # == SUMMARY ==
  message("\n=== Setup Summary ===")
  if (status$rbioformats_available) {
    message("✓ RBioFormats ready for OME-TIFF processing")
    message("  Java version: ", status$java_version)
    message("  Java heap: ", round(status$java_heap_gb, 1), " GB")
    message("  RBioFormats version: ", status$rbioformats_version)
    message("  JAR files: ", status$jar_count)
  } else if (!status$java_available) {
    message("✗ Java not available - RBioFormats setup skipped")
  } else {
    message("✗ RBioFormats not available (see errors above)")
  }

  return(invisible(status))
}

# If sourced directly (not as a function), run the setup
if (!exists("PROG_DIR")) {
  message("Note: Running as standalone script")
  message("For use in workflows, source this file and call: setup_java_rbioformats()")
}

# Auto-run if not being sourced as a library
if (sys.nframe() == 0) {
  setup_java_rbioformats()
}
