# ==============================================================================
# .Rprofile - Project-specific R startup configuration
# ==============================================================================
# This file runs BEFORE .Renviron and BEFORE renv is loaded
# Critical for Java configuration with renv + RBioFormats

# Set Java options BEFORE any package loading (including renv)
# This ensures Java settings are available during package installation testing
options(java.parameters = "-Xmx768g")

# Set R memory limits for large image analysis (SJ_04.2 project)
# Increase expression limit to prevent protection stack overflow
options(expressions = 500000) # Default ~5000, increase to 500k

# Source renv activation script (if it exists)
source("renv/activate.R")

# Inform user
if (interactive()) {
  message("Project .Rprofile loaded")
  message("  Java heap space: 768 GB")
  message("  Expression limit: ", format(getOption("expressions"), big.mark = ","))

  # Check if R was started with --max-ppsize
  message("\nNote: For large cell quantification (86k+ cells):")
  message("  Start R with: R --max-ppsize=500000")
  message("  Or set in Cursor: .vscode/settings.json (already configured)")
}
