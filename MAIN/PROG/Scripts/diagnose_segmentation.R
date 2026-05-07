# Quick diagnostic to check segmentation results
# Run this after segmentation to see cell count changes

diagnose_segmentation <- function(masks_file) {
  library(EBImage)
  
  # Load masks
  if (!file.exists(masks_file)) {
    stop("Masks file not found: ", masks_file)
  }
  
  masks <- readRDS(masks_file)
  mask_matrix <- as.matrix(masks[[1]])
  
  # Count cells
  n_cells <- length(unique(as.vector(mask_matrix))) - 1
  
  # Get object sizes
  object_sizes <- table(as.vector(mask_matrix))
  object_sizes <- object_sizes[names(object_sizes) != "0"]
  
  cat("\n=== Segmentation Diagnostics ===\n")
  cat("Total cells: ", n_cells, "\n")
  cat("\nSize Distribution (pixels):\n")
  cat("  Min:    ", min(object_sizes), "\n")
  cat("  Q1:     ", quantile(object_sizes, 0.25), "\n")
  cat("  Median: ", quantile(object_sizes, 0.50), "\n")
  cat("  Mean:   ", round(mean(object_sizes), 1), "\n")
  cat("  Q3:     ", quantile(object_sizes, 0.75), "\n")
  cat("  Max:    ", max(object_sizes), "\n")
  
  # Size bins
  cat("\nSize Distribution by Bins:\n")
  bins <- cut(object_sizes, 
              breaks = c(0, 10, 20, 30, 50, 100, 200, 500, 1000, Inf),
              labels = c("<10", "10-20", "20-30", "30-50", "50-100", "100-200", "200-500", "500-1000", ">1000"))
  print(table(bins))
  
  # Very small objects (possible oversegmentation)
  very_small <- sum(object_sizes < 15)
  cat("\nVery small objects (<15 pixels): ", very_small, " (", 
      round(100 * very_small / n_cells, 1), "%)\n")
  
  # Very large objects (possible undersegmentation)
  very_large <- sum(object_sizes > 100)
  cat("Very large objects (>100 pixels): ", very_large, " (",
      round(100 * very_large / n_cells, 1), "%)\n")
  
  return(invisible(list(
    n_cells = n_cells,
    sizes = as.numeric(object_sizes)
  )))
}
