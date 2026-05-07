#' Import Python Segmentation & Quantification Results into SpatialExperiment
#' Author: Dr. Michael Zaiken
#'
#' Reads the CSV output from segment_and_quantify.py and constructs a
#' SpatialExperiment object compatible with SpicyFlow downstream analysis.
#'
#' @param image_name Character, e.g. "CAR_TREG_2"
#' @param prog_dir Character, path to MAIN/PROG directory containing the CSV
#' @return SpatialExperiment object
import_python_results <- function(image_name, prog_dir) {
    csv_file <- file.path(prog_dir, paste0(image_name, "_cells.csv"))
    if (!file.exists(csv_file)) {
        stop(
            "Python results not found: ", csv_file,
            "\nRun segment_and_quantify.py first."
        )
    }

    cat(sprintf(
        "[%s] Loading Python results: %s\n",
        format(Sys.time(), "%H:%M:%S"), basename(csv_file)
    ))

    # Read CSV
    cells_df <- data.table::fread(csv_file)
    cat(sprintf(
        "  %s cells × %d columns\n",
        format(nrow(cells_df), big.mark = ","), ncol(cells_df)
    ))

    # Extract channel columns (prefixed with "ch_")
    channel_cols <- grep("^ch_", colnames(cells_df), value = TRUE)
    cat(sprintf("  %d channels detected\n", length(channel_cols)))

    # Clean channel names (remove "ch_" prefix and numeric prefix like "01_")
    clean_names <- sub("^ch_", "", channel_cols)
    # Also remove the numeric prefix if present: "01_CD19-PerCP-Cy5_5" → "CD19-PerCP-Cy5_5"
    clean_names <- sub("^\\d+_", "", clean_names)

    # Build expression matrix (features × cells, as required by SpatialExperiment)
    expr_matrix <- t(as.matrix(cells_df[, channel_cols, with = FALSE]))
    rownames(expr_matrix) <- clean_names
    colnames(expr_matrix) <- paste0("cell_", cells_df$cell_id)

    # Build spatial coordinates
    spatial_coords <- as.matrix(cells_df[, c("x", "y")])
    rownames(spatial_coords) <- paste0("cell_", cells_df$cell_id)

    # Create SpatialExperiment
    library(SpatialExperiment)

    spe <- SpatialExperiment(
        assays = list(counts = expr_matrix),
        spatialCoords = spatial_coords
    )

    # Add cell metadata
    spe$image_name <- image_name
    spe$cell_id <- cells_df$cell_id
    spe$area <- cells_df$area

    # Set coordinate names
    spatialCoordsNames(spe) <- c("x", "y")

    cat(sprintf(
        "  ✓ SpatialExperiment: %s cells, %d features\n",
        format(ncol(spe), big.mark = ","), nrow(spe)
    ))

    # Save as .qs for compatibility with existing pipeline
    spe_file <- file.path(prog_dir, paste0(image_name, "_cells.qs"))
    qs::qsave(spe, spe_file, preset = "fast")
    cat(sprintf("  💾 Saved: %s\n", basename(spe_file)))

    return(spe)
}
