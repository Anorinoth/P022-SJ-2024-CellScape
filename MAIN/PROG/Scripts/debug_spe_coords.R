# Minimal reproduction script for SpatialExperiment coordinate modification

if (!requireNamespace("SpatialExperiment", quietly = TRUE)) {
    stop("SpatialExperiment not installed")
}

# Create a dummy SpatialExperiment
coords <- matrix(1:10, ncol = 2)
colnames(coords) <- c("x", "y")
spe <- SpatialExperiment::SpatialExperiment(spatialCoords = coords)

message("Original coords:")
print(SpatialExperiment::spatialCoords(spe))

# Attempt 1: In-place modification (Current 'buggy' approach)
msg <- "Attempt 1: In-place modification"
message("\n", msg)
tryCatch(
    {
        SpatialExperiment::spatialCoords(spe)[, 1] <- SpatialExperiment::spatialCoords(spe)[, 2]
        print(SpatialExperiment::spatialCoords(spe))
    },
    error = function(e) message("Error: ", conditionMessage(e))
)


# Attempt 2: Explicit replacement (My fix)
msg <- "Attempt 2: Explicit replacement"
message("\n", msg)
coords_new <- SpatialExperiment::spatialCoords(spe)
coords_new[, 1] <- coords_new[, 1] + 100
SpatialExperiment::spatialCoords(spe) <- coords_new

message("Coords after explicit replacement:")
print(SpatialExperiment::spatialCoords(spe))

# Attempt 3: Direct int_colData manipulation (Hypothesis)
msg <- "Attempt 3: Check internal storage"
message("\n", msg)
# spatialCoords are usually stored in int_colData under 'spatialCoords'
internal <- SingleCellExperiment::int_colData(spe)$spatialCoords
message("Internal storage:")
print(internal)
