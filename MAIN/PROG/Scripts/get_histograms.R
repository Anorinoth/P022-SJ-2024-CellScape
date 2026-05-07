library(SpatialExperiment)
library(deldir)
library(igraph)
library(dplyr)
library(ggplot2)

spe <- readRDS("/mnt/b/Projects/Image_Analysis_Projects/SJ_04.2/MAIN/PROG/Checkpoints/combined_spe_after_section7.rds")
df_all <- data.frame(
    id = colnames(spe),
    x = spatialCoords(spe)[, 1],
    y = spatialCoords(spe)[, 2],
    imageID = colData(spe)$imageID,
    domain = colData(spe)$SpatialDomain,
    stringsAsFactors = FALSE
)

# Store results
patch_areas_list <- list()
boundary_lengths_list <- list()

images <- unique(df_all$imageID)
for (img in images) {
    message("Processing ", img)
    df <- df_all[df_all$imageID == img, ]
    dup_idx <- duplicated(df[, c("x", "y")])
    if (any(dup_idx)) df <- df[!dup_idx, ]
    if (nrow(df) < 10) next
    
    vt <- tryCatch(deldir(df$x, df$y, suppressMsge = TRUE), error = function(e) NULL)
    if (is.null(vt)) next
    
    cell_areas <- vt$summary$dir.area
    delsgs <- vt$delsgs
    
    adj_cells <- data.frame(
        ind1 = delsgs$ind1,
        ind2 = delsgs$ind2,
        dom1 = df$domain[delsgs$ind1],
        dom2 = df$domain[delsgs$ind2],
        stringsAsFactors = FALSE
    )
    
    # Patch areas for B cell follicle and Treg Enriched Zone
    same_dom_edges <- adj_cells %>% filter(dom1 == dom2)
    nodes <- data.frame(id = 1:nrow(df), domain = df$domain, area = cell_areas, stringsAsFactors = FALSE)
    edges_graph <- data.frame(from = same_dom_edges$ind1, to = same_dom_edges$ind2, stringsAsFactors = FALSE)
    
    g <- graph_from_data_frame(edges_graph, directed = FALSE, vertices = nodes)
    comp <- components(g)
    nodes$patch_id <- paste0(nodes$domain, "_", img, "_", comp$membership)
    
    patch_areas <- nodes %>%
        group_by(patch_id, domain) %>%
        summarise(total_area = sum(area), .groups = "drop") %>%
        filter(domain %in% c("B_Cell_Zone", "T_Cell_Zone", "B Cell Zone", "T Cell Zone"))
    
    patch_areas$imageID <- img
    patch_areas_list[[img]] <- patch_areas
    
    # Boundary lengths between B cell follicle and Treg Enriched Zone
    # Need edge lengths. delsgs has x1, y1, x2, y2 for Voronoi edges in dirsgs?
    # Wait, vt$dirsgs contains the boundaries!
    dirsgs <- vt$dirsgs
    bound_df <- data.frame(
        ind1 = dirsgs$ind1,
        ind2 = dirsgs$ind2,
        dom1 = df$domain[dirsgs$ind1],
        dom2 = df$domain[dirsgs$ind2],
        length = sqrt((dirsgs$x1 - dirsgs$x2)^2 + (dirsgs$y1 - dirsgs$y2)^2),
        stringsAsFactors = FALSE
    )
    
    target_doms <- c("B_Cell_Zone", "T_Cell_Zone", "B Cell Zone", "T Cell Zone")
    # Sort dom1 and dom2 alphabetically so we can easily filter
    bounds <- bound_df %>% filter(dom1 %in% target_doms & dom2 %in% target_doms & dom1 != dom2)
    if(nrow(bounds) > 0) {
        bounds$imageID <- img
        boundary_lengths_list[[img]] <- bounds
    }
}

all_areas <- bind_rows(patch_areas_list)
all_lengths <- bind_rows(boundary_lengths_list)

# Plot Areas
p1 <- ggplot(all_areas, aes(x = total_area, fill = domain)) +
    geom_histogram(bins = 50, alpha = 0.6, position="identity") +
    scale_x_log10() +
    facet_wrap(~imageID) +
    theme_minimal() +
    labs(title = "Distribution of Patch Areas", x = "Log10(Patch Area)", y = "Count")
ggsave("/mnt/b/Projects/Image_Analysis_Projects/SJ_04.2/MAIN/PROG/Plots/Patch_Areas_Histogram.png", p1, width = 10, height = 8)

# Plot Lengths
p2 <- ggplot(all_lengths, aes(x = length)) +
    geom_histogram(bins = 50, alpha = 0.6, fill = "blue") +
    scale_x_log10() +
    facet_wrap(~imageID) +
    theme_minimal() +
    labs(title = "Distribution of Shared Boundary Lengths (T Cell Zone vs B Cell Zone)", x = "Log10(Boundary Edge Length)", y = "Count")
ggsave("/mnt/b/Projects/Image_Analysis_Projects/SJ_04.2/MAIN/PROG/Plots/Boundary_Lengths_Histogram.png", p2, width = 10, height = 8)

cat("Success\n")
