library(igraph)

set.seed(42)
N <- 150000
nodes <- data.frame(id = 1:N, domain = sample(c("A","B","C"), N, replace=TRUE))
# Generate some random edges
from <- sample(1:N, 400000, replace=TRUE)
to <- sample(1:N, 400000, replace=TRUE)
edges <- data.frame(from=from, to=to)

t0 <- Sys.time()
g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
print(paste("Graph creation:", Sys.time() - t0))

t1 <- Sys.time()
comp <- igraph::components(g)
print(paste("Components:", Sys.time() - t1))

nodes$patch_id <- paste0(nodes$domain, "_", comp$membership)

print("Done")
