# ============================================================
# Gene Ontology (GO) Term Network Visualization
# Based on Shared Genes (Jaccard Similarity)
# ============================================================
# This script generates a network visualization of enriched GO terms
# based on shared genes (Jaccard index). Outputs a JPG figure and
# CSV with cluster assignments.
# ============================================================


# ------------------------------------------------------------
# 1. Load required libraries
# ------------------------------------------------------------
library(igraph)
library(ggraph)
library(tidygraph)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(ggforce)
library(dplyr)


# ------------------------------------------------------------
# 2. Read input (adjust path)
# ------------------------------------------------------------
# Expected columns: term, over_represented_FDR, gene_ids (comma-separated), numDEInCat
df <- read.csv("your_path_here/GO_enrichment_results.csv", sep="\t", header=TRUE)
df <- df[!df$term %in% c("biological_process", "cellular_component", "molecular_function"), ]


# ------------------------------------------------------------
# 3. Filter significant terms and select top 100
# ------------------------------------------------------------
df <- df[df$over_represented_FDR < 0.05, ]
df <- head(df[order(df$over_represented_FDR), ], 100)


# ------------------------------------------------------------
# 4. Build co-occurrence matrix (Jaccard index)
# ------------------------------------------------------------
gene_list <- strsplit(df$gene_ids, ", ")
names(gene_list) <- df$term

similarity_matrix <- matrix(0, nrow = length(gene_list), ncol = length(gene_list))
rownames(similarity_matrix) <- colnames(similarity_matrix) <- names(gene_list)

for (i in 1:length(gene_list)) {
  for (j in 1:length(gene_list)) {
    genes_i <- gene_list[[i]]
    genes_j <- gene_list[[j]]
    similarity_matrix[i, j] <- length(intersect(genes_i, genes_j)) / length(union(genes_i, genes_j))
  }
}


# ------------------------------------------------------------
# 5. Build graph object from similarity matrix
# ------------------------------------------------------------
g <- graph_from_adjacency_matrix(similarity_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)

# Filter weak connections (optional threshold)
E(g)$weight[E(g)$weight < 0.1] <- 0
g <- delete.edges(g, which(E(g)$weight == 0))

# Remove low-degree nodes to reduce noise
g <- delete_vertices(g, which(degree(g) < 2))


# ------------------------------------------------------------
# 6. Detect communities (clusters) and add attributes
# ------------------------------------------------------------
clusters <- cluster_louvain(g)
V(g)$group <- as.factor(membership(clusters))
V(g)$size <- df$numDEInCat[match(V(g)$name, df$term)] * 2
V(g)$color <- -log10(df$over_represented_FDR[match(V(g)$name, df$term)])

# Export cluster assignments
cluster_membership <- data.frame(
  term = V(g)$name,
  cluster = as.integer(V(g)$group)
)
cluster_membership <- cluster_membership[order(cluster_membership$cluster), ]
write.csv(cluster_membership, "GO_clusters.csv", row.names = FALSE)


# ------------------------------------------------------------
# 7. Generate the network plot
# ------------------------------------------------------------
set.seed(456)
layout_fr <- create_layout(g, layout = "fr")

graph_plot <- ggraph(layout_fr) +
  geom_edge_link(aes(edge_alpha = 4), color = "#000714") +
  geom_mark_ellipse(aes(x = x, y = y, group = group, fill = group), alpha = 0.3) +
  geom_node_point(aes(size = size, color = color, stroke = 0.5)) +
  scale_size_continuous(range = c(10, 20)) +
  geom_label_repel(
    data = layout_fr %>% filter(size > 50),
    aes(x = x, y = y, label = name),
    size = 8,
    box.padding = 0.4,
    point.padding = 0.8,
    segment.colour = "grey10",
    segment.size = 1.2,
    fill = NA,
    label.size = 0,
    color = "black",
    max.overlaps = 0,
    family = "Verdana"
  ) +
  scale_color_gradientn(colors = c("blue", "purple", "orange", "red")) +
  scale_fill_manual(values = c("#F17F29", "#00A6ED", "#97cc04", "#F4E409")) +
  theme_minimal() +
  theme(
    legend.position = "left",
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 36),
    legend.key.size = unit(2.5, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


# ------------------------------------------------------------
# 8. Save output
# ------------------------------------------------------------
ggsave(
  "GO_network_filtered_label.jpg",
  plot = graph_plot,
  width = 30,
  height = 30,
  units = "in",
  dpi = 300
)


# ------------------------------------------------------------
# 9. Display plot
# ------------------------------------------------------------
print(graph_plot)


# ============================================================
# Notes:
# - Adjust input path to your GO enrichment results file.
# - Jaccard similarity threshold (< 0.1) is customizable.
# ============================================================
