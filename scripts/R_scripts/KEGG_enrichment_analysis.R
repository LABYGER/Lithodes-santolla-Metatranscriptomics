# ----------------------------
# Install required libraries (if not installed)
# ----------------------------

required_packages <- c("clusterProfiler", "dplyr", "tidyr", "ggplot2", "readr", "patchwork")

new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load libraries
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(patchwork)

# ------------------------
# Load required libraries
# ------------------------
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(patchwork)


# -----------------------------------
# Load DEG results and annotations
# -----------------------------------

# Differential expression results
degs <- read.delim("path/to/DESeq2_results.txt", header = TRUE, sep = "\t")
colnames(degs)[1] <- "Gene"

# KEGG annotations (e.g., from EggNOG)
annotations <- read.delim("path/to/annotations.txt", header = TRUE, sep = "\t", quote = "")
annotations$query <- gsub("\\.p\\d+", "", annotations$query)

kegg_anno <- annotations %>%
  select(Gene = query, KO = KEGG_ko) %>%
  filter(!is.na(KO), KO != "", KO != "-") %>%
  mutate(KO = gsub("ko:", "", KO))


# ------------------------------
# Filter DEGs by threshold
# ------------------------------
degs_filtered <- degs %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2)

# Separate upregulated genes by condition
upregulated_A <- degs_filtered %>% filter(log2FoldChange > 1)
upregulated_B <- degs_filtered %>% filter(log2FoldChange < -1)


# -------------------------------------------
# KEGG enrichment setup
# -------------------------------------------
degs_A <- inner_join(upregulated_A, kegg_anno, by = "Gene")
degs_B <- inner_join(upregulated_B, kegg_anno, by = "Gene")

ko_list_A <- degs_A$KO
ko_list_B <- degs_B$KO

organism_code <- "ko"  # Non-model organisms
universe_genes <- unique(kegg_anno$KO)


# -------------------------------------------
# Define KEGG category colors for plotting
# -------------------------------------------
category_colors <- c(
  "Cellular Processes" = "#FF715B",
  "Environmental Information Processing" = "#8ac926",
  "Genetic Information Processing" = "#F6AE2D",
  "Metabolism" = "#05A8AA",
  "Organismal Systems" = "#8963BA"
)


# -------------------------------------------
# Dataset 1 enrichment
# -------------------------------------------
upregulated_dataset1 <- read_delim("path/to/dataset1_genes.txt", col_names = TRUE)
colnames(upregulated_dataset1)[1] <- "Gene"

degs_dataset1 <- inner_join(upregulated_dataset1, kegg_anno, by = "Gene")
ko_list_dataset1 <- degs_dataset1$KO

kegg_dataset1 <- enrichKEGG(gene = ko_list_dataset1, organism = organism_code, keyType = "kegg", pvalueCutoff = 0.05)
write.csv(as.data.frame(kegg_dataset1@result), "KEGG_enrichment_dataset1.csv", row.names = FALSE)


# -------------------------------------------
# Dataset 2 enrichment
# -------------------------------------------
upregulated_dataset2 <- read_delim("path/to/dataset2_genes.txt", col_names = TRUE)
colnames(upregulated_dataset2)[1] <- "Gene"

degs_dataset2 <- inner_join(upregulated_dataset2, kegg_anno, by = "Gene")
ko_list_dataset2 <- degs_dataset2$KO

kegg_dataset2 <- enrichKEGG(gene = ko_list_dataset2, organism = organism_code, keyType = "kegg", pvalueCutoff = 0.05)
write.csv(as.data.frame(kegg_dataset2@result), "KEGG_enrichment_dataset2.csv", row.names = FALSE)


# ---------------------------------------------------
# PLOTS: Dataset 1 vs. Dataset 2 (Bubble plots)
# ---------------------------------------------------

filtered_dataset1 <- kegg_dataset1@result %>%
  filter(p.adjust <= 0.05) %>%
  mutate(Group = "Dataset_1") %>%
  arrange(p.adjust) %>%
  slice_head(n = 20)

filtered_dataset2 <- kegg_dataset2@result %>%
  filter(p.adjust <= 0.05) %>%
  mutate(Group = "Dataset_2") %>%
  arrange(p.adjust) %>%
  slice_head(n = 20)

filtered_data <- bind_rows(filtered_dataset1, filtered_dataset2) %>%
  filter(category != "Human Diseases") %>%
  mutate(
    Description = factor(Description, levels = rev(unique(Description))),
    logFDR = -log10(p.adjust),
    logFDR = scales::rescale(logFDR, to = c(0.5, 5)),
    logFDR = ifelse(logFDR > 10, 10, logFDR)
  )


p1 <- ggplot(filtered_data %>% filter(Group == "Dataset_1"),
             aes(x = Group, y = Description, color = category, size = logFDR)) +
  geom_point() +
  scale_color_manual(values = category_colors) +
  scale_size_continuous(name = "-log10(FDR)", range = c(5, 10)) +
  labs(title = "A) KEGG Dataset 1", x = "", y = "KEGG Pathways", color = "KEGG Category") +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black", face = "bold"),
    legend.text = element_text(size = 16, color = "black", face = "bold"),
    legend.title = element_text(size = 16, color = "black", face = "bold"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )


p2 <- ggplot(filtered_data %>% filter(Group == "Dataset_2"),
             aes(x = Group, y = Description, color = category, size = logFDR)) +
  geom_point() +
  scale_color_manual(values = category_colors) +
  scale_size_continuous(name = "-log10(FDR)", range = c(5, 10)) +
  labs(title = "B) KEGG Dataset 2", x = "", y = "KEGG Pathways", color = "KEGG Category") +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    axis.text = element_text(size = 16, color = "black", face = "bold"),
    legend.text = element_text(size = 16, color = "black", face = "bold"),
    legend.title = element_text(size = 16, color = "black", face = "bold"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )


# Combine plots
final_plot <- (p1 / p2) + plot_layout(guides = "collect") & 
  theme(legend.position = "right", legend.box = "vertical")


# Export final figure
ggsave("KEGG_Dataset1_Dataset2.png", 
       plot = final_plot, 
       width = 12, 
       height = 14, 
       dpi = 600)
