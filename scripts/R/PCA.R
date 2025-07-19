# ============================================================
# PCA with automatic variance explained calculation (PC1, PC2)
# ============================================================
# This script performs PCA on a numeric matrix, extracts the 
# percentage of variance explained by PC1 and PC2 automatically,
# and generates a PCA plot with sample labels and group ellipses.
# ============================================================


# ------------------------------------------------------------
# 1. Load required libraries
# ------------------------------------------------------------
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2")
if (!requireNamespace("ggpubr", quietly = TRUE)) stop("Please install ggpubr")
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Please install ggrepel")

library(ggplot2)
library(ggpubr)
library(ggrepel)


# ------------------------------------------------------------
# 2. Input data: Replace with your real numeric data
# ------------------------------------------------------------
# Example data: 10 samples, 5 variables (numeric)
set.seed(123)
data_matrix <- matrix(rnorm(50), nrow = 10)  # 10 samples x 5 variables
rownames(data_matrix) <- paste0("Sample_", 1:10)

# Optional: Create group labels
sample_groups <- rep(c("GroupA", "GroupB"), each = 5)


# ------------------------------------------------------------
# 3. PCA calculation
# ------------------------------------------------------------
pca_result <- prcomp(data_matrix, center = TRUE, scale. = TRUE)

# Extract % variance explained
explained_variance <- (pca_result$sdev)^2 / sum(pca_result$sdev^2) * 100
pc1_var <- round(explained_variance[1], 2)
pc2_var <- round(explained_variance[2], 2)

cat("PC1 explains:", pc1_var, "%\n")
cat("PC2 explains:", pc2_var, "%\n")


# ------------------------------------------------------------
# 4. Prepare PCA scores for plotting
# ------------------------------------------------------------
pca_data <- as.data.frame(pca_result$x)
pca_data$Sample <- rownames(pca_data)
pca_data$Group <- sample_groups


# ------------------------------------------------------------
# 5. Define colors and shapes
# ------------------------------------------------------------
point_colors <- c("GroupA" = "#1f78b4", "GroupB" = "#33a02c")
ellipse_colors <- point_colors
point_shapes <- c("GroupA" = 16, "GroupB" = 17)


# ------------------------------------------------------------
# 6. Generate PCA plot
# ------------------------------------------------------------
p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.1, level = 0.95) +
  geom_text_repel(aes(label = Sample), size = 4) +
  scale_color_manual(values = point_colors) +
  scale_fill_manual(values = ellipse_colors) +
  scale_shape_manual(values = point_shapes) +
  theme_minimal() +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  ) +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 18, face = "bold")
  )


# ------------------------------------------------------------
# 7. Display plot
# ------------------------------------------------------------
print(p)


# ------------------------------------------------------------
# 8. Save as high-resolution image
# ------------------------------------------------------------
ggsave("PCA_plot_auto_variance.jpg", plot = p, width = 8, height = 6, dpi = 300)


# ============================================================
# Notes:
# - Replace `data_matrix` with your real numeric dataset.
# - Group labels are arbitrary in this example (GroupA / GroupB).
# - PCA % variance is calculated directly from prcomp output.
# ============================================================
