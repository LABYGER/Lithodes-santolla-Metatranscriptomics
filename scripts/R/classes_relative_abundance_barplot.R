# ============================================================
# Stacked Barplots of Relative Abundance by Group and Site
# ============================================================
# This script generates barplots of relative abundance (%) for 
# different taxonomic Classes, grouped by a higher taxonomic level 
# (e.g., Phylum), and separated by Site. 
# Output is saved as a high-resolution image.
# ============================================================


# ------------------------------------------------------------
# 1. Load required libraries
# ------------------------------------------------------------
library(tidyverse)
library(paletteer)
library(patchwork)


# ------------------------------------------------------------
# 2. Read input table (replace with your file path)
# ------------------------------------------------------------
# Expected input format:
# - Columns: 'Group' (e.g., Phylum), 'Class', 'Choiseul', 'Ballena'
# - Values: Raw counts for each Class at each Site

df <- read.table("your_path_here/class.txt", header = TRUE, sep = "\t")


# ------------------------------------------------------------
# 3. Prepare data
# ------------------------------------------------------------
# Keep the original order of Class for plotting consistency
ordered_classes <- df$Class[!duplicated(df$Class)]
df$Class <- factor(df$Class, levels = ordered_classes)

# Pivot wider to longer format and calculate relative abundance (%)
df_long <- df %>%
  pivot_longer(cols = c(Choiseul, Ballena), names_to = "Site", values_to = "RawCount") %>%
  group_by(Group, Site) %>%
  mutate(Total = sum(RawCount),
         Percent = RawCount / Total * 100) %>%
  ungroup()


# ------------------------------------------------------------
# 4. Color palette
# ------------------------------------------------------------
my_colors <- c(
  paletteer_d("ggthemes::calc")[1:11],
  paletteer_d("basetheme::clean")[1:5],
  paletteer_d("miscpalettes::pastel")[1:8],
  paletteer_d("rcartocolor::Vivid")[1:5]
)


# ------------------------------------------------------------
# 5. Generate barplots for each 'Group'
# ------------------------------------------------------------
groups <- unique(df_long$Group)

plots <- lapply(groups, function(g) {
  df_group <- df_long %>% filter(Group == g)
  names(my_colors) <- levels(df_group$Class)
  
  ggplot(df_group, aes(x = Site, y = Percent, fill = Class)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = my_colors) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(
      title = g,
      x = NULL,
      y = "Relative abundance (%)",
      fill = "Class"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      text = element_text(face = "bold", size = 20, color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 22),  
      legend.text = element_text(face = "bold", size = 20, color = "black"),
      axis.text = element_text(face = "bold", size = 20, color = "black"),
      legend.title = element_text(face = "bold", size = 20, color = "black"),
      panel.border = element_rect(color = "black", size = 1.5, fill = NA),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
})


# ------------------------------------------------------------
# 6. Combine plots into a single figure
# ------------------------------------------------------------
final_plot <- wrap_plots(plots, ncol = 2)


# ------------------------------------------------------------
# 7. Save output as high-resolution JPG
# ------------------------------------------------------------
ggsave(
  filename = "relative_abundance_by_group.jpg",
  plot = final_plot,
  width = 14,
  height = 14,
  units = "in",
  dpi = 300
)


# ============================================================
# Notes:
# - Adjust file path for input table (line 26).
# - Sites (columns) can be changed to your dataset.
# - This outputs a publication-quality figure.
# ============================================================
