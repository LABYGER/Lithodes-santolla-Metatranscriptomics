# ----------------------------
# Install required libraries (if not installed)
# ----------------------------
required_packages <- c("tidyverse", "dplyr", "ggplot2", "paletteer", "patchwork", "tidyr")

new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(paletteer)
library(patchwork)
library(tidyr)


# -----------------------------------------
# Load and transform data for plotting
# -----------------------------------------

# Replace with the path to your data
df <- read.table("path/to/Percentagephyla.txt", header = TRUE, sep = "\t")

# Order phyla for consistent plotting
ordered_phyla <- df$Phylum[!duplicated(df$Phylum)]
df$Phylum <- factor(df$Phylum, levels = ordered_phyla)

# Reshape data: from wide (columns = sites) to long format
df_long <- df %>%
  pivot_longer(cols = c(Site_A, Site_B), names_to = "Site", values_to = "RawCount") %>%
  group_by(Group, Site) %>%
  mutate(
    Total = sum(RawCount),
    Percent = RawCount / Total * 100
  ) %>%
  ungroup()


# -----------------------------------------
# Define color palette
# -----------------------------------------
# Using paletteer for diverse color palettes
my_colors <- c(
  paletteer_d("miscpalettes::pastel")[1:8],
  paletteer_d("rcartocolor::Bold")[1:5],
  paletteer_d("basetheme::clean")[1:5],
  paletteer_d("rcartocolor::Vivid")[1:5],
  paletteer_d("ggthemes::calc")[1:11]
)


# -----------------------------------------
# Create individual barplots for each Group
# -----------------------------------------
groups <- unique(df_long$Group)

plots <- lapply(groups, function(g) {
  df_group <- df_long %>% filter(Group == g)
  n_phyla <- length(unique(df_group$Phylum))
  
  # Assign colors dynamically to Phylum factor
  names(my_colors) <- levels(df_group$Phylum)
  
  ggplot(df_group, aes(x = Site, y = Percent, fill = Phylum)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = my_colors) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(
      title = g,
      x = NULL,
      y = "Relative abundance (%)",
      fill = "Phylum"
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


# -----------------------------------------
# Combine plots and export
# -----------------------------------------
final_plot <- wrap_plots(plots, ncol = 2)

ggsave(
  filename = "phyla_relative_abundance.jpg",
  plot = final_plot,
  width = 14.5,
  height = 14,
  units = "in",
  dpi = 300
)
