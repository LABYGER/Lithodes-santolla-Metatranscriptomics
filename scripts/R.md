# README

## R scripts included in this repository

This repository contains a collection of R scripts for data analysis and visualization associated with the analyses described. Each script addresses a specific figure or dataset visualization step.

### List of scripts

**1. `PCA_from_data_with_variance_explained.R`**  
Performs Principal Component Analysis (PCA) from a numeric dataset. Automatically calculates the percentage of variance explained by PC1 and PC2 from the PCA output and generates a scatter plot including ellipses and sample labels.

**2. `PCA_with_example_data_auto_variance.R`**  
Example script using simulated data to demonstrate the PCA workflow, from numeric matrix to figure, including variance explained and sample grouping.

**3. `barplots_relative_abundance_by_group.R`**  
Generates stacked bar plots of relative abundance at the Class level, grouped by higher-level taxonomic categories and across sampling sites.

**4. `chord_diagram_GO_categories.R`**  
Creates chord diagrams to illustrate the distribution of categories (e.g., Gene Ontology domains) between conditions or compartments. Includes export of figures in both PNG and PDF formats.

**5. `GO_network_cooccurrence_visualization.R`**  
Constructs a network of Gene Ontology terms based on shared gene membership using Jaccard similarity, applies community detection, and visualizes the network with node attributes reflecting significance and gene counts.

---

## Notes on figures included in the repository

Figures that are not explicitly associated with a dedicated script within this repository were generated using outputs obtained directly from **Trinity** or **VENNY 2.1** and do not require additional R code for their generation. The scripts included here correspond exclusively to analyses and visualizations developed in R as detailed in the Methods section of the manuscript.

---

## Required R packages

```r
# Core visualization and data manipulation packages
ggplot2
dplyr
patchwork
ggpubr
ggrepel
tidyr
paletteer

# Network visualization packages
igraph
ggraph
tidygraph
ggforce
reshape2
RColorBrewer

# Chord diagram package
circlize
