# scRNA Marker Identification

A Shiny app for comparing marker genes between single-cell RNA-seq clusters to validate cell type annotations.

## Features
- Compare marker genes between two datasets
- Visualize overlap with interactive heatmaps
- Export results for publication

## Quick Start

```r
# Install dependencies
install.packages(c("shiny", "tidyverse", "reshape2", "DT"))

# Run the app
shiny::runApp()
