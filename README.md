# scRNA Marker Identification

A Shiny app for comparing marker genes between single-cell RNA-seq clusters to validate cell type annotations. 

## Features
- **Compare marker genes** between two datasets
- **Interactive heatmaps** for visualizing marker overlap
- **Machine Learning predictions** using a Random Forest classifier
- **Confidence scoring** for cell type assignments
- **Export results**

## Quick Start

```r
# Install dependencies
install.packages(c("shiny", "tidyverse", "reshape2", "DT"))

# Run the app from the project root
shiny::runApp()
```


### Input 
The app expects a CSV file with the following columns:
- cluster: Cluster identifier (integer or string)
- gene: Gene symbol
- avg_log2FC (optional): Average log fold change
- p_val_adj (optional): Adjusted p-value

### Basic Overlap Analysis
1. Upload two CSV files with marker genes
2. Click "Analyze Overlap" 
3. View heatmap (red = high overlap, white = no overlap)
4. Download results as CSV

### Machine Learning Prediction
1. Upload reference dataset with known cell types
2. Click "Train Model" to train Random Forest classifier
3. View predictions with confidence scores
4. Export predicted cell type annotations

## Methods
- **Overlap Analysis**: Jaccard index and overlap coefficient
- **Machine Learning**: Random Forest with feature engineering including:
   - Gene set enrichment scores
   - Pathway activation patterns
   - Marker specificity metrics
   - Expression level features

