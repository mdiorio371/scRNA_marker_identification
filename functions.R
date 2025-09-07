# Marker gene overlap functions
library(tidyverse)
library(reshape2)

computeMarkersOverlap <- function(sample_1, sample_2, denom = "Sample_1") {
  # Your existing function code here (simplified version)
  sample_1 <- split(sample_1, sample_1$cluster)
  sample_2 <- split(sample_2, sample_2$cluster)
  
  overlap_matrix <- matrix(0, 
                           nrow = length(sample_1), 
                           ncol = length(sample_2))
  
  for (i in seq_along(sample_1)) {
    for (j in seq_along(sample_2)) {
      genes_1 <- unique(sample_1[[i]]$gene)
      genes_2 <- unique(sample_2[[j]]$gene)
      
      shared <- length(intersect(genes_1, genes_2))
      
      if (denom == "Sample_1") {
        total <- length(genes_1)
      } else {
        total <- length(genes_2)
      }
      
      overlap_matrix[i, j] <- if(total > 0) shared/total else 0
    }
  }
  
  rownames(overlap_matrix) <- names(sample_1)
  colnames(overlap_matrix) <- names(sample_2)
  
  melt(overlap_matrix)
}

heatmapMarkersOverlap <- function(markers_overlap) {
  ggplot(markers_overlap, aes(Var2, Var1)) +
    geom_tile(aes(fill = value)) +
    geom_text(aes(label = round(value, 2))) +
    scale_fill_gradient(low = "white", high = "red", limits = c(0,1)) +
    labs(x = "Dataset 2 Clusters", y = "Dataset 1 Clusters",
         title = "Marker Gene Overlap") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}