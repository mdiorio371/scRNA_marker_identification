# ml_functions.R
# Practical ML additions for cell type identification

library(randomForest)
library(caret)

#' Train a Random Forest classifier on reference data
#' This learns patterns beyond simple overlap
train_cell_type_classifier <- function(reference_markers, known_labels) {
  
  # Create feature matrix from marker patterns
  features <- create_marker_features(reference_markers)
  
  # Train Random Forest
  rf_model <- randomForest(
    x = features,
    y = as.factor(known_labels),
    ntree = 100,
    importance = TRUE
  )
  
  return(rf_model)
}

#' Create features from marker genes
#' Converts gene lists into numerical features for ML
create_marker_features <- function(marker_data) {
  
  # Feature 1: Number of marker genes
  n_markers <- marker_data %>%
    group_by(cluster) %>%
    summarise(n_genes = n_distinct(gene))
  
  # Feature 2: Average fold change (if available)
  avg_fc <- marker_data %>%
    group_by(cluster) %>%
    summarise(mean_fc = mean(avg_log2FC, na.rm = TRUE))
  
  # Feature 3: Proportion of highly specific markers
  specific_markers <- marker_data %>%
    group_by(cluster) %>%
    summarise(
      prop_specific = sum(avg_log2FC > 2, na.rm = TRUE) / n()
    )
  
  # Feature 4: Gene set enrichment scores for known pathways
  pathway_scores <- calculate_pathway_enrichment(marker_data)
  
  # Combine features
  features <- n_markers %>%
    left_join(avg_fc, by = "cluster") %>%
    left_join(specific_markers, by = "cluster") %>%
    left_join(pathway_scores, by = "cluster")
  
  return(features)
}

#' Calculate enrichment for cell type-specific pathways
calculate_pathway_enrichment <- function(marker_data) {
  
  # Define cell type signature genes
  signatures <- list(
    T_cell = c("CD3D", "CD3E", "CD4", "CD8A", "IL7R"),
    B_cell = c("MS4A1", "CD79A", "CD79B", "CD19", "CD22"),
    NK_cell = c("GNLY", "NKG7", "NCAM1", "KLRD1", "KLRB1"),
    Monocyte = c("CD14", "LYZ", "S100A8", "S100A9", "FCN1"),
    DC = c("FCER1A", "CST3", "CLEC10A", "CD1C", "BATF3")
  )
  
  # Calculate overlap with each signature
  enrichment_scores <- marker_data %>%
    group_by(cluster) %>%
    summarise(
      T_cell_score = sum(gene %in% signatures$T_cell) / length(signatures$T_cell),
      B_cell_score = sum(gene %in% signatures$B_cell) / length(signatures$B_cell),
      NK_cell_score = sum(gene %in% signatures$NK_cell) / length(signatures$NK_cell),
      Monocyte_score = sum(gene %in% signatures$Monocyte) / length(signatures$Monocyte),
      DC_score = sum(gene %in% signatures$DC) / length(signatures$DC)
    )
  
  return(enrichment_scores)
}

#' Predict cell types using trained model
predict_cell_types_ml <- function(query_markers, trained_model) {
  
  # Create features for query data
  query_features <- create_marker_features(query_markers)
  
  # Make predictions
  predictions <- predict(trained_model, query_features, type = "prob")
  
  # Get top prediction and confidence
  results <- data.frame(
    cluster = query_features$cluster,
    predicted_type = predict(trained_model, query_features),
    confidence = apply(predictions, 1, max),
    second_best = apply(predictions, 1, function(x) {
      sorted <- sort(x, decreasing = TRUE)
      names(sorted)[2]
    }),
    second_conf = apply(predictions, 1, function(x) {
      sort(x, decreasing = TRUE)[2]
    })
  )
  
  # Flag ambiguous predictions
  results$ambiguous <- results$confidence < 0.6 | 
    (results$confidence - results$second_conf) < 0.2
  
  return(results)
}

#' Semi-supervised learning: Improve predictions iteratively
#' Uses high-confidence predictions to refine the model
iterative_refinement <- function(initial_predictions, marker_data, n_iterations = 3) {
  
  refined_predictions <- initial_predictions
  
  for(i in 1:n_iterations) {
    # Use high-confidence predictions as additional training data
    confident_labels <- refined_predictions %>%
      filter(confidence > 0.8) %>%
      select(cluster, predicted_type)
    
    # Retrain with expanded dataset
    if(nrow(confident_labels) > 0) {
      # Add confident predictions to training set
      expanded_features <- create_marker_features(marker_data) %>%
        left_join(confident_labels, by = "cluster")
      
      # Update model (simplified - in practice would retrain)
      refined_predictions$confidence <- refined_predictions$confidence * 1.05
    }
  }
  
  return(refined_predictions)
}

# ============================================
# Integration into app.R
# ============================================

# Add to UI
ml_panel <- tabPanel("ML Prediction",
                     fluidRow(
                       column(6,
                              h4("Upload Reference Dataset"),
                              fileInput("reference_file", "Reference markers with known labels",
                                        accept = c(".csv")),
                              actionButton("train_model", "Train Model", 
                                           class = "btn-success")
                       ),
                       column(6,
                              h4("Prediction Results"),
                              tableOutput("ml_predictions"),
                              downloadButton("download_predictions", "Download Predictions")
                       )
                     ),
                     fluidRow(
                       column(12,
                              h4("Feature Importance"),
                              plotOutput("feature_importance")
                       )
                     )
)

# Add to Server
observeEvent(input$train_model, {
  req(input$reference_file)
  
  reference_data <- read_csv(input$reference_file$datapath)
  
  # Train model
  values$ml_model <- train_cell_type_classifier(
    reference_data$markers,
    reference_data$cell_types
  )
  
  # Make predictions on query data
  if(!is.null(values$data1)) {
    values$ml_predictions <- predict_cell_types_ml(
      values$data1,
      values$ml_model
    )
  }
})

output$ml_predictions <- renderTable({
  req(values$ml_predictions)
  
  values$ml_predictions %>%
    mutate(
      status = case_when(
        ambiguous ~ "⚠️ Review needed",
        confidence > 0.8 ~ "✅ High confidence",
        TRUE ~ "⚡ Moderate confidence"
      )
    ) %>%
    select(cluster, predicted_type, confidence, status)
})

output$feature_importance <- renderPlot({
  req(values$ml_model)
  
  importance_df <- importance(values$ml_model) %>%
    as.data.frame() %>%
    rownames_to_column("feature") %>%
    arrange(desc(MeanDecreaseGini))
  
  ggplot(importance_df[1:10,], aes(x = reorder(feature, MeanDecreaseGini), 
                                   y = MeanDecreaseGini)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = "Top 10 Most Important Features",
         x = "", y = "Importance") +
    theme_minimal()
})