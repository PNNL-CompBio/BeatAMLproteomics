library(dplyr)
library(ggplot2)
library(MSnSet.utils)
library(limma)
library(pheatmap)


## Runs differential expression comparing mutated vs non mutated samples
## in the given datasets. Limma is used for this part. The significant
## features are collected. Output is a small table containing the summary
## of the number of significant features found in each dataset. Meant
## to be used in lapply to run many mutations at once.

check_diff_exp_individual <- function(mutation, mutation_data, datasets){
  
  if (is.null(names(datasets))){
    stop("Please name the datasets list.")
  }
  samples <- mutation_data %>%
    filter(Gene == mutation) %>%
    pull(Barcode.ID) %>% 
    unique()
  
  differentially_expressed <- data.frame(mutation = mutation, 
                                         'Number mutated' = length(samples))
  all_results <- data.frame(feature = c(), logFC = c(), adj.P.Val = c(),
                            dataset = c())
  
  for (i in 1:length(datasets)){
    pData(datasets[[i]]) <- pData(datasets[[i]]) %>%
      mutate(mutation = case_when(rownames(.) %in% samples ~ "TRUE",
                                  TRUE ~ "FALSE"))
    diff_exp <- limma_gen(datasets[[i]], "~mutation", "mutation") %>%
      arrange(P.Value) %>%
      mutate(feature = rownames(.), datatype = names(datasets)[[i]]) %>%
      select(feature, logFC, P.Value, adj.P.Val, datatype) %>%
      ## If the mutated group happens to have missing data in a feature, 
      ## that row gives NA from limma. So we filter those out.
      filter(!is.na(P.Value))
    sig_features <- diff_exp %>%
      filter(adj.P.Val < 0.05)
    col_name = paste0("Differentially expressed in ", names(datasets)[[i]])
    differentially_expressed[[col_name]] <- nrow(sig_features)
    all_results <- rbind(all_results, diff_exp)
  }

  top_features <- all_results %>%
    filter(adj.P.Val < 0.05)
  
  return(list("Counts" = differentially_expressed, "Top features" = top_features))
}


## This is meant as a complementary part of the function above. Passing 
## the mutation and mutations_data determines the labeling of the columns
## on the heatmap, as we want to highlight the comparison between
## mutated and non mutated samples. We pass the features
## to plot using the vector plot_features. Finally, suffix should be
## the data type we are using. This is added to the title and file name
## of the heatmap.


plot_diffexp_heatmap <- function(mutation, mutation_data, datasets, 
                                 plot_features, suffix, ...){
  samples <- mutation_data %>%
    filter(Gene == mutation) %>%
    pull(Barcode.ID) %>% 
    unique()
  
  for (i in 1:length(datasets)){
    pData(datasets[[i]]) <- pData(datasets[[i]]) %>%
      mutate(mutation = case_when(rownames(.) %in% samples ~ "TRUE",
                                  TRUE ~ "FALSE"))
  }
  
  heatmap_mat <- lapply(1:length(datasets), function(i){
    mat <- exprs(datasets[[i]])
  }) %>% do.call("rbind", .)
  heatmap_mat <- heatmap_mat[plot_features, ]
  row_means <- apply(heatmap_mat, 1, mean, na.rm = T)
  row_sds <- apply(heatmap_mat, 1, sd, na.rm = T)
  heatmap_mat <- sweep(heatmap_mat, 1, row_means, FUN = '-')
  heatmap_mat <- sweep(heatmap_mat, 1, row_sds, FUN = '/')
    
  colors_mutation <- c("TRUE" = "Black", "FALSE" = "White")
  colors_other <- c("TRUE" = "Black", "FALSE" = "White")
  
  annotation_df <- pData(datasets[[1]]) %>%
    select(InitialAMLDiagnosis, PostChemotherapy, mutation)
  
  annotation_colors <- lapply(colnames(annotation_df), function(cat){
    if (cat == "mutation"){
      colors_mutation
    } else {
      colors_other
    }
  })
  names(annotation_colors) <- colnames(annotation_df)
  
  plot_title <- paste(mutation, suffix)
  path <- paste0("./differentially_expressed_features_in_", mutation, 
                 "_mutated_samples_", suffix, ".png")
  pheatmap(heatmap_mat, annotation_col = annotation_df, file = path, main = plot_title,
           annotation_colors = annotation_colors, show_colnames = FALSE, ...)
}











