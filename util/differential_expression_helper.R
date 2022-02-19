library(dplyr)
library(ggplot2)
library(MSnSet.utils)
library(limma)
library(pheatmap)
library(gridExtra)
library(stringr)


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
      dplyr::select(feature, logFC, P.Value, adj.P.Val, datatype) %>%
      ## If the mutated group happens to have missing data in a feature, 
      ## that row gives NA from limma. So we filter those out.
      filter(!is.na(P.Value))
    sig_features <- diff_exp %>%
      filter(adj.P.Val < 0.05)
    col_name = paste0("Differentially expressed in ", names(datasets)[[i]])
    differentially_expressed[[col_name]] <- nrow(sig_features)
    all_results <- rbind(all_results, diff_exp)
  }

  top_features <- all_results
  
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
                                 plot_features, folder, suffix, ...){
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
    dplyr::select(InitialAMLDiagnosis, PostChemotherapy, mutation)
  
  annotation_colors <- lapply(colnames(annotation_df), function(cat){
    if (cat == "mutation"){
      colors_mutation
    } else {
      colors_other
    }
  })
  names(annotation_colors) <- colnames(annotation_df)
  
  plot_title <- paste(mutation, suffix)
  path <- paste0(folder, "differentially_expressed_features_in_", mutation, 
                 "_mutated_samples_", suffix, ".png")
  pheatmap(heatmap_mat, annotation_col = annotation_df, file = path, main = plot_title,
           annotation_colors = annotation_colors, show_colnames = FALSE, ...)
}


plot_enrichment_result <- function(enrichment_result, enrichment_label, 
                                   enrich_width = 9, sig_width = 4,
                                   enrichment_title = "Enrichment"){
  
  #' Used to make reversed logarithmic scales
  #' @import scales
  reverselog_trans <- function(base = exp(1)) {
    library(scales)
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    scales::trans_new(paste0("reverselog-", format(base)), trans, inv,
                      scales::log_breaks(base = base),
                      domain = c(1e-100, Inf))
  }
  
  
  plot_df <- enrichment_result %>%
    dplyr::select(pathway, adj_p_val, enrichment) %>%
    mutate(sign = case_when(enrichment > 0 ~ "Up",
                            TRUE ~ "Down"))
  
  p_enrichment <- ggplot(plot_df, aes(x = enrichment, y = reorder(pathway, enrichment))) +
    geom_bar(stat = 'identity', aes(fill = sign)) +
    scale_fill_manual(values = c("Down" = "dodgerblue3", "Up" = "firebrick2", "Not significant" = "black")) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 55)) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
          legend.position = "none",
          axis.title.x = element_text(size=16),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 13),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(x = enrichment_label) +
    ggtitle(enrichment_title)
  
  
  p_sig <- ggplot(plot_df, aes(x = adj_p_val, y = reorder(pathway, enrichment))) +
    geom_bar(stat='identity') +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
          legend.position="none",
          axis.title.x = element_text(size=16),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 14),
          axis.ticks.y = element_blank(),
          axis.line.y = element_line(color = "black"),
          axis.text.y = element_blank()) +
    scale_x_continuous(trans = reverselog_trans(10)) +
    labs(x = "Adjusted p-value") +
    ggtitle("Significance")
  
  arrange_matrix <- as.matrix(c(rep(1, enrich_width), rep(2, sig_width))) %>% t()
  p_both <- grid.arrange(p_enrichment, p_sig, layout_matrix = arrange_matrix)
}








