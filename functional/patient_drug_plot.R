
# Collection of functions using pheatmap and ggplot to
# make nice drug-specific heatmaps and PCA.

#' @param m (MSnSet) msnset for plotting
#' @param funcData (data.frame) patient response data
#' @param data_type (character) either `"proteinLevels"` or `"Phosphosite"`
#' @param drug (character) drug/inhibitor
#' @param model (character) either LASSO or Logistic
#' @param filename (character) filename to save
#' @param ... additional arguments to pheatmap
#' 
#' @return (pheatmap) heatmap for plotting
library(MSnbase)
library(pheatmap)
plot_patient_drug_heatmap <- function(m, funcData, data_type, drug, model, filename=NA, ...) {
  
  m <- attach_msnset_with_drug_data(m, funcData, data_type, drug,model)
  x <- exprs(m)
  ann.dat <- get_annotation_data(m)

  
  cellheight = max(10, 40 - 3*dim(x)[1]) # makes cells larger
  # when there are less of them.
  if (!is.na(filename)) message("Saving ", filename)
  
  
  try(pheatmap(x,
               cellwidth = 10,
               cellheight = cellheight,
               annotation_col = ann.dat,
               # only cluster if there is >1 data points
               cluster_rows = (dim(x)[1] > 1),
               cluster_cols = (dim(x)[2] > 1),
               clustering_distance_cols = 'correlation',
               clustering_method = 'ward.D2',
               filename=filename, ...),
      silent= TRUE)
}

#' @param m (MSnSet) msnset for plotting
#' @param funcData (data.frame) patient response data
#' @param data_type (character) either `"proteinLevels"` or `"Phosphosite"`
#' @param drug (character) drug/inhibitor
#' @param model (character) either LASSO or Logistic
#' @param filename (character) filename to save
#' @param ... additional arguments to pheatmap
#' 
#' @return None
library(MSnbase)
library(ggplot2)
library(RColorBrewer)
plot_patient_drug_pca <- function(m, funcData, data_type, drug, model, filename=NA, ...) {
  
  m <- attach_msnset_with_drug_data(m, funcData, data_type, drug,model)
  x <- exprs(m)
  ann.dat <- get_annotation_data(m)
  
  
  
  z <- t(x)
  z <- sweep(z, 1, rowMeans(z), FUN = "-")
  z <- sweep(z, 1, apply(z, 1, sd), FUN = "/")
  
  print(dim(z))
  pca1 <- prcomp(z, scale. = F)
  scores <- as.data.frame(pca1$x)
  exp_var <- 100 * summary(pca1)$importance[2, ][c(1, 2)]
  
  axes <- paste0("PC", c(1, 2))
  axes <- paste0(axes, " (", round(exp_var, 2), "%)")
  
  # needs to contain AUC, FLT3 variant, Resistance
  ann.dat <- cbind(ann.dat, scores)
  ggdata <- ann.dat
  
  p <- ggplot(ann.dat)
  if ("FLT3 Variant" %in% names(ann.dat)) {
    p <- p + 
      geom_point(aes(x = PC1, y = PC2, color = AUC, shape=`FLT3 Variant`),
                 size=2)
  } else {
    p <- p + 
      geom_point(aes(x = PC1, y = PC2, color = AUC),
                 size=2)
  }
  p <- p +
    scale_colour_gradientn(colours = colorRampPalette(rev(brewer.pal(11, "Spectral")))(8)) +
    stat_ellipse(aes(x = PC1, y = PC2, fill = Resistance), 
                 geom = "polygon", type = "norm", level = 0.5, 
                        alpha = 0.1, show.legend = TRUE) #+ guides(color = guide_legend(phenotype_str), 
                                                       #           fill = guide_legend(phenotype_str))
  if (!is.na(filename)) {
    ggsave(filename)
  }
  return(p)
}






#' @param feature_table (data.frame) contains output of previous models
#' @param output_dir (character) output path for figures
#' @param suffix (character) file name description
#' @param fn (function) plotting function to call
#' @param (...) additional arguments to `fn`
#' 
#' @return None
make_plots <- function(feature_table,
                       output_dir,
                       suffix, fn, ...) {
  filename.l <- file.path(output_dir,
                          paste(feature_table$Drug,
                                feature_table$Model,
                                "predictors",
                                paste0(suffix, ".pdf"),
                                sep="_"))
  .l <- list(drug = feature_table$Drug,
             features = feature_table$features,
             model = feature_table$Model,
             filename = filename.l)
  purrr::pmap(.l = .l,
              .f = fn,
              ...)
}