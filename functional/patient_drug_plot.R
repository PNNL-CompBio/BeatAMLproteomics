
# Script to make patient/drug heatmaps and PCA plots in the 
# chosen output directory.
output_dir <- "./patient_drug_plots/"
if (!dir.exists(output_dir)) dir.create(output_dir)


library(dplyr)
library(ggplot2)
library(ggrepel)
library(MSnbase)
library(purrr)
library(readxl)
library(pheatmap)

source("../load_beatAML_data.R")
source("patient_drug_helper_functions.R")


msnset_gl <- load_beataml_global_msnset()
msnset_ph <- load_beataml_phospho_msnset()
funcData <- querySynapseTable("syn25830473")
feature_table <- read_xlsx("supplemental_tables.xlsx",
                           sheet=2)

#' @param m (MSnSet) msnset for plotting
#' @param funcData (data.frame) patient response data
#' @param data_type (character) either `"proteinLevels"` or `"Phosphosite"`
#' @param drug (character) drug/inhibitor
#' @param model (character) either LASSO or Logistic
#' @param filename (character) filename to save
#' @param ... additional arguments to pheatmap
#' 
#' @return (pheatmap) heatmap for plotting
plot_patient_drug_heatmap <- function(m, funcData, data_type, drug, model, filename=NA, ...) {
  
  m <- attach_msnset_with_drug_data(m, funcData, data_type, drug,model)
  x <- exprs(m)
  ann.dat <- get_annotation_data(m)
  
  
  if (is.na(filename)) {
    filename <- file.path("patient_drug_heatmaps",
                          paste(drug_i, model_i,
                                "predictors_heatmap.pdf", sep="_"))
  }
  
  cellheight = max(10, 40 - 3*dim(x)[1]) # makes cells larger
  # when there are less of them.
  
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
  message("Saving ", filename)
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
  
  
  pdf(filename)
  message("Saving ", filename)
  p <- ggplot(ann.dat) + geom_point(aes(x = PC1, y = PC2, color = AUC, shape=`FLT3 Variant`), size=2) +
    scale_colour_gradientn(colours = colorRampPalette(rev(brewer.pal(11, "Spectral")))(8))
  print(p)
  dev.off()
  
  
  
  
  
  
  
  x <- feature_table
  
  purrr::pmap(list(drug=x$Drug,
                   features=x$features,
                   model=x$Model,
                   filename=file.path(output_dir,
                                      paste(x$Drug, x$Model,
                                            "proteinLevels_predictors_heatmap.pdf", sep="_"))),
              plot_patient_drug_heatmap,
              m = msnset_gl,
              data_type="proteinLevels")
  
  purrr::pmap(list(drug=x$Drug,
                   features=x$features,
                   model=x$Model,
                   filename=file.path(output_dir,
                                      paste(x$Drug, x$Model,
                                            "phosphosite_predictors_heatmap.pdf", sep="_"))),
              plot_patient_drug_heatmap,
              m = msnset_ph,
              data_type="Phosphosite")
  
  
  
  
  purrr::pmap(list(drug=x$Drug,
                   features=x$features,
                   model=x$Model,
                   filename=file.path(output_dir,
                                      paste(x$Drug, x$Model,
                                            "proteinLevels_predictors_pca.pdf", sep="_"))),
              function(...) try(plot_patient_drug_pca(...)),
              m = msnset_gl,
              data_type="proteinLevels")
  
  purrr::pmap(list(drug=x$Drug,
                   features=x$features,
                   model=x$Model,
                   filename=file.path(output_dir,
                                      paste(x$Drug, x$Model,
                                            "phosphosite_predictors_pca.pdf", sep="_"))),
              function(...) try(plot_patient_drug_pca(...)),
              m = msnset_ph,
              data_type="Phosphosite")