
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


plot_patient_drug_heatmap <- function(m, data_type, drug, model, filename=NA, ...) {
  
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


x <- feature_table



# Plot heatmaps
purrr::pmap(list(drug=x$Drug,
                 features=x$features,
                 model=x$Model,
                 filename=file.path("patient_drug_heatmaps",
                                    paste(x$Drug, x$Model,
                                          "proteinLevels_predictors_heatmap.pdf", sep="_"))),
            plot_patient_drug_heatmap,
            m = msnset_gl,
            data_type="proteinLevels")

purrr::pmap(list(drug=x$Drug,
                 features=x$features,
                 model=x$Model,
                 filename=file.path("patient_drug_heatmaps",
                                    paste(x$Drug, x$Model,
                                          "phosphosite_predictors_heatmap.pdf", sep="_"))),
            plot_patient_drug_heatmap,
            m = msnset_ph,
            data_type="Phosphosite")



