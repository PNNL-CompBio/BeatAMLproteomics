
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
source("patient_drug_plot.R")

msnset_gl <- load_beataml_global_msnset()
msnset_ph <- load_beataml_phospho_msnset()
funcData <- querySynapseTable("syn25830473")
feature_table <- read_xlsx("supplemental_tables.xlsx",
                           sheet=2)


x <- feature_table

purrr::pmap(list(drug=x$Drug,
                 features=x$features,
                 model=x$Model,
                 filename=file.path(output_dir,
                                    paste(x$Drug, x$Model,
                                          "proteinLevels_predictors_heatmap.pdf", sep="_"))),
            plot_patient_drug_heatmap,
            m = msnset_gl,
            funcData=funcData,
            data_type="proteinLevels")

purrr::pmap(list(drug=x$Drug,
                 features=x$features,
                 model=x$Model,
                 filename=file.path(output_dir,
                                    paste(x$Drug, x$Model,
                                          "phosphosite_predictors_heatmap.pdf", sep="_"))),
            plot_patient_drug_heatmap,
            m = msnset_ph,
            funcData=funcData,
            data_type="Phosphosite")




purrr::pmap(list(drug=x$Drug,
                 features=x$features,
                 model=x$Model,
                 filename=file.path(output_dir,
                                    paste(x$Drug, x$Model,
                                          "proteinLevels_predictors_pca.pdf", sep="_"))),
            function(...) try(plot_patient_drug_pca(...)),
            m = msnset_gl,
            funcData=funcData,
            data_type="proteinLevels")

purrr::pmap(list(drug=x$Drug,
                 features=x$features,
                 model=x$Model,
                 filename=file.path(output_dir,
                                    paste(x$Drug, x$Model,
                                          "phosphosite_predictors_pca.pdf", sep="_"))),
            function(...) try(plot_patient_drug_pca(...)),
            m = msnset_ph,
            funcData=funcData,
            data_type="Phosphosite")
