
# Script to make patient/drug heatmaps and PCA plots in the 
# chosen output directory.

# Prepare output folder
output_dir <- "./patient_drug_plots/"

if (dir.exists(output_dir)) {
  do.call(file.remove, list(list.files(output_dir, full.names = TRUE)))
} else {
  dir.create(output_dir)
}

# Load data

{
library(purrr)
library(readxl)
library(amlresistancenetworks)
  
source("../load_beatAML_data.R")
source("patient_drug_helper_functions.R")
source("patient_drug_plot.R")

msnset_gl <- load_beataml_global_msnset()
msnset_ph <- load_beataml_phospho_msnset()
funcData <- querySynapseTable("syn25830473")
feature_table <- read_xlsx("supplemental_tables.xlsx",
                           sheet=2) %>%
  mutate(Drug = if_else(Drug == "Gilteritinib (ASP-2215)",
                        "Gilteritinib",
                        Drug))
}

# Main loop

for (plot_type in c("heatmap", "pca")) {
  for (data_type in c("proteinLevels", "Phosphosite")) {
    if (data_type == "proteinLevels") {
      m <- msnset_gl
    } else if (data_type == "Phosphosite") {
      m <- msnset_ph
    }
    if (plot_type=="heatmap") {
      fn = plot_patient_drug_heatmap
    } else if (plot_type == "pca") {
      fn = plot_patient_drug_pca
    }
    make_plots(funcData=funcData,
               feature_table=feature_table,
               m = m,
               data_type = data_type,
               output_dir=output_dir,
               suffix=paste(data_type, plot_type,sep="_"),
               fn=function(...) try(fn(...)))
    make_plots(funcData=funcData,
               feature_table=feature_table,
               m = m[,m$FLT3.ITD],
               data_type = data_type,
               output_dir=output_dir,
               suffix=paste("FLT3_ITD",data_type, plot_type, sep="_"),
               fn=function(...) try(fn(...)))
    make_plots(funcData=funcData,
               feature_table=feature_table,
               m = m[,!m$FLT3.ITD],
               data_type = data_type,
               output_dir=output_dir,
               suffix=paste("FLT3_WT",data_type, plot_type,sep="_"),
               fn=function(...) try(fn(...)))
  }
}

