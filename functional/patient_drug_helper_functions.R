
# Collection of helper functions for patient drug plots.


#' @param funcData (data.frame) patient functional response data
#' @param drug (character) name of drug/inhibitor
#' 
#' @return single plate and multiplate AUC data for the given
#' drugs, with lab IDs as rownames
library(dplyr)
fetch_AUC_drug_data <- function(funcData, drug) {
  drug.dat <- funcData %>%
    dplyr::select(lab_id, inhibitor, auc, type) %>%
    dplyr::mutate(Resistance = if_else(auc >= 100,
                                       "Resistant",
                                       "Sensitive")) %>%
    # Select only single and multi plate
    dplyr::filter(type != "qc_salvage_data",
                  inhibitor == drug)
  
  # Check if the mapping was successful
  if (any(duplicated(drug.dat$lab_id))) {
    warning("Some duplicate lab IDs.")
  }
  if (any(is.na(drug.dat$auc))) {
    warning("Some AUC values missing.")
  }
  rownames(drug.dat) <- drug.dat$lab_id
  return(drug.dat)
}


#' @param funcData (data.frame) patient functional response data
#' @param drug (character) name of drug/inhibitor
#' 
#' @return (character) vector of feature names
library(dplyr)
fetch_features <- function(feature_table, data_type, drug, model) {
  x <- feature_table %>%
    filter(Drug == drug,
           Model == model,
           dataType == data_type)
  if (nrow(x) > 1) {
    stop("More than one row after filtering.")
  }
  features <- unlist(strsplit(x$features, ","))
  
}

#' @param m (MSnSet) contains the crosstab and patient metadata
#' @param funcData (data.frame) patient response data
#' @param data_type (character) either `"proteinLevels"` or `"Phosphosite"`
#' @param drug (character) name of drug/inhibitor
#' @param model (character) either `"LASSO"` or `"Logistic"`
#' 
#' @return (MSnSet) another MSnSet object with AUC data
library(MSnbase)
library(dplyr)
attach_msnset_with_drug_data <- function(m, funcData,
                                         data_type,
                                         drug,
                                         model) {
  drug.data <- fetch_AUC_drug_data(funcData, drug=drug)
  
  features <- fetch_features(feature_table, data_type=data_type, drug=drug, model=model)
  
  m <- m[complete.cases(exprs(m)),]
  m <- m[intersect(features, featureNames(m)),]
  m <- m[,intersect(rownames(drug.data), sampleNames(m))]
  
  p <- inner_join(pData(m), drug.data,
                  by=c("Barcode.ID"="lab_id"))
  rownames(p) <- p$Barcode.ID
  pData(m) <- p[sampleNames(m),]
  
  
  return(m)
}

#' @param m (MSnSet)
#' 
#' @return (data.frame) annotation data for plotting
library(MSnbase)
library(dplyr)
get_annotation_data <- function(m) {
  ann.dat <- pData(m) %>%
    dplyr::select(auc, Resistance, FLT3.ITD) %>%
    mutate(Resistance = as.factor(Resistance)) %>%
    mutate(FLT3.ITD = as.factor(if_else(FLT3.ITD,
                                        "ITD Mutant",
                                        "Wild Type"))) %>%
    dplyr::rename(`FLT3 Variant` = FLT3.ITD) %>%
    dplyr::rename(AUC = auc)
  # Extra step: If the dataset contains only one
  # FLT3 variant, drop that column
  if (length(unique(ann.dat$`FLT3 Variant`)) == 1) {
    ann.dat <- ann.dat %>% select(-`FLT3 Variant`)
  }
  return(ann.dat)
}
