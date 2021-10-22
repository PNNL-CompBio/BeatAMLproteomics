#' query synapse table
#' This is how you get data from the project
#' @param tableid
#' @export
querySynapseTable<-function(tableid){
  res = read.csv(tableid)
  if('Gene'%in%names(res))
    res$Gene<-unlist(res$Gene)
  if('site'%in%names(res))
    res$site<-unlist(res$site)
  return(res)
}

## Loads phospho data for the 210 patient experiment.
load.phospho.data <- function(filePath_phospho){
  require(dplyr)
  data <- querySynapseTable(filePath_phospho)
  
  meta.cols <- c("Sample", "SampleID.full", "Barcode.ID", 
                 "Plex", "Channel", "Loading.Mass", 
                 "specimen.type", "specimen.location", 
                 "Specimen.access.group.concatenated", 
                 "InitialAMLDiagnosis", "PostChemotherapy", 
                 "FLT3.ITD")
  
  meta <- unique(data[, meta.cols])  %>%
    dplyr::rename(SpecimenType = specimen.type)
  rownames(meta) <- meta$Barcode.ID
  meta <- meta[order(meta$Sample), ]
  
  data <- data %>%
    select(Gene, SiteID, Sample, 'AML Sample' = Barcode.ID, LogRatio)
  
  return(list("Long-form phospho" = data, "Metadata" = meta))
}


## Loads global data for the 210 patient experiment.
load.global.data <- function(filePath_global){
  require(dplyr)
  data <- querySynapseTable(filePath_global)
  
  meta.cols <- c("Sample", "SampleID.full", "Barcode.ID", 
                 "Plex", "Channel", "Loading.Mass", 
                 "specimen.type", "specimen.location", 
                 "Specimen.access.group.concatenated", 
                 "InitialAMLDiagnosis", "PostChemotherapy", 
                 "FLT3.ITD")
  
  meta <- unique(data[, meta.cols]) %>%
    dplyr::rename(SpecimenType = specimen.type)
  rownames(meta) <- meta$Barcode.ID
  meta <- meta[order(meta$Sample), ]
  
  data <- data %>%
    select(Gene, Sample, 'AML Sample' = Barcode.ID, LogRatio)
  
  return(list("Long-form global" = data, "Metadata" = meta))
}


## Loads functional data for the 210 patient experiment.
load.functional.data <- function(filePath_functional){
  require(dplyr)
  data <- querySynapseTable(filePath_functional)
  
  meta.cols <- c("proteomic_lab_id", "patient_id", 
                 "paired_sample", "qc_salvage")
  
  meta <- unique(data[, meta.cols]) %>%
    dplyr::rename(Barcode.ID = proteomic_lab_id,
           Paired = paired_sample)
  rownames(meta) <- meta$Barcode.ID
  
  ### AUC is the average of auc, done if multiple auc values for a given
  ### patient + inhibitor combination. percAUC tells us how large the AUC
  ### value is compared to the median AUC for that patient, measured in 
  ### percentage points.
  data <- data %>%
    select(proteomic_lab_id, inhibitor, converged, deviance, auc) %>%
    dplyr::rename('AML Sample' = proteomic_lab_id,
           Condition = inhibitor) %>%
    group_by('AML Sample', Condition) %>%
    mutate(AUC = mean(auc, na.rm = T)) %>%
    ungroup(Condition) %>%
    mutate(medAUC = median(AUC, na.rm = T)) %>%
    mutate(percAUC = AUC/medAUC*100)
    
    
  return(list("Long-form functional" = data, "Metadata" = meta))
}

