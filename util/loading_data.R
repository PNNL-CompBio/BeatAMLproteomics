## Loads phospho data for the 210 patient experiment.
load.phospho.data <- function(){
  require(dplyr)
  data <- querySynapseTable("syn25808662")
  
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
    select(Gene, SiteID, Sample, Barcode.ID, LogRatio)
  
  return(list("Long-form phospho" = data, "Metadata" = meta))
}


## Loads global data for the 210 patient experiment.
load.global.data <- function(){
  require(dplyr)
  data <- querySynapseTable("syn25808020")
  
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
    select(Gene, Sample, Barcode.ID, LogRatio)
  
  return(list("Long-form global" = data, "Metadata" = meta))
}


## Loads functional data for the 210 patient experiment.
load.functional.data <- function(){
  require(dplyr)
  data <- querySynapseTable("syn25830473")
  
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
    dplyr::rename(Barcode.ID = proteomic_lab_id,
           Inhibitor = inhibitor) %>%
    group_by(Barcode.ID, Inhibitor) %>%
    mutate(AUC = mean(auc, na.rm = T)) %>%
    ungroup(Inhibitor) %>%
    mutate(medAUC = median(AUC, na.rm = T)) %>%
    mutate(percAUC = AUC/medAUC*100)
    
    
  return(list("Long-form functional" = data, "Metadata" = meta))
}

