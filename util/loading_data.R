## Loads phospho data for the 210 patient experiment.
load.phospho.data <- function(){
  require(dplyr)
  
  data <- querySynapseTable("syn25808662") %>% 
    mutate(Specimen.access.group.concatenated = 
             unlist(Specimen.access.group.concatenated))
  
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
    select(Gene, SiteID, Sample, Barcode.ID, LogRatio)
  
  return(list("Long-form phospho" = data, "Metadata" = meta))
}


## Loads global data for the 210 patient experiment.
load.global.data <- function(){
  require(dplyr)
  
  data <- querySynapseTable("syn25808020")  %>%
    mutate(Specimen.access.group.concatenated = 
             unlist(Specimen.access.group.concatenated))
  
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
  
  data <- querySynapseTable("syn25830473") %>%
    mutate(probit_qc_error = unlist(probit_qc_error)) %>%
    ungroup()
  
  meta.cols <- c("proteomic_lab_id", "lab_id", "patient_id",
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
    mutate(percAUC = AUC/medAUC*100) %>%
    ungroup(Barcode.ID)
    
  return(list("Long-form functional" = data, "Metadata" = meta))
}


load.drugfam.data <- function(){
  require(dplyr)
  ## WARNING!! Having some issue with strings which appear to be equal
  ## But have slight encoding differences, leading to possible issues when filtering.
  ## FOR EXAMPLE: Row 99 contains inhibitor "Lestaurtinib (CEP−701)", as does
  ## row 224. However these two strings ARE NOT EXACTLY EQUAL. Row 99 uses a "hyphen minus"
  ## while row 224 uses the "minus sign", which has unicode \u2212.
  
  data <- querySynapseTable("syn22156956") %>%
    mutate(inhibitor = gsub("\u2212", "-", inhibitor))
  
  return(data)
}


load.WES.data <- function(){
  require(dplyr)
  
  data <- querySynapseTable("syn26428827") %>%
    dplyr::rename(Barcode.ID = labId) %>%
    mutate(refseq = unlist(refseq),
           hgvsp = unlist(hgvsp),
           Specimen.access.group.concatenated = 
             unlist(Specimen.access.group.concatenated))
  
  meta.cols <- c("Barcode.ID", "alt_ID", "Plex", "Channel", "Loading.Mass", 
                 "SpecimenType", "specimen.location", "Specimen.access.group.concatenated",
                 "InitialAMLDiagnosis", "PostChemotherapy", "FLT3.ITD")
  
  meta <- data[, meta.cols] %>%
    unique()
  rownames(meta) <- meta$Barcode.ID
  
  data <- data %>%
    select(Barcode.ID, alt_ID, gene, symbol, refseq, t_vaf, hgvsc, hgvsp)
  
  return(list("Long-form WES" = data, "Metadata" = meta))
}


load.RNA.data <- function(){
  require(dplyr)
  
  data <- querySynapseTable("syn26428813") %>%
    dplyr::rename(Barcode.ID = labId) %>%
    mutate(description = unlist(description),
           Specimen.access.group.concatenated = 
             unlist(Specimen.access.group.concatenated),
           Channel = unlist(Channel))
  
  meta.cols <- c("Barcode.ID", "alt_ID", "Plex", "Channel", "Loading.Mass", 
                 "SpecimenType", "specimen.location", "Specimen.access.group.concatenated",
                 "InitialAMLDiagnosis", "PostChemotherapy", "FLT3.ITD")
  
  meta <- data[, meta.cols] %>%
    unique()
  rownames(meta) <- meta$Barcode.ID 
  
  data <- data %>%
    select(Barcode.ID, alt_ID, stable_id, display_label, description, biotype, `RNA counts`)
  
  return(list("Long-form RNA" = data, "Metadata" = meta))
}











