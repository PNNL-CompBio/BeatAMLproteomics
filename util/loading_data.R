library(MSnbase)
library(tidyr)
library(tibble)

## Loads latest corrected phospho data for the 210 patient experiment.
load.phospho.data <- function(type = "Corrected"){
  require(dplyr)
  
  if (type == "Corrected") {
    syn <- "syn26477193"
  } else if (type == "Uncorrected") {
    syn <- "syn25808685"
  } else if (type == "Old Corrected") {
    syn <- "syn25808662"
  } else {
    stop("Type not recognized. Available types are 'Corrected', 
         'Uncorrected', and 'Old Corrected'")
  }
  
  data <- querySynapseTable(syn) %>% 
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
    select(Gene, SiteID, Barcode.ID, LogRatio)
  
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
    select(Gene, Barcode.ID, LogRatio)
  
  return(list("Long-form global" = data, "Metadata" = meta))
}


## Loads functional data for the 210 patient experiment.
load.functional.data <- function(threshold = 0.05, fam.threshold = 0.05){
  require(dplyr)
  
  syn <- synapseLogin()
  
  summary.syn <- "syn25796769"
  workbook.syn <- "syn26427390"
  summary.table <- readxl::read_xlsx(syn$get(summary.syn)$path) %>%
    select(labId, overallSurvival) %>%
    dplyr::rename(Barcode.ID = labId)
  
  data <- querySynapseTable("syn25830473") %>%
    mutate(probit_qc_error = unlist(probit_qc_error))
  
  meta.cols <- c("proteomic_lab_id", "lab_id", "patient_id",
                 "paired_sample", "qc_salvage")
  
  meta <- unique(data[, meta.cols]) %>%
    dplyr::rename(Barcode.ID = proteomic_lab_id,
           Paired = paired_sample)
  rownames(meta) <- meta$Barcode.ID
  
  ## When both labID and alt_ID exist, both the WES and RNA datasets contain 
  ## The alternate ID's instead of labId!
  sample.workbook <- readxl::read_xlsx(syn$get(workbook.syn)$path) %>%
    select(labId, patientId, `Has WES +/- 14 days`, `Has RNAseq same day sample`, Comments) %>%
    mutate(alt_ID = case_when(grepl("Alt", Comments) ~ sub("Alt -[ ]*", "", Comments),
                              T ~ labId))
  
  ## Getting drug family data
  drug.fam <- load.drugfam.data() %>%
    dplyr::rename(Inhibitor = inhibitor)
  
  ### AUC is the average of auc, done if multiple auc values for a given
  ### patient + inhibitor combination. percAUC tells us how large the AUC
  ### value is compared to the median AUC for that patient, measured in 
  ### percentage points.
  data <- data %>%
    select(proteomic_lab_id, inhibitor, converged, deviance, auc) %>%
    filter(!is.na(auc)) %>%
    dplyr::rename(Barcode.ID = proteomic_lab_id,
           Inhibitor = inhibitor) %>%
    group_by(Barcode.ID, Inhibitor) %>%
    mutate(AUC = mean(auc, na.rm = T)) %>%
    slice(1) %>%
    ungroup(Inhibitor) %>%
    mutate(medAUC = median(AUC, na.rm = T)) %>%
    mutate(percAUC = AUC/medAUC*100) %>%
    ungroup(Barcode.ID)
  
  data <- left_join(data, drug.fam, by = "Inhibitor") %>%
    left_join(summary.table, by = "Barcode.ID")
  
  ## Counting sensitivity
  numSens<-data%>%
    group_by(Inhibitor) %>%
    ## Inhibitors can have multiple families. In such cases, counting rows
    ## to get number of sensitive samples is over counting
    subset(AUC<100)%>%summarize(numSens=length(unique(Barcode.ID)))
  
  ##counting sensitivity by family
  numSensFam <-data%>%
    subset(!is.na(family))%>%
    group_by(family)%>%
    subset(AUC<100)%>%summarize(numSensFam=n())
  
  ## sensitivity frequency
  fracSens<-data%>%group_by(Inhibitor)%>%
    summarize(nSamps=length(unique(Barcode.ID)))%>%
    left_join(numSens)%>%mutate(fracSens=numSens/nSamps)
  
  ## frequency by family
  fracSensFam<-data%>%
    subset(!is.na(family))%>%
    group_by(family)%>%
    summarize(nSampsFam=n())%>%
    left_join(numSensFam)%>%mutate(fracSens=numSensFam/nSampsFam)
  
  ## subset data to inhibitors passing threshold
  withSens=subset(fracSens,fracSens>threshold)%>%
    subset(fracSens<(1-threshold))%>%
    subset(numSens>1)
  include<-withSens$Inhibitor
  data<-subset(data, Inhibitor %in% include)
  
  # drug.combos<-unique(auc.dat$Condition[grep(" - |\\+",auc.dat$Condition)])
  # print("Removing drug combinations")
  # auc.dat<<-subset(auc.dat,!Condition%in%drug.combos)
  
  ## also subseting by family sensitivity
  withSensFam=subset(fracSensFam,fracSens>fam.threshold)%>%
    subset(fracSens<(1-fam.threshold))%>%
    subset(numSensFam>1)
  
  data.fam <- subset(data, family %in% withSensFam$family) %>%
    group_by(Barcode.ID, Inhibitor) %>%
    ungroup(Barcode.ID) %>%
    ungroup(Inhibitor) %>%
    filter(!is.na(AUC))
    
  return(list("Long-form functional" = data, 
              "Long-form functional sensitive family" = data.fam,
              "Metadata" = meta,
              "Sensitivity by Inhibitor" = fracSens,
              "Sensitivity by Family" = fracSensFam))
}


load.drugfam.data <- function(){
  require(dplyr)
  ## WARNING!! Having some issue with strings which appear to be equal
  ## But have slight encoding differences, leading to possible issues when filtering.
  ## FOR EXAMPLE: Row 99 contains inhibitor "Lestaurtinib (CEP−701)", as does
  ## row 224. However these two strings ARE NOT EXACTLY EQUAL. Row 99 uses a "hyphen minus"
  ## while row 224 uses the "minus sign", which has unicode \u2212.
  
  data <- querySynapseTable("syn22156956") %>%
    mutate(inhibitor = gsub("\u2212", "-", inhibitor)) %>%
    ## Fixing typos in the drug names
    mutate(inhibitor = case_when(inhibitor == "Sorafinib" ~ "Sorafenib",
                                 inhibitor == "Gilterinitib  (ASP−2215)" ~ "Gilteritinib  (ASP−2215)",
                                 TRUE ~ inhibitor))
  
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
    select(Barcode.ID, alt_ID, stable_id, display_label, description, biotype, `RNA counts`) %>%
    dplyr::rename(Gene = display_label)
  
  return(list("Long-form RNA" = data, "Metadata" = meta))
}


load.combined.data <- function(){
  
  global.data.list <- load.global.data()
  meta <<- global.data.list$Metadata
  global.data <<- global.data.list$`Long-form global`
  
  phospho.data <<- load.phospho.data()$`Long-form phospho`
  
  functional.data.list <- load.functional.data()
  functional.data <<- functional.data.list$`Long-form functional`
  functional.data.sensitive.family <<- functional.data.list$`Long-form functional sensitive family`
  
  RNA.data <<- load.RNA.data()$`Long-form RNA` %>%
    select(Barcode.ID, Gene, `RNA counts`)
  
  WES.data <<- load.WES.data()$`Long-form WES` %>%
    select(Barcode.ID, symbol, t_vaf) %>%
    dplyr::rename(Gene = symbol) %>%
    mutate(wesMutation = t_vaf > 0)
  
  combined <<- full_join(global.data, RNA.data, by = c("Barcode.ID", "Gene")) %>%
    select(Barcode.ID, Gene, LogRatio, `RNA counts`) %>%
    full_join(WES.data, by = c("Barcode.ID", "Gene")) %>%
    dplyr::rename(globalLogRatio = LogRatio) %>%
    full_join(phospho.data, by = c("Barcode.ID", "Gene")) %>%
    dplyr::rename(phosphoLogRatio = LogRatio) %>%
    mutate(wesMutation = case_when(!is.na(wesMutation) & wesMutation ~ 1,
                                       TRUE ~ 0)) %>%
    mutate(wesMutation = as.factor(wesMutation)) %>%
    select(Barcode.ID, Gene, SiteID, globalLogRatio, 
           phosphoLogRatio, `RNA counts`, wesMutation, t_vaf)
  
  samples.g <- unique(global.data$Barcode.ID)
  samples.p <- unique(phospho.data$Barcode.ID)
  samples.r <- unique(RNA.data$Barcode.ID)
  samples.w <- unique(WES.data$Barcode.ID)

  sample.summary <<- meta %>%
    select(Barcode.ID) %>%
    mutate(Global = Barcode.ID %in% samples.g,
           Phospho = Barcode.ID %in% samples.p,
           RNA = Barcode.ID %in% samples.r,
           WES = Barcode.ID %in% samples.w)
}


make.msnset <- function(data, feature.col, sample.col = "Barcode.ID",
                        value.col = "LogRatio", metadata) {
  mat <- data %>%
    select(sym(sample.col), sym(feature.col), sym(value.col)) %>%
    pivot_wider(names_from = sym(sample.col), values_from = sym(value.col)) %>%
    column_to_rownames(feature.col) %>%
    as.matrix()
  
  m <- MSnSet(exprs = mat)
  pData(m) <- metadata[sampleNames(m), ]
  return(m)
}



