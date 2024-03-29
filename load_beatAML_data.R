
load_beataml_functional_data <- function(){
  library(amlresistancenetworks)
  syn <- reticulate::import("synapseclient")
  syn$login()
  querySynapseTable("syn25830473")
  
}

load_beataml_genomic_data <- function() {
  library(amlresistancenetworks)
  syn <- reticulate::import('synapseclient')
  syn$login()
  amlresistancenetworks::querySynapseTable('syn26428827')
}

load_beat_aml_transcript_data <- function(){
  
  library(amlresistancenetworks)
  syn <- reticulate::import('synapseclient')
  syn$login()
  amlresistancenetworks::querySynapseTable('syn26428813')
  
}


load_beataml_global_data <- function() {
  library(amlresistancenetworks)
  syn <- reticulate::import("synapseclient")
  syn$login()
  amlresistancenetworks::querySynapseTable("syn25808020")
}

load_beataml_phospho_data <- function() {
  library(amlresistancenetworks)
  syn <- reticulate::import("synapseclient")
  syn$login()
  amlresistancenetworks::querySynapseTable("syn25808662")
}


load_beataml_global_msnset <- function(global_data = NULL) {
  library(dplyr)
  library(tidyr)
  library(MSnbase)
  if (is.null(global_data)) {
    global_data <- load_beataml_global_data() 
  }
  
  phenoData <- global_data %>%
    dplyr::select(-Gene, -LogRatio) %>%
    distinct()
  rownames(phenoData) <- phenoData$Barcode.ID
  
  exprsData <- global_data %>%
    dplyr::select(Gene, LogRatio, Barcode.ID) %>%
    tidyr::pivot_wider(id_cols = "Gene",
                       names_from = "Barcode.ID",
                       values_from = "LogRatio",
                       values_fill = NA_real_) %>%
    as.data.frame()
  rownames(exprsData) <- exprsData$Gene
  exprsData <- exprsData %>% select(-Gene)
  
  m <- MSnbase::MSnSet(as.matrix(exprsData))
  pData(m) <- phenoData
  return(m)
}

load_beataml_phospho_msnset <- function(phospho_data = NULL) {
  library(dplyr)
  library(tidyr)
  library(MSnbase)
  if (is.null(phospho_data)) {
    phospho_data <- load_beataml_phospho_data() 
  }
  library(dplyr)
  phenoData <- phospho_data %>%
    dplyr::select(-Gene, -SiteID, -LogRatio) %>%
    distinct()
  rownames(phenoData) <- phenoData$Barcode.ID
  
  exprsData <- phospho_data %>%
    dplyr::select(SiteID, LogRatio, Barcode.ID) %>%
    tidyr::pivot_wider(id_cols="SiteID",
                       names_from="Barcode.ID",
                       values_from="LogRatio",
                       values_fill=NA_real_) %>%
    as.data.frame()
  rownames(exprsData) <- exprsData$SiteID
  exprsData <- exprsData %>% select(-SiteID)
  
  m <- MSnSet(as.matrix(exprsData))
  pData(m) <- phenoData
  return(m)
}