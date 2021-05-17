#'
#' Tools for post-processing global and phospho-proteomics data
#' Author: Michael Nestor `(michael.nestor@pnnl.gov)`
#' 
#' @importFrom dplyr select inner_join group_by_at summarize filter pull
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather
#' @importFrom MSnbase exprs pData
#' 
#' 

library(dplyr)
library(tibble)
library(tidyr)

#' @param msnset (MSnSet object) 
normalizeBySampleMedians <- function(msnset) {
  global_coeffs <- apply(exprs(msnset), 2, median, na.rm=T)
  exprs(msnset) <- sweep(exprs(msnset), 2, global_coeffs, "-")
  return(msnset)
}

#' @param msnset (MSnSet object) 
#' @param msnset_global (MSnSet object) 
normalizeByGlobalSampleMedians <- function(msnset, msnset_global) {
  global_coeffs <- apply(exprs(msnset_global), 2, median, na.rm=T)
  exprs(msnset) <- sweep(exprs(msnset), 2, global_coeffs, "-")
  return(msnset)
}


#' @param msnset (MSnSet object) 
normalizeByMedpolish <- function(msnset, ...) {
  out <- medpolish(exprs(msnset), eps = .Machine$double.eps,
                 maxiter = 100, na.rm = T, trace.iter = F, ...)
  exprs(msnset) <- out$residuals
  return(msnset)
}


#' @param msnset (MSnSet object) 
#' @param least_proportion_threshold (numeric)
filterByProportionMissingValues <- function(msnset, least_proportion_threshold = 0.5) {
  sufficiently_present_features <- which(apply(!is.na(exprs(msnset)), 1, mean) >= least_proportion_threshold)
  msnset <- msnset[sufficiently_present_features,]
}

#' @param msnset (MSnSet object) 
#' @param batch_name (character) 
#' @param least_count_threshold (integer)
filterByMissingPerBatch <- function(msnset, batch_name, least_count_threshold=1L) {
  batch_to_sample <- pData(msnset) %>%
    select(batch_name) %>%
    rownames_to_column("sample_name")
  
  sufficiently_present_features <- exprs(msnset) %>%
    as.data.frame() %>% 
    rownames_to_column("feature_name") %>%
    tidyr::gather(sample_name, abundance, -feature_name) %>%
    inner_join(batch_to_sample, by = "sample_name") %>%
    group_by_at(c(batch_name, "feature_name")) %>% 
    summarize(cnt = sum(!is.na(abundance))) %>%
    group_by_at("feature_name") %>% 
    summarize(min_cnt = min(cnt)) %>% 
    filter(min_cnt >= least_count_threshold) %>%
    pull(feature_name)
  
  msnset <- msnset[sufficiently_present_features, ]
}



#' @param msnset (MSnSet object) 
#' @param removed_cov_name (character) 
#' @param retained_cov_name (character)
correct_batch_effect_empiricalBayesLM <- function (msnset, removed_cov_name, retained_cov_name = NULL, ...) 
{
  x <- exprs(msnset)
  x <- as.data.frame(t(x))
  
  if (!all(c(removed_cov_name, retained_cov_name) %in% names(pData(msnset)))) {
    stop("The covariates are not recognized")
  }
  
  if (is.null(retained_cov_name)) {
    soln <- WGCNA::empiricalBayesLM(data = x,
                                    removedCovariates = pData(msnset)[, removed_cov_name])
  }
  else {
    soln <- WGCNA::empiricalBayesLM(data = x,
                                    removedCovariates = pData(msnset)[, removed_cov_name],
                                    retainedCovariates = pData(msnset)[, retained_cov_name])
  }
  
  x <- soln[["adjustedData"]]
  exprs(msnset) <- t(as.matrix(x))
  return(msnset)
}