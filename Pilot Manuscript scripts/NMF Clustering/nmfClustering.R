library(NMF)
library(dplyr)
library(tidyr)
library(tibble)
library(DreamAI)
library(impute)


cophenetic.correlation <- function(data.distance, nmf.H, method = "ward.D") {
  cop.distance <- dist(t(nmf.H)) %>%
    hclust(., method) %>%
    cophenetic() %>%
    as.matrix()
  idx <- lower.tri(data.distance)
  nmf.distance <- cop.distance[idx]
  data.distance <- data.distance[idx]
  
  out <- cor(nmf.distance, data.distance, 
             use = "everything", method = "pearson")
  return(out)
}


nmf.single.run <- function(n.clusters, mat, ...) {
  start <- Sys.time()
  xx <- nmf(mat, n.clusters, method = "brunet", ...)
  nmf.H <- xx@fit@H
  print(paste(n.clusters, Sys.time() - start))
  return(nmf.H)
}


nmf.clustering <- function(mat, n.clusters = c(3,4,5), N.trials = 50, 
                           clust = NULL, ...){
  set.seed(117117)
  clusters <- rep(n.clusters, N.trials)
  if (is.null(clust)){
    results <- lapply(clusters, nmf.single.run, mat)
  } else {
    clusterExport(clust, c("clusters", "mat", "nmf.single.run", "nmf"), envir = environment())
    results <- parLapply(clust, clusters, nmf.single.run, mat)
  }
  return(results)
}


prepare.nmf.mat <- function(datasets) {
  samples <- Reduce(intersect, lapply(datasets, colnames))
  datasets <- lapply(datasets, function(dataset){
    ## KNN imputation
    if (any(rowSums(is.na(dataset)) > 0)){
      dataset <- DreamAI(dataset, k=10, maxiter_MF = 10, ntree = 100,
                         maxnodes = NULL, maxiter_ADMIN = 30, tol = 10^(-2),
                         gamma_ADMIN = 0, gamma = 50, CV = FALSE, fillmethod = "row_mean",
                         maxiter_RegImpute = 10,conv_nrmse = 1e-6, iter_SpectroFM = 40,
                         method = c("KNN"), out="Ensemble")$Ensemble
    } else {
      dataset
    }
  })
  mat <- Reduce(rbind, lapply(datasets, function(dataset){dataset[, samples]}))
  means <- apply(mat, 1, mean)
  mat <- sweep(mat, 1, means, FUN = '-')
  sds <- apply(mat, 1, sd)
  mat <- sweep(mat, 1, sds, FUN = '/')
  
  mat.plus <- mat
  mat.plus[mat.plus < 0] <- 0
  mat.minus <- mat
  mat.minus[mat.minus > 0] <- 0
  mat.minus <- -mat.minus
  mat.nmf <- rbind(mat.plus, mat.minus)
  
  return(mat.nmf)
}









