library(NMF)
library(dplyr)
library(tidyr)
library(tibble)
library(DreamAI)
library(impute)
library(data.table)


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


### Runs nmf for each k value (n.clusters) provided, a total of 50 times each by default.
### The raw output is saved to an RDS file.
nmf.clustering <- function(mat, n.clusters = c(3,4,5), N.trials = 50, 
                           prefix){
  results <- list()
  for (k in n.clusters){
    start <- Sys.time()
    message(paste0("Working on k = ", k))
    message(start)
    name <- paste(prefix, "NMF k =", k, "using", N.trials, "runs.RDS")
    
    nmf.result <- nmf(mat, k, nrun = N.trials, seed = 117117, 
                      .options = list(keep.all = TRUE, parallel = TRUE))
    saveRDS(nmf.result, name)
    results <- append(results, list(nmf.result))
  }
  
  return(results)
}


### combines a list of datasets into a single matrix, ready to be factored
### using NMF. Note that missing values are imputed using KNN imputation.
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


### Given a matrix H from the NMF output, this creates the connectivity matrix, or 
### adjacency matrix, associated with H.
connectivity.matrix <- function(mat){
  samples <- colnames(mat)
  out <- matrix(0, nrow = length(samples), ncol = length(samples))
  cluster <- apply(mat, 2, which.max)
  N <- length(unique(cluster))
  members <- lapply(1:N, function(i){
    which(cluster == i) %>%
      names()
  })
  out <- lapply(samples, function(sample){
    xx = members[[cluster[[sample]]]]
    xx = samples %in% xx %>%
      as.numeric()
  }) %>% do.call(rbind, .)
  
  rownames(out) <- samples
  colnames(out) <- samples
  return(out)
}


### Given the output for a single k value from nmf.clustering, this outputs
### the cluster each sample belongs to. Cluster naming is arbitrary.
get.clusters <- function(results){
  nmf.C <- connectivity(results)
  samples <- colnames(nmf.C)
  
  clusters <- lapply(samples, function(x){
    values = nmf.C[, x]
    names(which(values == 1))
  }) %>%
    unique()
  
  out <- lapply(1:length(clusters), function(i){
    members = data.frame(Cluster = i, Sample = clusters[[i]])
  }) %>%
    rbindlist() %>%
    as.data.frame()
  rownames(out) <- out$Sample
  
  return(out)
}

# results.3.C <- lapply(results.3, function(x){
#   connectivity.matrix(x@fit@H)
# })
# 
# C <- matrix(0, ncol = 159, nrow = 159)
# rownames(C) <- rownames(xx)
# colnames(C) <- colnames(xx)
# lapply(results.3.C, function(x){
#   C <<- C + x
# })
# C <- C/50
# 
# results.3.C.og <- lapply(results.3, connectivity)
# diff <- lapply(1:50, function(i){results.3.C.og[[i]] - results.3.C[[i]]})
# sums <- lapply(1:50, function(i){
#   total.diff <- apply(diff[[i]], 2, sum) %>%
#     sum()
# })







