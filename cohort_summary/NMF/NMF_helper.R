library(NMF)
library(dplyr)
library(tidyr)
library(tibble)
library(DreamAI)
library(impute)
library(data.table)


### Given a matrix H from the NMF output, this creates the connectivity matrix, or 
### adjacency matrix, associated with H. The NMF package uses the function 'connectivity', and they 
### are equivalent. This function is really a sanity check to make sure I understand the clustering process.
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


### Computes the cophenetic correlation. Use the NMF function cophcor for this instead. 
### This is here as a sanity check, to check my understanding of each step.
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
### Note the NMF functionality for this exists, the issue the implementation
### leads to RAM issues when running multiples k values. So we split just
### use a simple for loop to run through the k values. See the following link for some discussion.
### https://github.com/renozao/NMF/issues/14
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
### This is the same processing performed prior to NMF in this paper
### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7373300/
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


### Given the output for a single k value from nmf.clustering, this outputs
### the cluster each sample belongs to. Cluster naming is arbitrary.
get.clusters.individual <- function(results){
  ### 'what = chc' derives the clusters from the consensus matrix, averaging the runs
  ### If not using 'chc', then it outputs the connectivity of a SINGLE nmf run, namely
  ### the one with smallest residual score.
  if (length(results) == 1) {
    nmf.C <- connectivity(results)
  } else {
    cat("Using chc (all results) to infer clusters\n")
    nmf.C <- connectivity(results, what = "chc")
  }
  samples <- colnames(nmf.C)
  
  clusters <- lapply(samples, function(x){
    values = nmf.C[, x]
    names(which(values == 1))
  }) %>%
    unique()
  
  out <- lapply(1:length(clusters), function(i){
    members = data.frame(Cluster = i, Barcode.ID = clusters[[i]])
  }) %>%
    rbindlist() %>%
    as.data.frame()
  rownames(out) <- out$Barcode.ID
  
  return(out)
}


get.clusters <- function(results.list, type = "Cluster Membership"){
  if (type == "Cluster Membership") {
    out.list <- lapply(results.list, function(result){
      k = dim(result[[1]]@fit@H)[1]
      xx <- get.clusters.individual(result) %>%
        select(Cluster) %>%
        plyr::rename(replace = c("Cluster" = paste0("k=", k))) %>%
        as.data.frame() %>%
        mutate(Barcode.ID = rownames(.))
      return(xx)
    })
  }
  
  out <- purrr::reduce(out.list, left_join, by = "Barcode.ID")
  return(out)
}

### sample.categories should have a column with sample names (Barcode.ID)
### computes enrichment using fisher test for each cluster in each category from
### sample.categories. Outputs a mtrix containing the significance of enrichment,
### which can be used as heatmap column/row labels
get.enrichment <- function(results, sample.categories, log.scale = FALSE) {
  cluster.membership <- get.clusters.individual(results)
  k = length(unique(cluster.membership$Cluster))
  
  sample.categories <- sample.categories %>%
    filter(rownames(.) %in% rownames(cluster.membership)) %>%
    left_join(cluster.membership, by = "Barcode.ID") %>%
    column_to_rownames("Barcode.ID") %>% 
    as.data.frame()
  
  categories <- sample.categories %>%
    select(-Cluster) %>%
    colnames()
  
  out <- lapply(categories, function(category){
    
    p.values <- sapply(1:k, function(i){
      cluster.df <- sample.categories %>%
        mutate(Cluster = case_when(Cluster == i ~ TRUE,
                                   TRUE ~ FALSE)) %>%
        mutate(Cluster = as.factor(Cluster))
      dat <- table(cluster.df[[category]], cluster.df$Cluster)
      fisher.test(dat, workspace = 2e8, alternative = "greater")$p.value
    })
    
    if (log.scale) {
      p.values <- -log10(p.values)
    } else {
      p.values <- p.values
    }
    
  }) %>% do.call(cbind, .) %>% as.data.frame()
  
  colnames(out) <- paste0("Enrichment in ", categories)
  out <- out %>%
    mutate(Cluster = 1:k) %>%
    select(Cluster, everything())
  
}


### result is the NMF output from a single k. Outputs the purity (as defined by ...) of the clustering
### for each category given.
clustering.purity <- function(result, sample.categories, normalized = FALSE) {
  clusters <- get.clusters.individual(result) %>%
    mutate(Cluster = as.factor(Cluster))
  sample.categories <- sample.categories[clusters$Barcode.ID, ] %>%
    mutate_all( ~ as.factor(.))
  
  out <- sapply(colnames(sample.categories), function(id){
    xx <- clusters$Cluster
    yy <- sample.categories[[id]]
    lower.limit <- max(table(yy))/length(yy)
    if (normalized) {
      (purity(xx, yy) - lower.limit)/(1-lower.limit)
    } else {
      purity(xx,yy)
    }
  })
}



### result is the output from NMF for a SINGLE k value. result consisting of multiple runs is allowed.
### A single element from 'result' is being used to extract the features for each cluster.
### This element is chosen as the member most closely approximating the cluster structure
### induced by the consensus matrix using ALL runs within the result.
extract.features <- function(result){
  ### Using WHOLE result to infer clusters (adjacency matrix)
  aggregate.connectivity <- connectivity(result, what = 'chc')
  
  ### Choosing an individual run which most closely resembles aggregate clusters.
  index <- lapply(1:length(result), function(i){
    sum(abs(connectivity(result[[i]]) - aggregate.connectivity))}
    ) %>%
    which.min()
  
  chosen.result <- result[[index]]
  W <- chosen.result@fit@W
  
  chosen.sample.assignments <- apply(chosen.result@fit@H, 2, which.max) 
  combined.assignment <- data.frame(chosenCluster = chosen.sample.assignments, 
                                    Barcode.ID = names(chosen.sample.assignments))
  
  
  combined.assignment <- get.clusters.individual(result) %>%
    full_join(combined.assignment, by = "Barcode.ID")
  
  cluster.mapping <- table(combined.assignment$Cluster, combined.assignment$chosenCluster) %>%
    apply(1, which.max)
  W_out <- W[, cluster.mapping]
  colnames(W_out) <- paste("Cluster", 1:ncol(W_out))
  
  if (length(cluster.mapping) != length(unique(cluster.mapping))){
    stop("Unable to unambiguously infer cluster labels. Please inspect further.")
  }
  
  extracted.features.index <- extractFeatures(chosen.result)[cluster.mapping]
  
  extracted.features <- extracted.features.index %>%
    lapply(function(indeces){
      rownames(W)[indeces]
    })
  
  feature_score <- featureScore(chosen.result)
  feature_score <- data.frame(NMF_index = 1:length(feature_score), feature_score = feature_score)
  
  out.df <- lapply(1:length(extracted.features), function(i){
    features <- extracted.features[[i]]
    indeces <- extracted.features.index[[i]]
    if (all(is.na(features))){
      data.frame(Cluster = i, Feature_NMF_label = NA, NMF_index = NA, feature_score = NA)
    } else {
      yy <- data.frame(Cluster = i, Feature_NMF_label = features, NMF_index = indeces) %>%
        cbind(W_out[indeces, ]) %>%
        left_join(feature_score, by = "NMF_index") %>%
        select(Feature_NMF_label, Cluster, NMF_index, feature_score, everything())
    }
  }) %>% do.call(rbind, .)
  
  return(out.df)
}




















