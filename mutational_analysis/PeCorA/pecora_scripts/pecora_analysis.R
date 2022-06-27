library(PeCorA)
library(grid)

## PeCorA fits linear models to assess whether a peptide’s change across treatment 
## groups differs from all other peptides assigned to the same protein. 


#' Performs the PeCorA analysis described in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8592057/.
#' Take a look at the linked paper to get an idea of the procedure.
#'
#' Takes an msnset, a list of proteins, and a treatment variable. msnset is assumed to be normalized, for instance it can be the output of
#' pecora_preprocessing, which normalizes as in the original PeCorA paper.
#' The PeCorA function is then used for the analysis. Returns a table containing Peptide, Protein, and significance.
#'
#' @param m msnset object containing a pData AND fData table. fData must have a Protein and Peptide column.
#' @param treatment_string the name of the column in pData containing treatment information. The data can be numerical of categorical.
#' @param proteins character vector indicating which proteins to analyze.
#' @param median_mod logical for whether to use the original PeCorA (median_mod = FALSE), or the median version (median_mod = TRUE). See the help of PeCorA_mod for more details.
#'
#' @return PeCorA results table for the supplied proteins. Saves plots of the significant peptides to the given folder.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom grid grobTree textGrob gpar
#'
#' @export pecora_analysis
#'
#'
pecora_analysis <- function(m, treatment_string, proteins, 
                            median_mod = FALSE, m_protein = NULL){
  
  mat <- exprs(m)
  peptide_mapping <- fData(m)
  protein_mod <- FALSE
  
  if (!("Protein" %in% colnames(peptide_mapping))){
    stop("fData should have 'Protein' and 'Peptide' columns.\n")
  }
  
  if (!is.null(m_protein)){
    proteins <- intersect(proteins, rownames(exprs(m_protein)))
  }
  
  peptide_mapping <- peptide_mapping %>%
    mutate(Protein = as.character(Protein)) %>%
    select(Protein) %>%
    filter(Protein %in% proteins) %>%
    mutate(Peptide = rownames(.))
  
  feature_names <- rownames(peptide_mapping)
  
  metadata <- pData(m) %>%
    select(sym(treatment_string))
  metadata$Sample <- rownames(metadata)
  
  mat <- mat[feature_names, ]
  
  if (typeof(metadata[[1]]) == "character") {
    metadata[[1]] <- as.factor(metadata[[1]])
  }
  
  PeCorA_input <- mat %>%
    as.data.frame() %>%
    mutate(Peptide = rownames(.)) %>%
    pivot_longer(cols = -Peptide, names_to = "Sample", values_to = "LogRatio") %>%
    mutate(modpep_z = Peptide, ms1adj = LogRatio) %>%
    merge(peptide_mapping, by = "Peptide") %>%
    merge(metadata, by = "Sample") %>%
    dplyr::rename(Condition = sym(treatment_string))
  
  if (!is.null(m_protein)){
    protein_input <- exprs(m_protein)[proteins, ] %>%
      as.data.frame() %>%
      mutate(Protein = rownames(.)) %>%
      pivot_longer(cols = -Protein, names_to = "Sample", values_to = "LogRatio") %>%
      mutate(Peptide = paste(Protein, "@DUMMY_PEPTIDE"), 
             modpep_z = paste(Protein, "@DUMMY_PEPTIDE"), 
             ms1adj = LogRatio) %>%
      merge(metadata, by = "Sample") %>%
      dplyr::rename(Condition = sym(treatment_string)) %>%
      select(Sample, Peptide, LogRatio, modpep_z, ms1adj, Protein, Condition)
    
    PeCorA_input <- PeCorA_input %>%
      rbind(protein_input)
    protein_mod <- TRUE
  }
  
  pval_test_label <- paste0("pvalue_", treatment_string)
  
  PeCorA_result <- PeCorA_mod(PeCorA_input, median_mod, protein_mod) %>%
    dplyr::rename(Peptide = peptide,
                  Protein = protein) %>%
    group_by(Protein) %>%
    ## Adjust raw pvalues using only pvalues from the same protein (group by protein)
    mutate(adj_pval_protein_wise = p.adjust(pvalue, method = "BH"))
  
  pval_index <- which(colnames(PeCorA_result) == "pvalue")[[1]]
  colnames(PeCorA_result)[[pval_index]] <- pval_test_label
  
  return(PeCorA_result)
}











#' Plots a particular peptide using the PeCorA input and output tables.
#'
#' @param m msnset object containing a pData AND fData table. fData must have a Protein and Peptide column.
#' @param pecora_results PeCorA results table from pecora_analysis.
#' @param chosen_protein The protein of interest.
#' @param chosen_peptide The peptide of interest. Should map to the supplied protein.
#'
#' @return Plot of the specified peptide + protein.
#'
#' @import ggplot2
#' @importFrom grid grobTree textGrob gpar
#'
#' @export pecora_plot
#'
#'
pecora_plot <- function(m, pecora_results,
                        chosen_protein, chosen_peptide,
                        median_mod = FALSE, m_protein = NULL) {
  
  index <- which(grepl("pvalue_", colnames(pecora_results)))
  
  if (length(index) > 1) {
    stop("Found multiple treatment variables in the given results.\n")
  }
  
  p.value <- pecora_results %>%
    filter(Peptide == chosen_peptide) %>%
    pull(adj_pval) %>% 
    format(digits = 4)
  
  # chosen_feature <- pecora_results %>%
  #   filter(Peptide == chosen_peptide) %>%
  #   pull(feature)
  
  if (length(p.value) == 0){
    message <- paste0(chosen_protein, " + ", chosen_peptide, " not found within the results.\n")
    stop(message)
  }
  
  treatment_string <- colnames(pecora_results)[[index]] %>%
    sub("pvalue_", "", .)
  
  features <- fData(m) %>%
    filter(Protein == chosen_protein) %>%
    mutate(feature = rownames(.))
  
  metadata <- pData(m) %>%
    select(sym(treatment_string)) %>%
    mutate(Sample = rownames(.))
  
  if (typeof(metadata[[1]]) == "character") {
    metadata[[1]] <- as.factor(metadata[[1]])
  }
  
  plot_df <- exprs(m)[features$feature, ] %>%
    as.data.frame() %>%
    mutate(feature = rownames(.)) %>%
    pivot_longer(cols = -feature, names_to = "Sample", values_to = "ms1adj") %>%
    mutate(modpep_z = feature) %>%
    merge(features, by = "feature") %>%
    merge(metadata, by = "Sample") %>%
    dplyr::rename(Condition = sym(treatment_string)) %>%
    mutate(peptide_group = case_when(feature == chosen_peptide ~ chosen_peptide,
                                     TRUE ~ "All other peptides")) %>%
    mutate(peptide_group = factor(peptide_group,
                                  levels = c("All other peptides", chosen_peptide)))
  
  if (median_mod) {
    plot_df <- plot_df %>%
      group_by(peptide_group, Sample) %>%
      mutate(ms1adj = median(ms1adj, na.rm = T)) %>%
      select(Sample, ms1adj, Condition, peptide_group) %>%
      unique()
  } else if (!is.null(m_protein)) {
    plot_prot_df <- data.frame(ms1adj = exprs(m_protein)[chosen_protein, ], 
                               Sample = colnames(exprs(m))) %>%
      mutate(feature = paste(chosen_protein, "@DUMMY_PEPTIDE")) %>%
      mutate(modpep_z = feature) %>%
      merge(metadata, by = "Sample") %>%
      dplyr::rename(Condition = sym(treatment_string)) %>%
      mutate(peptide_group =  "Protein")
    plot_df <- plot_df %>%
      select(Sample, feature, ms1adj, modpep_z, Condition, peptide_group) %>% 
      rbind(plot_prot_df) %>%
      filter(peptide_group != "All other peptides") %>%
      mutate(peptide_group = factor(peptide_group,
                                    levels = c("Protein", chosen_peptide))) %>%
      unique()
  }
  
  ## Boxplot if condition is a factor. Otherwise, scatterplots when using numerical data.
  if (is.factor(plot_df$Condition)) {
    
    grob <- grobTree(textGrob(paste0("Adj_pval = ", p.value),
                              x = 0.1,  y = 0.97, hjust = 0,
                              gp = gpar(col = "black", fontsize = 11)))
    
    p <- ggplot(plot_df, aes(x = Condition, y = ms1adj, fill = peptide_group)) +
      geom_boxplot(notch = TRUE, outlier.shape = NA) +
      ggtitle(chosen_peptide) + annotation_custom(grob) +
      xlab(treatment_string) + ylab("Log Intensity") +
      geom_point(position = position_jitterdodge(), alpha = 0.1) +
      theme(plot.title = element_text(hjust = 0.5))
    
    
  } else {
    
    plot_df <- plot_df %>%
      mutate(alpha = case_when(peptide_group == "All other peptides" ~ 0.005,
                               TRUE ~ 0.25))
    
    lm_df <- plot_df %>%
      filter(peptide_group == "All other peptides")
    allothers_lm <- lm(lm_df$ms1adj ~ lm_df$Condition)
    
    lm_df <- plot_df %>%
      filter(peptide_group == chosen_peptide)
    chosen_lm <- lm(lm_df$ms1adj ~ lm_df$Condition)
    
    grob <- grobTree(textGrob(paste0("Adj_pval = ", p.value),
                              x = 0.1,  y = 0.97, hjust = 0,
                              gp = gpar(col = "black", fontsize = 11)))
    
    p <- ggplot(plot_df, aes(x = Condition, y = ms1adj, color = peptide_group, alpha = alpha)) + geom_point() +
      ggtitle(chosen_peptide) + annotation_custom(grob) + xlab(treatment_string) +
      ylab("Log Intensity") + guides(alpha = "none") +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_abline(intercept = allothers_lm$coefficients[[1]],
                  slope = allothers_lm$coefficients[[2]],
                  color = "red", size = 1) +
      geom_abline(intercept = chosen_lm$coefficients[[1]],
                  slope = chosen_lm$coefficients[[2]],
                  color = scales::hue_pal()(2)[[2]], size = 1)
  }
  
  return(p)
}




#' Performs the PeCorA analysis described in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8592057/.
#' Take a look at the linked paper to get an idea of the procedure.
#'
#' Pass an msnset here to preprocess the data. The result is again an msnset. Can perform two normalization
#' steps, also carried out in the PeCorA paper. Use the parameters peptide_standardize and sample_standardize
#' to select which steps should be performed.
#'
#' If it's desired that both normalization steps be skipped (both set to FALSE), then pecora_preprocess does nothing
#' at all to the data, and this step can be safely skipped.
#'
#' @param m msnset object containing a pData AND fData table. fData must have a Protein and Peptide column.
#' @param treatment_string the name of the column in pData containing treatment information. The data can be numerical of categorical.
#' @param sample_standardize logical for whether exprs data should be standardized so that sample wise we have mean zero and variance one.
#' @param peptide_standardize logical for whether exprs data should be normalized so that logratio data is relative to the control group. When using numerical data, this determines whether to zero center each peptide.
#'
#' @return PeCorA results table for the supplied proteins. Saves plots of the significant peptides to the given folder.
#'
#' @export pecora_preprocess
#'
#'
pecora_preprocess <- function(m, treatment_string, control_group = NULL,
                              sample_standardize = TRUE, peptide_standardize = TRUE){
  
  mat <- exprs(m)
  
  if (sample_standardize) {
    cat("Standardizing by sample.\n")
    sample_means <- apply(mat, 2, mean, na.rm = T)
    sample_sds <- apply(mat, 2, sd, na.rm = T)
    mat <- sweep(mat, 2, sample_means, FUN = '-')
    mat <- sweep(mat, 2, sample_sds, FUN = '/')
  }
  
  metadata <- pData(m) %>%
    select(sym(treatment_string))
  metadata$Sample <- rownames(metadata)
  
  if (typeof(metadata[[1]]) == "character") {
    metadata[[1]] <- as.factor(metadata[[1]])
  }
  
  if (peptide_standardize & is.factor(metadata[[1]])) {
    cat("Standardizing to control group.\n")
    if (is.null(control_group)) {
      stop("Must define 'control_group' in order to normalize relative to control\n")
    }
    
    samples_control <- metadata %>%
      filter(!!sym(treatment_string) == control_group) %>%
      rownames()
    
    if (length(samples_control) == 0){
      stop("No samples in control group. Please check value of 'control_group' and pData\n")
    }
    
    ## In case samples_control consists of a single sample we need as.matrix() in the apply call.
    peptide_means <- apply(as.matrix(mat[, samples_control]), 1, mean, na.rm = T)
    mat <- sweep(mat, 1, peptide_means, FUN = '-')
    
  } else if (peptide_standardize) {
    
    cat("Mean centering each peptide.\n")
    peptide_means <- apply(mat, 1, mean, na.rm = T)
    mat <- sweep(mat, 1, peptide_means, FUN = '-')
    
  }
  
  m.processed <- MSnSet(exprs = mat, fData = fData(m), pData = pData(m))
  pecora_design <- paste0("pecora_design_", treatment_string)
  pData(m.processed)[[pecora_design]] <- metadata[[1]]
  
  return(m.processed)
}




#' This is a modified version of the PeCorA function.
#' For the original paper see here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8592057/.
#'
#' The modification involves aggregating the intensities of "all other peptides" within
#' the test described in the PeCorA paper. We aggregate "all other peptides" by sample using the median, this
#' way each peptide mapping to a particular protein is not treated as an 'independent' observation.
#' The effect is that p-values are much higher (less significant) than in the original version of PeCorA.
#'
#' This median modification is turned OFF by default, and only affects the p-values reported by PeCorA.
#'
#'
#' @param t PeCorA input table. See the original PeCorA function for details.
#' @param median_mod logical for whether to use the median version of PeCorA (median_mod = TRUE), or the original (median_mod = FALSE)
#'
#' @return PeCorA results table.
#'
#' @importFrom stats lm
#'
#' @export PeCorA_mod
#'
PeCorA_mod <- function (t, median_mod = FALSE, protein_mod = FALSE) {
  
  if (median_mod & protein_mod){
    stop("The protein pecora modification is incompatible with the median modification.")
  }
  print("checking which proteins still have at least 2 peptides")
  pgs <- levels(as.factor(t$Protein))
  pgs_morethan2 <- c()
  # Extremely slow
  # for (x in pgs) {
  #   if (length(unique(t[t$Protein %in% x, "modpep_z"])) >
  #       1) {
  #     pgs_morethan2 <- c(pgs_morethan2, x)
  #   }
  # }
  pgs_morethan2 <- t %>%
    select(Protein, modpep_z) %>%
    unique() %>%
    group_by(Protein) %>%
    summarize(n_peptides = n()) %>%
    filter(n_peptides > 1) %>%
    pull(Protein)
  
  allp <- list()
  j = 1
  t0 <- Sys.time()
  print("computing the interaction p-values")
  pb <- txtProgressBar(min = 0, max = length(pgs_morethan2),
                       style = 3)
  for (x in pgs_morethan2) {
    tmpdf <- t[which(t$Protein == x), ]
    tmpdf["allothers"] <- rep("allothers", times = nrow(tmpdf))
    pvalues <- c(rep(0, length(unique(tmpdf$modpep_z))))
    i = 1
    for (y in setdiff(unique(tmpdf$modpep_z), paste(x, "@DUMMY_PEPTIDE"))) {
      subtmpdf <- tmpdf
      subtmpdf[which(tmpdf$modpep_z == y), "allothers"] <- y
      if (median_mod){
        subtmpdf <- subtmpdf %>%
          group_by(Sample, allothers) %>%
          mutate(ms1adj = median(ms1adj, na.rm = T)) %>%
          ungroup() %>%
          select(Sample, ms1adj, Condition, allothers) %>%
          unique()
      } else if (protein_mod){
        subtmpdf <- subtmpdf %>%
          filter(allothers != "allothers" | grepl("@DUMMY_PEPTIDE", Peptide)) %>%
          select(Sample, ms1adj, Condition, allothers) %>%
          unique()
      }
      tmplm <- lm(subtmpdf$ms1adj ~ subtmpdf$Condition *
                    subtmpdf$allothers)
      tmpanova <- car::Anova(tmplm)
      pvalues[i] <- tmpanova$`Pr(>F)`[3]
      i = i + 1
    }
    allp[[x]] <- pvalues
    setTxtProgressBar(pb, j)
    j = j + 1
  }
  print(" ")
  print(paste("PeCorA finished in ", round(Sys.time() - t0,
                                           2), " minutes", sep = ""))
  print(paste("number of proteins tested =", length(allp),
              sep = " "))
  print(paste("number of peptides tested =", length(unlist(allp)),
              sep = " "))
  print("started making data table")
  alldf = data.frame()
  x <- names(allp)[1]
  for (x in names(allp)) {
    tmpdf <- t[which(t$Protein == x), ]
    tmp_peps <- unique(tmpdf$modpep_z)
    if (length(tmp_peps) > 0) {
      tmp_pval <- allp[[x]]
      tmpout = cbind.data.frame(protein = rep(x, length(allp[[x]])),
                                tmp_peps, tmp_pval = as.numeric(tmp_pval))
      alldf = rbind(alldf, tmpout)
    }
  }
  print("correcting p-values")
  alldf$adj_pval <- p.adjust(alldf$tmp_pval, method = "BH")
  alldf_ordered <- alldf[order(alldf$adj_pval), ]
  print(paste("number of uncorrelated peptides =", nrow(alldf[alldf$adj_pval <=
                                                                0.01, ]), sep = " "))
  print(paste("number of proteins with uncorrelated peptides =",
              length(unique(alldf[alldf$adj_pval <= 0.01, ]$protein)),
              sep = " "))
  colnames(alldf_ordered)[2] <- "peptide"
  colnames(alldf_ordered)[3] <- "pvalue"
  alldf_ordered
}












