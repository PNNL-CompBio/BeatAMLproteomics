# Title     : Python Data Class example
# Objective : Demonstrate how to use the python class in R
# Created by: pinojc
# Created on: 12/14/2021
library(reticulate)
library(DreamAI)
library(NMF)

reticulate::use_condaenv("beatAML_env_37", required = TRUE)
reticulate::source_python("load_data.py")
# load data class
data <- AMLData()

# can access functions with $
mat <- data$all_data

# All Camilos code below, can be refactored
mat <- DreamAI(mat, k = 10, maxiter_MF = 10, ntree = 100,
               maxnodes = NULL, maxiter_ADMIN = 30, tol = 10^(-2),
               gamma_ADMIN = 0, gamma = 50, CV = FALSE, fillmethod = "row_mean",
               maxiter_RegImpute = 10, conv_nrmse = 1e-6, iter_SpectroFM = 40,
               method = c("KNN"), out = "Ensemble")$Ensemble


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
print("Made it to NMF")
nmf_result <- nmf(t(mat.nmf), c(2, 4, 8, 12), "lee", nrun = 5, seed = 117117,
                  .options = list(keep.all = TRUE, parallel = 12, verbose = TRUE))
print("Done with NMF")
summary(nmf_result)
png('test.png')
NMF::plot(nmf_result)

consensusmap(nmf_result,
             tracks = "consensus:",
             main = "Consensus matrix", info = FALSE,)
dev.off()