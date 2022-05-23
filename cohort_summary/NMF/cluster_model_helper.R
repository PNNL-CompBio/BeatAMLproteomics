library(glmnet)
library(groupdata2)
library(dplyr)
library(DreamAI)
library(impute)

set.seed(1)

source("../../../util/synapseUtil.R")
source("../../../util/loading_data.R")
source("../../../util/mutational_analysis_helper.R")
source("../../../util/make_plots_util.R")


## This script loads all the data needed for the modeling scripts. We also
## normalize the train and test matrices of for global, phospho and RNA data.
## This way they are ready for the follow up analysis scripts.

syn <- synapseLogin()

# load.combined.data()
load("../../../Misc/load.combined.data 3-09-2022.RData")

cluster.syn <- "syn26642544"
cluster_assignments <- read.table(syn$get(cluster.syn)$path, sep = "\t")
summary.syn <- "syn26642974"
summary.table <- read.table(syn$get(summary.syn)$path) %>%
  as.data.frame()

rownames(summary.table) <- summary.table$labId
summary.table <- summary.table[rownames(meta), ]
meta$vitalStatus <- summary.table$vitalStatus
meta <- meta %>%
  mutate(vitalStatus = case_when(vitalStatus == "Dead" ~ 2,
                                 TRUE ~ 1))
cluster_samples <- cluster_assignments$Barcode.ID

## filtering for complete features
global_features <- global.data %>%
  group_by(Gene) %>%
  summarize(total = n()) %>%
  filter(total == 210) %>%
  pull(Gene)

phospho_features <- phospho.data %>%
  group_by(SiteID) %>%
  summarize(total = n()) %>%
  filter(total == 210) %>%
  pull(SiteID)

RNA_features <- RNA.data %>%
  group_by(Gene) %>%
  summarize(total = n()) %>%
  filter(total == 159) %>%
  pull(Gene)

## Splitting data into train and test.
## Global
global_mat_train <- global.data %>%
  filter(Barcode.ID %in% cluster_samples) %>%
  dplyr::rename(feature = Gene) %>%
  filter(feature %in% global_features) %>%
  pivot_wider(names_from = Barcode.ID, 
              values_from = LogRatio) %>% 
  as.data.frame()
rownames(global_mat_train) <- global_mat_train$feature
global_mat_train <- global_mat_train[, -1] %>% as.matrix()

global_mat_test <- global.data %>%
  filter(!(Barcode.ID %in% cluster_samples)) %>%
  dplyr::rename(feature = Gene) %>%
  filter(feature %in% global_features) %>%
  pivot_wider(names_from = Barcode.ID, 
              values_from = LogRatio) %>% 
  as.data.frame()
rownames(global_mat_test) <- global_mat_test$feature
global_mat_test <- global_mat_test[, -1] %>% as.matrix()

## Phospho
phospho_mat_train <- phospho.data %>%
  filter(Barcode.ID %in% cluster_samples) %>%
  select(Barcode.ID, SiteID, LogRatio) %>%
  dplyr::rename(feature = SiteID) %>%
  filter(feature %in% phospho_features) %>%
  pivot_wider(names_from = Barcode.ID, 
              values_from = LogRatio) %>% 
  as.data.frame()
rownames(phospho_mat_train) <- phospho_mat_train$feature
phospho_mat_train <- phospho_mat_train[, -1] %>% as.matrix()

phospho_mat_test <- phospho.data %>%
  filter(!(Barcode.ID %in% cluster_samples)) %>%
  select(Barcode.ID, SiteID, LogRatio) %>%
  dplyr::rename(feature = SiteID) %>%
  filter(feature %in% phospho_features) %>%
  pivot_wider(names_from = Barcode.ID, 
              values_from = LogRatio) %>% 
  as.data.frame()
rownames(phospho_mat_test) <- phospho_mat_test$feature
phospho_mat_test <- phospho_mat_test[, -1] %>% as.matrix()

set.seed(42)
#### Imputing phospho data, lots of missing values.
phospho_imp <- phospho.data %>%
  select(Barcode.ID, SiteID, LogRatio) %>%
  dplyr::rename(feature = SiteID) %>%
  pivot_wider(names_from = Barcode.ID, 
              values_from = LogRatio) %>% 
  as.data.frame()
rownames(phospho_imp) <- phospho_imp$feature
phospho_imp <- phospho_imp[, -1] %>% as.matrix()
phospho_imp <- DreamAI(phospho_imp, k=10, maxiter_MF = 10, ntree = 100,
                       maxnodes = NULL, maxiter_ADMIN = 30, tol = 10^(-2),
                       gamma_ADMIN = 0, gamma = 50, CV = FALSE, fillmethod = "row_mean",
                       maxiter_RegImpute = 10,conv_nrmse = 1e-6, iter_SpectroFM = 40,
                       method = c("KNN"), out="Ensemble")$Ensemble
# ## Writing to table for easier use.
# write.table(phospho_imp, "phospho_imputed.txt", sep = "\t")

phospho_mat_train_imp <- phospho_imp[, colnames(global_mat_train)]
phospho_mat_test_imp <- phospho_imp[, colnames(global_mat_test)]

## RNA
RNA_mat_train <- RNA.data %>%
  filter(Barcode.ID %in% cluster_samples) %>%
  dplyr::rename(feature = Gene,
                LogRatio = `RNA counts`) %>%
  filter(feature %in% RNA_features)  %>%
  pivot_wider(names_from = Barcode.ID, 
              values_from = LogRatio) %>% 
  as.data.frame()
rownames(RNA_mat_train) <- RNA_mat_train$feature
RNA_mat_train <- RNA_mat_train[, -1] %>% as.matrix()

## Standardizing
global_mat_train <- sweep(global_mat_train, 1, apply(global_mat_train, 1, mean), FUN = '-')
global_mat_train <- sweep(global_mat_train, 1, apply(global_mat_train, 1, sd), FUN = '/')
global_mat_test <- sweep(global_mat_test, 1, apply(global_mat_test, 1, mean), FUN = '-')
global_mat_test <- sweep(global_mat_test, 1, apply(global_mat_test, 1, sd), FUN = '/')

phospho_mat_train <- sweep(phospho_mat_train, 1, apply(phospho_mat_train, 1, mean), FUN = '-')
phospho_mat_train <- sweep(phospho_mat_train, 1, apply(phospho_mat_train, 1, sd), FUN = '/')
phospho_mat_test <- sweep(phospho_mat_test, 1, apply(phospho_mat_test, 1, mean), FUN = '-')
phospho_mat_test <- sweep(phospho_mat_test, 1, apply(phospho_mat_test, 1, sd), FUN = '/')

phospho_mat_train_imp <- sweep(phospho_mat_train_imp, 1, apply(phospho_mat_train_imp, 1, mean), FUN = '-')
phospho_mat_train_imp <- sweep(phospho_mat_train_imp, 1, apply(phospho_mat_train_imp, 1, sd), FUN = '/')
phospho_mat_test_imp <- sweep(phospho_mat_test_imp, 1, apply(phospho_mat_test_imp, 1, mean), FUN = '-')
phospho_mat_test_imp <- sweep(phospho_mat_test_imp, 1, apply(phospho_mat_test_imp, 1, sd), FUN = '/')

RNA_mat_train <- sweep(RNA_mat_train, 1, apply(RNA_mat_train, 1, mean), FUN = '-')
RNA_mat_train <- sweep(RNA_mat_train, 1, apply(RNA_mat_train, 1, sd), FUN = '/')

combined_mat_train <- rbind(global_mat_train, phospho_mat_train_imp)
combined_mat_test <- rbind(global_mat_test, phospho_mat_test_imp)

## Setting up metadata for training models
enet_meta <- left_join(meta, cluster_assignments, by = "Barcode.ID") %>%
  filter(!is.na(k.5)) %>%
  mutate(overallSurvival = as.numeric(overallSurvival),
         k.5 = as.character(k.5),
         k.4 = as.character(k.4)) %>%
  select(Barcode.ID, overallSurvival, vitalStatus, k.5, k.4)
rownames(enet_meta) <- enet_meta$Barcode.ID


