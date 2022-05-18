library(glmnet)
library(groupdata2)
library(dplyr)


source("../util/synapseUtil.R")
source("../util/loading_data.R")
source("../util/mutational_analysis_helper.R")
source("../util/make_plots_util.R")


## This script loads all the data needed for the modeling scripts. We also
## normalize the train and test matrices of for global, phospho and RNA data.
## This way they are ready for the follow up analysis scripts.

syn <- synapseLogin()

# load.combined.data()
load("../Misc/load.combined.data 3-09-2022.RData")

all_samples <- global.data %>%
  pull(Barcode.ID) %>% unique()
dummy_mutation <- data.frame(Barcode.ID = all_samples, 
                             Gene = "dummy")

mutation_df <- load_mutational_sample_data()
mutation_mat <- mutation_df %>%
  rbind(dummy_mutation) %>%
  unique() %>%
  group_by(Gene) %>%
  mutate(total = n()) %>%
  ungroup(Gene) %>%
  filter(total > 3,
         Gene != "NPM1") %>%
  mutate(dummy = "TRUE") %>%
  pivot_wider(-total, names_from = "Gene", 
              values_from = "dummy", values_fill = "FALSE") %>%
  select(-dummy) %>%
  as.data.frame()
  

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
global_mat <- global.data %>%
  dplyr::rename(feature = Gene) %>%
  filter(feature %in% global_features) %>%
  pivot_wider(names_from = Barcode.ID, 
              values_from = LogRatio) %>% 
  as.data.frame()
rownames(global_mat) <- global_mat$feature
global_mat <- global_mat[, -1] %>% as.matrix()

## Phospho
phospho_mat <- phospho.data %>%
  select(Barcode.ID, SiteID, LogRatio) %>%
  dplyr::rename(feature = SiteID) %>%
  filter(feature %in% phospho_features) %>%
  pivot_wider(names_from = Barcode.ID, 
              values_from = LogRatio) %>% 
  as.data.frame()
rownames(phospho_mat) <- phospho_mat$feature
phospho_mat <- phospho_mat[, -1] %>% as.matrix()

## RNA
RNA_mat <- RNA.data %>%
  dplyr::rename(feature = Gene,
                LogRatio = `RNA counts`) %>%
  filter(feature %in% RNA_features)  %>%
  pivot_wider(names_from = Barcode.ID, 
              values_from = LogRatio) %>% 
  as.data.frame()
rownames(RNA_mat) <- RNA_mat$feature
RNA_mat <- RNA_mat[, -1] %>% as.matrix()

## Standardizing
global_mat <- sweep(global_mat, 1, apply(global_mat, 1, mean), FUN = '-')
global_mat <- sweep(global_mat, 1, apply(global_mat, 1, sd), FUN = '/')


phospho_mat <- sweep(phospho_mat, 1, apply(phospho_mat, 1, mean), FUN = '-')
phospho_mat <- sweep(phospho_mat, 1, apply(phospho_mat, 1, sd), FUN = '/')


RNA_mat <- sweep(RNA_mat, 1, apply(RNA_mat, 1, mean), FUN = '-')
RNA_mat <- sweep(RNA_mat, 1, apply(RNA_mat, 1, sd), FUN = '/')


