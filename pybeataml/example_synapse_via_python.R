# Title     : Synapse functions in R via python
# Objective : Demonstrate how to use the python class in R
# Created by: pinojc
# Created on: 12/14/2021
library(reticulate)


reticulate::source_python("load_data_from_synpase.py")
# 4 main functions
# load_excel, load_table, load_file

# NMF meta genes
nmf_meta <- load_file('syn26718015')

# pilot project model performance
model_performance <- load_table('syn26469964')


# load meta info
meta <-  load_excel('syn26532699')