##get phospho substrate sequences

#here we want to get phospho substrate sequences for processing from the cantley lab

library(dplyr)
source("../util/synapseUtil.R")

syn<-synapseLogin()
tab<-read.table(syn$get('syn25714902')$path)

write.table(rownames(tab),'phosphoSites.txt',quote=F,row.names = F)