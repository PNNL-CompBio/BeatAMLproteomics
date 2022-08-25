##get data

source("../util/synapseUtil.R")
library(dplyr)

syn <- synapseLogin()

prot <- syn$tableQuery("select Gene,LogRatio,\"Barcode.ID\" from syn25808020")$asDataFrame()%>%
  tidyr::pivot_wider(names_from='Barcode.ID',values_from='LogRatio')%>%
  tibble::column_to_rownames('Gene')
rna <- syn$tableQuery("select labId, display_label, `RNA counts` from syn26545877")$asDataFrame()%>%
  tidyr::pivot_wider(names_from='labId',values_from='RNA counts')%>%
  tibble::column_to_rownames('display_label')

write.table(prot,file="./AML-tumor-prot-raw.tsv",sep='\t',row.names=T,quote=F)

write.table(rna,file="./AML-tumor-mrna-raw.tsv",sep='\t',row.names=T,quote=F)
