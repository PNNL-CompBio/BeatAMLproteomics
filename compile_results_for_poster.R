# Compile Beat AML model results w acetyl, met, lipidomics
# Author: Belinda B. Garana
# Date created: 2024-02-20
# Last edit: 2024-02-20

library(plyr); library(dplyr); library(ggplot2); library(ggvenn)
library(data.table); library(reshape2)
# load theme for plots
bg.theme2 <- ggplot2::theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.border = element_rect(fill = NA),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"), 
  axis.title.x = element_text(size = 20, colour = "black"),
  axis.title.y = element_text(size = 20, colour = "black"),
  axis.text.x = element_text(colour = "black", size =16), 
  axis.text.y = element_text(colour = "black", size = 16),
  axis.ticks.x = element_line(colour = "black"), 
  axis.ticks.y = element_line(colour = "black"),
  legend.title = element_blank(), legend.background = element_rect(), 
  legend.position = "top",
  legend.text = element_text(size = 14), legend.key = element_blank(),
  plot.title = element_text(lineheight = .8, face = "bold", size = 36)
)

base.path <- "~/OneDrive - PNNL/Documents/GitHub/BeatAMLproteomics/pybeataml"
setwd(base.path)
setwd("data")

#### look at features ####
data.files <- list.files(pattern = ".*.csv", full.names = TRUE)

non.omics.files <- c("./meta_ids.csv", "./cluster_pred.csv", "./meta_labels.csv", "./drug_response.csv")
all.files <- sort(data.files[!(data.files %in% non.omics.files)])
data.sources <- c("Acetyl", "Global", "Lipid", "Met_HILIC", "Met_RP", "Phospho", "RNA", "WES")

# load data
data.list <- list()
setwd(base.path)
setwd("data")
for (i in 1:length(all.files)) {
  data.list[[data.sources[i]]] <- read.csv(all.files[i]) 
}
all.data <- data.table::rbindlist(data.list, use.names = TRUE, fill = TRUE, idcol = "plotLabel")

# count number of features
n.features <- plyr::ddply(all.data, .(plotLabel), summarize,
                          N_features = length(unique(na.omit(label))))
n.features <- n.features[order(n.features$N_features), ]

dir.create("analysis_BG")
setwd("analysis_BG")

# bar plot of N features
bp <- barplot(n.features$N_features, ylab = "# of Features", xlab = "Omics Type", 
        names.arg = n.features$plotLabel, horiz = FALSE, 
        #las = 1, 
        cex.names = 1,
        cex.lab = 1.5)
text(x=bp, y = ifelse(n.features$N_features < 20000, n.features$N_features + 2000, 21000), labels = n.features$N_features, cex=1.5)

# venn diagram of shared gene_symbols
venn.list <- list("Acetyl" = unique(na.omit(all.data[all.data$Source == "Acetyl", ]$gene_symbol)),
                  "Global" = unique(na.omit(all.data[all.data$Source == "Global", ]$gene_symbol)),
                  #"Phospho" = unique(na.omit(all.data[all.data$Source == "Acetyl", ]$gene_symbol)),
                  "RNA" = unique(na.omit(all.data[all.data$Source == "RNA", ]$gene_symbol)),
                  "WES" = unique(na.omit(all.data[all.data$Source == "WES", ]$gene_symbol)))
ggvenn::ggvenn(venn.list)
ggplot2::ggsave("Gene_symbol_venn_diagram.pdf")

venn.list <- list("Acetyl" = unique(na.omit(all.data[all.data$Source == "Acetyl", ]$gene_symbol)),
                  "Global" = unique(na.omit(all.data[all.data$Source == "Global", ]$gene_symbol)),
                  "Phospho" = unique(na.omit(all.data[all.data$Source == "Phospho", ]$gene_symbol)),
                  "RNA" = unique(na.omit(all.data[all.data$Source == "RNA", ]$gene_symbol))
                  #,
                  #"WES" = unique(na.omit(all.data[all.data$Source == "WES", ]$gene_symbol))
                  )
ggvenn::ggvenn(venn.list)
ggplot2::ggsave("Gene_symbol_venn_diagram_w_phospho.pdf")
ggvenn::ggvenn(venn.list, show_percentage = FALSE)
ggplot2::ggsave("Gene_symbol_venn_diagram_w_phospho_wo_percentages.pdf")
ggvenn::ggvenn(venn.list, show_percentage = FALSE, set_name_size = 5)
ggplot2::ggsave("Gene_symbol_venn_diagram_w_phospho_wo_percentages_smallerFont.pdf")
ggvenn::ggvenn(venn.list, show_percentage = FALSE, set_name_size = 5.5, text_size = 5.5)
ggplot2::ggsave("Gene_symbol_venn_diagram_w_phospho_wo_percentages_biggerFont.pdf")

#### compile newest results ####
setwd("model_results")
result.files <- list.files(pattern = ".*.csv", full.names = TRUE)
result.list <- list()
for (i in 1:length(result.files)) {
  result.list[[i]] <- read.csv(result.files[i])
}
new.results <- data.table::rbindlist(result.list, use.names = TRUE, fill = TRUE)
new.results$X <- NULL
new.results.kmax15 <- new.results[new.results$k < 15, ]
write.csv(new.results.kmax15, "Results_kmax15_2024-02-21.csv", row.names = FALSE)
new.results.kmax15 <- read.csv("Results_kmax15_2024-02-21.csv")
data.types <- unique(new.results.kmax15$data_type)

### save pearson results (max, avg, median, sd)
pearson.agg <- plyr::ddply(new.results.kmax15, .(drug_name, data_type), summarize,
                           mean_pearson = mean(pearson, na.rm = TRUE),
                           sd_pearson = sd(pearson, na.rm = TRUE),
                           max_pearson = max(pearson, na.rm = TRUE),
                           min_pearson = min(pearson, na.rm = TRUE),
                           median_pearson = median(pearson, na.rm = TRUE),
                           feature_importance_most_to_least = paste0(sort(table(
                             stringr::str_split(feature_names,"[|]")[[1]]), decreasing = TRUE), collapse = "|"),
                           features_most_to_least_important = paste0(names(sort(table(
                             stringr::str_split(feature_names,"[|]")[[1]]), decreasing = TRUE)), collapse = "|"))
pearson.agg$minusLogSD_pearson <- -log(pearson.agg$sd_pearson, 10)
write.csv(pearson.agg, paste0("pearson_kMax15_", Sys.Date(), ".csv"), row.names = FALSE)
pearson.median <- reshape2::dcast(pearson.agg, drug_name ~ data_type, value.var = "median_pearson", fill = NA)
write.csv(pearson.median, paste0("median_pearson_kMax15_", Sys.Date(), ".csv"), row.names = FALSE)
pearson.minusLogSD <- reshape2::dcast(pearson.agg, drug_name ~ data_type, value.var = "minusLogSD_pearson", fill = NA)
write.csv(pearson.minusLogSD, paste0("minusLogSD_pearson_kMax15_", Sys.Date(), ".csv"), row.names = FALSE)

# save results in all data_types
setwd(base.path)
setwd("data")
setwd("analysis_BG")
write.csv(pearson.median[, colSums(is.na(pearson.median)) == 0], paste0("median_pearson_kMax15_noMissing_", Sys.Date(), ".csv"), row.names = FALSE)
write.csv(pearson.minusLogSD[, colSums(is.na(pearson.minusLogSD)) == 0], paste0("minusLogSD_pearson_kMax15_noMissing_", Sys.Date(), ".csv"), row.names = FALSE)

### look at # of drugs for which each source is the best predictor
ind_data_types <- c("proteomics", "rna_seq", "phospho", "acetyl", "metabolomics_HILIC", "metabolomics_RP", "lipidomics", "wes")
#individual.results <- na.omit(new.results.kmax15[new.results.kmax15$data_type %in% ind_data_types, ])
ind.pearson.median <- pearson.median[ , c("drug_name", ind_data_types)]
ind.pearson.sd <- pearson.minusLogSD[ , c("drug_name", ind_data_types)]
write.csv(ind.pearson.median, paste0("median_pearson_kMax15_individualSourcesOnly_", Sys.Date(), ".csv"), row.names = FALSE)
write.csv(ind.pearson.sd, paste0("minusLogSD_pearson_kMax15_individualSourcesOnly_", Sys.Date(), ".csv"), row.names = FALSE)

best.df <- data.frame(ind_data_types)
best.df$best_drugs_by_median_pearson <- NA
best.df$N_best_drugs_by_median_pearson <- NA
pearson.agg.long <- reshape2::melt(pearson.agg, id.vars = c("drug_name", "data_type"))
pearson.median.long <- pearson.agg.long[pearson.agg.long$variable == "median_pearson", ]
pearson.median.ind.long <- pearson.median.long[pearson.median.long$data_type %in% ind_data_types, ]
# identify best data source for each drug
drug <- unique(pearson.median.long$drug_name)
best.df <- data.frame(drug)
best.df$best_source <- NA
best.df$best_ind_source <- NA
best.df$best_median_pearson <- NA
best.df$best_ind_median_pearson <- NA

best.results <- list()
best.ind.results <- list()
all.best.sources <- c()
all.ind.best.sources <- c()

all.sources <- na.omit(unique(new.results.kmax15$data_type))
n.best <- data.frame(all.sources)
n.best$N_drugs <- 0
n.best$drugs <- ""
for (i in 1:length(drug)) {
  # isolate drug data
  drug.data <- pearson.median.long[pearson.median.long$drug_name == drug[i], ]
  drug.ind.data <- pearson.median.ind.long[pearson.median.ind.long$drug_name == drug[i], ]
  drug.results <- new.results.kmax15[new.results.kmax15$drug_name == drug[i], ]
  
  # identify best median_pearson
  best_median_pearson <- unique(drug.data[which.max(drug.data$value), ]$value)
  best_ind_median_pearson <- unique(drug.ind.data[which.max(drug.ind.data$value), ]$value)
  best.df$best_median_pearson[i] <- best_median_pearson
  best.df$best_ind_median_pearson[i] <- best_ind_median_pearson
  
  # identify sources with best median_pearson
  best.sources <- na.omit(unique(drug.data[drug.data$value == best_median_pearson, ]$data_type))
  best.ind.sources <- na.omit(unique(drug.ind.data[drug.ind.data$value == best_ind_median_pearson, ]$data_type))
  all.best.sources <- unique(c(all.best.sources, best.sources))
  all.ind.best.sources <- unique(c(all.ind.best.sources, best.ind.sources))
  best.df$best_source[i] <- paste0(best.sources, collapse = "|")
  best.df$best_ind_source[i] <- paste0(best.ind.sources, collapse = "|")
  
  # store results with best median_pearson
  best.results[[i]] <- drug.results[drug.results$data_type %in% best.sources, ]
  best.ind.results[[i]] <- drug.results[drug.results$data_type %in% best.ind.sources, ]
  
  # keep track of drugs for which each source is the best
  complete.best.sources <- unique(c(best.sources, best.ind.sources))
  for (j in complete.best.sources) {
    n.best[n.best$all.sources == j, ]$N_drugs <- 
      n.best[n.best$all.sources == j, ]$N_drugs + 1
    n.best[n.best$all.sources == j, ]$drugs <- 
      paste0(c(n.best[n.best$all.sources == j, ]$drugs, drug[i]), collapse="|") 
  }
}
write.csv(best.df, paste0("Best_sources_for_each_drug_", Sys.Date() ,".csv"), row.names = FALSE)
best.results.df <- data.table::rbindlist(best.results, use.names = TRUE, fill = TRUE)
best.ind.results.df <- data.table::rbindlist(best.ind.results, use.names = TRUE, fill = TRUE)
write.csv(best.results.df, paste0("Best_source_results_for_each_drug_", Sys.Date() ,".csv"), row.names = FALSE)
write.csv(best.ind.results.df, paste0("Best_individual_source_results_for_each_drug_", Sys.Date() ,".csv"), row.names = FALSE)
write.csv(n.best, paste0("Drugs_best_predicted_by_each_source_", Sys.Date() ,".csv"), row.names = FALSE)

best.df <- read.csv("Best_sources_for_each_drug_2024-02-21.csv")
best.results.df <- read.csv("Best_source_results_for_each_drug_2024-02-21.csv")
best.ind.results.df <- read.csv("Best_individual_source_results_for_each_drug_2024-02-21.csv")
n.best <- read.csv("Drugs_best_predicted_by_each_source_2024-02-21.csv")

# # create pie chart for best individual sources
# order n.best.ind by alphabetical order of omics type so colors match with violin plot
n.best.ind <- n.best[n.best$all.sources %in% ind_data_types, ]
n.best.ind <- n.best.ind[order(n.best.ind$all.sources), ]
n.best.ind.nonZero <- n.best.ind[n.best.ind$N_drugs > 0, ]
library(RColorBrewer)
my.labels <- paste0(n.best.ind.nonZero$all.sources, " (", n.best.ind.nonZero$N_drugs, ")")
pie(n.best.ind.nonZero$N_drugs, labels = my.labels, border = "white",
    init.angle=90,
    cex=1.5,
    col = brewer.pal(length(my.labels), "Set2"))

# 
# # create pie chart for best sources
# library(RColorBrewer)
# my.labels <- paste0(n.best$best.sources, " (", n.best$N_drugs, ")")
my.pal <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))
# pie(n.best$N_drugs, labels = my.labels, border = "white",
#     col = my.pal(length(my.labels)), radius = 1, cex = 0.5)

# create bar chart of # drugs best predicted for each source
n.best.filtered <- n.best[n.best$N_drugs > 0, ]
n.best.filtered <- n.best.filtered[order(n.best.filtered$N_drugs), ]
par(mar=c(20.4,4,1,1))
bp <- barplot(n.best.filtered$N_drugs, ylab = "# of Drugs Best Predicted", xlab = "", 
              names.arg = n.best.filtered$all.sources, 
              #horiz = TRUE, 
              las = 2,
              srt=45,
              cex.names = 0.57
              )

# 
# par(mar=c(4,20.4,0,0))
# bp <- barplot(n.best.filtered$N_drugs, xlab = "# of Drugs Best Predicted", ylab = "", 
#               names.arg = n.best.filtered$all.sources, 
#               horiz = TRUE, 
#               las = 1,
#               cex.names = 0.5
# )
#text(x=bp, y = ifelse(n.features$N_features < 20000, n.features$N_features + 1000, 21500), labels = n.features$N_features)

### look at how many drugs each source is the best for
new.results.kmax15 <- read.csv("Results_kmax15_2024-02-21.csv")


# sort by median_pearson
drug_order <- unique(best.df[order(best.df$best_median_pearson), ]$drug)
ind_drug_order <- unique(best.df[order(best.df$best_ind_median_pearson), ]$drug)

## violin plot: best sources
my.pal <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))
my.pal2 <- my.pal(length(unique(best.results.df$data_type)))
best.violin <- ggplot2::ggplot(best.results.df, 
                               aes(fill = data_type, x=drug_name, y=pearson)) + 
  geom_violin(position=position_dodge(width=0.4), alpha=0.5) + scale_x_discrete(limits=drug_order) +
  # guides(fill=guide_legend(nrow=30, byrow=TRUE)) +
  # guides(shape = guide_legend(override.aes = list(size = 0.1))) +
  # guides(color = guide_legend(override.aes = list(size = 0.1))) +
  # theme(legend.title = element_text(size = 0.1), 
  #       legend.text = element_text(size = 0.1)) + 
  scale_fill_manual(values = my.pal(91)) + 
  geom_boxplot(width=0.1, position = position_dodge(width=0.4), alpha=0.5) + bg.theme2 +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + 
  #ggtitle("Individual Omics Most Predictive of Each Drug") +
  xlab("Drug") + ylab("Pearson Estimate") + theme(legend.position = "none")
# legend not useful because it's hard to tell which data source best predicted each drug's sensitivity
ggsave("Best_omics_combos_for_each_drugs.pdf")

## violin plot: best individual sources
best.violin <- ggplot2::ggplot(best.ind.results.df, 
                               aes(fill = data_type, x=drug_name, y=pearson)) + 
  geom_violin(position="dodge", alpha = 0.5) + scale_x_discrete(limits=ind_drug_order) +
  geom_boxplot(width=0.1) + bg.theme2 + 
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + 
  #ggtitle("Individual Omics Most Predictive of Each Drug") +
  xlab("Drug") + ylab("Pearson Estimate")
ggsave("Best_individual_omics_for_each_drug_v11.pdf")

#+
  # theme(legend.position = "bottom") + 
  # guides(fill=guide_legend(nrow=10, byrow=TRUE)) +
  # guides(shape = guide_legend(override.aes = list(size = 0.5))) +
  # guides(color = guide_legend(override.aes = list(size = 0.5))) +
  # theme(legend.title = element_text(size = 3), 
  #       legend.text = element_text(size = 3))

# not sure why drug.results above is only getting results for drug_name[1]
# try to compile it right
best.ind.results <- list()
drugs <- unique(new.results.kmax15$drug_name)
for (i in 1:length(drugs)) {
  best.ind.results[[i]] <- new.results.kmax15[new.results.kmax15$drug_name == drugs[i] &
                                                new.results.kmax15$data_type == best.df$best_ind_source[i], ]
}

### look at effect of acetyl

## find data combos w acetyl

## find equivalent wo acetyl

### look at effect of # of patients
setwd(file.path(base.path, "data", "model_results", "wo_lipid_met_for_210_patients"))
# load results w 210 patients
all.files210 <- list.files(pattern = ".*.csv", full.names = TRUE)
data.list210 <- list()
for (i in 1:length(all.files210)) {
  data.list210[[i]] <- read.csv(all.files210[i]) 
}
all.data210 <- data.table::rbindlist(data.list210, use.names = TRUE, fill = TRUE)

# gather venetoclax data
all.data210$N_patients <- "Up to 210 patients"
new.results.kmax15$N_patients <- "Up to 81 patients"
all.data210$X <- NULL
ven.all.data210 <- all.data210[all.data210$drug_name == "Venetoclax", ]
data.types210 <- unique(all.data210$data_type)
ven.new <- new.results.kmax15[new.results.kmax15$drug_name == "Venetoclax" & 
                                new.results.kmax15$data_type %in% data.types210, ]
ven.data <- rbind(ven.all.data210, ven.new)
ven.agg <- plyr::ddply(ven.data, .(drug_name, data_type), summarize,
                                      mean_pearson = mean(pearson, na.rm = TRUE),
                                      sd_pearson = sd(pearson, na.rm = TRUE),
                                      max_pearson = max(pearson, na.rm = TRUE),
                                      min_pearson = min(pearson, na.rm = TRUE),
                                      median_pearson = median(pearson, na.rm = TRUE))
ven.agg <- ven.agg[order(ven.agg$median_pearson), ]
ven.order <- ven.agg$data_type

ven.violin <- ggplot2::ggplot(ven.data, 
                               aes(fill = N_patients, x=data_type, y=pearson)) + 
  geom_violin(position=position_dodge(width=0.4), alpha = 0.5) + scale_x_discrete(limits=ven.order) +
  geom_boxplot(width=0.1, position = position_dodge(width=0.4), alpha = 0.5) + bg.theme2 + 
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + 
  #ggtitle("Individual Omics Most Predictive of Each Drug") +
  xlab("Drug") + ylab("Pearson Estimate")
ven.violin <- ggplot2::ggplot(ven.data, 
                              aes(fill = N_patients, x=data_type, y=pearson)) + 
  geom_violin(position=position_dodge(width=0.4), alpha = 0.5) + scale_x_discrete(limits=ven.order) +
  geom_boxplot(width=0.1, position = position_dodge(width=0.4), alpha = 0.5) + bg.theme2 + 
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + 
  #ggtitle("Individual Omics Most Predictive of Each Drug") + 
  xlab("Omics Data Type") + ylab("Pearson Estimate") + theme(axis.text=element_text(size=3))
ggsave("Venetoclax_210_vs_81_patients_v3.pdf")

# gather pano data
pano.data <- all.data210[all.data210$drug_name == "Panobinostat", ] # pano is not in the newer data w 81 patients max
pano.agg <- plyr::ddply(pano.data, .(drug_name, data_type), summarize,
                       mean_pearson = mean(pearson, na.rm = TRUE),
                       sd_pearson = sd(pearson, na.rm = TRUE),
                       max_pearson = max(pearson, na.rm = TRUE),
                       min_pearson = min(pearson, na.rm = TRUE),
                       median_pearson = median(pearson, na.rm = TRUE))
pano.agg <- pano.agg[order(pano.agg$median_pearson), ]
pano.order <- pano.agg$data_type

pano.violin <- ggplot2::ggplot(pano.data, 
                              aes(x=data_type, y=pearson)) + 
  geom_violin(position=position_dodge(width=0.4), alpha = 0.5) + scale_x_discrete(limits=pano.order) +
  geom_boxplot(width=0.1, position = position_dodge(width=0.4), alpha = 0.5) + bg.theme2 + 
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) + 
  #ggtitle("Individual Omics Most Predictive of Each Drug") +
  xlab("Drug") + ylab("Pearson Estimate")
ggsave("Panobinostat_210patients.pdf")

# make a table of the top 5 features for each individual omics type's best drug predictions
best.feat <- merge(best.ind.results.df, pearson.agg)
best.feat <- plyr::ddply(best.feat, .(data_type), summarize,
                         feature_importance_most_to_least = paste0(sort(table(
                           stringr::str_split(feature_names,"[|]")[[1]]), decreasing = TRUE), collapse = "|"),
                         features_most_to_least_important = paste0(names(sort(table(
                           stringr::str_split(feature_names,"[|]")[[1]]), decreasing = TRUE)), collapse = "|"),
                         top_5_features = paste0(names(sort(table(
                           stringr::str_split(feature_names,"[|]")[[1]]), decreasing = TRUE)[1:5]), collapse = "|"))

# # remove data_type tags
# tag.map <- list("acetyl" = "",
#                 "lipidomics" = "_lip",
#                 "phospho" = "",
#                 "proteomics" = "_prot",
#                 "rna_seq" = "_rna",
#                 "metabolomics_HILIC" = "_met_HILIC",
#                 "metabolomics_RP" = "_met_RP")
# for (i in 1:nrow(best.feat)) {
#   temp.omics <- best.feat$data_type[i]
#   temp.tag <- tag.map[[temp.omics]]
#   temp.features <- stringr::str_split(best.feat$top_5_features, "[|]")
#   best.feat$top_5_features[i] <- paste0(sub(temp.tag, "", temp.features), collapse = ", ")
# }
# #best.feat$top_5_features2 <- gsub("^.*_", "", best.feat$top_5_features) # remove data_type tag
# #best.feat[grepl("met", best.feat$top_5_features)]
write.csv(best.feat, "Best_features_for_individual_omics_for_best_drug_predictions.csv", row.names = FALSE)
