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
setwd("z_scaled_cols")
met.files <- list.files(pattern = ".*.csv", full.names = TRUE)

non.omics.files <- c("./meta_ids.csv", "./cluster_pred.csv", "./meta_labels.csv", "./drug_response.csv")
all.files <- sort(c(data.files[!(data.files %in% non.omics.files)], met.files))
data.sources <- c("Acetyl", "Global", "Lipid", "Met_HILIC", "Met_RP", "Phospho", "RNA", "WES")
met.sources <- c("Lipid", "Met_HILIC", "Met_RP")

# load data
data.list <- list()
setwd(base.path)
setwd("data")
for (i in 1:length(all.files)) {
  if (data.sources[i] %in% met.sources) {
    setwd("z_scaled_cols")
    data.list[[data.sources[i]]] <- read.csv(all.files[i])
    setwd(file.path(base.path, "data"))
  } else {
    data.list[[data.sources[i]]] <- read.csv(all.files[i]) 
  }
}
all.data <- data.table::rbindlist(data.list, use.names = TRUE, fill = TRUE, idcol = "Source")

# count number of features
Data_source <- data.sources
n.features <- data.frame(Data_source)
n.features$N_features <- NA
for (i in 1:length(data.sources)) {
  n.features$N_features[i] <- length(unique(na.omit(data.list[[data.sources[i]]]$label)))
}
n.features <- n.features %>% arrange(desc(N_features))
n.features$Thousand_features <- n.features$N_features/1000
n.features$Log_features <- log(n.features$N_features, 10)
n.features <- n.features[order(n.features$N_features), ]

# bar plot of N features
bp <- barplot(n.features$N_features, ylab = "# of Features", xlab = "Data Type", 
        names.arg = n.features$Data_source, horiz = FALSE, 
        #las = 1, 
        cex.names = 1)
text(x=bp, y = ifelse(n.features$N_features < 20000, n.features$N_features + 1000, 21500), labels = n.features$N_features)

#bar.feat <- ggplot2::ggplot(all.data, aes(x=Source)) + ggplot2::geom_bar() + bg.theme2

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
data.types <- unique(new.results.kmax15$data_type)

### save pearson results (max, avg, median, sd)
pearson.agg <- plyr::ddply(new.results.kmax15, .(drug_name, data_type), summarize,
                           mean_pearson = mean(pearson, na.rm = TRUE),
                           sd_pearson = sd(pearson, na.rm = TRUE),
                           max_pearson = max(pearson, na.rm = TRUE),
                           min_pearson = min(pearson, na.rm = TRUE),
                           median_pearson = median(pearson, na.rm = TRUE))
pearson.agg$minusLogSD_pearson <- -log(pearson.agg$sd_pearson, 10)
write.csv(pearson.agg, "pearson_kMax15.csv", row.names = FALSE)
pearson.median <- reshape2::dcast(pearson.agg, drug_name ~ data_type, value.var = "median_pearson", fill = NA)
write.csv(pearson.median, "median_pearson_kMax15.csv", row.names = FALSE)
pearson.minusLogSD <- reshape2::dcast(pearson.agg, drug_name ~ data_type, value.var = "minusLogSD_pearson", fill = NA)
write.csv(pearson.minusLogSD, "minusLogSD_pearson_kMax15.csv", row.names = FALSE)

# save results in all data_types
setwd(base.path)
setwd("data")
setwd("analysis_BG")
write.csv(pearson.median[, colSums(is.na(pearson.median)) == 0], "median_pearson_kMax15_noMissing.csv", row.names = FALSE)
write.csv(pearson.minusLogSD[, colSums(is.na(pearson.minusLogSD)) == 0], "minusLogSD_pearson_kMax15_noMissing.csv", row.names = FALSE)

### look at # of drugs for which each source is the best predictor
ind_data_types <- c("proteomics", "rna_seq", "phospho", "acetyl", "metabolomics_HILIC", "metabolomics_RP", "lipidomics", "wes")
#individual.results <- na.omit(new.results.kmax15[new.results.kmax15$data_type %in% ind_data_types, ])
ind.pearson.median <- pearson.median[ , c("drug_name", ind_data_types)]
ind.pearson.sd <- pearson.minusLogSD[ , c("drug_name", ind_data_types)]
write.csv(ind.pearson.median, "median_pearson_kMax15_individualSourcesOnly.csv", row.names = FALSE)
write.csv(ind.pearson.sd, "minusLogSD_pearson_kMax15_individualSourcesOnly.csv", row.names = FALSE)

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

# # create pie chart for best individual sources
library(RColorBrewer)
my.labels <- paste0(n.best.ind$best.ind.sources, " (", n.best.ind$N_drugs, ")")
pie(n.best.ind$N_drugs, labels = my.labels, border = "white",
    init.angle=90,
    col = brewer.pal(5, "Set2"))

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
drug_order <- unique(best.df[order(best.df$best_median_pearson), ]$drug_name)
ind_drug_order <- unique(best.df[order(best.df$best_ind_median_pearson), ]$drug_name)

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
ggsave("Best_individual_omics_for_each_drug_v5.pdf")

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
setwd(file.path(base.path, "data", "model_results", "actyl_results"))
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
