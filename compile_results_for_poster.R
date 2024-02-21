# Compile Beat AML model results w acetyl, met, lipidomics
# Author: Belinda B. Garana
# Date created: 2024-02-20
# Last edit: 2024-02-20

library(plyr); library(dplyr); library(ggplot2); library(ggvenn)

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
drug_name <- unique(pearson.median.long$drug_name)
best.df <- data.frame(drug_name)
best.df$best_source <- NA
best.df$best_ind_source <- NA
for (i in 1:length(drug_name)) {
  drug.data <- pearson.median.long[pearson.median.long$drug_name == drug_name[i], ]
  drug.ind.data <- pearson.median.ind.long[pearson.median.ind.long$drug_name == drug_name[i], ]
  best.df$best_source[i] <- drug.data[which.max(drug.data$value), ]$data_type
  best.df$best_ind_source[i] <- drug.ind.data[which.max(drug.ind.data$value), ]$data_type
}
write.csv(best.df, "Best_sources_for_each_drug.csv", row.names = FALSE)

best.sources <- unique(best.df$best_source)
best.ind.sources <- unique(best.df$best_ind_source)

n.best <- data.frame(best.sources)
n.best$N_drugs <- NA
n.best$drugs <- NA
for (i in 1:length(best.sources)) {
  n.best$N_drugs[i] <- nrow(best.df[best.df$best_source == best.sources[i], ])
  n.best$drugs[i] <- paste0(best.df[best.df$best_source == best.sources[i], ]$drug_name, collapse = "|")
}

n.best.ind <- data.frame(best.ind.sources)
n.best.ind$N_drugs <- NA
n.best.ind$drugs <- NA
for (i in 1:length(best.ind.sources)) {
  n.best.ind$N_drugs[i] <- nrow(best.df[best.df$best_ind_source == best.ind.sources[i], ])
  n.best.ind$drugs[i] <- paste0(best.df[best.df$best_ind_source == best.ind.sources[i], ]$drug_name, collapse = "|")
}

# create pie chart for best individual sources
library(RColorBrewer)
my.labels <- paste0(n.best.ind$best.ind.sources, " (", n.best.ind$N_drugs, ")")
pie(n.best.ind$N_drugs, labels = my.labels, border = "white",
    col = brewer.pal(5, "Set2"))

for (i in ind_data_types) {
  source.results <- pearson.agg.long[pearson.agg.long$data_type == i,]
}

### look at effect of acetyl

## find data combos w acetyl

## find equivalent wo acetyl


