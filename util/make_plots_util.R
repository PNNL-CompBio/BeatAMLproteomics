## This script contain wrappers around pheatmap and others
## in order to save the picture locally, and have the choice to upload it
## to a synapse folder as a pdf or png.


library(pheatmap)
library(nationalparkcolors)
library(RColorBrewer)


data_type_colors <<- park_palette("Badlands", 4) %>% as.character()
subtype_colors <<- brewer.pal(8,'Dark2')
## Rearranging a bit.
subtype_colors <<- c(subtype_colors[[2]], subtype_colors[[3]],
                     subtype_colors[[5]], subtype_colors[[6]], 
                     subtype_colors[[1]], subtype_colors[[4]],
                     subtype_colors[[7]], subtype_colors[[8]])



make.pheatmap <- function(mat, filename, width = 7, height = 7, 
                          format = "pdf", parentId = NULL, ...) {
  filename = paste(filename, format, sep = ".")
  pheatmap(mat, filename = filename, height = height, 
           width = width, ...)
  
  if (!is.null(parentId)) {
    synapseStore(filename, parentId)
  }
}

### Used for plots which do not have a properly flexible wrapper. 
### For example, when uploading ggplot plots
upload.plot <- function(filename, parentId) {
  synapseStore(filename, parentId)
}

