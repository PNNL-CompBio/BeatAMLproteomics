## This script contain wrappers around pheatmap and others
## in order to save the picture locally, and have the choice to upload it
## to a synapse folder as a pdf or png.


library(pheatmap)


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

