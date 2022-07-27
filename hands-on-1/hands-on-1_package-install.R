#' Packages should be pre-installed in the project.
#' If for some reason they are not, run the code below to install all
#' necessary dependencies for the first Hands-on Session.

## install cran packages
install.packages('ggplot2')
install.packages('ggpubr')
install.packages('Seurat')
install.packages('cluster')
install.packages('factoextra')
install.packages('dplyr')
install.packages("BiocManager")
## install bioconductor packages
BiocManager::install("celldex")
BiocManager::install("pheatmap")
BiocManager::install("SingleR")
## install umap dependency w/ python
library(reticulate)
virtualenv_create("hands_on_1")use_virtualenv("hands_on_1", required = TRUE)
py_install("umap-learn")