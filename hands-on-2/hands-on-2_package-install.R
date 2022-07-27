#' Packages should be pre-installed in the project.
#' If for some reason they are not, run the code below to install all
#' necessary dependencies for the first Hands-on Session.

## install cran packages
install.packages("BiocManager")
install.packages('devtools')
install.packages('DirichletReg')
install.packages('ggplot2')
install.packages('RSpectra')
install.packages('uwot')
## install bioconductor packages
BiocManager::install("biomaRt")
BiocManager::install("celldex")
BiocManager::install("ComplexHeatmap")
BiocManager::install("densvis")
BiocManager::install("SingleR")
## install pisces from github
devtools::install_github("califano-lab/PISCES")
