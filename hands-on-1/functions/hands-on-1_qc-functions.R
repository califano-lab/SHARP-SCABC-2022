#' Generates basic quality control plots from raw gene expression data.
#' 
#' @param raw.mat Matrix of raw gene expression data (genes X samples).
#' @param mt.genes Character vector of MT genes
QCPlots <- function(raw.mat, mt.genes) {
  # packages
  require(ggplot2)
  require(ggpubr)
  ## sequencing depth plot
  p1.dat <- data.frame('Depth' = colSums(raw.mat), 'Sample' = as.factor(rep("", ncol(raw.mat))))
  p1 <- ggplot(p1.dat, aes(x=Sample, y=Depth)) + geom_violin(color = '#F8766D', fill = '#F8766D') + 
    theme_bw()+ theme(text = element_text(size = 20))+xlab("")
  ## detected gene plot
  p2.dat <- data.frame('dgenes' = colSums(raw.mat > 0), 'Sample' = as.factor(rep("", ncol(raw.mat))))
  p2 <- ggplot(p2.dat, aes(x=Sample, y=dgenes)) + geom_violin(color = '#00BA38', fill = '#00BA38') + 
    ylab('Datected Genes') + theme_bw()+ theme(text = element_text(size = 20))+xlab("")
  ## mt percentage plot
  mt.perc <- MTPercent(raw.mat, mt.genes)*100
  p3.dat <- data.frame('mt' = mt.perc, 'Sample' = as.factor(rep("", length(mt.perc))))
  p3 <- ggplot(p3.dat, aes(x=Sample, y=mt)) + geom_violin(color = '#619CFF', fill = '#619CFF') +
    ylab('MT%') + theme_bw()+ theme(text = element_text(size = 20))+xlab("")
  ggarrange(plotlist = list(p1, p2, p3), ncol = 3)
}


#' Filters out genes with no expression and low quality cells.
#' 
#' @param raw.mat Matrix of raw gene expression data (genes X samples).
#' @param minCount Minimum number of UMIs in a cell. Default of 1000.
#' @param maxCount Maximum number of UMIs in a cell. Default of 100000.
#' @param minGeneReads Minimum number of reads for a gene to be kept. Default of 1 (any gene with no reads will be removed).
#' @return Quality filtered matrix.
QCTransform <- function(raw.mat, minCount = 1000, maxCount = 100000, minGeneReads = 1) {
  filt.mat <- raw.mat[, colSums(raw.mat) > minCount & colSums(raw.mat) < maxCount]
  filt.mat <- filt.mat[ rowSums(filt.mat) >= minGeneReads ,]
  rem.genes <- nrow(raw.mat) - nrow(filt.mat); rem.cells <- ncol(raw.mat) - ncol(filt.mat)
  print(paste('Removed ', rem.cells, ' cells.', sep =''))
  return(filt.mat)
}


#' Returns vector of mitochondrial percentages for the given samples.
#'
#' @param dat.mat Matrix of raw gene expression data (genes X samples).
#' @param mt.genes List of mitochondrial genes
#' @retun Returns named vector of mitochondrial gene percentage.
MTPercent <- function(dat.mat, mt.genes) {
  mt.count <- colSums(dat.mat[ intersect(rownames(dat.mat), mt.genes) ,])
  total.count <- colSums(dat.mat)
  mt.percent <- mt.count / total.count
  head(mt.percent)
  return( mt.percent )
}

#' Filters data based on percentage of mitochondrial gens.
#'
#' @param raw.mat Matrix of raw gene expression data (genes X samples)
#' @param mt.genes Character vector containing the names of the MT genes
#' @param mt.thresh Threshold above which cells will be removed. Default of 0.15
#' @return Quality filtered matrix.
MTFilter <- function(dat.mat, mt.genes, mt.thresh = 0.1) {
  ## find mt percentages
  mt.perc <- MTPercent(dat.mat, mt.genes)
  ## filter matrix
  thresh.cells <- names(mt.perc)[which(mt.perc < mt.thresh)]
  rem.cells <- ncol(dat.mat) - length(thresh.cells)
  print(paste('Removed', rem.cells, 'cells with too many MT reads', sep = ' '))
  dat.mat <- dat.mat[, thresh.cells ]
  return(dat.mat)
}