#Function to draw heatmap of genes in seurat object dat, with cluster labels clust
#seurat object must have labels "tissue" and "patient" 
gExprHeatmap=function(genes,expression_to_plot,Clusters,colors){
  require(pheatmap)
  annotation_columns<-data.frame(Clusters)
  rownames(annotation_columns)<-colnames(expression_to_plot)
  colnames(annotation_columns)<-c("clusters")
  
  anno_colors <- list(clusters = colors)
  names(anno_colors[[1]]) <- levels(annotation_columns$cluster)
  
  cellType_Order<-unique(Clusters)
  idx<-c()
  for (i in 1:length(cellType_Order)){
    idx<-c(idx, which(annotation_columns[,1]==cellType_Order[i]))
  }
  
  expression_to_plot=expression_to_plot[,idx]
  annotation_columns=data.frame(annotation_columns[idx,])
  rownames(annotation_columns)<-colnames(expression_to_plot)
  colnames(annotation_columns)=c("clusters")
  
  quantile_breaks <- function(xs, n = 10) {
    breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
    breaks[!duplicated(breaks)]
  }
  mat_breaks <- quantile_breaks(as.matrix(apply(expression_to_plot,1,function(x){(x-mean(x))/sd(x)})), n = 30)
  
  library(pheatmap)
  file=c("Gene_Expression_Heatmap")
  p<-pheatmap(expression_to_plot, cluster_rows=FALSE, 
              show_rownames=TRUE, cluster_cols=FALSE, 
              annotation_col=annotation_columns, breaks=mat_breaks, 
              color = colorRampPalette(colors = c('blue', 'white', 'red'))(length(mat_breaks)),
              fontsize_row = 8, show_colnames = FALSE, annotation_colors = anno_colors,
              scale="row")
  pdf(file= paste(file, "pdf", sep="."), height = 10, width=15)
  print(p)
  no_show<-dev.off()
  
}




