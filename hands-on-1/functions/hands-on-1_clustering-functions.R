#' Application of the PAM algorithm for various values of the k parameter and 
#' computation of the Silhouette score
#' 
#' @param dist.mat Dissimilarity matrix
#' @param kmin Minimum value for the k parameter
#' @param kmax Maximum value for the k parameter
pam_k <- function(dist.mat, kmin = 2, kmax = 5) {
  require(cluster)
  
  # generate clustering for each value of k
  cat('Generating clusterings...\n')
  clustering.objects <- lapply(kmin:kmax, function(k) {pam(dist.mat, k, diss = TRUE)})
  # generate silhouette score for each cluster
  cat('Generating silhouette scores...\n')
  sil.scores <- sapply(clustering.objects, function(clust.obj) { mean(silhouette(clust.obj, dist.mat)[,3]) })
  # name the lists
  names(clustering.objects) <- paste('k', kmin:kmax, sep = '')
  names(sil.scores) <- paste('k', kmin:kmax, sep = '')
  # identify the optimal clustering
  opt.clust <- clustering.objects[[which.max(sil.scores)]]$clustering
  # return objects
  ret.list <- list('opt.clust' = opt.clust,
                   'clustering.objs' = clustering.objects,
                   'sil.scores' = sil.scores)
  return(ret.list)
}


#' Silhouette Scores computation for the different Louvain clustering solutions
#' 
#' @param Seurat_obj Seurat object containing the different Louvain clustering solutions
#' @param dist.mat Dissimilarity matrix

Louvain_SilhouetteScore <- function(Seurat_obj, dist.mat) {
  require(cluster)
  
  clusters<-Seurat_obj@meta.data[,which(grepl("SCT_snn_res.",colnames(Seurat_obj@meta.data)))]
  sil.scores <- apply(clusters, 2, function(x) { mean(silhouette(as.numeric(x), dist.mat)[,3]) })
  return(sil.scores)
}