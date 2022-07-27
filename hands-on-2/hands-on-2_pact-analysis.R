#' Load packages and Data
set.seed(343)
library(celldex)
library(ggplot2)
library(PISCES)
library(Seurat)
library(SingleR)
blueprint.encode <- BlueprintEncodeData()
# pbmc data
pbmc.raw <- readRDS('data/pbmc_raw-counts.rds')
# ARACNe networks
lymphoid.net <- readRDS('data/lymphoid_net.rds')
myeloid.net <- readRDS('data/myeloid_net.rds')
# ADT markers
adt.counts <- readRDS('data/adt_cpm-counts.rds')
adt.set.df <- readRDS('data/adt_marker-df.rds')

#' Step 1 - Initial Gene Expression Clustering
#' Here, we utilize the Seurat pipelien to perform initial clustering analysis.
#' This is identical to what you did in the first hands-on session.
#' Our goal here is to identify coarse groupings for network generation.
#' Note that you can apply many different approaches here.
pbmc.filt <- qc_filt(pbmc.raw, max.depth = 25000, max.mt = 0.25)
pbmc.seurat <- CreateSeuratObject(counts = pbmc.filt, assay = 'RNA', project = 'pbmc')
pbmc.seurat <- PercentageFeatureSet(pbmc.seurat, pattern = "^MT-", col.name = "percent.mt")
pbmc.seurat <- SCTransform(pbmc.seurat, vars.to.regress = "percent.mt", verbose = FALSE)
pbmc.seurat <- RunPCA(pbmc.seurat, verbose = FALSE)
pbmc.seurat <- RunUMAP(pbmc.seurat, dims = 1:30, verbose = FALSE)
pbmc.seurat <- FindNeighbors(pbmc.seurat, dims = 1:30, verbose = FALSE)
pbmc.seurat <- FindClusters(pbmc.seurat, verbose = FALSE)
DimPlot(pbmc.seurat)
## we'll also run singleR as before
bp.singleR <- SingleR(test = pbmc.filt, ref = blueprint.encode, 
                      labels = blueprint.encode$label.main, assay.type.test = 1)
lymphoid.cell.types <- c('B-cells', 'CD4+ T-cells', 'CD8+ T-cells', 'NK cells')
myeloid.cell.types <- c('Monocytes')
lymphoid.cells <- colnames(pbmc.filt)[which(bp.singleR$first.labels %in% lymphoid.cell.types)]
myeloid.cells <- colnames(pbmc.filt)[which(bp.singleR$first.labels %in% myeloid.cell.types)]

#' Step 2 - Prepare matrices for network generation.
#' Below is an example of how you can use PISCES to generate matrices for network generation w/ ARACNe3.
#' Running ARACNe3 for single-cells requires an HPC and is not feasible in a bootcamp setting.
#' However, we've pre-generated the networks you need for this session.
#' Additionally, you can find the script we use on our HPC in the Github.
## construct a distance matrix
pca.dist <- dist(pbmc.seurat@reductions$pca@cell.embeddings[,1:30])
## generate meta cells; NOTE that we've lowered the min.samps for purposes of demonstration ONLY
meta.cells <- make_metacells(pbmc.filt, pca.dist, pbmc.seurat$seurat_clusters,
                             num.neighbors = 5, min.samps = 100)

#' Step 3 - Generate an internal Gene Expression Signature.
#' For this exercise, we'll focus on the identification of populations or subclusters in this data.
#' Here, we use the Seurat SCT transformation to generate a signature.
#' Note that there are many other options here, as we discussed in the lecture.
#' Explore the PISCES codebase further for implementations of some of these other options.
ges.mat <- as.matrix(pbmc.seurat@assays$SCT@scale.data)
net.list <- list('lymphoid' = lymphoid.net, 'myeloid' = myeloid.net)
## alternatively, you can use PISCES' built in method (OPTIONAL)
internal.ges.mat <- internal_ges(pbmc.filt, norm.method = 'pflpf', est.method = 'map')

#' Step 4 - Protein Activity Signature Generation
#' We now use lineage specific networks and the NaRnEA algorithm to generate protein activity signatures.
#' Here, we use a matched network approach; the lymphoid network will be used for lymphoid cells, while
#' the myeloid network will be used for myeloid cells. MetaNaRnEA will integrate all networks for unlabeled cells.
#' You can use MetaNaRnEA in situations where lineage relationships are less evident.
## create network matching vector, labeling lymphoid cells, myeloid cells, and unmatched cells
net.match.vec <- rep(NA, ncol(ges.mat)); names(net.match.vec) <- colnames(ges.mat)
net.match.vec[myeloid.cells] <- 'myeloid'
net.match.vec[lymphoid.cells] <- 'lymphoid'
## run narnea
narnea.obj <- network_match_narnea(ges.mat, net.list, net.match.vec)

#' Step 5 - Clustering Analysis
#' We now use the Louvain algorithm to identify clusters within our protein activity signature.
#' For our distance metric, we use Euclidean distance calculated from the first 10 principle components (PCs).
#' Other algorithms or distance metrics are available in the PISCES package and should be explored for other datasets.
#' Finally, we visualize our clustering results using both UMAP and MDS.
## create a distance matrix
narnea.pca <- fast_pca(narnea.obj$PES, num.pcs = 10)
narnea.dist <- dist(narnea.pca$x)
## run silhouette optimized louvain clustering
narnea.clust <- louvain_k(narnea.dist, kmin = 5, kmax = 50)
## generate umap and mds
narnea.umap <- make_umap(narnea.obj$PES, dense = TRUE)
narnea.mds <- make_mds(narnea.dist)
## plot the clusters
plot.df <- data.frame('UMAP1' = narnea.umap[,1], 'UMAP2' = narnea.umap[,2],
                      'MDS1' = narnea.mds[,1], 'MDS2' = narnea.mds[,2],
                      'Cluster' = as.factor(narnea.clust$opt.clust))
ggplot(plot.df, aes(UMAP1, UMAP2)) + geom_point(aes(color = Cluster)) +
  ggtitle('PISCES Clustering - UMAP')
ggplot(plot.df, aes(MDS1, MDS2)) + geom_point(aes(color = Cluster)) +
  ggtitle('PISCES Clustering - MDS')

#' Step 6 - Master Regulator Analysis
#' From our protein activity signature, we identify candiate cluster master regulators using a Kruskal-Wallis test.
#' We can then visualize these results split into TFs - which likely drive the biological differences of these cell types - 
#' and markers - which present potential differentiating surface proteins that are not necessarily regulators.
#' We can visualize these results in a heatmap.
narnea.mrs <- kw_cluster_mrs(narnea.obj, narnea.clust$opt.clust)
## plot the top regulator results
cluster_mr_heatmap(narnea.obj$PES, dat.type = 'pact', clust.vec = narnea.clust$opt.clust, 
                   mr.list = narnea.mrs, reg.class = 'regulator', num.mrs = 5)
## plot the top marker results
cluster_mr_heatmap(narnea.obj$PES, dat.type = 'pact', clust.vec = narnea.clust$opt.clust, 
                   mr.list = narnea.mrs, reg.class = 'marker', num.mrs = 5)
## you can also try a different MR detection function such as Cohen's Kappa
## OPTIONAL - this may be a little slow
narnea.kappa.mrs <- kappa_cluster_mrs(narnea.obj, narnea.clust$opt.clust)

#' Addendum - CITE-Seq validation
#' These PBMC data were generated using CITE-Seq, which provide orthogonal antibody staining via sequencing.
#' We can demonstrate the similarity - and differences! - of our protein activity clustering to known biology
#' using a heatmap of these ADT genes.
adt.set <- rownames(adt.counts)
cluster_mr_heatmap(adt.counts, dat.type = 'gexp', clust.vec = narnea.clust$opt.clust, 
                   marker.set = adt.set.df, clust.rows = FALSE, scale.rows = TRUE, group.means = TRUE)

#' We can do a quick analysis of the correlation between our PES values and the ADT markers
#' and compare that to the correlation between GES values and ADT markers.
shared.genes <- intersect(rownames(adt.counts), intersect(rownames(narnea.obj$PES), rownames(ges.mat)))
ges.adt.cor <- lapply(shared.genes, function(x) {cor(adt.counts[x,], ges.mat[x,], method = 'spearman')})
pes.adt.cor <- lapply(shared.genes, function(x) {cor(adt.counts[x,], narnea.obj$PES[x,], method = 'spearman')})
