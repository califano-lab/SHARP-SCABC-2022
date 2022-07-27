## load packages
library(celldex)
library(cluster)
library(dplyr)
library(factoextra)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(Seurat)
library(SingleR)
## source functions
source('functions/hands-on-1_clustering-functions.R')
source('functions/hands-on-1_heatmap-function.R')
source('functions/hands-on-1_qc-functions.R')
## load data
pbmc_raw_counts <- readRDS('data/pbmc_raw-counts.rds')


#### This file is the UMI-count matrix with genes in the rows and cells in the columns
#### Check the initial number of cells with the following command:
ncol(pbmc_raw_counts)


###### QUALITY CONTROL FILTERING STEP ################################################
#### First thing to do is to have a sense of the quality of the data
#### To this purpose we look at 3 indicators:
### 1) Distribution of number of UMIs per cell (Sequencing Depth)
### 2) Distribution of number of detected genes per cell
### 3) Distribution of the proportion of UMIs per cell coming from the Mitochondrial genes

#### In the file called QC_Filtering_Functions.R there is the implementation of some functions 
#### to generate violin plots and to filter out the low-quality cells

#### By using the 'source' command you automatically load into your RStudio environment 
#### the functions implemented in the QC_Filtering_Functions.R the script, 
#### ready to be used for your analysis

#### MT gene names are save in an external .csv file that has to be loaded
mt.genes <- read.table('data/mt-genes.csv', header = TRUE, sep = ',', stringsAsFactors = FALSE)
hum.mt <- mt.genes$hum.symb

QCPlots(pbmc_raw_counts, hum.mt)

#### Manually compute the average number of UMIs per cell and 
#### the average number of detected genes per cell
mean(colSums(pbmc_raw_counts))
mean(apply(pbmc_raw_counts,2,function(x){sum(x!=0)}))


#### Use the functions 'MTFilter' and 'QCTransform' to remove from the dataset the low-quality cells.
#### As a rule of thumb, it is recommended to use 1000 as a lower threshold for UMIs/cell
### and 25% as upper threshold for MT genes
pbmc_MTFiltered_counts <- MTFilter(pbmc_raw_counts, hum.mt, mt.thresh = 0.25)
pbmc_QCFiltered_counts<-QCTransform(pbmc_MTFiltered_counts, minCount = 1000, minGeneReads = 1)

#### After QC filtering, cehck the number of cells, and again the average number of UMIs per cell and 
#### the average number of detected genes per cell
ncol(pbmc_QCFiltered_counts)
mean(colSums(pbmc_QCFiltered_counts))
mean(apply(pbmc_QCFiltered_counts,2,function(x){sum(x!=0)}))

#### Remove from the workspace useless variables to save memory and save the final QC filtered UMI-count matrix
rm(mt.genes, hum.mt, pbmc_raw_counts, pbmc_MTFiltered_counts)
saveRDS(pbmc_QCFiltered_counts, "pbmc_QCFiltered_counts.rds")

######################################################################################


###### DIMENSIONALITY REDUCTION ################################################
#### Now we want to plot the cells in two dimensions, and to this purpose we will  
#### use two dimensionlity reduction techniques: PCA and UMAP.

### We will use some functions from the Seurat pipeline.
### First we need to create a 'Seurat-object' which is a wrapper of our UMI-count matrix.
pbmc_QCFiltered_Seurat <- CreateSeuratObject(counts = pbmc_QCFiltered_counts, project = "PBMC", assay = "RNA")

### Now that the Seurat object is created we can use other functions of the Seurat pipeline
### for the downstream analysis and all the results will be included and added in the current Seurat-object

### We will first normalize the data and generate the gene expression signature using the 
### SCTransform function
pbmc_QCFiltered_Seurat <- SCTransform(pbmc_QCFiltered_Seurat, conserve.memory = T, verbose = FALSE)

#### Compute the Principal components to generate a 2D plot using the RunPCA function
pbmc_QCFiltered_Seurat <- RunPCA(pbmc_QCFiltered_Seurat, features = VariableFeatures(object = pbmc_QCFiltered_Seurat))
DimPlot(pbmc_QCFiltered_Seurat, reduction = 'pca', dims = c(1,2), 
        cols = 'black', pt.size = 1.5) + NoLegend()

#### Check the standard deviation explained by each principal component
ElbowPlot(pbmc_QCFiltered_Seurat)

#### Plot the cells using two higher principal components (#15 and #16 as an example) 
DimPlot(pbmc_QCFiltered_Seurat, reduction = 'pca', dims = c(15,16), 
        cols = 'black', pt.size = 1.5) + NoLegend()

#### Perform UMAP dimensionality reduction using the RunUMAP function
pbmc_QCFiltered_Seurat <- RunUMAP(pbmc_QCFiltered_Seurat, dims = 1:10, verbose = FALSE, 
                                  metric="correlation")
DimPlot(pbmc_QCFiltered_Seurat, reduction = 'umap', cols = 'black', pt.size = 1.5) + NoLegend()
##################################################################################################

###### Unsupevised Clustering Analysis ################################################

### Louvain Algorithm
### Use the function FindNeighbors, of the Seurat pipeline, to generate the 
### Nearest-Neighbor graph from the data. Then with function FindClusters run the Louvain algorithm 
### varying the resolution paramter from 0.1 to 1
pbmc_QCFiltered_Seurat <- FindNeighbors(pbmc_QCFiltered_Seurat, dims = 1:10)
pbmc_QCFiltered_Seurat <- FindClusters(pbmc_QCFiltered_Seurat, verbose = TRUE, algorithm=1, 
                                       resolution=seq(0.1,1,by=0.05))

#### Distance metric: Euclidean distance computed on the top-10 principal components
top10_PCs<-as.data.frame(pbmc_QCFiltered_Seurat$pca@cell.embeddings[,1:10])
euc_distance<-dist(top10_PCs, method="euclidean")

Louvain_SilScores<-Louvain_SilhouetteScore(pbmc_QCFiltered_Seurat, euc_distance)

#### Plot the average Silhouette Score values 
plot(seq(0.1,1,by=0.05), Louvain_SilScores, ylab="Slihouette Scores", xlab="Resolution Parameter", pch=18, 
     cex=1.5, main = "Louvain - Average Silhouette Scores", xaxt='n ', ylim = c(-1,1))
axis(1, at = seq(0.1,1,by=0.05), las=2)

#### For the optimal clustering solution, visualize the silhouette scores of each single-cell
fviz_silhouette(silhouette(as.numeric(pbmc_QCFiltered_Seurat$SCT_snn_res.0.1), euc_distance),
                palette=c("red", "blue", "green", "purple"))

#### save the Seurat object
saveRDS(pbmc_QCFiltered_Seurat, "pbmc_QCFiltered_Seurat.rds")

### Use ggplot to color the cells in the PCA and UMAP plots according
### to the clustering solution
UMAP_dataframe<-data.frame(pbmc_QCFiltered_Seurat@reductions[["umap"]]@cell.embeddings)
PCA_dataframe<-data.frame(pbmc_QCFiltered_Seurat@reductions[["pca"]]@cell.embeddings)
Clusters<-as.factor(as.numeric(pbmc_QCFiltered_Seurat$SCT_snn_res.0.1))

ggplot(PCA_dataframe,aes(x = PC_1, y = PC_2, color=Clusters)) + theme_bw()+
  geom_point(size=2)+ggtitle("PCA Plot - Louvain Clusters")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  scale_color_manual(values=c("red", "blue", "green", "purple"))

ggplot(UMAP_dataframe,aes(x = UMAP_1, y = UMAP_2, color=Clusters)) + theme_bw()+
  geom_point(size=2)+ggtitle("UMAP Plot - Louvain Clusters")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  scale_color_manual(values=c("red", "blue", "green", "purple"))

#######################################################################################################

##### Semi-supervised Cell Type Inference ######################

#### Run SingleR to infer cell identity
#### Reference paper of the BlueprintEncode https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3789449/
blueprint.encode <- BlueprintEncodeData()
SingleR_Results <- SingleR(test = pbmc_QCFiltered_Seurat@assays$SCT@counts, ref = blueprint.encode, 
                           labels = blueprint.encode$label.main, assay.type.test = 1)


### Pull-out from the SingleR object the top SingleR prediction for each single-cell
### and plot them in PCA and UMAP coordinates
SingleR_Labels <- SingleR_Results$first.labels
SingleR_Labels[which(SingleR_Labels %in% names(which(table(SingleR_Labels)<10)))] <- "Others"

ggplot(PCA_dataframe,aes(x = PC_1, y = PC_2, color=SingleR_Labels)) + theme_bw()+
  geom_point(size=2)+ggtitle("PCA Plot - SingleR Results")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  scale_color_manual(values=c("pink", "red", "blue", "orange", "green", "grey"))

ggplot(UMAP_dataframe,aes(x = UMAP_1, y = UMAP_2, color=SingleR_Labels)) + theme_bw()+
  geom_point(size=2)+ggtitle("UMAP Plot - SingleR Results")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  scale_color_manual(values=c("pink", "red", "blue", "orange", "green", "grey"))

###########################


######## Differential Gene Expression Analysis #################

#### Use the FindAllMarkers function of the Seurat pipeline to 
#### identify the most differnetially expressed genes 
Idents(pbmc_QCFiltered_Seurat)<-"SCT_snn_res.0.1"
pbmc_markers <- FindAllMarkers(pbmc_QCFiltered_Seurat, test.use = "wilcox", 
                               only.pos = TRUE, min.pct = 0.25)
saveRDS(pbmc_markers, "pbmc_markers.rds")

#### Select the top 10 genes per cluster 
top10 <- pbmc_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

#### Plot the expression of the top10 genes in a heatmap.
### In the script Heatmap_Function.R there is the implementation of a function
### to generate and save in a .pdf file the figure of the heatmap

### Prepare the expression matrix to plot the genes
genes<-top10$gene
expression_to_plot=pbmc_QCFiltered_Seurat[["SCT"]]@data[genes,]
Clusters<-as.factor(as.numeric(pbmc_QCFiltered_Seurat$SCT_snn_res.0.1))
colors<-c("red", "blue", "green", "purple")

gExprHeatmap(genes, expression_to_plot, Clusters, colors)




######## Supplemental Section: PAM Algorithm   ###################
#### PAM Algorithm
#### We will run the PAM algorithm varying the parameter k from 2 to 20
#### in order to determine the optimal clustering solution based on the Silhouette score

#### Distance metric: Euclidean distance on the top-10 principal components
top10_PCs<-as.data.frame(pbmc_QCFiltered_Seurat$pca@cell.embeddings[,1:10])
euc_distance<-dist(top10_PCs, method="euclidean")

#### Load the script 'Clustering_Functions.R' where there is the implementation
#### of the functions for the clustering analysis and for the computation of the 
#### Silhouette score
pam_clusters<-pam_k(euc_distance, 2, 20)

#### Plot the average values of the Silhouette scores
plot(2:20, pam_clusters$sil.scores, ylab="Slihouette Scores", xlab="k", pch=18, 
     cex=1.5, main = "PAM - Average Silhouette Scores")
axis(1, at = seq(2, 20, by = 1))

### The optimal clustering solution is for k=2
### Use ggplot to color the cells in the PCA and UMAP plots according
### to the clustering solution
UMAP_dataframe<-data.frame(pbmc_QCFiltered_Seurat@reductions[["umap"]]@cell.embeddings)
PCA_dataframe<-data.frame(pbmc_QCFiltered_Seurat@reductions[["pca"]]@cell.embeddings)
Clusters<-as.factor(pam_clusters$opt.clust)

ggplot(PCA_dataframe,aes(x = PC_1, y = PC_2, color=Clusters)) + theme_bw()+
  geom_point(size=2)+ggtitle("PCA Plot - PAM Clusters")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  scale_color_manual(values=c("red", "blue"))

ggplot(UMAP_dataframe,aes(x = UMAP_1, y = UMAP_2, color=Clusters)) + theme_bw()+
  geom_point(size=2)+ggtitle("UMAP Plot - PAM Clusters")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  scale_color_manual(values=c("red", "blue"))


##################################################################