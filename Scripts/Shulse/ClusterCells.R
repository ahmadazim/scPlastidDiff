## Functions

# Prepare a matrix for Seurat - remove genes' column and add rownames
seuratPrep <- function(d) {
  rownames(d) <- toupper(d[,1])
  d[,1] <- NULL
  as.matrix(d)
} 

normalizeData <- function(d, genes) {
  rownames(d) <- toupper(d[,1])
  d[,1] <- NULL
  
  genes$id <- substr(genes$Geneid, 1, 9)
  d <- d[is.element(row.names(d), genes$id),]
  genes <- genes[is.element(genes$id, row.names(d)),]
  d <- d[genes$id, ]
  
  x <- prop.table(as.matrix(d), 2)*1000000  ## reads per million, rpm
  1000*x/genes$Length  ## reads per kb of gene length, rpkb

}


## GEO Omnibus accession: GSE116614
## Read Data Shulse et al data files 
d1 <- read.delim("Data/Shulse/GSM3243802_star_gene_exon_tagged_at_2183cells_DDCS040517A.clean.dge.txt")
d2 <- read.delim("Data/Shulse/GSM3243803_star_gene_exon_tagged_at_1860cells_DDCS040517B.clean.dge.txt")
genes <- read.delim("Data/GSM3243804_LIB1_counts.txt", skip = 1)

## Get gene info:
## ==============
# Format of chromosomal based nomenclature
#   AT (Arabidopsis thaliana)
#   1,2,3,4,5 (chromosome number) or M for mitochondrial or C for chloroplast.
#   G (gene), other letters possible for repeats etc.)
#   12300 (five-digit code, numbered from top/north to bottom/south of chromosome)

x <- substr(genes$Geneid, 1, 9)
table(substr(x, 1, 2))  # all AT, ie, all are Arabidopsis thaliana.
table(substr(x, 3, 3)) # 5 chromosomes, 88 chloroplast and 122 mitochondrial genes

#    1    2    3    4    5    C    M 
# 7078 4245 5437 4128 6318   88  122 

table(substr(d$GENE, 3, 3)) # different numbers with data
#    1    2    3    4    5    C    M 
# 6139 3540 4612 3517 5277   91   41 


# Normalize data - use Seurat's functions instead!
#d.norm <- normalizeData(d, genes)  
#dim(d.norm) # [1] 21378  2183 - 21,378 genes and 2,183 cells

library(Seurat)
library(dplyr)
library(Matrix) # probably loaded by default
library(gdata) # # probably loaded by default
library(seriation)
d1.seurat <- seuratPrep(d1) # pre matrix for Seurat
d2.seurat <- seuratPrep(d2) 

common.genes <- intersect(row.names(d1.seurat), row.names(d2.seurat))
d1.seurat <- d1.seurat[is.element(row.names(d1.seurat), common.genes),]
d2.seurat <- d2.seurat[is.element(row.names(d2.seurat), common.genes),]
d.seurat <- cbind(d1.seurat, d2.seurat)

object.size(d.seurat)
# 358747936 bytes
# The following cleaning up is only precautionary 
counts_per_cell <- Matrix::colSums(d.seurat) # same as apply(d.seurat, 2, sum)
counts_per_gene <- Matrix::rowSums(d.seurat) # apply(d.seurat, 1, sum)
# count gene only if it has non-zero reads mapped.
genes_per_cell <- Matrix::colSums(d.seurat > 0) # same as apply(d.seurat, 2, function(x) length(x[x>0]))
                                         # Note that (d.seurat > 0) is just a matrix of TRUE & FALSE
                                        # that adds 1 when an element is > 0.
cells_per_gene <- Matrix::rowSums(d.seurat > 0) 
# Log10 is an easy way to look at large numbers - log10(1e3) = 3 and log10(1e4) = 4, etc.
hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')
plot(counts_per_cell, genes_per_cell, log='xy', col='wheat') # both axes are logarithmic
title('counts vs genes per cell')
# Plot cells ranked by their number of detected genes
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)') # y axis is logarithmic

# Create Seurat Object
# include only genes expressed in 3 or more cells and cells with complexity of 350 genes or more
d.obj <- CreateSeuratObject(counts = d.seurat, project = "plantRNA", min.cells = 3, min.features = 350)
d.obj

# The previous line filtered out 2,728 genes expressed in less than 3 cells, see below
raw.data.in.object <- GetAssayData(d.obj, slots = "counts")
dim(raw.data.in.object) # 20489  2183 which is 2728 rows below dim(d.seurat)
length(cells_per_gene[cells_per_gene < 3]) # 2728

# Seurat object contains “slots” (designated by object@slotname) that store not only the raw count data, but also 
# the results from various computations below. This has the advantage that we do not need to keep track of inidividual 
# variables of interest - they can all be collapsed into a single object as long as these slots are pre-defined.
# d.obj@raw.data is a slot that stores the original gene count matrix. We can view the first 10 rows (genes) and 
# the first 10 columns (cells) by

dim(GetAssayData(object = d.obj, slot = "counts"))  # get raw data

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# Add a feature for percent mitochondria and plastids of each cell
d.obj[["percent.mt"]] <- PercentageFeatureSet(d.obj, pattern = "M")
d.obj[["percent.pl"]] <- PercentageFeatureSet(d.obj, pattern = "C")

head(d.obj[["percent.mt"]]); summary(d.obj[["percent.mt"]])
head(d.obj[["percent.pl"]]); summary(d.obj[["percent.pl"]])
head(d.obj[["nCount_RNA"]]); summary(d.obj[["nCount_RNA"]])
head(d.obj[["nFeature_RNA"]]); summary(d.obj[["nFeature_RNA"]])

# Visualize QC metrics as a violin plot
VlnPlot(d.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.pl"), ncol = 2)

# Can also visualize genes as features - 
VlnPlot(d.obj, features = c("AT1G01090", "AT1G01100", "AT1G01050"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(d.obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(d.obj, feature1 = "nCount_RNA", feature2 = "percent.pl")
plot3 <- FeatureScatter(d.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2, plot3), ncol=1)

# Per vln and scatter plots above, filter out cells that have 
#     unique feature counts > 8,500
#     percent mt > 5
#     percent plastids > 5
d.obj <- subset(d.obj, subset = nFeature_RNA < 8500 & percent.mt < 5 & percent.pl < 5)

# Normalize data. Normalized values are stored 
#     in: d.obj[["RNA"]]@data  
#     or: GetAssayData(d.obj, slot = "counts)
d.obj <- NormalizeData(d.obj, normalization.method = "LogNormalize", scale.factor = 10000)
# Can just use NormalizeData(d.obj) as all features above are the defaults
# normalization.method = "CR" and scale.factor = 1e6 gives CPM (counts per million) w/t log

# Identify highly variable genes (feature selection)
# https://www.nature.com/articles/nmeth.2645 (Nature Methods vol 10, 1093–95. 2013)
# Focusing on highly variable genes in downstream analysis helps to highlight biological signal 
# in single-cell datasets
d.obj <- FindVariableFeatures(d.obj, selection.method = "vst", nfeatures = 4000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(d.obj), 10)

# plot variable features with and without labels
# To plot, run plot1 (with no labels) or plot2 (with labels for to 10 genes)
plot1 <- VariableFeaturePlot(d.obj) 
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
CombinePlots(plots = list(plot1, plot2))

# Scaling the data
# ScaleData function is used to shift gene exp mean to 0 and gene exp var to 1.
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# The results are stored
    # in: d.obj[["RNA"]]@scale.data
    # or: GetAssayData(d.obj, slot = "scale.data")
all.genes <- rownames(d.obj)
d.obj <- ScaleData(d.obj, features = all.genes)
# Look at function SCTransform() - at Satija website -- 

# May also scale using variable genes if genes are in hundreds of thousands mabye?
# d.obj <- ScaleData(d.obj, features = VariableFeatures(d.obj))

# Perform linear dimensional reduction (both PCA and tSNE)
d.obj <- RunPCA(d.obj, features = VariableFeatures(object = d.obj))
# d.obj <- RunTSNE(d.obj, features = VariableFeatures(object = d.obj))

# Seurat provides several useful ways of visualizing both cells and features that define the PCA, 
# including VizDimReduction,  DimPlot, and DimHeatmap

# # Examine and visualize PCA results a few different ways
print(d.obj[["pca"]], dims = 1:5, nfeatures = 7)

# Viz to 20 genes of the first 2 PC's
VizDimLoadings(d.obj, dims = 1:2, reduction = "pca", nfeatures=20)
DimPlot(d.obj, reduction = "pca")

# Both cells and features are ordered according to their PCA scores. 
# Setting cells to a number plots the ‘extreme’ cells on both ends of the spectrum, 
# which dramatically speeds plotting for large datasets. 
# Though clearly a supervised analysis, we find this to be a valuable tool for exploring correlated feature sets.
DimHeatmap(d.obj, dims = 1, cells = 300, balanced = TRUE)
# DimPlot(d.obj, reduction = "tsne") ## try tsne

# Check the first 15 PC's
DimHeatmap(d.obj, dims = 1:3, cells = 300, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
# To overcome the extensive technical noise in any single feature for scRNA-seq data, 
# Seurat clusters cells based on their PCA scores, with each PC essentially representing 
# a ‘metafeature’ that combines information across a correlated feature set. 
# The top principal components therefore represent a robust compression of the dataset. 
# However, how many componenets should we choose to include? 10? 20? 100?
  
# Decide how many components to use. Define the significant ones.
d.obj <- JackStraw(d.obj, prop.freq = 0.01, num.replicate = 500, dims=20)
d.obj <- ScoreJackStraw(d.obj, dims = 1:20)
JackStrawPlot(d.obj, dims = 1:20)  # looks all first 20 PC's are highly significant

# Elbow plot is approximate and faster for large data sets
ElbowPlot(d.obj, ndims = 50)  # usually shows a difference b/n highly important and not so important PC's

## The advice here is to err on the higher side when choosing the number
## of PC's to use in downstream analysis

# We encourage users to repeat downstream analyses with a different number of PCs (10, 15, or even 50!).
# As you will observe, the results often do not differ dramatically.

# Clustering - 
# First construct a KNN graph based on the euclidean distance in PCA space, and 
# refine the edge weights between any two cells based on the shared overlap in their 
# local neighborhoods (Jaccard similarity). This step is performed using the 
# FindNeighbors function, and takes as input the previously defined dimensionality of the 
# dataset (first 20 PCs).
d.obj <- FindNeighbors(d.obj, dims = 1:20)

# To cluster the cells, we next apply modularity optimization techniques such as
# the Louvain algorithm (default) or SLM [SLM, Blondel et al., Journal of Statistical Mechanics],
# to iteratively group cells together, with the goal of optimizing the standard modularity function. 
# The FindClusters function implements this procedure, and contains a resolution parameter 
# that sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater 
# number of clusters. We find that setting this parameter between 0.4 - 1.2 typically returns good 
# results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets. 
d.obj <- FindClusters(d.obj, resolution = 0.6) # resolution set to 0.6 as in Shulse's paper

# The clusters can be found using the Idents function.
head(Idents(d.obj), 5)  # # Look at cluster IDs of the first 5 cells 
## =================================================
## Run non-linear dimensional reduction (UMAP/tSNE)
## =================================================

# The following creates a python environment inside rstudio
    conda_create("r-reticulate") # this is the environment
    reticulate::py_install(packages = 'umap-learn') # install python package umap-learn
      
# As input to the UMAP and tSNE, we suggest using the same PCs as input to the clustering analysis.
d.obj <- RunUMAP(d.obj, dims = 1:20)

# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(d.obj, reduction = "umap", label=TRUE)

# You can save the object at this point so that it can easily be loaded back in without having to rerun the 
# computationally intensive steps performed above, or easily shared with collaborators.
saveRDS(d.obj, file = "./output/d.obj_20dims.rds")  # read with readRDS()!
d.obj <- readRDS("./output/d.obj_20dims.rds")

## ==========================================
## Finding differentially expressed features
##         -- Cluster Biomarkers --
## ==========================================
# By default, Seurat identifes positive and negative markers of a single cluster 
# (specified in ident.1), compared to all other cells. FindAllMarkers automates 
# this process for all clusters, but you can also test groups of clusters vs. each other, 
# or against all cells.

# Find markers of single clusters -- place in 'all.markers' data frame
# Keep only markers with p values < 0.01 
all.markers <- FindAllMarkers(d.obj, min.pct = 0.25)
all.markers <- all.markers[all.markers$p_val_adj < 0.01,]

# Print cluster stats:
# ====================
cat("Number of significant genes in cluster\n")
table(all.markers$cluster)

cat("Maximum p-value of top 100 genes in cluster\n")
tapply(all.markers$p_val_adj, all.markers$cluster, function(x) max(x[1:100]))

# Let's compare our top markers to Shulse's top 50 markers
shulse50 <- read.csv("./Data/Shulse/Top50.csv") 
i <- 0 # Endodermis I (38 Endodermis I and 19 Endodermis II)
i <- 1 # Non-hair cells III (49 Non-hair cells III, 17 Columella)
i <- 2 # Columella 
i <- 3 # Cortex I (49 Cortex I, 16 Cortex II)
i <- 4 # Hair cells I (50 Hair cells I, 17 Hair cells II)
i <- 5 # Cortex IV (50 Cortex IV, 13 Cortex II)
i <- 6 # Cortex II (56 Cortex II, 13 Cortex IV)
i <- 7 # Non-hair cells II 
i <- 8 # Hair cells II (23 Hair cells I, 50 Hair cells II)
i <- 9 # Cortex II (60 Cortex II, 13 Cortex I)
i <- 10 # Stele
i <- 11 # Endodermis I (40 Endodermis I, 17 Endodermis II)
i <- 12 # Non-hair cells I

ext <- shulse50[is.element(shulse50$Gene.ID, all.markers[all.markers$cluster == i, 'gene'][1:100] ),]
table(ext$Cell.Type)

# Average expression of all cells within a cluster?
cluster.averages <- AverageExpression(d.obj, use.scale=T)
head(cluster.averages[["RNA"]][, 1:5])
# For each cluster, get the average expression of the 4000 higly variable genes
cluster.averages.hv <- cluster.averages[["RNA"]][VariableFeatures(d.obj),]

library(seriation)
newOrder <- seriate(t(cluster.averages.hv)) 
get_order(newOrder) - 1

# Define an order of cluster identities
newLevels <- get_order(newOrder) - 1
# 5  6  1  8  9  3  2 11  0  4  7 10 12

# Rename identities with relative values to average gene expression
obj <- d.obj
levels(Idents(obj))[newLevels+1] <- 0:12
DimPlot(obj, reduction = "umap", label=TRUE, pt.size = 0.5) + NoLegend()
# note 0 -> 1, 4 -> 5, and 7 -> 8


# Find all markers of the group 7, 8, 12 and keep the significant ones
cluster7812.markers <- FindMarkers(d.obj, ident.1 = c(7,8,12), min.pct = 0.25)
head(cluster7812.markers)
dim(cluster7812.markers)  # [1] 1579    5
dim(cluster7812.markers[cluster7812.markers$p_val_adj < 0.01,])  # [1] 1297    5
cluster7812.markers <- cluster7812.markers[cluster7812.markers$p_val_adj < 0.01,]

par(mfrow=c(2,1))
hist(cluster7812.markers$p_val_adj, main = "Adjusted p-values", xlab="Clusters: [7, 8, 12]")
hist(log10(cluster7812.markers$p_val_adj), main = expression("Log"[10]*" of adjusted p-values"),xlab="Clusters: [7, 8, 12]") 

# Find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(d.obj, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
d.obj.markers <- FindAllMarkers(d.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# You need dplyr for the following command to group by clusters and print to 2 genes of each
d.obj.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)

# Get the 'classification power' for each gene in cluster 1 where 0 is random 
# and 1 is perfect 
cluster1.markers <- FindMarkers(d.obj, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster1.markers, n = 5) # top one is AT3G22600
tail(cluster1.markers, n = 5) # bottom one is AT3G23690

VlnPlot(d.obj, features = c("AT3G22600", "AT3G23690"))
# Note that AT3G22600 is highly expressed in both 0 and 11 which are likely the 
# same cluster of endodermis I and II -