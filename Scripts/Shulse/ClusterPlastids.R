
# Data --
d1 <- read.delim("Data/Shulse/GSM3243802_star_gene_exon_tagged_at_2183cells_DDCS040517A.clean.dge.txt")
d2 <- read.delim("Data/Shulse/GSM3243803_star_gene_exon_tagged_at_1860cells_DDCS040517B.clean.dge.txt")

## http://ppdb.tc.cornell.edu/dbsearch/subproteome.aspx
## Read plastid genes - first column of gene names 
pl <- unique(substr(as.character(read.csv("Data/plastidgenes.csv")$Accession), 1, 9))

## Get gene info:
table(substr(pl, 1, 2))  # 1733 Arabidopsis thaliana plastid genes
table(substr(pl, 3, 3)) 

# Use only nucleus genome, for now
pl <- pl[substr(pl, 3, 3) != 'C']
table(substr(pl, 3, 3)) 
#   1   2   3   4   5  
# 429 242 324 264 386  

dp1 <- d1[is.element(d1$GENE, pl),]
dim(dp1) 

dp2 <- d2[is.element(d2$GENE, pl),]
dim(dp2) # 1656 plastid genes found in data

table(substr(dp1$GENE, 3, 3)) # different numbers with data

table(substr(dp2$GENE, 3, 3)) # different numbers with data

library(Seurat)
library(dplyr)
library(Matrix) # probably loaded by default
library(gdata) # # probably loaded by default
library(seriation)

dp1.seurat <- seuratPrep(dp1) # pre matrix for Seurat
dp2.seurat <- seuratPrep(dp2) # pre matrix for Seurat
common.genes <- intersect(row.names(dp1.seurat), row.names(dp2.seurat))
dp1.seurat <- dp1.seurat[is.element(row.names(dp1.seurat), common.genes),]
dp2.seurat <- dp2.seurat[is.element(row.names(dp2.seurat), common.genes),]
dp.seurat <- cbind(dp1.seurat, dp2.seurat)

# Clean-up is only precautionary 
counts_per_cell <- Matrix::colSums(dp.seurat) # same as apply(d1.seurat, 2, sum)
counts_per_gene <- Matrix::rowSums(dp.seurat) # apply(d1.seurat, 1, sum)
# count gene only if it has non-zero reads mapped.
genes_per_cell <- Matrix::colSums(dp.seurat > 0) # same as apply(d1.seurat, 2, function(x) length(x[x>0]))
# Note that (dp.seurat > 0) is just a matrix of TRUE & FALSE
# that adds 1 when an element is > 0.
cells_per_gene <- Matrix::rowSums(dp.seurat > 0) 
# Log10 is an easy way to look at large numbers - log10(1e3) = 3 and log10(1e4) = 4, etc.
hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')
plot(counts_per_cell, genes_per_cell, log='xy', col='wheat') # both axes are logarithmic
title('counts vs genes per cell')
# Plot cells ranked by their number of detected genes
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)') # y axis is logarithmic

# Create Seurat Object
dp.obj <- CreateSeuratObject(counts = dp.seurat, project = "plastidRNA", min.cells = 3, min.features = 100)
dp.obj

# The previous line filtered out 70 genes expressed in less than 3 cells
dim(GetAssayData(dp.obj, slots = "counts")) # 1586 2183
dim(dp.seurat) - dim(GetAssayData(dp.obj, slots = "counts")) # 89 1463 : genes and cells filtered

# Add a feature for percent plastids of each cell
# dp.obj[["percent.pl"]] <- PercentageFeatureSet(dp.obj, pattern = "C")

#head(dp.obj[["percent.pl"]]); summary(dp.obj[["percent.pl"]])
head(dp.obj[["nCount_RNA"]]); summary(dp.obj[["nCount_RNA"]])
head(dp.obj[["nFeature_RNA"]]); summary(dp.obj[["nFeature_RNA"]])

# Visualize QC metrics as a violin plot
VlnPlot(dp.obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

plot2 <- FeatureScatter(dp.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2), ncol=1)

dp.obj <- subset(dp.obj, subset = nFeature_RNA < 500 & nCount_RNA < 4000)

# Normalize data. Normalized values are stored 
#     in: dp.obj[["RNA"]]@data  
#     or: GetAssayData(dp.obj, slot = "counts)
dp.obj <- NormalizeData(dp.obj, normalization.method = "LogNormalize", scale.factor = 10000)
# Can just use NormalizeData(dp.obj) as all features above are the defaults
# normalization.method = "CR" and scale.factor = 1e6 gives CPM (counts per million) w/t log

# Identify highly variable genes (feature selection)
# https://www.nature.com/articles/nmeth.2645 (Nature Methods vol 10, 1093–95. 2013)
# Focusing on highly variable genes in downstream analysis helps to highlight biological signal 
# in single-cell datasets
dp.obj <- FindVariableFeatures(dp.obj, selection.method = "vst", nfeatures = 500)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(dp.obj), 10)

# plot variable features with and without labels
# To plot, run plot1 (with no labels) or plot2 (with labels for to 10 genes)
plot1 <- VariableFeaturePlot(dp.obj) 
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
CombinePlots(plots = list(plot1, plot2))

# Scaling the data
# ScaleData function is used to shift gene exp mean to 0 and gene exp var to 1.
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# The results are stored
# in: dp.obj[["RNA"]]@scale.data
# or: GetAssayData(dp.obj, slot = "scale.data")
all.genes <- rownames(dp.obj)
dp.obj <- ScaleData(dp.obj, features = all.genes)

# May also scale using variable genes if genes are in hundreds of thousands mabye?
# d1.obj <- ScaleData(d1.obj, features = VariableFeatures(d1.obj))

# Perform linear dimensional reduction (both PCA and tSNE)
dp.obj <- RunPCA(dp.obj, features = VariableFeatures(object = dp.obj))
# d1.obj <- RunTSNE(d1.obj, features = VariableFeatures(object = d1.obj))

# Seurat provides several useful ways of visualizing both cells and features that define the PCA, 
# including VizDimReduction,  DimPlot, and DimHeatmap

# # Examine and visualize PCA results a few different ways
print(dp.obj[["pca"]], dims = 1:5, nfeatures = 7)

# Viz to 20 genes of the first 2 PC's
VizDimLoadings(dp.obj, dims = 1:2, reduction = "pca", nfeatures=20)
DimPlot(dp.obj, reduction = "pca")

# Both cells and features are ordered according to their PCA scores. 
# Setting cells to a number plots the ‘extreme’ cells on both ends of the spectrum, 
# which dramatically speeds plotting for large datasets. 
# Though clearly a supervised analysis, we find this to be a valuable tool for exploring correlated feature sets.
DimHeatmap(dp.obj, dims = 1, cells = 200, balanced = TRUE)

# Check the first 6 PC's
DimHeatmap(d1.obj, dims = 1:6, cells = 200, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
# To overcome the extensive technical noise in any single feature for scRNA-seq data, 
# Seurat clusters cells based on their PCA scores, with each PC essentially representing 
# a ‘metafeature’ that combines information across a correlated feature set. 
# The top principal components therefore represent a robust compression of the dataset. 
# However, how many componenets should we choose to include? 10? 20? 100?

# Decide how many components to use. Define the significant ones.
dp.obj <- JackStraw(dp.obj, prop.freq = 0.01, num.replicate = 500, dims=20)
dp.obj <- ScoreJackStraw(dp.obj, dims = 1:20)
JackStrawPlot(dp.obj, dims = 1:20)  # First 10 were all significant but let's go with 20!

# Elbow plot is approximate and faster for large data sets
ElbowPlot(dp.obj, ndims = 20)  # usually shows a difference b/n highly important and not so important PC's

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
dp.obj <- FindNeighbors(dp.obj, dims = 1:20)

# To cluster the cells, we next apply modularity optimization techniques such as
# the Louvain algorithm (default) or SLM [SLM, Blondel et al., Journal of Statistical Mechanics],
# to iteratively group cells together, with the goal of optimizing the standard modularity function. 
# The FindClusters function implements this procedure, and contains a resolution parameter 
# that sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater 
# number of clusters. We find that setting this parameter between 0.4 - 1.2 typically returns good 
# results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets. 
dp.obj <- FindClusters(dp.obj, resolution = .8) # higher resultion = more clusters

# The clusters can be found using the Idents function.
head(Idents(dp.obj), 5)  # # Look at cluster IDs of the first 5 cells 
## =================================================
## Run non-linear dimensional reduction (UMAP/tSNE)
## =================================================

# The following creates a python environment inside rstudio
conda_create("r-reticulate") # this is the environment
reticulate::py_install(packages = 'umap-learn') # install python package umap-learn

# As input to the UMAP and tSNE, we suggest using the same PCs as input to the clustering analysis.
dp.obj <- RunUMAP(dp.obj, dims = 1:20)

# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(dp.obj, reduction = "umap", label=TRUE)

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
all.pmarkers <- FindAllMarkers(dp.obj, min.pct = 0.25)
all.pmarkers <- all.pmarkers[all.pmarkers$p_val_adj < 0.01,]

# Print cluster stats:
# ====================
cat("Number of significant genes in cluster\n")
table(all.pmarkers$cluster)

cat("Maximum p-value of top 50 genes in cluster\n")
tapply(all.pmarkers$p_val_adj, all.pmarkers$cluster, function(x) max(x[1:20]))

# Average expression of all cells within a cluster?
cluster.averages <- AverageExpression(dp.obj, use.scale=T)
head(cluster.averages[["RNA"]][, 1:5])
# For each cluster, get the average expression of the 4000 higly variable genes
cluster.averages.hv <- cluster.averages[["RNA"]][VariableFeatures(dp.obj),]

library(seriation)
newOrder <- seriate(t(cluster.averages.hv)) 
get_order(newOrder) - 1

# Define an order of cluster identities
newLevels <- get_order(newOrder) - 1
# 6 1 3 7 2 0 5 8 4

# Rename identities so that clusters with similar expression get close numbers
pobj <- dp.obj
levels(Idents(pobj))[newLevels+1] <- 0:(length(newLevels)-1)
DimPlot(pobj, reduction = "umap", label=TRUE, pt.size = 0.5) + NoLegend()

# Get the 'classification power' 
cluster6.markers <- FindMarkers(dp.obj, ident.1 = 6, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster6.markers, n = 5) # AT1G79040
tail(cluster6.markers, n = 5) # AT3G10420

VlnPlot(dp.obj, features = c("AT1G79040", "AT3G10420"))
