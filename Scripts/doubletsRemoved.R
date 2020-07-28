##======================================================================================
## Removing doublets from pWTdata and pWTdata.psub
## Making separate data sets NDpWTdata and NDpWTdata.psub
##======================================================================================
NDpWTdata <- merge(WT1.obj, y = c(WT2.obj, WT3.obj), add.cell.ids = c("R1", "R2", "R3"), project = "plantRNA")
nonPgene <- row.names(NDpWTdata)[which(row.names(NDpWTdata) %in% gene.uncur == F)]
# Also removing mitochondrial genes
nonPgene <- nonPgene[!(nonPgene %in% as.character(genes$V2[27207:27416]))]
NDpWTdata <- subset(NDpWTdata, features = nonPgene, cells = Cells(pWTdata)[!(Cells(pWTdata) %in% doublets)])

NDpWTdata[["percent.mt"]] <- PercentageFeatureSet(NDpWTdata, pattern = "^MT-")
VlnPlot(NDpWTdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
NDpWTdata <- subset(NDpWTdata, subset = nFeature_RNA > 1000 & nFeature_RNA < 10500 & percent.mt < 0.01 & nCount_RNA < 200000)
library(sctransform)
NDpWTdata <- SCTransform(NDpWTdata, verbose = TRUE)
NDpWTdata <- RunPCA(NDpWTdata, verbose = TRUE)
NDpWTdata <- RunTSNE(NDpWTdata, verbose = TRUE)
library(reticulate)
NDpWTdata <- RunUMAP(NDpWTdata, dims = 1:30, verbose = TRUE)
NDpWTdata <- FindNeighbors(NDpWTdata, dims = 1:30, verbose = TRUE)
NDpWTdata <- FindClusters(NDpWTdata, verbose = TRUE, resolution = 1)   #clustering with high resolution and will regroup manually if necessary

DimPlot(NDpWTdata, reduction = "umap")
NDpWTdata.markers <- FindAllMarkers(NDpWTdata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "bimod")
NDpWTdata.markers <- NDpWTdata.markers[NDpWTdata.markers$pct.1 >= 0.1 & NDpWTdata.markers$pct.2 <= 0.1,]

png("./NDcellClust.png", units = "in", width = 10, height = 8, res = 600)
DimPlot(NDpWTdata, reduction = "umap", label = TRUE, pt.size = 0.5, label.size =4.5)
dev.off()

## Cluster ID's Found...
# Cluster0 = meristem
# Cluster1 = non hair cells
# Cluster2 = root cap cells
# Cluster3 = stele
# Cluster4 = phloem
# Cluster5 = hair cells
# Cluster6 = stele
# Cluster7 = endodermis
# Cluster8 = meristem
# Cluster9 = hair cells
# Cluster10 = cortex
# Cluster11 = hair cells
# Cluster12 = hair cells
# Cluster13 = non hair cells
# Cluster14 = non hair Cells
# Cluster15 = xylem
# Cluster16 = cortex
# Cluster17 = root cap cells
# Cluster18 = stele
# Cluster19 = protophloem
# Cluster20 = root cap cells
# Cluster21 = meristem
# Cluster22 = non hair cells

newClusterIDs <- c("Meristem", "Non Hair Cells", "Root Cap Cells", "Stele", "Phloem",
    "Hair Cells", "Stele", "Endodermis", "Meristem", "Hair Cells", "Cortex", "Hair Cells",
    "Hair Cells", "Non Hair Cells", "Non Hair Cells", "Xylem", "Cortex",
    "Root Cap Cells", "Stele", "Protophloem", "Root Cap Cells", "Meristem", "Non Hair Cells")
names(newClusterIDs) <- levels(NDpWTdata)
NDpWTdata <- RenameIdents(NDpWTdata, newClusterIDs)

png("./NDcellLabeled.png", units = "in", width = 10, height = 8, res = 600)
DimPlot(NDpWTdata, reduction = "umap", label = TRUE, pt.size = 0.5, label.size =4.5)
dev.off()

png("./NDcellUnlabeled.png", units = "in", width = 10, height = 8, res = 600)
DimPlot(NDpWTdata, reduction = "umap", pt.size = 0.5)
dev.off()



##=========================================================================
## Work with plastid clustering
##=========================================================================
NDpWTdata.psub <- merge(WT1.obj, y = c(WT2.obj, WT3.obj), add.cell.ids = c("R1", "R2", "R3"), project = "plantRNA")
NDpWTdata.psub <- subset(NDpWTdata.psub, features= gene.uncur, cells = Cells(pWTdata)[!(Cells(pWTdata) %in% doublets)])

VlnPlot(NDpWTdata.psub, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
NDpWTdata.psub <- subset(NDpWTdata.psub, subset = nFeature_RNA < 700 & nCount_RNA < 5500)
NDpWTdata.psub <- SCTransform(NDpWTdata.psub, verbose = TRUE)
NDpWTdata.psub <- RunPCA(NDpWTdata.psub, verbose = TRUE)
NDpWTdata.psub <- RunTSNE(NDpWTdata.psub, verbose = TRUE)
ElbowPlot(NDpWTdata.psub, ndims = 50)
NDpWTdata.psub <- RunUMAP(NDpWTdata.psub, dims = 1:30, verbose = TRUE)
NDpWTdata.psub <- FindNeighbors(NDpWTdata.psub, dims = 1:30, verbose = TRUE)
NDpWTdata.psub <- FindClusters(NDpWTdata.psub, verbose = TRUE, resolution = 1)

DimPlot(NDpWTdata.psub, reduction = "umap", label = TRUE, pt.size = 0.5, label.size =4.5)
NDpWTdata.psub.markers <- FindAllMarkers(NDpWTdata.psub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "bimod")

pClustIDs <- c("2", "3", "4", "5", "1", "1", "1", "6", "7", "8", "9", "10", "11", "1",
                "12", "1", "13")
names(pClustIDs) <- levels(NDpWTdata.psub)
NDpWTdata.psub <- RenameIdents(NDpWTdata.psub, pClustIDs)

png("./NDplastClust.png", units = "in", width = 11, height = 8, res = 600)
DimPlot(NDpWTdata.psub, reduction = "umap", label = TRUE, label.size = 7, pt.size = 0.5)
dev.off()


# looking more into cluster 7
FindMarkers(NDpWTdata.psub, ident.1 = "7", ident.2 = c("4", "12")) -> sevenMark
row.names(sevenMark) -> seven
which(is.na(seven %>% convertAT)) -> x
seven[x] <- substr(seven[x], 1, nchar(seven[x])-1)
convertAT(seven)




##=========================================================================
## Matching cell and plastid type clustering
## Working with nonnegative matrix factorization to match "clusterings"
##=========================================================================
x <- as.data.frame(Idents(NDpWTdata)); names(x) <- "cells"
x[,1] <- as.numeric(x[,1])
y <- as.data.frame(Idents(NDpWTdata.psub)); names(y) <- "plastids"
x$genes <- row.names(x); y$genes <- row.names(y)
z <- merge(x, y, by = 0)
for(i in c(2,4)) z[,i] <- as.numeric(as.character(z[,i]))
a.mat <- as.matrix(z[,c(2,4)])

# We then use NMF to extract pattern information:
# install.packages("NMF")
library(NMF)

fct <- nmf(a.mat+1, rank=2, .options='v')
# get factor matrics w and h where a = w %*% h
w <-  basis(fct)
h <- coef(fct)

# Reconstruct based on plastid types in column 1 of w and row 1 of h
x <- round(w[order(w[,1]),] %*% h[,order(h[,1])])

x <- data.frame(ids = names(Idents(NDpWTdata)), ct = as.character(Idents(NDpWTdata)))
x[x$ct == "Protophloem", "ct"] <- "Stele"
x[x$ct == "Phloem", "ct"] <- "Stele"
x[x$ct == "Xylem", "ct"] <- "Stele"
y <- data.frame(ids = names(Idents(NDpWTdata.psub)), pt = as.character(Idents(NDpWTdata.psub)))

a.mat <- merge(x, y)

## Identify high frequency matches:
matchFreq <- as.data.frame(table(paste(a.mat[,2], a.mat[,3], sep='-')))
names(matchFreq) <- c("clusterMatch", "freq")

matchFreq$clusterMatch <- as.character(matchFreq$clusterMatch)
matchFreq$Cell <- ''
matchFreq$Plastid <- ''
for(i in 1:nrow(matchFreq)) {
  matchFreq$Cell[i] <- strsplit(matchFreq$clusterMatch, '-')[[i]][1]
  matchFreq$Plastid[i] <- strsplit(matchFreq$clusterMatch, '-')[[i]][2]
}

cellTypes <- unique(matchFreq$Cell)

plastidTypes <- unique(matchFreq$Plastid)
plastidTypes <- plastidTypes[c(1, 4:5, 11, 6:7, 12, 8, 13, 2:3, 9:10)]
plastidTypes <- as.character(plastidTypes)

## When matching plastid classes to cell type, must account for size of plastid class
matchFreqtoc <- matchFreq
pclustFreq <- as.numeric(table(NDpWTdata.psub$seurat_clusters))
pclustFreq <- c(sum(pclustFreq[c(5, 6, 7, 14, 16)]), pclustFreq[c(1:4, 8:13, 15, 17)])
for(i in 1:13){
  matchFreqtoc$freq[matchFreqtoc$Plastid == plastidTypes[i]] <- (matchFreqtoc$freq[matchFreqtoc$Plastid == plastidTypes[i]])/pclustFreq[i]
}

png("./NDmatchTOC.png", units = "in", height = 9.5, width = 8, res = 600)
par(mfrow=c(3,3))
for(i in cellTypes){
  barplot((matchFreqtoc[matchFreqtoc$Cell==i,"freq"])*100,
          names.arg=matchFreqtoc[matchFreqtoc$Cell==i,"Plastid"], main=i,
          ylim = c(0,100))
}
dev.off()

##==================================================================
## Making the figures for cluster matching
##==================================================================
png("./matchCortex.png", units = "in", width = 8, height = 6, res = 600)
DimPlot(NDpWTdata.psub, cells.highlight = names(NDpWTdata@active.ident[NDpWTdata@active.ident == "Cortex"]), cols.highlight = "red")
dev.off()
png("./matchEndo.png", units = "in", width = 8, height = 6, res = 600)
DimPlot(NDpWTdata.psub, cells.highlight = names(NDpWTdata@active.ident[NDpWTdata@active.ident == "Endodermis"]), cols.highlight = "red")
dev.off()
png("./matchHair.png", units = "in", width = 8, height = 6, res = 600)
DimPlot(NDpWTdata.psub, cells.highlight = names(NDpWTdata@active.ident[NDpWTdata@active.ident == "Hair Cells"]), cols.highlight = "red")
dev.off()
png("./matchMeristem.png", units = "in", width = 8, height = 6, res = 600)
DimPlot(NDpWTdata.psub, cells.highlight = names(NDpWTdata@active.ident[NDpWTdata@active.ident == "Meristem"]), cols.highlight = "red")
dev.off()
png("./matchNHair.png", units = "in", width = 8, height = 6, res = 600)
DimPlot(NDpWTdata.psub, cells.highlight = names(NDpWTdata@active.ident[NDpWTdata@active.ident == "Non Hair Cells"]), cols.highlight = "red")
dev.off()
png("./matchRC.png", units = "in", width = 8, height = 6, res = 600)
DimPlot(NDpWTdata.psub, cells.highlight = names(NDpWTdata@active.ident[NDpWTdata@active.ident == "Root Cap Cells"]), cols.highlight = "red")
dev.off()
png("./matchStele.png", units = "in", width = 8, height = 6, res = 600)
DimPlot(NDpWTdata.psub, cells.highlight = names(NDpWTdata@active.ident[NDpWTdata@active.ident == "Stele" | NDpWTdata@active.ident == "Phloem" | NDpWTdata@active.ident == "Xylem" | NDpWTdata@active.ident == "Protophloem"]), cols.highlight = "red")
dev.off()



##====================================================================================
## Use TensorFlow to classify plastid gene expression profiles into cell type classes
## Read/Visualize/Section data into training, validation and test sets
##====================================================================================
t(NDpWTdata.psub@assays$SCT@scale.data) -> plastGE
x <- data.frame(ids = names(Idents(NDpWTdata)), ct = as.character(Idents(NDpWTdata)))
x[x$ct == "Protophloem", "ct"] <- "Phloem"
plastGE[row.names(plastGE) %in% x[,1],] -> plastGE
x[x$ids %in% row.names(plastGE),] -> x
plastGE <- as.data.frame(plastGE)
plastGE$ids <- row.names(plastGE)

plastML <- merge(x, plastGE)
plastML$ids <- NULL
write.csv(plastML, file = "./NDplastML.csv")

library(data.table)
plastData <- fread('./NDplastML.csv', header = T, sep = ',')
plastData[sample(1:nrow(plastData), 10), 1:5]  # look at 10 random rows
plastData$V1 <- NULL

## data set ratios
trnRatio <- 0.7; valRatio <- 0.2; tstRatio <- 0.1; N <- nrow(plastData)
trnN <- round(trnRatio * N); valN <- round(valRatio * N); tstN = N - trnN - valN;

rows <- 1:N
trnRows <- sample(rows, size=trnN); rows <- rows[-trnRows] # values and vector indexes are the same
valRows <- sample(rows, size=valN); rows <- rows[!is.element(rows, valRows)] # values are not the indexes
tstRows <- rows # whatever remains is used as the test dataset --

write.csv(plastData[trnRows,], file = "./trnData.csv", row.names=F)
write.csv(plastData[valRows,], file = "./valData.csv", row.names=F)
write.csv(plastData[tstRows,], file = "./tstData.csv", row.names=F)






##======================================================================================
## Idea: aberrant plastid-cell types are intermediate types that are "mixed up"
## Run pseudotime --> aberrant-plastid cells should be in "gap"
##======================================================================================
library(Seurat)
library(dyno)
library(dplyr)
library(tidyverse)
nonPlast <- row.names(NDpWTdata)[!is.element(row.names(NDpWTdata), gene.uncur)]
# Also removing mitochondrial genes
nonPlast <- nonPlast[!(nonPlast %in% as.character(genes$V2[27207:27416]))]
meriRC <-  names(NDpWTdata@active.ident[NDpWTdata@active.ident == "Root Cap Cells" | NDpWTdata@active.ident == "Meristem"])
ptimeNoPlast.meriRC <- subset(NDpWTdata, features = nonPlast, cells = meriRC)

ptimeNoPlast.meriRC <- FindVariableFeatures(ptimeNoPlast.meriRC)
varGenes.meriRC <- VariableFeatures(ptimeNoPlast.meriRC)
allExp.meriRC <-  Matrix::t(ptimeNoPlast.meriRC@assays$SCT@scale.data)
allCount.meriRC <- Matrix::t(ptimeNoPlast.meriRC$SCT@counts)
varExp.meriRC <- allExp.meriRC[,colnames(allExp.meriRC) %in% varGenes.meriRC]
varCount.meriRC <-  varExp.meriRC[,colnames(varExp.meriRC) %in% varGenes.meriRC]

meriRCNoPlast.dyn <- wrap_expression(
  expression = varExp.meriRC,
  counts = varCount.meriRC
)
meriRCNoPlast.dyn <- add_grouping(meriRCNoPlast.dyn, ptimeNoPlast.meriRC@active.ident)
initials.meriRC <- names(ptimeNoPlast.meriRC@active.ident[ptimeNoPlast.meriRC@active.ident == "Meristem"])
meriRCNoPlast.dyn <- add_prior_information(meriRCNoPlast.dyn, start_id = initials.meriRC)

ptimemeriRC <- infer_trajectory(meriRCNoPlast.dyn, ti_slingshot(), verbose = TRUE)
ptimemeriRC.simp <- simplify_trajectory(ptimemeriRC)

png("./NDmeriRCTI.png", units = "in", width = 10, height = 9, res = 600)
plot_dimred(ptimemeriRC, label_milestones = TRUE, color_cells = "grouping", size_cells = 1.5, grouping = ptimeNoPlast.meriRC@active.ident) + ggtitle("Slingshot TI Method used with Mersitem and Root Cap Cells")
dev.off()

png("./NDmeriRCTISimp.png", units = "in", width = 10, height = 9, res = 600)
plot_dimred(ptimemeriRC.simp, label_milestones = TRUE, color_cells = "grouping", size_cells = 1.9, grouping = ptimeNoPlast.meriRC@active.ident) + ggtitle("Simplified Slingshot TI Method used with Meristem and Root Cap Cells")
dev.off()


# Trying Meristem with Endodermis
meriEndo <-  names(NDpWTdata@active.ident[NDpWTdata@active.ident == "Endodermis" | NDpWTdata@active.ident == "Meristem"])
ptimeNoPlast.meriEndo <- subset(NDpWTdata, features = nonPlast, cells = meriEndo)

ptimeNoPlast.meriEndo <- FindVariableFeatures(ptimeNoPlast.meriEndo)
varGenes.meriEndo <- VariableFeatures(ptimeNoPlast.meriEndo)
allExp.meriEndo <-  Matrix::t(ptimeNoPlast.meriEndo@assays$SCT@scale.data)
allCount.meriEndo <- Matrix::t(ptimeNoPlast.meriEndo$SCT@counts)
varExp.meriEndo <- allExp.meriEndo[,colnames(allExp.meriEndo) %in% varGenes.meriEndo]
varCount.meriEndo <-  varExp.meriEndo[,colnames(varExp.meriEndo) %in% varGenes.meriEndo]

meriEndoNoPlast.dyn <- wrap_expression(
  expression = varExp.meriEndo,
  counts = varCount.meriEndo
)
meriEndoNoPlast.dyn <- add_grouping(meriEndoNoPlast.dyn, ptimeNoPlast.meriEndo@active.ident)
initials.meriEndo <- names(ptimeNoPlast.meriEndo@active.ident[ptimeNoPlast.meriEndo@active.ident == "Meristem"])
meriEndoNoPlast.dyn <- add_prior_information(meriEndoNoPlast.dyn, start_id = initials.meriEndo)

ptimemeriEndo <- infer_trajectory(meriEndoNoPlast.dyn, ti_slingshot(), verbose = TRUE)
ptimemeriEndo.simp <- simplify_trajectory(ptimemeriEndo)

png("./NDmeriEndoTI.png", units = "in", width = 10, height = 9, res = 600)
plot_dimred(ptimemeriEndo, label_milestones = TRUE, color_cells = "grouping", size_cells = 1.5, grouping = ptimeNoPlast.meriEndo@active.ident) + ggtitle("Slingshot TI Method used with Mersitem and Endodermis")
dev.off()

png("./NDmeriEndoTISimp.png", units = "in", width = 10, height = 9, res = 600)
plot_dimred(ptimemeriEndo.simp, label_milestones = TRUE, color_cells = "grouping", size_cells = 1.9, grouping = ptimeNoPlast.meriEndo@active.ident) + ggtitle("Simplified Slingshot TI Method used with Meristem and Endodermis")
dev.off()





##==============================================================================
## Normalizing separate dataset before pseudotime
##==============================================================================
meriEndo.names <-  names(NDpWTdata@active.ident[NDpWTdata@active.ident == "Endodermis" | NDpWTdata@active.ident == "Meristem"])
meriEndo.obj <- merge(WT1.obj, y = c(WT2.obj, WT3.obj), add.cell.ids = c("R1", "R2", "R3"), project = "plantRNA")
nonPgene <- row.names(meriEndo.obj)[which(row.names(meriEndo.obj) %in% gene.uncur == F)]
nonPgene <- nonPgene[!(nonPgene %in% as.character(genes$V2[27207:27416]))]
meriEndo.obj <- subset(meriEndo.obj, features = nonPgene, cells = meriEndo.names)

meriEndo.obj[["percent.mt"]] <- PercentageFeatureSet(meriEndo.obj, pattern = "^MT-")
VlnPlot(meriEndo.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
meriEndo.obj <- subset(meriEndo.obj, subset = nFeature_RNA > 2000 & nFeature_RNA < 8000 & percent.mt < 0.0075 & nCount_RNA < 75000)
meriEndo.obj <- SCTransform(meriEndo.obj, verbose = TRUE)
meriEndo.obj <- RunPCA(meriEndo.obj, verbose = TRUE)
meriEndo.obj <- RunTSNE(meriEndo.obj, verbose = TRUE)
meriEndo.obj <- RunUMAP(meriEndo.obj, dims = 1:30, verbose = TRUE)
meriEndo.obj <- FindNeighbors(meriEndo.obj, dims = 1:30, verbose = TRUE)
meriEndo.obj <- FindClusters(meriEndo.obj, verbose = TRUE, resolution = 0.1)

newClusterIDs <- c("Meristem", "Endodermis")
names(newClusterIDs) <- levels(meriEndo.obj)
meriEndo.obj <- RenameIdents(meriEndo.obj, newClusterIDs)

png("meriEndoUMAP.png", units = "in", width = 10, height = 8, res = 600)
DimPlot(meriEndo.obj, reduction = "umap", label = TRUE, label.size = 7, pt.size = 0.5)
dev.off()

varGenes.x <- VariableFeatures(meriEndo.obj)
allExp.x <-  Matrix::t(meriEndo.obj@assays$SCT@scale.data)
allCount.x <- Matrix::t(meriEndo.obj$SCT@counts)
varExp.x <- allExp.x[,colnames(allExp.x) %in% varGenes.x]
varCount.x <-  varExp.x[,colnames(varExp.x) %in% varGenes.x]

xNoPlast.dyn <- wrap_expression(
    expression = varExp.x,
    counts = varCount.x
)

xNoPlast.dyn <- add_grouping(xNoPlast.dyn, meriEndo.obj@active.ident)
initials.x <- names(meriEndo.obj@active.ident[meriEndo.obj@active.ident == "Meristem"])
xNoPlast.dyn <- add_prior_information(xNoPlast.dyn, start_id = initials.x)
ptimemeriEndo.norm <- infer_trajectory(xNoPlast.dyn, ti_slingshot(), verbose = TRUE)
ptimemeriEndo.norm.simp <- simplify_trajectory(ptimemeriEndo.norm)
slingProg <- ptimemeriEndo.norm.simp$progressions
slingProg <- slingProg[order(slingProg$percentage),]

png("meriEndoNormSimp.png", units = "in", width = 10, height = 8, res = 600)
plot_dimred(ptimemeriEndo.norm.simp, label_milestones = TRUE, color_cells = "grouping", size_cells = 1.5, grouping = meriEndo.obj@active.ident) + ggtitle("Slingshot TI Method used with Re-Normalized Meristem and Endodermis")
dev.off()

# Trying scorpius ti method
meriEndoPtime.scorpius <- infer_trajectory(xNoPlast.dyn, ti_scorpius(), verbose = TRUE)
png("meriEndoNormSCORPIUS.png", units = "in", width = 10, height = 8, res = 600)
plot_dimred(meriEndoPtime.scorpius, label_milestones = TRUE, color_cells = "grouping", size_cells = 1.5, grouping = meriEndo.obj@active.ident) + ggtitle("SCORPIUS TI Method used with Re-Normalized Meristem and Endodermis")
dev.off()
scorpProg <- meriEndoPtime.scorpius$progressions
scorpProg <- scorpProg[order(scorpProg$percentage),]


# Trying embeddr ti method
meriEndoPtime.embeddr <- infer_trajectory(xNoPlast.dyn, ti_embeddr(), verbose = TRUE)
png("meriEndoNormEmbeddr.png", units = "in", width = 10, height = 8, res = 600)
plot_dimred(meriEndoPtime.embeddr, label_milestones = TRUE, color_cells = "grouping", size_cells = 1.5, grouping = meriEndo.obj@active.ident) + ggtitle("Embeddr TI Method used with Re-Normalized Meristem and Endodermis")
dev.off()
embProg <- meriEndoPtime.embeddr$progressions
embProg <- embProg[order(embProg$percentage),]


# Comparing scorpius to slingshot to embeddr
png("scorpVSslingVSembeddr.png", units = "in", width = 10, height = 8, res = 600)
plot(slingProg$percentage, col = "red", ylab = "Percentage Along Pseudotime Axis", xlab = "Cell Index", main = "Comparison of Pseudotime Progression with SCORPIUS versus Slingshot")
points(embProg$percentage, col = "navyblue")
points(scorpProg$percentage, col = "black")
legend(1050, 0.20, col = c("red", "navyblue", "black"), legend = c("Slingshot", "Embeddr", "SCORPIUS"), pch = c(16, 16, 16))
dev.off()
# CHOSE SCORPIUS AS TI INFERENCE METHOD --> most stable increasing along ptime



# Making meriEndo.plast to score proplastid genes in separate plastid object
meriEndoPlast <- merge(WT1.obj, y = c(WT2.obj, WT3.obj), add.cell.ids = c("R1", "R2", "R3"), project = "plantRNA")
meriEndoPlast <- subset(meriEndoPlast, features = gene.uncur, cells = meriEndo.names)
meriEndoPlast[["percent.mt"]] <- PercentageFeatureSet(meriEndoPlast, pattern = "^MT-")
meriEndoPlast <- SCTransform(meriEndoPlast, verbose = TRUE)
meriEndoPlast <- RunPCA(meriEndoPlast, verbose = TRUE)
meriEndoPlast <- RunTSNE(meriEndoPlast, verbose = TRUE)
meriEndoPlast <- RunUMAP(meriEndoPlast, dims = 1:30, verbose = TRUE)
meriEndoPlast <- FindNeighbors(meriEndoPlast, dims = 1:30, verbose = TRUE)
meriEndoPlast <- FindClusters(meriEndoPlast, verbose = TRUE, resolution = 1)

png("meriEndoPlastUMAP.png", units = "in", width = 10, height = 8, res = 600)
DimPlot(meriEndoPlast, reduction = "umap", label = TRUE, label.size = 7, pt.size = 0.5)
dev.off()

proplastMarkers <- c("LPA3","PAC")

meriEndoPlast <- CellCycleScoring(meriEndoPlast, s.features = proplastMarkers, g2m.features = amyloplastMarkers$gene, set.ident = TRUE)
names(meriEndoPlast@meta.data)[names(meriEndoPlast@meta.data) == "S.Score"] <- "proplastScore"
names(meriEndoPlast@meta.data)[names(meriEndoPlast@meta.data) == "G2M.Score"] <- "amyloplastScore"

proplastScore <- as.data.frame(meriEndoPlast@meta.data[,"proplastScore"])
row.names(proplastScore) <- row.names(meriEndoPlast@meta.data)
proplastScore$extra <- 0

scorpProg <- meriEndoPtime.scorpius$progressions
scorpProg <- scorpProg[order(scorpProg$percentage),]
row.names(scorpProg) <- scorpProg$cell_id

proplastScore <- proplastScore[row.names(proplastScore) %in% row.names(scorpProg),]
proplastScore$extra <- NULL
scorpProg <- scorpProg[row.names(scorpProg) %in% row.names(proplastScore),]
scorpProg$cell_id <- NULL; scorpProg$from <- NULL; scorpProg$to <- NULL

proplastScore$id <- row.names(proplastScore); row.names(proplastScore) <- NULL
scorpProg$id <- row.names(scorpProg); row.names(scorpProg) <- NULL
linkage <- merge(scorpProg, proplastScore, by = "id")
colnames(linkage) <- c("id", "Pseudotime", "ProplastidScore")

plot(linkage$Pseudotime, linkage$ProplastidScore)

# for some reason, proplastid expression is increasinjg along pseudotime....
# Ideas: maybe check where paseudotime is actually starting and make sure it
#        is at meristem. Also, find better marker genes as these could just be
#        trashy genes for proplastid expression. Maybe look at plastid active
#        differentiation genes.


## Using dynfeature to see cell-type features that change anywhere along the trajectory
overall_feature_importances <- dynfeature::calculate_overall_feature_importance(meriEndoPtime.scorpius, expression_source= t(meriEndo.obj@assays$SCT@scale.data))
features <- overall_feature_importances %>% top_n(40, importance) %>% pull(feature_id)
  #plot_dimred(meriEndoPtime.scorpius)
png("TIheatmap.png", units = "in", width = 16, height = 10, res = 700)
dynplot::plot_heatmap(meriEndoPtime.scorpius, expression_source= t(meriEndo.obj@assays$SCT@scale.data))  #red represents high expression of gene while blue represents low expression
dev.off()

# After looking at top 20 important genes...
#   meristem (beiginning) shows high expression of genes involved in differentiation and formation of cell parts
#   endodermis (end) shows high expression of genes chracteristic to endodermis, especially casparian strip
#   Next, see if I can do this same thing but based on plastid gene expression profiles

## Using dynfeature to see plastid-type features that change anywhere along the trajectory
pExp <- t(meriEndoPlast@assays$SCT@scale.data)
pExp <- pExp[row.names(pExp) %in% meriEndoPtime.scorpius$cell_ids,]
overall_pfeature_importances <- dynfeature::calculate_overall_feature_importance(meriEndoPtime.scorpius, expression_source= pExp)
pfeatures <- overall_feature_importances %>% top_n(40, importance) %>% pull(feature_id)
  #plot_dimred(meriEndoPtime.scorpius)
png("TIheatmapPlastid.png", units = "in", width = 16, height = 10, res = 700)
dynplot::plot_heatmap(meriEndoPtime.scorpius, expression_source= pExp)  #red represents high expression of gene while blue represents low expression
dev.off()
