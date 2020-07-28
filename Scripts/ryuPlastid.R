##=========================================================================
## Post RSI work with plastid clustering
##=========================================================================
pWTdata.psub <- merge(WT1.obj, y = c(WT2.obj, WT3.obj), add.cell.ids = c("R1", "R2", "R3"), project = "plantRNA")
pWTdata.psub <- subset(pWTdata.psub, features= gene.uncur)

VlnPlot(pWTdata.psub, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
pWTdata.psub <- subset(pWTdata.psub, subset = nFeature_RNA < 700 & nCount_RNA < 6000)
pWTdata.psub <- SCTransform(pWTdata.psub, verbose = TRUE)

pWTdata.psub <- RunPCA(pWTdata.psub, verbose = TRUE)
pWTdata.psub <- RunTSNE(pWTdata.psub, verbose = TRUE)

ElbowPlot(pWTdata.psub, ndims = 50)

library(reticulate)
pWTdata.psub <- RunUMAP(pWTdata.psub, dims = 1:30, verbose = TRUE)

pWTdata.psub <- FindNeighbors(pWTdata.psub, dims = 1:30, verbose = TRUE)
pWTdata.psub <- FindClusters(pWTdata.psub, verbose = TRUE, resolution = 1)
pWTdata.psub.markers <- FindAllMarkers(pWTdata.psub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "bimod")

pClustIDs <- c("2", "3", "4", "1", "5", "1", "1", "6", "7",
                "8", "9", "10", "11", "1", "12", "1", "1", "13")
names(pClustIDs) <- levels(pWTdata.psub)
pWTdata.psub <- RenameIdents(pWTdata.psub, pClustIDs)

png("./plastClustPostRSI.png", units = "in", width = 11, height = 8, res = 600)
DimPlot(pWTdata.psub, reduction = "umap", label = TRUE, label.size = 7, pt.size = 0.5)
dev.off()



##=========================================================================
## Working with nonnegative matrix factorization to match "clusterings"
##=========================================================================

x <- as.data.frame(Idents(pWTdata)); names(x) <- "cells"
x[,1] <- as.numeric(x[,1])
y <- as.data.frame(Idents(pWTdata.psub)); names(y) <- "plastids"
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
# This should be an âorderedâ reconstructed 'a' matrix that aligns plastid clusters
# nicely against cell types.

# library(gplots)
#
# # Run heat map on matched clusters without ordering
# png("./HMno", units = "in", width = 20, height = 10, res = 800)
# heatmap.2(a.mat,
#           density.info="none",
#           trace="none",
#           margins =c(12,9),
#           #col=my_palette,
#           #breaks=col_breaks,
#           dendrogram='none',
#           Rowv=FALSE,
#           Colv=FALSE)
# dev.off()
#
# # Run heat map on reconstructed matrix 'x'
# # -- results are diff from(heatmap(a)
# png("./HMre", units = "in", width = 20, height = 10, res = 800)
# heatmap.2(x,
#           density.info="none",
#           trace="none",
#           margins =c(12,9),
#           #col=my_palette,
#           #breaks=col_breaks,
#           dendrogram='none',
#           Rowv=FALSE,
#           Colv=FALSE)
# dev.off()


x <- data.frame(ids = names(Idents(pWTdata)), ct = as.character(Idents(pWTdata)))
x[x$ct == "Protoxylem", "ct"] <- "Stele"
x[x$ct == "Protophloem", "ct"] <- "Stele"
x[x$ct == "Phloem", "ct"] <- "Stele"
x[x$ct == "Xylem", "ct"] <- "Stele"
y <- data.frame(ids = names(Idents(pWTdata.psub)), pt = as.character(Idents(pWTdata.psub)))

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
plastidTypes <- plastidTypes[c(1, 4:9, 13, 10, 11, 2, 3, 12)]
plastidTypes <- as.character(plastidTypes)

## When matching plastid classes to cell type, must account for size of plastid class
matchFreqtoc <- matchFreq
pclustFreq <- as.numeric(table(pWTdata.psub$seurat_clusters))
pclustFreq <- c(sum(pclustFreq[c(4, 6, 7, 14, 16, 17)]), pclustFreq[c(1, 2, 3, 5, 8:13, 15, 18)])
for(i in 1:13){
  matchFreqtoc$freq[matchFreqtoc$Plastid == plastidTypes[i]] <- (matchFreqtoc$freq[matchFreqtoc$Plastid == plastidTypes[i]])/pclustFreq[i]
}

png("./matchPostRSItoc.png", units = "in", height = 12, width = 16, res = 600)
par(mfrow=c(3,3))
for(i in cellTypes){
  barplot((matchFreqtoc[matchFreqtoc$Cell==i,"freq"])*100,
          names.arg=matchFreqtoc[matchFreqtoc$Cell==i,"Plastid"], main=i,
          xlab = "Plastid Clusters", ylab = "Percent of Plastid Cluster",
          ylim = c(0,100))
}
dev.off()

## When matching cell types to plastid classes, must account for size of cell cluster
matchFreqtop <- matchFreq
table(pWTdata@active.ident)
cclustFreq <- data.frame("Cortex" = 550, "Endodermis" = 399, "Hair Cells" = 1550, "Meristem" = 1030,
                         "Non Hair Cells" = 1104, "Root Cap Cells" = 785,
                         "Stele" = (1113 + 474 + 307 + 154 + 46))
cclustFreq <- as.numeric(cclustFreq)
for(i in 1:7){
  matchFreqtop$freq[matchFreqtop$Cell == cellTypes[i]] <- (matchFreqtop$freq[matchFreqtop$Cell == cellTypes[i]])/cclustFreq[i]
}
png("./matchPostRSItop.png", units = "in", height = 24, width = 16, res = 600)
par(mfrow=c(6,3))
for(i in plastidTypes){
  barplot(matchFreq[matchFreq$Plastid==i,"freq"],
          names.arg=matchFreq[matchFreq$Plastid==i,"Cell"], main=i)
}
dev.off()



##==================================================================
## Making the figures for cluster matching
##==================================================================
png("./matchMeristem", units = "in", width = 8, height = 6, res = 600)
DimPlot(pWTdata.psub, cells.highlight = names(pWTdata@active.ident[pWTdata@active.ident == "Cortex"]), cols.highlight = "red")
DimPlot(pWTdata.psub, cells.highlight = names(pWTdata@active.ident[pWTdata@active.ident == "Endodermis"]), cols.highlight = "red")
DimPlot(pWTdata.psub, cells.highlight = names(pWTdata@active.ident[pWTdata@active.ident == "Hair Cells"]), cols.highlight = "red")
DimPlot(pWTdata.psub, cells.highlight = names(pWTdata@active.ident[pWTdata@active.ident == "Meristem"]), cols.highlight = "red")
DimPlot(pWTdata.psub, cells.highlight = names(pWTdata@active.ident[pWTdata@active.ident == "Non Hair Cells"]), cols.highlight = "red")
DimPlot(pWTdata.psub, cells.highlight = names(pWTdata@active.ident[pWTdata@active.ident == "Root Cap Cells"]), cols.highlight = "red")
DimPlot(pWTdata.psub, cells.highlight = names(pWTdata@active.ident[pWTdata@active.ident == "Stele" | pWTdata@active.ident == "Phloem" | pWTdata@active.ident == "Xylem" | pWTdata@active.ident == "Protoxylem" | pWTdata@active.ident == "Protophloem"]), cols.highlight = "red")



##=======================================================================================
## Looking into cells that did not cluster to same plastid cluster
## Will try to reconstruct matchFreq matrix to include UMI's and investigate those cells
##=======================================================================================
cType <- data.frame(ids = names(Idents(pWTdata)), ct = as.character(Idents(pWTdata)))
cType[cType$ct == "Protoxylem", "ct"] <- "Stele"
cType[cType$ct == "Protophloem", "ct"] <- "Stele"
cType[cType$ct == "Phloem", "ct"] <- "Stele"
cType[cType$ct == "Xylem", "ct"] <- "Stele"
pType <- data.frame(ids = names(Idents(pWTdata.psub)), pt = as.character(Idents(pWTdata.psub)))
cpType <- merge(cType, pType)

# Root Cap (pClust = 3)
rcHet <- cpType$ids[cpType$ct == "Root Cap Cells" & cpType$pt != 3]
rc <- names(pWTdata@active.ident[pWTdata@active.ident == "Root Cap Cells"])

rcData <- merge(WT1.obj, y = c(WT2.obj, WT3.obj), add.cell.ids = c("R1", "R2", "R3"), project = "plantRNA")
rcData <- subset(rcData, features= gene.uncur, cells = rc)

VlnPlot(rcData, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
rcData <- subset(rcData, subset = nFeature_RNA < 500 & nCount_RNA < 1500)
rcData <- SCTransform(rcData, verbose = TRUE)
rcData <- RunPCA(rcData, verbose = TRUE)
rcData <- RunTSNE(rcData, verbose = TRUE)
ElbowPlot(rcData, ndims = 50)
rcData <- RunUMAP(rcData, dims = 1:20, verbose = TRUE)
rcData <- FindNeighbors(rcData, dims = 1:20, verbose = TRUE)
rcData <- FindClusters(rcData, verbose = TRUE, resolution = 1)
rcData.markers <- FindAllMarkers(rcData, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "bimod")
DimPlot(rcData, cells.highlight = as.character(rcHet), cols.highlight = "red")

rcData@assays$SCT@scale.data[, colnames(rcData@assays$SCT@scale.data) == rcHet]
rcHetGE <- rcData@assays$SCT@scale.data[, colnames(rcData@assays$SCT@scale.data) %in% rcHet]
rcSameGE <- (rcData@assays$SCT@scale.data)[,!(is.element(colnames(rcData@assays$SCT@scale.data), as.character(rcHet)))]
