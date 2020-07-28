##=============================================================================
## Working with data from the following paper by Ryu:
## http://www.plantphysiol.org/content/179/4/1444
## GEO Accession number: Series "GSE123013"
##
## I downloaded the "GSE123013_RAW_matrices.tar.gz" supplementary file, which
## downloaded 5 samples, 3 of which were wild-type replicates. The other two
## were mutated samples.
##
## I un-".tar"-ed, un-".gz"-ed, etc. and imported the three folders (3 WT
## samples) into my GCP environment. To use the "Read10X()" function from
## Seurat the folder must have 3 files (named correctly) in a directory...
##=============================================================================

sort(sapply(ls(),function(x){object.size(get(x))}))

##=========================================================================
## Post RSI work with pWTdata
##=========================================================================
library(Seurat)
library(reticulate)
list.files("/home/ahmad/Ryu/WT1")
list.files("/home/ahmad/Ryu/WT2")   ##for these three lines, output must exactly be...
list.files("/home/ahmad/Ryu/WT3")   # [1] "barcodes.tsv" "genes.tsv"    "matrix.mtx"

WT1 <-  Read10X(data.dir = "/home/ahmad/Ryu/WT1")
WT2 <-  Read10X(data.dir = "/home/ahmad/Ryu/WT2")
WT3 <-  Read10X(data.dir = "/home/ahmad/Ryu/WT3")

genes <- read.csv("genes.csv", header = F)

#Creating the Seurat object
WT1.obj <- CreateSeuratObject(counts = WT1)
WT2.obj <- CreateSeuratObject(counts = WT2)
WT3.obj <- CreateSeuratObject(counts = WT3)

pWTdata <- merge(WT1.obj, y = c(WT2.obj, WT3.obj), add.cell.ids = c("R1", "R2", "R3"), project = "plantRNA")
nonPgene <- row.names(pWTdata)[which(row.names(pWTdata) %in% gene.uncur == F)]
pWTdata <- subset(pWTdata, features = nonPgene)

pWTdata[["percent.mt"]] <- PercentageFeatureSet(pWTdata, pattern = "^MT-")
VlnPlot(pWTdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pWTdata <- subset(pWTdata, subset = nFeature_RNA > 1000 & nFeature_RNA < 10500 & percent.mt < 0.015)

# Making sure no bad quality cells, namely cells that express less than 200 genes (>3 counts to be "expressed")
a <- 0
for(i in 1:7512){
  x <- pWTdata@assays$RNA@counts[,i]
  a <- c(a, length(x[x > 3]))
  cat(i, " : ", length(x[x > 3]), "\n")
}


library(sctransform)
pWTdata <- SCTransform(pWTdata, verbose = TRUE)

pWTdata <- RunPCA(pWTdata, verbose = TRUE)
pWTdata <- RunTSNE(pWTdata, verbose = TRUE)

library(reticulate)
pWTdata <- RunUMAP(pWTdata, dims = 1:30, verbose = TRUE)

pWTdata <- FindNeighbors(pWTdata, dims = 1:30, verbose = TRUE)
pWTdata <- FindClusters(pWTdata, verbose = TRUE, resolution = 1)   #clustering with high resolution and will regroup manually if necessary

#Vizulizing clusters
DimPlot(pWTdata, reduction = "pca")
DimPlot(pWTdata, reduction = "tsne")
DimPlot(pWTdata, reduction = "umap")

#Find markers for each cluster
pWTdata.markers <- FindAllMarkers(pWTdata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "bimod")
#Only keep markers that are expressed more than 0.1 in their cluster and less than 0.1 in other clusters
pWTdata.markers <- pWTdata.markers[pWTdata.markers$pct.1 >= 0.1 & pWTdata.markers$pct.2 <= 0.1,]

mark <- matrix(0, nrow = 300, ncol = 23)
colnames(mark) <- c("Cluster0", "Cluster1","Cluster2","Cluster3","Cluster4","Cluster5","Cluster6","Cluster7","Cluster8","Cluster9","Cluster10","Cluster11","Cluster12","Cluster13","Cluster14","Cluster15","Cluster16","Cluster17","Cluster18","Cluster19","Cluster20","Cluster21","Cluster22")
mark <- as.data.frame(mark)
for(i in 1:23){
    mark[,i] <- c(pWTdata.markers$gene[pWTdata.markers$cluster == (i-1)], rep(0, length = (300 - (length(pWTdata.markers$gene[pWTdata.markers$cluster == (i-1)])))))
}

#convert genes in mark to AT name
for(j in 1:23){
  for(i in 1:300){
    if(mark[i,j] == 0){mark[i,j] <- mark[i,j]}
    if(mark[i,j] != 0){
      x <- as.character(genes$V1[genes$V2 == mark[i,j]])
      mark[i,j] <- x[1]
    }
    cat("row", i, "of column", j, "done.", "\n")
  }
}
write.csv(mark,file="markersPostRSI.csv")


# Finding differentially expressed genes within each cluster
pWTdataMark0 <- FindMarkers(pWTdata, ident.1 = 0, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
pWTdataMark0 <- pWTdataMark0[pWTdataMark0$pct.1 >= 0.1 & pWTdataMark0$pct.2 <= 0.4,]
write.csv(pWTdataMark0, file = "./postRSImarkers/pWTdataMark0.csv")

pWTdataMark1 <- FindMarkers(pWTdata, ident.1 = 1, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
pWTdataMark1 <- pWTdataMark1[pWTdataMark1$pct.1 >= 0.1 & pWTdataMark1$pct.2 <= 0.4,]
write.csv(pWTdataMark1, file = "./postRSImarkers/pWTdataMark1.csv")

pWTdataMark2 <- FindMarkers(pWTdata, ident.1 = 2, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
pWTdataMark2 <- pWTdataMark2[pWTdataMark2$pct.1 >= 0.1 & pWTdataMark2$pct.2 <= 0.4,]
write.csv(pWTdataMark2, file = "./postRSImarkers/pWTdataMark2.csv")

pWTdataMark3 <- FindMarkers(pWTdata, ident.1 = 3, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
pWTdataMark3 <- pWTdataMark3[pWTdataMark3$pct.1 >= 0.1 & pWTdataMark3$pct.2 <= 0.4,]
write.csv(pWTdataMark3, file = "./postRSImarkers/pWTdataMark3.csv")

pWTdataMark4 <- FindMarkers(pWTdata, ident.1 = 4, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
pWTdataMark4 <- pWTdataMark4[pWTdataMark4$pct.1 >= 0.1 & pWTdataMark4$pct.2 <= 0.4,]
write.csv(pWTdataMark4, file = "./postRSImarkers/pWTdataMark4.csv")

pWTdataMark5 <- FindMarkers(pWTdata, ident.1 = 5, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
pWTdataMark5 <- pWTdataMark5[pWTdataMark5$pct.1 >= 0.1 & pWTdataMark5$pct.2 <= 0.4,]
write.csv(pWTdataMark5, file = "./postRSImarkers/pWTdataMark5.csv")

pWTdataMark6 <- FindMarkers(pWTdata, ident.1 = 6, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
pWTdataMark6 <- pWTdataMark6[pWTdataMark6$pct.1 >= 0.1 & pWTdataMark6$pct.2 <= 0.4,]
write.csv(pWTdataMark6, file = "./postRSImarkers/pWTdataMark6.csv")

pWTdataMark7 <- FindMarkers(pWTdata, ident.1 = 7, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
pWTdataMark7 <- pWTdataMark7[pWTdataMark7$pct.1 >= 0.1 & pWTdataMark7$pct.2 <= 0.4,]
write.csv(pWTdataMark7, file = "./postRSImarkers/pWTdataMark7.csv")

pWTdataMark8 <- FindMarkers(pWTdata, ident.1 = 8, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
pWTdataMark8 <- pWTdataMark8[pWTdataMark8$pct.1 >= 0.1 & pWTdataMark8$pct.2 <= 0.4,]
write.csv(pWTdataMark8, file = "./postRSImarkers/pWTdataMark8.csv")

pWTdataMark9 <- FindMarkers(pWTdata, ident.1 = 9, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
pWTdataMark9 <- pWTdataMark9[pWTdataMark9$pct.1 >= 0.1 & pWTdataMark9$pct.2 <= 0.4,]
write.csv(pWTdataMark9, file = "./postRSImarkers/pWTdataMark9.csv")

pWTdataMark10 <- FindMarkers(pWTdata, ident.1 = 10, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
pWTdataMark10 <- pWTdataMark10[pWTdataMark10$pct.1 >= 0.1 & pWTdataMark10$pct.2 <= 0.4,]
write.csv(pWTdataMark10, file = "./postRSImarkers/pWTdataMark10.csv")

pWTdataMark11 <- FindMarkers(pWTdata, ident.1 = 11, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
pWTdataMark11 <- pWTdataMark11[pWTdataMark11$pct.1 >= 0.1 & pWTdataMark11$pct.2 <= 0.4,]
write.csv(pWTdataMark11, file = "./postRSImarkers/pWTdataMark11.csv")

pWTdataMark12 <- FindMarkers(pWTdata, ident.1 = 12, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
pWTdataMark12 <- pWTdataMark12[pWTdataMark12$pct.1 >= 0.1 & pWTdataMark12$pct.2 <= 0.4,]
write.csv(pWTdataMark12, file = "./postRSImarkers/pWTdataMark12.csv")

pWTdataMark13 <- FindMarkers(pWTdata, ident.1 = 13, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
pWTdataMark13 <- pWTdataMark13[pWTdataMark13$pct.1 >= 0.1 & pWTdataMark13$pct.2 <= 0.4,]
write.csv(pWTdataMark13, file = "./postRSImarkers/pWTdataMark13.csv")

pWTdataMark14 <- FindMarkers(pWTdata, ident.1 = 14, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
pWTdataMark14 <- pWTdataMark14[pWTdataMark14$pct.1 >= 0.1 & pWTdataMark14$pct.2 <= 0.4,]
write.csv(pWTdataMark14, file = "./postRSImarkers/pWTdataMark14.csv")

pWTdataMark15 <- FindMarkers(pWTdata, ident.1 = 15, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
pWTdataMark15 <- pWTdataMark15[pWTdataMark15$pct.1 >= 0.1 & pWTdataMark15$pct.2 <= 0.4,]
write.csv(pWTdataMark15, file = "./postRSImarkers/pWTdataMark15.csv")

pWTdataMark16 <- FindMarkers(pWTdata, ident.1 = 16, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
pWTdataMark16 <- pWTdataMark16[pWTdataMark16$pct.1 >= 0.1 & pWTdataMark16$pct.2 <= 0.4,]
write.csv(pWTdataMark16, file = "./postRSImarkers/pWTdataMark16.csv")

pWTdataMark17 <- FindMarkers(pWTdata, ident.1 = 17, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
pWTdataMark17 <- pWTdataMark17[pWTdataMark17$pct.1 >= 0.1 & pWTdataMark17$pct.2 <= 0.4,]
write.csv(pWTdataMark17, file = "./postRSImarkers/pWTdataMark17.csv")

pWTdataMark18 <- FindMarkers(pWTdata, ident.1 = 18, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
pWTdataMark18 <- pWTdataMark18[pWTdataMark18$pct.1 >= 0.1 & pWTdataMark18$pct.2 <= 0.4,]
write.csv(pWTdataMark18, file = "./postRSImarkers/pWTdataMark18.csv")

pWTdataMark19 <- FindMarkers(pWTdata, ident.1 = 19, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
pWTdataMark19 <- pWTdataMark19[pWTdataMark19$pct.1 >= 0.1 & pWTdataMark19$pct.2 <= 0.4,]
write.csv(pWTdataMark19, file = "./postRSImarkers/pWTdataMark19.csv")

pWTdataMark20 <- FindMarkers(pWTdata, ident.1 = 20, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
pWTdataMark20 <- pWTdataMark20[pWTdataMark20$pct.1 >= 0.1 & pWTdataMark20$pct.2 <= 0.4,]
write.csv(pWTdataMark20, file = "./postRSImarkers/pWTdataMark20.csv")

pWTdataMark21 <- FindMarkers(pWTdata, ident.1 = 21, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
pWTdataMark21 <- pWTdataMark21[pWTdataMark21$pct.1 >= 0.1 & pWTdataMark21$pct.2 <= 0.4,]
write.csv(pWTdataMark21, file = "./postRSImarkers/pWTdataMark21.csv")

pWTdataMark22 <- FindMarkers(pWTdata, ident.1 = 22, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
pWTdataMark22 <- pWTdataMark22[pWTdataMark22$pct.1 >= 0.1 & pWTdataMark22$pct.2 <= 0.4,]
write.csv(pWTdataMark22, file = "./postRSImarkers/pWTdataMark22.csv")


## Cluster ID's Found...
# Cluster0 = non hair cells (intermediate)
# Cluster1 = root cap cells
# Cluster2 = mersitem (endo/meristem)
# Cluster3 = hair cells (differentiating)
# Cluster4 = stele
# Cluster5 = phloem
# Cluster6 = stele
# Cluster7 = endodermis
# Cluster8 = hair cell (more mature)
# Cluster9 = meristem
# Cluster10 = hair cells
# Cluster11 = cortex
# Cluster12 = xylem
# Cluster13 = hair cells
# Cluster14 = non hair Cells
# Cluster15 = cortex
# Cluster16 = non hair cells (differentiating)
# Cluster17 = root cap cells
# Cluster18 = protophloem
# Cluster19 = stele
# Cluster20 = meristem
# Cluster21 = protoxylem
# Cluster22 = non hair cells

convertAT <- function(x){
  for(i in 1:length(x)){
    a <- as.character(genes$V1[genes$V2 == x[i]])
    x[i] <- a[1]
  }
  x
}

convertSym <- function(x){
  for(i in 1:length(x)){
    a <- as.character(genes$V2[genes$V1 == x[i]])
    x[i] <- a[1]
  }
  x
}


## Verifying clusters are distinct by comparing average gene expression profiles
pWTdataSCTmat <- pWTdata@assays$SCT@scale.data
clustAvgExp <- matrix(0, nrow= 3000, ncol = 23)
row.names(clustAvgExp) <- row.names(pWTdataSCTmat)
colnames(clustAvgExp) <- c("Cluster0", "Cluster1","Cluster2","Cluster3","Cluster4","Cluster5","Cluster6","Cluster7","Cluster8","Cluster9","Cluster10","Cluster11","Cluster12","Cluster13","Cluster14","Cluster15","Cluster16","Cluster17","Cluster18","Cluster19","Cluster20","Cluster21","Cluster22")
for(i in 1:23){
  x <- rowMeans(pWTdataSCTmat[,which(colnames(pWTdataSCTmat) %in% names(pWTdata$seurat_clusters[pWTdata$seurat_clusters == (i-1)]))])
  clustAvgExp[,i] <- x
}
clustAvgExp <- as.data.frame(clustAvgExp)



newClusterIDs <- c("Non Hair Cells", "Root Cap Cells", "Meristem", "Hair Cells", "Stele",
    "Phloem", "Stele", "Endodermis", "Hair Cells", "Meristem", "Hair Cells", "Cortex", "Xylem", "Hair Cells",
    "Non Hair Cells", "Cortex", "Non Hair Cells", "Root Cap Cells", "Protophloem", "Stele", "Meristem", "Protoxylem", "Non Hair Cells")
names(newClusterIDs) <- levels(pWTdata)
pWTdata <- RenameIdents(pWTdata, newClusterIDs)

png("./postRSIcellIdent.png", units = "in", width = 10, height = 8, res = 600)
DimPlot(pWTdata, reduction = "umap", label = TRUE, pt.size = 0.5, label.size =4.5)
dev.off()
