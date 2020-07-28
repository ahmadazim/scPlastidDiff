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

library(Seurat)
library(reticulate)
list.files("/home/ahmad/Ryu/WT1")
list.files("/home/ahmad/Ryu/WT2")   ##for these three lines, output must exactly be...
list.files("/home/ahmad/Ryu/WT3")   # [1] "barcodes.tsv" "genes.tsv"    "matrix.mtx"

WT1 <-  Read10X(data.dir = "/home/ahmad/Ryu/WT1")
WT2 <-  Read10X(data.dir = "/home/ahmad/Ryu/WT2")
WT3 <-  Read10X(data.dir = "/home/ahmad/Ryu/WT3")

#Creating the Seurat object
WT1.obj <- CreateSeuratObject(counts = WT1)
WT2.obj <- CreateSeuratObject(counts = WT2)
WT3.obj <- CreateSeuratObject(counts = WT3)

#Normalizing the data
WT1.obj <- NormalizeData(WT1.obj, normalization.method = "LogNormalize", scale.factor = 10000)
WT2.obj <- NormalizeData(WT2.obj, normalization.method = "LogNormalize", scale.factor = 10000)
WT3.obj <- NormalizeData(WT3.obj, normalization.method = "LogNormalize", scale.factor = 10000)

#Merging the data, while maintaining the normalization (merge.data = T)
WTdata <- merge(WT1.obj, y = c(WT2.obj, WT3.obj), add.cell.ids = c("R1", "R2", "R3"), project = "plantRNA", merge.data = TRUE)

#Examining nFeature_RNA, nCount_RNA, and percent.mt for quality control
#Note that high percentage of mitochondrial DNA indicates cell lysis during sequencing protocol
WTdata[["percent.mt"]] <- PercentageFeatureSet(WTdata, pattern = "^MT-")
VlnPlot(WTdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
WTdata <- subset(WTdata, subset = nFeature_RNA > 1000 & nFeature_RNA < 11500 & percent.mt < 0.015)

#Indentifying highly variable genes
WTdata <- FindVariableFeatures(WTdata, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(WTdata), 10)
LabelPoints(VariableFeaturePlot(WTdata), points = top10, repel = TRUE)

#Scaling/Centering the data with sctransform()
#Removes confounding sources of variation such as mitochondrial mapping percentage
#For more information, visit: https://satijalab.org/seurat/v3.0/sctransform_vignette.html
#Note that I ran lines 52-65 in the background because these functions take a long time to run (~2hrs)
library(sctransform)
WTdata <- SCTransform(WTdata, vars.to.regress = "percent.mt", verbose = TRUE)

WTdata <- RunPCA(WTdata, verbose = TRUE)
WTdata <- RunTSNE(WTdata, verbose = TRUE)

library(reticulate)
WTdata <- RunUMAP(WTdata, dims = 1:30, verbose = TRUE)

WTdata <- FindNeighbors(WTdata, dims = 1:30, verbose = TRUE)
WTdata <- FindClusters(WTdata, verbose = TRUE, resolution = 0.8)   #clustering with high resolution and will regroup manually if necessary

#Find markers for each cluster
WTdata.markers <- FindAllMarkers(WTdata, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "bimod")
#Only keep markers that are expressed more than 0.1 in their cluster and less than 0.1 in other clusters
WTdata.markers <- WTdata.markers[WTdata.markers$pct.1 >= 0.1 & WTdata.markers$pct.2 <= 0.1,]

#Vizulizing clusters
DimPlot(WTdata, reduction = "pca")
DimPlot(WTdata, reduction = "tsne")
DimPlot(WTdata, reduction = "umap")

#ID'ing clusters according to Ryu's Table S3
#upload clusterID.csv file
clustID <- read.csv("/home/ahmad/scRNA-seqAnalysis/clusterID.csv", header = T)
clustID$expID[clustID$expID == "Stele(YoungVascular"] <- "Stele(YoungVascularTissue)"

png("./Stele.png", units = "in", width = 20, height = 16, res = 300)
VlnPlot(WTdata, features = (clustID[clustID$expID == "Stele",])$gene)
dev.off()
png("./RootCap.png", units = "in", width = 20, height = 16, res = 300)
VlnPlot(WTdata, features = (clustID[clustID$expID == "RootCap(Columella)",])$gene)
dev.off()
png("./HairCells.png", units = "in", width = 20, height = 16, res = 300)
VlnPlot(WTdata, features = (clustID[clustID$expID == "Epidermis-RootHairCells(Trichoblasts)",])$gene)
dev.off()
png("./Cortex.png", units = "in", width = 20, height = 16, res = 300)
VlnPlot(WTdata, features = (clustID[clustID$expID == "Cortex",])$gene)
dev.off()
png("./NonHairCells.png", units = "in", width = 20, height = 16, res = 300)
VlnPlot(WTdata, features = (clustID[clustID$expID == "Epidermis-Non-hairCells(Atrichoblasts)",])$gene)
dev.off()
png("./YVStele.png", units = "in", width = 20, height = 16, res = 300)
VlnPlot(WTdata, features = (clustID[clustID$expID == "Stele(YoungVascularTissue)",])$gene)
dev.off()
png("./EpidermisLatRootCap.png", units = "in", width = 20, height = 16, res = 300)
VlnPlot(WTdata, features = (clustID[clustID$expID == "EpidermisandLateralRootCap",])$gene)
dev.off()
png("./Meristem.png", units = "in", width = 20, height = 16, res = 300)
VlnPlot(WTdata, features = (clustID[clustID$expID == "Meristem(Young)Endodermis/Cortex",])$gene)
dev.off()
png("./Endodermis.png", units = "in", width = 20, height = 16, res = 300)
VlnPlot(WTdata, features = (clustID[clustID$expID == "Endodermis",])$gene)
dev.off()
png("./Phloem.png", units = "in", width = 20, height = 16, res = 300)
VlnPlot(WTdata, features = (clustID[clustID$expID == "Phloem",])$gene)
dev.off()
png("./ProtophloemSieveElements.png", units = "in", width = 20, height = 16, res = 300)
VlnPlot(WTdata, features = (clustID[clustID$expID == "ProtophloemSieveElements",])$gene)
dev.off()
png("./Protoxylem.png", units = "in", width = 20, height = 16, res = 300)
VlnPlot(WTdata, features = (clustID[clustID$expID == "Protoxylem",])$gene)
dev.off()
png("./PhloemCompCells.png", units = "in", width = 20, height = 16, res = 300)
VlnPlot(WTdata, features = (clustID[clustID$expID == "PhloemCompanionCells",])$gene)
dev.off()
png("./Xylem.png", units = "in", width = 20, height = 16, res = 300)
VlnPlot(WTdata, features = (clustID[clustID$expID == "Xylem",])$gene)
dev.off()
png("./QC1.png", units = "in", width = 20, height = 16, res = 300)
VlnPlot(WTdata, features = (clustID[clustID$expID == "QC",])$gene[1:15])
dev.off()
png("./QC2.png", units = "in", width = 20, height = 16, res = 300)
VlnPlot(WTdata, features = (clustID[clustID$expID == "QC",])$gene[16:30])
dev.off()
png("./QC3.png", units = "in", width = 20, height = 16, res = 300)
VlnPlot(WTdata, features = (clustID[clustID$expID == "QC",])$gene[31:45])
dev.off()
png("./QC4.png", units = "in", width = 20, height = 16, res = 300)
VlnPlot(WTdata, features = (clustID[clustID$expID == "QC",])$gene[46:52])
dev.off()

## Cluster ID's Found...
#Cluster0 = stele
#Cluster1 = non hair cell
#Cluster2 = QC/columella
#Cluster3 = endodermis
#Cluster4 = hair cells
#Cluster5 = cortex
#Cluster6 = hair cells
#Cluster7 = stele
#Cluster8 = meristem
#Cluster9 = endodermis
#Cluster10 = hair cells
#Cluster11 = stele
#Cluster12 = lateral root cap
#Cluster13 = non hair cells
#Cluster14 = QC/columella
#Cluster15 = phloem
#Cluster16 = stele
#Cluster17 = lateral root cap
#Cluster18 = protoxylem
#Cluster19 = cortex
newClusterIDs <- c("Stele", "Non Hair Cells", "QC/Columella", "Endodermis", "Hair Cells", "Cortex",
    "Hair Cells", "Stele", "Meristem", "Endodermis", "Hair Cells", "Stele", "Lateral Root Cap",
    "Non Hair Cells", "QC/Columella", "Phloem", "Stele", "Lateral Root Cap", "Protoxylem", "Cortex")
names(newClusterIDs) <- levels(WTdata)
WTdata <- RenameIdents(WTdata, newClusterIDs)

png("./Clustering.png", units = "in", width = 10, height = 8, res = 500)
DimPlot(WTdata, reduction = "umap", label = TRUE, pt.size = 0.5, label.size =5.5)
dev.off()

png("./ClusteringNotIdent.png", units = "in", width = 10, height = 8, res = 500)
DimPlot(WTdata, reduction = "umap", label = TRUE, pt.size = 0.5, label.size =7, group.by = "seurat_clusters")
dev.off()


#Separating Columella from rest of Non-hair Cells
nonHairClust <- subset(WTdata, idents = "Non Hair Cells")

nonHairClust <- RunUMAP(nonHairClust, dims = 1:30, verbose = TRUE)
nonHairClust <- FindNeighbors(nonHairClust, dims = 1:30, verbose = TRUE)
nonHairClust <- FindClusters(nonHairClust, verbose = TRUE, resolution = 0.8)

nonHairClust.markers <- FindAllMarkers(nonHairClust, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "bimod")
nonHairClust.markers <- nonHairClust.markers[nonHairClust.markers$pct.1 >= 0.1 & nonHairClust.markers$pct.2 <= 0.1,]
DimPlot(nonHairClust, reduction = "umap")
VlnPlot(nonHairClust, features = c("RGF3", "YUC3", "MT2B", "DI21", "EMB2171","NDPK1",
"UBQ1","ATML3","RPS5A","SAC52","WIH1","DI19","TCTP","DJ-1a","EMB2386","LTI45","YCF1.2","AtANN2"))
VlnPlot(nonHairClust, features = c("PRH-26","GSTU25","CEP2","UBQ10","OPR1","NSP5","GSTU7","SAG24","SRC2","ATAF2",
"UBQ14","GSTU22","SERAT2","GLT1","TR-BAMY","GSTU19","ATAGP2","ATPIS1","NHL3","SAG21","LSR1","RCI2A","GDH2", "IAA7"))

        #Cluster 17 was renamed columella above...

top10.clust <- WTdata.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
png("./clustHeatMap.png", units = "in", width = 20, height = 10, res= 500)
DoHeatmap(WTdata, features = top10.clust$gene)
dev.off()

PlotClusterTree(BuildClusterTree(WTdata))



# #========================================================================================
# # Examine expression of plastid genes in each cluster
# #========================================================================================
# plastGenes <- read.delim("/home/ahmad/scRNA-seqAnalysis/plastidGenes.txt", header = T)
# plastGenes <- plastGenes[,c(2,5,6)]  #create data frame with just locus,  otherName, and description
# colnames(plastGenes) <- c("Locus","Description","otherName")
# plastGenes.AT <- grep("AT", plastGenes$Locus, value = T)
# notAT <- plastGenes$Locus[!is.element(plastGenes$Locus, plastGenes.AT)]
#     #Output: [1] SCAFFOLD_303808.1              SCAFFOLD_201109.1
#     #        [3] FGENESH1_PM.C_SCAFFOLD_6000658 SCAFFOLD_201022.1      ##ALL PLASTID GENES!!!!
#     #        [5] AL_SCAFFOLD_0002_979
# plastGenes <- plastGenes[unique(plastGenes$Locus),]
# write.csv(plastGenes, "./plastGenes.csv")
#
# plastGenes.AT <- plastGenes.AT[is.element(plastGenes.AT, row.names(WTdata)),]  #197
#
# #Convert locus to gene names
# wtn <- gsub(row.names(WTdata), pattern = " ", replace = "")
# plastGenes.u <- plastGenes
# plastGenes.u$Locus <- as.character(plastGenes.u$Locus)
# plastGenes.u <- plastGenes.u[!duplicated(plastGenes.u$Locus), ]
#
# plastGenes.u$nameInData <- '0'
# for(i in 1:nrow(plastGenes.u)) {
#   x <- c(as.character(plastGenes.u$Locus[i]), (strsplit(gsub(as.character(plastGenes.u[i,3]), pattern=" ", replace = ""), ';'))[[1]])
#   if(sum(is.element(x, wtn)) == 1) {
#     plastGenes.u[i, 'nameInData'] <- x[is.element(x, wtn)]
#     #print(as.character(plastGenes.u$Locus[i]))
#   }
# }
# myPlastids <- plastGenes.u[plastGenes.u$nameInData != '0', ]
#
#
# ATCgenes <- c("ATCG00170","ATCG00800", "ATCG01010","ATCG01250","ATCG00160","ATCG00150","ATCG01070","ATCG01060","ATCG01050","ATCG01280","ATCG00580","ATCG00340","ATCG00570","ATCG00330","ATCG00560","ATCG00790","ATCG00550","ATCG00380",
# "ATCG00130","ATCG00360","ATCG00120",'ATCG00350','ATCG00500','ATCG00740',"ATCG00730","ATCG00720","ATCG00700","ATCG00710","ATCG00300","ATCG00540","ATCG00780","ATCG00530","ATCG00770","ATCG00520","ATCG00760","ATCG00750","ATCG00065","ATCG00280","ATCG01120","ATCG00270","ATCG00080","ATCG00070",
# "ATCG00220","ATCG00690","ATCG00210","ATCG00680","ATCG00440","ATCG00430","ATCG00670","ATCG00020","ATCG01110","ATCG00490","ATCG01100","ATCG00480","ATCG00470","ATCG00860","ATCG00600","ATCG00840","ATCG00660","ATCG00650","ATCG00890","ATCG00630","ATCG01090","ATCG01080","ATCG00810","ATCG00820")
# FeaturePlot(WTdata, features = sample(myPlastids$nameInData, 4, replace = FALSE))
#
# # Differential Expression of Plastid Genes between Cell Types
# WTdata.plast <- subset(WTdata, features= myPlastids$nameInData)
# varblPlast <- FindVariableFeatures(WTdata.plast, selection.method = "vst", nfeatures = 100)
# png("./VarGenes.png", units = "in", width = 12, height = 9, res= 300)
# LabelPoints(plot = VariableFeaturePlot(WTdata.plast), points = head(VariableFeatures(WTdata.plast), 10), repel = TRUE)
# dev.off()
#
# # Clustering with all plastid genes
# WTdata.plast <- RunPCA(WTdata.plast, verbose = TRUE)
# WTdata.plast <- RunTSNE(WTdata.plast, verbose = TRUE)
#
# library(reticulate)
# WTdata.plast <- RunUMAP(WTdata.plast, dims = 1:30, verbose = TRUE)
#
# WTdata.plast <- FindNeighbors(WTdata.plast, dims = 1:30, verbose = TRUE)
# WTdata.plast <- FindClusters(WTdata.plast, verbose = TRUE, resolution = 0.8)
#
# WTdata.plast.markers <- FindAllMarkers(WTdata.plast, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "bimod")
# WTdata.plast.markers <- WTdata.plast.markers[WTdata.plast.markers$pct.1 >= 0.1 & WTdata.plast.markers$pct.2 <= 0.1,]
# png("./plastClustAll.png", units = "in", width = 12, height = 9, res= 300)
# DimPlot(WTdata.plast, reduction= "umap")
# dev.off()
#
# # Clustering only with significant PCs
# WTdata.plast <- JackStraw(WTdata.plast, num.replicate = 100)
# WTdata.plast <- ScoreJackStraw(WTdata.plast, dims = 1:20)
# JackStrawPlot(WTdata.plast, dims = 1:20)
# png("./JSplastClustAll.png", units = "in", width = 12, height = 9, res= 300)
# JackStrawPlot(WTdata.plast, dims = 1:20)
# dev.off()
#
# WTdata.plast.dim <- RunUMAP(WTdata.plast, dims = 1:2, verbose = TRUE)
# WTdata.plast.dim <- FindNeighbors(WTdata.plast.dim, dims = 1:2, verbose = TRUE)
# WTdata.plast.dim <- FindClusters(WTdata.plast.dim, verbose = TRUE, resolution = 0.8)
#
# WTdata.plast.dim.markers <- FindAllMarkers(WTdata.plast.dim, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "bimod")
# WTdata.plast.dim.markers <- WTdata.plast.dim.markers[WTdata.plast.dim.markers$pct.1 >= 0.1 & WTdata.plast.dim.markers$pct.2 <= 0.1,]
# png("./plastClustAll.dim.png", units = "in", width = 12, height = 9, res= 300)
# DimPlot(WTdata.plast.dim, reduction= "umap")
# dev.off()
#
# WTdata.plast.dim <- RunTSNE(WTdata.plast.dim, verbose = TRUE, initial_dims = 3)
#
#
# # Clustering with 3 PCs instead of only 2 PCs
# WTdata.plast.3dim <- RunUMAP(WTdata.plast, dims = 1:3, verbose = TRUE)
# WTdata.plast.3dim <- FindNeighbors(WTdata.plast.3dim, dims = 1:3, verbose = TRUE)
# WTdata.plast.3dim <- FindClusters(WTdata.plast.3dim, verbose = TRUE, resolution = 0.8)
# WTdata.plast.3dim.markers <- FindAllMarkers(WTdata.plast.3dim, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "bimod")
# WTdata.plast.3dim.markers <- WTdata.plast.3dim.markers[WTdata.plast.3dim.markers$pct.1 >= 0.1 & WTdata.plast.3dim.markers$pct.2 <= 0.1,]
# DimPlot(WTdata.plast.3dim, reduction= "umap") ## VERY WIERD PLOT?
#
#
# # Adjusting clustering resolution to identify distinct clusters
# WTdata.plast.dim <- FindClusters(WTdata.plast.dim, verbose = TRUE, resolution = 0.02)
# WTdata.plast.dim.markers <- FindAllMarkers(WTdata.plast.dim, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "bimod")
# WTdata.plast.dim.markers <- WTdata.plast.dim.markers[WTdata.plast.dim.markers$pct.1 >= 0.1 & WTdata.plast.dim.markers$pct.2 <= 0.1,]
#
# DimPlot(WTdata.plast.dim, reduction = "umap")
#
# # Examining PC1 and PC2
# png("./PCHeatMap.png", units = "in", width = 12, height = 16, res= 300)
# DimHeatmap(WTdata.plast.dim, dims = 1:15, cells = 500, balanced = TRUE)
# dev.off()
# png("./PCAclust.png", units = "in", width = 12, height = 9, res= 300)
# DimPlot(WTdata.plast.dim, reduction= "pca")
# dev.off()
#
# # Seeing which cell types are in each plastid cluster
# clustIdent <- WTdata.plast.dim
# WTdata.plast.dim@active.ident <- WTdata@active.ident
#
# prop.table(table(clustIdent@active.ident))   # must adjsut for cluster counts by dividing
#   #         0           1
#   # 0.6745176   0.3254824
#
# x <- which(WTdata.plast.dim@active.ident == "Stele")
# prop.table(table(clustIdent@active.ident[x])/c(0.6745176, 0.3254824))  #stele cells in 0 and 1 clusters
# #         0         1
# # 0.4805179 0.5194821
# x <- which(WTdata.plast.dim@active.ident == "Non Hair Cells")
# prop.table(table(clustIdent@active.ident[x])/c(0.6745176, 0.3254824))  #non hair cells in 0 and 1 clusters
# #         0         1
# # 0.5127443 0.4872557
# x <- which(WTdata.plast.dim@active.ident == "Root Cap")
# prop.table(table(clustIdent@active.ident[x])/c(0.6745176, 0.3254824))  #root cap cells in 0 and 1 clusters
# #         0         1
# # 0.5483137 0.4516863
# x <- which(WTdata.plast.dim@active.ident == "Endodermis")
# prop.table(table(clustIdent@active.ident[x])/c(0.6745176, 0.3254824))  #endodermis cells in 0 and 1 clusters
# #         0         1
# # 0.5213602 0.4786398
# x <- which(WTdata.plast.dim@active.ident == "Hair Cells")
# prop.table(table(clustIdent@active.ident[x])/c(0.6745176, 0.3254824))  #hair cells in 0 and 1 clusters
# #         0         1
# # 0.4984066 0.5015934
# x <- which(WTdata.plast.dim@active.ident == "Cortex")
# prop.table(table(clustIdent@active.ident[x])/c(0.6745176, 0.3254824))  #cortex cells in 0 and 1 clusters
# #         0         1
# # 0.5338738 0.4661262
# x <- which(WTdata.plast.dim@active.ident == "Meristem")
# prop.table(table(clustIdent@active.ident[x])/c(0.6745176, 0.3254824))  #meristem cells in 0 and 1 clusters
# #         0         1
# # 0.4885168 0.5114832
# x <- which(WTdata.plast.dim@active.ident == "Lateral Root Cap")
# prop.table(table(clustIdent@active.ident[x])/c(0.6745176, 0.3254824))  #lateral root cap cells in 0 and 1 clusters
# #         0         1
# # 0.4088005 0.5911995
# x <- which(WTdata.plast.dim@active.ident == "Phloem")
# prop.table(table(clustIdent@active.ident[x])/c(0.6745176, 0.3254824))  #phloem cells in 0 and 1 clusters
# #         0         1
# # 0.4790471 0.5209529
# x <- which(WTdata.plast.dim@active.ident == "Columella")
# prop.table(table(clustIdent@active.ident[x])/c(0.6745176, 0.3254824))  #columella cells in 0 and 1 clusters
# #        0        1
# # 0.333822 0.666178
# x <- which(WTdata.plast.dim@active.ident == "Protoxylem")
# prop.table(table(clustIdent@active.ident[x])/c(0.6745176, 0.3254824))  #protoxylem cells in 0 and 1 clusters
# #         0         1
# # 0.5775599 0.4224401
# propClust <- as.data.frame(matrix(c(0.4805179,0.5194821,0.5127443,0.4872557,0.5483137,0.4516863,0.5213602, 0.4786398,0.4984066,0.5015934,0.5338738,0.4661262,0.4885168,0.5114832, 0.4088005,0.5911995,0.4790471, 0.5209529,0.333822, 0.666178,0.5775599,0.4224401), nrow = 11, ncol = 2, byrow=T))
# row.names(propClust) <- c("Stele", "Non Hair Cells","Root Cap","Endodermis","Hair Cells","Cortex","Meristem","Lateral Root Cap","Phloem","Columella","Protoxylem")
# colnames(propClust) <- c("clust0","clust1")
# barplot(c(propClust$clust1, propClust$clust0), col=c("darkblue","red"))
#
# #Heat map of plastid gene expression
# DoHeatmap(WTdata.plast.dim, features = row.names(WTdata.plast.dim))
#
# #Subsetting clusters 1 and o
# plastClust1 <- subset(clustIdent, idents = 1)
# plastClust0 <- subset(clustIdent, idents = 0)
#
#
# # Examining PC1 and PC2 highly expressed genes
# as.matrix(WTdata.plast.dim@reductions$pca[,1][order(WTdata.plast.dim@reductions$pca[,1]),])  ## ordering genes by PC1 estimates
#                   [,1]
# AT5G44020 -0.178520045
# AT5G23750 -0.115888998
# FLA13     -0.096120091
# AT5G46870 -0.095294399
# CYSD2     -0.079554495
# GASA14    -0.044028832
# STR16     -0.041843548
# AT5G25460 -0.034982685
# LTI78     -0.034286388
# BGAL4     -0.033261005
# GDPD1     -0.031592071
# AT5G55050 -0.024115569
# AT5G16010 -0.023211682
# UXS3      -0.010056024
# LHT1      -0.000124398
# FKBP65     0.003488096
# AT5G62350  0.004203450
# KIN2       0.004682596  **
# AT5G22580  0.022336290
# HIPL2      0.023634865
# BCB        0.043772612  **
# PDCB1      0.055948302
# FLA1       0.073230636
# SRK2H      0.074042664  **
# HTB4       0.122542227  **
# HTA7       0.131221951  **
# AT5G20160  0.145027613
# TIM8       0.160139155
# AT5G19510  0.184953931
# AT5G59970  0.213855385  **
# H2B        0.227091766  **
# HDT2       0.233439973
# RPP1C      0.247118262
# AT5G47210  0.248381590
# RPS24B     0.282099599
# HTA6       0.288382286  **
# RPS10B     0.297463644
# AT5G59690  0.342742277  **
# RPL5B      0.383736305
#
# PC1.genes <- c("FKBP65","AT5G62350","KIN2","AT5G22580","HIPL2","BCB","PDCB1","FLA1","SRK2H","HTB4","HTA7","AT5G20160","TIM8","AT5G19510","AT5G59970","H2B","HDT2","RPP1C","AT5G47210","RPS24B","HTA6","RPS10B","AT5G59690","RPL5B")
#
# as.matrix(WTdata.plast.dim@reductions$pca[,2][order(WTdata.plast.dim@reductions$pca[,2]),])  ## ordering genes by PC2 estimates
#                   [,2]
# AT5G59690 -0.454807564
# HTA6      -0.391548249
# H2B       -0.322306508
# AT5G59970 -0.281079708
# BCB       -0.212281340
# HTA7      -0.201846065
# HTB4      -0.167677308
# AT5G44020 -0.074879014
# CYSD2     -0.059160829
# GDPD1     -0.057984744
# KIN2      -0.048856315
# LHT1      -0.041966737
# AT5G25460 -0.028282856
# GASA14    -0.027394756
# AT5G23750 -0.024664105
# SRK2H     -0.022041885
# BGAL4     -0.015570584
# LTI78     -0.015480374
# AT5G46870 -0.014647036
# AT5G16010 -0.012572452
# FLA13     -0.012289466
# STR16     -0.010713180
# UXS3      -0.005129242
# FKBP65     0.003710199
# AT5G62350  0.004590530
# AT5G55050  0.007028956  **
# AT5G22580  0.008643310
# PDCB1      0.021183407
# FLA1       0.044683531
# HIPL2      0.053709181
# AT5G20160  0.070048750
# TIM8       0.123627145
# AT5G19510  0.136932287
# AT5G47210  0.181507771
# RPP1C      0.190895032
# HDT2       0.194403572
# RPS24B     0.210969944
# RPS10B     0.219811484
# RPL5B      0.279623330
#
# # Exploring plastid-encoded genes
# ATCgenes.plast <- ATCgenes[is.element(ATCgenes, plastGenes.u$Locus)]
# ATCgenes.plast <- (plastGenes.u[which(is.element(plastGenes.u$Locus, ATCgenes.plast)),])$nameInData
# ATCgenes.plast <- ATCgenes.plast[-(which(ATCgenes.plast == "0"))]
#
# VlnPlot(WTdata.plast.dim, features = ATCgenes.plast)
# FeaturePlot(WTdata.plast.dim, features = "RPS19")
# FeaturePlot(WTdata.plast.dim, features = "RPL16")
#
# # Subset Meristem and Endodermis, and examine expression of plastid genes along developmental trajectory
# # clust.EndoMeri <- subset(WTdata, idents= c("Endodermis", "Meristem"))
# # Find what genes distinguish endodermis from meristem cluster; those are important in trajectory analysis






#===============================================================================================
## DECIDED IT WAS UNECESSARY TO LOOK AT ONLY HV GENES SO COMMENTED OUT
# Clustering with 39 highly variable plastid genes
# WTdata.plastHV <- subset(WTdata.plast, features= VariableFeatures(WTdata.plast))
#
# WTdata.plastHV <- RunPCA(WTdata.plastHV, verbose = TRUE)
# WTdata.plastHV <- RunTSNE(WTdata.plastHV, verbose = TRUE)
#
# library(reticulate)
# WTdata.plastHV <- RunUMAP(WTdata.plastHV, dims = 1:30, verbose = TRUE)
#
# WTdata.plastHV <- FindNeighbors(WTdata.plastHV, dims = 1:30, verbose = TRUE)
# WTdata.plastHV <- FindClusters(WTdata.plastHV, verbose = TRUE, resolution = 0.8)
#
# WTdata.plastHV.markers <- FindAllMarkers(WTdata.plastHV, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "bimod")
# WTdata.plastHV.markers <- WTdata.plastHV.markers[WTdata.plastHV.markers$pct.1 >= 0.1 & WTdata.plastHV.markers$pct.2 <= 0.1,]
# png("./plastClustHV.png", units = "in", width = 12, height = 9, res= 300)
# DimPlot(WTdata.plastHV, reduction= "umap")
# dev.off()
#
# WTdata.plastHV <- JackStraw(WTdata.plastHV, num.replicate = 100)
# WTdata.plastHV <- ScoreJackStraw(WTdata.plastHV, dims = 1:20)
# JackStrawPlot(WTdata.plastHV, dims = 1:20)
# png("./JSplastClustHV.png", units = "in", width = 12, height = 9, res= 300)
# JackStrawPlot(WTdata.plastHV, dims = 1:20)
# dev.off()
#
# WTdata.plastHV.dim <- RunUMAP(WTdata.plastHV, dims = 1:2, verbose = TRUE)
# WTdata.plastHV.dim <- FindNeighbors(WTdata.plastHV.dim, dims = 1:2, verbose = TRUE)
# WTdata.plastHV.dim <- FindClusters(WTdata.plastHV.dim, verbose = TRUE, resolution = 0.8)
#
# WTdata.plastHV.dim.markers <- FindAllMarkers(WTdata.plastHV.dim, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "bimod")
# WTdata.plastHV.dim.markers <- WTdata.plastHV.dim.markers[WTdata.plastHV.dim.markers$pct.1 >= 0.1 & WTdata.plastHV.dim.markers$pct.2 <= 0.1,]
# png("./plastClustHV.dim.png", units = "in", width = 12, height = 9, res= 300)
# DimPlot(WTdata.plastHV.dim, reduction= "umap")
# dev.off()
#
#
# #Clustering with highly variable plastid genes minus 4 very HV genes
# WTdata.plastHV.wVery <- subset(WTdata.plastHV, features= tail(VariableFeatures(WTdata.plast), n=35))
#
# WTdata.plastHV.wVery <- RunPCA(WTdata.plastHV.wVery, verbose = TRUE)
# WTdata.plastHV.wVery <- RunTSNE(WTdata.plastHV.wVery, verbose = TRUE)
#
# library(reticulate)
# WTdata.plastHV.wVery <- RunUMAP(WTdata.plastHV.wVery, dims = 1:30, verbose = TRUE)
#
# WTdata.plastHV.wVery <- FindNeighbors(WTdata.plastHV.wVery, dims = 1:30, verbose = TRUE)
# WTdata.plastHV.wVery <- FindClusters(WTdata.plastHV.wVery, verbose = TRUE, resolution = 0.8)
#
# WTdata.plastHV.wVery.markers <- FindAllMarkers(WTdata.plastHV.wVery, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "bimod")
# WTdata.plastHV.wVery.markers <- WTdata.plastHV.wVery.markers[WTdata.plastHV.wVery.markers$pct.1 >= 0.1 & WTdata.plastHV.wVery.markers$pct.2 <= 0.1,]
# png("./plastClustHV.wVery.png", units = "in", width = 12, height = 9, res= 300)
# DimPlot(WTdata.plastHV.wVery, reduction= "umap")
# dev.off()
#
# WTdata.plastHV.wVery <- JackStraw(WTdata.plastHV.wVery, num.replicate = 100)
# WTdata.plastHV.wVery <- ScoreJackStraw(WTdata.plastHV.wVery, dims = 1:20)
# JackStrawPlot(WTdata.plastHV.wVery, dims = 1:20)
# png("./JSplastClustHV.wVery.png", units = "in", width = 12, height = 9, res= 300)
# JackStrawPlot(WTdata.plastHV.wVery, dims = 1:20)    #Only 1 PC significant?????
# dev.off()
#
# WTdata.plastHV.wVery.dim <- RunUMAP(WTdata.plastHV.wVery, dims = 1, verbose = TRUE)
# WTdata.plastHV.wVery.dim <- FindNeighbors(WTdata.plastHV.wVery.dim, dims = 1, verbose = TRUE)
# WTdata.plastHV.wVery.dim <- FindClusters(WTdata.plastHV.wVery.dim, verbose = TRUE, resolution = 0.8)
#
# WTdata.plastHV.wVery.dim.markers <- FindAllMarkers(WTdata.plastHV.wVery.dim, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "bimod")
# WTdata.plastHV.wVery.dim.markers <- WTdata.plastHV.wVery.dim.markers[WTdata.plastHV.wVery.dim.markers$pct.1 >= 0.1 & WTdata.plastHV.wVery.dim.markers$pct.2 <= 0.1,]
#=======================================================================================================================================================================


#===============================================================================
## Using Monocle for Pseudotime Analysis
#===============================================================================
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("monocle")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment'))
install.packages("reticulate")
reticulate::py_install("louvain")
devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)






##================================================================================
## New Plastid Gene List from PPDB (more and less curated)
##================================================================================
# Analysis for "raw" list - less curated
plastGenes.uncur <- read.csv("./plastGenes.uncur.csv", header=TRUE)
plastGenes.uncur$Accession <- substr(plastGenes.uncur$Accession, 1, 9)
pGenes.uncur <- unique(plastGenes.uncur$Accession)
# Convert locus names to gene names
geneNames <- read.table(file = "/home/ahmad/Ryu/WT1/genes.tsv", sep = '\t', header = F)
gene.uncur <- c("remove")
for(i in 1:length(pGenes.uncur)){
  index <- match(pGenes.uncur[i], geneNames$V1)
  gene.uncur <- c(gene.uncur, as.character(geneNames$V2[index]))
}
gene.uncur <- as.character(na.omit(gene.uncur[-1]))
plastidSubset.uncur <- subset(WTdata, features= gene.uncur)

# Analysis for more curated list - experimentally verified
plastGenes.cur <- read.csv("./plastGenes.cur.csv", header=TRUE)
plastGenes.cur$Accession <- substr(plastGenes.cur$Accession, 1, 9)
pGenes.cur <- unique(plastGenes.cur$Accession)
# Convert locus names to gene names
gene.cur <- c("remove")
for(i in 1:length(pGenes.cur)){
  index <- match(pGenes.cur[i], geneNames$V1)
  gene.cur <- c(gene.cur, as.character(geneNames$V2[index]))
}
gene.cur <- as.character(na.omit(gene.cur[-1]))
plastidSubset.cur <- subset(WTdata, features= gene.cur)

# Run FindAllMarkers() on subsetted plastid genes (no re-clustering)...
plastMarkersCell <- FindAllMarkers(plastidSubset.cur, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "bimod")
plastMarkersCell <- plastMarkersCell[plastMarkersCell$pct.1 >= 0.1 & plastMarkersCell$pct.2 <= 0.1,]

table(plastMarkersCell$cluster)/table(WTdata@active.ident)

png("./plastGenesClust1.png", units = "in", width = 50, height = 35, res= 600)
FeaturePlot(WTdata, feature= plastMarkersCell$gene[1:20])
dev.off()
png("./plastGenesClust2.png", units = "in", width = 50, height = 35, res= 600)
FeaturePlot(WTdata, feature= plastMarkersCell$gene[21:40])
dev.off()
png("./plastGenesClust3.png", units = "in", width = 50, height = 35, res= 600)
FeaturePlot(WTdata, feature= plastMarkersCell$gene[41:47])
dev.off()


# Expression levels for genes for different plastid types
#Amyloplast
png("./amyloplastFP.png", units = "in", width = 15, height = 9, res = 600)
FeaturePlot(WTdata, features = c("SBE2.1","SVR1","PAC","SS2","SBE3","SS4","ASP5","SGR9","SBE2.2","SS1","DPE1","SGR2","GRV2","SYP22"))
dev.off()
#Proplastid
png("./proplastidFP.png", units = "in", width = 10, height = 5, res = 600)
FeaturePlot(WTdata, features = c("PAC", "LPA3"))
dev.off()

avgGeneExp <- AverageExpression(WTdata)
head(avgGeneExp[["RNA"]][,1:5])

as.numeric(avgGeneExp[["RNA"]][(which(row.names(avgGeneExp[["RNA"]])=="LPA3")),])
as.numeric(avgGeneExp[["RNA"]][(which(row.names(avgGeneExp[["RNA"]])=="PAC")),])

# Looking for proplastid genes in meristem subcluster
# meristem.obj <- subset(WTdata, idents = "Meristem")


# Clustering plastid genes (with **re** normalization)
WTdata.psub <- merge(WT1.obj, y = c(WT2.obj, WT3.obj), add.cell.ids = c("R1", "R2", "R3"), project = "plantRNA")
WTdata.psub <- subset(WTdata.psub, features= gene.uncur)

WTdata.psub[["percent.mt"]] <- PercentageFeatureSet(WTdata.psub, pattern = "^MT-")
VlnPlot(WTdata.psub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
WTdata.psub <- SCTransform(WTdata.psub, verbose = TRUE)

WTdata.psub <- RunPCA(WTdata.psub, verbose = TRUE)
WTdata.psub <- RunTSNE(WTdata.psub, verbose = TRUE)
library(reticulate)
WTdata.psub <- RunUMAP(WTdata.psub, dims = 1:30, verbose = TRUE)

WTdata.psub <- FindNeighbors(WTdata.psub, dims = 1:30, verbose = TRUE)
WTdata.psub <- FindClusters(WTdata.psub, verbose = TRUE, resolution = 0.8)   #clustering with high resolution and will regroup manually if necessary
WTdata.psub.markers <- FindAllMarkers(WTdata.psub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "roc")

png("./plastClust.png", units = "in", width = 10, height = 8, res = 500)
DimPlot(WTdata.psub, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()


png("./amyloplastExp.png", units = "in", width = 20, height = 16, res = 600)
VlnPlot(WTdata.psub, features = c("SBE2.1","SVR1","PAC","SS2","SBE3","SS4","ASP5","SGR9","SBE2.2","SS1","DPE1","SGR2","GRV2","SYP22"))
dev.off()
png("./proplastidExp.png", units = "in", width = 20, height = 16, res = 600)
VlnPlot(WTdata.psub, features = c("LPA3", "PAC"))
dev.off()
png("./chromoplastExp.png", units = "in", width = 20, height = 16, res = 600)
VlnPlot(WTdata.psub, features = c("OASB", "ZKT", "ZDS1", "PAC", "AOX4", "AT5G61670"))
dev.off()
png("./etioplastExp.png", units = "in", width = 20, height = 16, res = 600)
VlnPlot(WTdata.psub, features = c("PAC", "OEP161", "GGPS3", "GGPPS1", "CRTISO", "FLN2", "TOC75-4"))
dev.off()

#Fairly certain that cluster 3 = etioplast
ind <- which(WTdata.psub@active.ident == 3)
etioplastCells <- subset(WTdata.psub, cells = colnames(WTdata)[ind])


# Looking at plastid marker genes in meristem.obj
meristemPlast <- subset(meristem.obj, features = gene.uncur)
meristemPlast <- FindVariableFeatures(meristemPlast, selection.method = "vst", nfeatures = 400)
LabelPoints(VariableFeaturePlot(meristemPlast), points = head(VariableFeatures(meristemPlast),20), repel = TRUE)


# Clustering plastid genes (with **original** normalization)
plastidSubset.uncur <- RunPCA(plastidSubset.uncur, verbose = T)
plastidSubset.uncur <- RunTSNE(plastidSubset.uncur, verbose = T)
plastidSubset.uncur <- RunUMAP(plastidSubset.uncur, dims = 1:30, verbose = TRUE)
plastidSubset.uncur <- FindNeighbors(plastidSubset.uncur, dims = 1:30, verbose = TRUE)
plastidSubset.uncur <- FindClusters(plastidSubset.uncur, verbose = TRUE, resolution = 0.8)
DimPlot(plastidSubset.uncur, reduction = "umap", label = TRUE, pt.size = 0.5)
plastidSubset.uncur.markers <- FindAllMarkers(plastidSubset.uncur, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "bimod")



##=================================================================================
## Working with nonnegative matrix factorization to allign "clusterings"
##=================================================================================
x <- as.data.frame(Idents(WTdata)); names(x) <- "cells"
x[,1] <- as.numeric(x[,1])
y <- as.data.frame(Idents(WTdata.psub)); names(y) <- "plastids"
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

## Get high frequency mnatches:
matchFreq <- as.data.frame(table(paste(x[,1], x[,2], sep='-')))
names(matchFreq) <- c("clusterMatch", "freq")
matchFreq <- matchFreq[matchFreq$freq >= 20,]
matchFreq$clusterMatch <- as.character(matchFreq$clusterMatch)
matchFreq$Cell <- 0
matchFreq$Plastid <- 0
for(i in 1:nrow(matchFreq)) {
  matchFreq$Cell[i] <- as.numeric(strsplit(matchFreq$clusterMatch, '-')[[i]][1])-1
  matchFreq$Plastid[i] <- as.numeric(strsplit(matchFreq$clusterMatch, '-')[[i]][2])-1
}

cellType <- c(
"Stele",
"Non Hair Cells",
"Columella Root Cap",
"Endodermis",
"Hair Cells",
"Cortex",
"Meristem",
"Lateral Root Cap",
"Quiescent Center",
"Stele",
"Stele"
)

matchFreq$cellType <- cellType[matchFreq$Cell]


x <- read.csv("./cell.Idents.Ryu.csv"); names(x) <- c("ids", "ct");
y <- read.csv("./plastid.Idents.Ryu.csv"); names(y) <- c("ids", "pt");
a.mat <- merge(x, y);

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
png("./matchFreq.png", units = "in", width = 8, height = 10, res = 600)
par(mfrow=c(3,3))
for(i in cellTypes)
  barplot(matchFreq[matchFreq$Cell==i,"freq"],
          names.arg=matchFreq[matchFreq$Cell==i,"Plastid"], main=i)
dev.off()


##================================================================================================================================================================================================
##================================================================================================================================================================================================
##================================================================================================================================================================================================
## Working with Dynverse and trajectory inference
## Installed docker on linux with (curl -fsSL https://get.docker.com | sh;)  and (sudo service docker start) and gave it permission using (sudo usermod -a -G docker $USER)
##================================================================================================================================================================================================
##================================================================================================================================================================================================
##================================================================================================================================================================================================
## FROM TERMINAL, RUN (sudo service docker start) and (sudo systemctl enable docker) TO START DOCKER
library(Seurat)
library(dyno)
library(dplyr)
library(tidyverse)
library(umap)

WTdata.dyn <- wrap_expression(
  expression = Matrix::t(WTdata@assays$SCT@scale.data),
  counts = Matrix::t(WTdata$SCT@counts)
)
WTdata.dyn <- add_grouping(WTdata.dyn, WTdata@active.ident)
WTdata.dyn <- add_prior_information(WTdata.dyn, start_id = names(WTdata@active.ident[WTdata@active.ident == "Meristem"]))

## Nature paper describing TI methods... "https://www.nature.com/articles/s41587-019-0071-9"  (A comparison of single-cell trajectory inference methods)
model <- infer_trajectory(WTdata.dyn, ti_paga())

dimred <- dyndimred::dimred_umap(WTdata.dyn$expression)
# Trajectory Inferenence (TI) Plotting based on 4 different approaches to coloring (cell ordering, grouping, feature expression, pseudotime)
patchwork::wrap_plots(
  plot_dimred(model) + ggtitle("Cell ordering"),
  plot_dimred(model, grouping = group_onto_nearest_milestones(model)) + ggtitle("Cell grouping"),
  plot_dimred(model, "pseudotime", pseudotime = calculate_pseudotime(model)) + ggtitle("Pseudotime")
)

##==================================================
## Pseudotime only with Highly Variable (HV) genes
##==================================================
varGenes <- VariableFeatures(WTdata)  ## 3000 genes
all_expression <-  Matrix::t(WTdata@assays$SCT@scale.data)
all_counts <- Matrix::t(WTdata$SCT@counts)
vargenes_expression <- all_expression[,colnames(all_expression) %in% varGenes]
vargenes_counts <-  all_counts[,colnames(all_counts) %in% varGenes]

WTdataHV.dyn <- wrap_expression(
  expression = vargenes_expression,
  counts = vargenes_counts
)
WTdataHV.dyn <- add_grouping(WTdataHV.dyn, WTdata@active.ident)
WTdataHV.dyn <- add_prior_information(WTdataHV.dyn, start_id = names(WTdata@active.ident[WTdata@active.ident == "Meristem"]))

modelHV <- infer_trajectory(WTdataHV.dyn, ti_paga(), verbose = TRUE)
png("./clustHVpaga.png", units = "in", width = 8, height = 8, res = 600)
plot_dimred(model, dimred = dyndimred::dimred_umap(WTdata.dyn$expression), color_cells = "grouping", size_cells = 1, grouping = WTdata@active.ident) + ggtitle("Cell grouping")
dev.off()

origTSNE <- WTdata@reductions$tsne@cell.embeddings
origUMAP <- WTdata@reductions$umap@cell.embeddings

png("./OrigClustHVpaga.png", units = "in", width = 8, height = 8, res = 600)
plot_dimred(model, dimred = origUMAP, color_cells = "grouping", size_cells = 1, grouping = WTdata@active.ident) + ggtitle("Cell grouping")
dev.off()
# ordering makes sense but a little to "vague"... will try slingshot method and then subset known differentiation pathways

##========================
## Slingshot model of TI
##========================
slingshot_model <- infer_trajectory(WTdataHV.dyn, ti_slingshot(), verbose=TRUE)
png("./origClustHVsling.png", units = "in", width = 8, height = 8, res = 600)
plot_dimred(slingshot_model, dimred = origUMAP, color_cells = "grouping",size_cells = 1, grouping = WTdata@active.ident)
dev.off()

##=========================================
## Using PAGA TI with 3 subsets...
##=========================================
# Meristem + Endodermis
TIsub1 <- subset(WTdata, idents = c("Meristem", "Endodermis"))
TIsub1$orig_ident <- TIsub1@active.ident
TIsub1 <- SCTransform(TIsub1, verbose = TRUE)
TIsub1 <- RunPCA(TIsub1, verbose = TRUE)
library(reticulate)
TIsub1 <- RunUMAP(TIsub1, dims = 1:15, verbose = TRUE)
TIsub1 <- FindNeighbors(TIsub1, dims = 1:15, verbose = TRUE)
TIsub1 <- FindClusters(TIsub1, verbose = TRUE, resolution = 0.8)

png("./sub1Recluster.newclust.png", units = "in", width = 10, height = 8, res = 500)
DimPlot(TIsub1, reduction = "umap", group.by= "ident") + ggtitle("Group by New Clusters")
dev.off()
png("./sub1Recluster.png", units = "in", width = 10, height = 8, res = 500)
DimPlot(TIsub1, reduction = "umap", group.by = "orig_ident") + ggtitle("Group by Putative Cell Identity")
dev.off()

            ## TI done only with 3000 HV genes ##
TIsub1 <- FindVariableFeatures(TIsub1, selection.method = "vst", nfeatures = 3000)
varGenes1 <- VariableFeatures(TIsub1)
all_expression1 <-  Matrix::t(TIsub1@assays$SCT@scale.data)
all_counts1 <- Matrix::t(TIsub1$SCT@counts)
vargenes_expression1 <- all_expression1[,colnames(all_expression1) %in% varGenes1]
vargenes_counts1 <-  all_counts1[,colnames(all_counts1) %in% varGenes1]

TIsub1.dyn <- wrap_expression(
  expression = vargenes_expression1,
  counts = vargenes_counts1
)
TIsub1.dyn <- add_grouping(TIsub1.dyn, TIsub1$orig_ident)
TIsub1.dyn <- add_prior_information(TIsub1.dyn, end_id = names(TIsub1$orig_ident[TIsub1$orig_ident == "Endodermis"]), start_id = names(TIsub1$orig_ident[TIsub1$orig_ident == "Meristem"]))
PAGAmodelSub1 <- infer_trajectory(TIsub1.dyn, ti_paga(), verbose = TRUE)
SLINGmodelSub1 <- infer_trajectory(TIsub1.dyn, ti_slingshot(), verbose=TRUE)

png("./sub1TImethodsPAGA.png", units = "in", width = 10, height = 8, res = 500)
plot_dimred(PAGAmodelSub1, color_cells = "grouping", size_cells = 1, grouping = TIsub1$orig_ident) + ggtitle("PAGA TI Method")
dev.off()
png("./sub1TImethodsSLING.png", units = "in", width = 10, height = 8, res = 500)
plot_dimred(SLINGmodelSub1, color_cells = "grouping", size_cells = 1, grouping = TIsub1$orig_ident) + ggtitle("Slingshot TI Method")
dev.off()

# QC + Epidermal cells + Meristem
TIsub2 <- subset(WTdata, idents = c("Quiescent Center", "Meristem", "Non Hair Cells", "Hair Cells", "Columella Root Cap", "Lateral Root Cap"))
TIsub2$orig_ident <- TIsub2@active.ident
TIsub2 <- SCTransform(TIsub2, verbose = TRUE)
TIsub2 <- RunPCA(TIsub2, verbose = TRUE)
library(reticulate)
TIsub2 <- RunUMAP(TIsub2, dims = 1:15, verbose = TRUE)
TIsub2 <- FindNeighbors(TIsub2, dims = 1:15, verbose = TRUE)
TIsub2 <- FindClusters(TIsub2, verbose = TRUE, resolution = 0.8)

png("./sub2Recluster.newclust.png", units = "in", width = 10, height = 8, res = 500)
DimPlot(TIsub2, reduction = "umap", group.by= "ident") + ggtitle("Group by New Clusters")
dev.off()
png("./sub2Recluster.png", units = "in", width = 10, height = 8, res = 500)
DimPlot(TIsub2, reduction = "umap", group.by = "orig_ident") + ggtitle("Group by Putative Cell Identity")
dev.off()

            ## TI done only with 3000 HV genes ##
TIsub2 <- FindVariableFeatures(TIsub2, selection.method = "vst", nfeatures = 3000)
varGenes2 <- VariableFeatures(TIsub2)
all_expression2 <-  Matrix::t(TIsub2@assays$SCT@scale.data)
all_counts2 <- Matrix::t(TIsub2$SCT@counts)
vargenes_expression2 <- all_expression2[,colnames(all_expression2) %in% varGenes2]
vargenes_counts2 <-  all_counts2[,colnames(all_counts2) %in% varGenes2]

TIsub2.dyn <- wrap_expression(
  expression = vargenes_expression2,
  counts = vargenes_counts2
)
TIsub2.dyn <- add_grouping(TIsub2.dyn, TIsub2$orig_ident)
TIsub2.dyn <- add_prior_information(TIsub2.dyn, start_id = names(TIsub2$orig_ident[TIsub2$orig_ident == "Quiescent Center"]))
PAGAmodelSub2 <- infer_trajectory(TIsub2.dyn, ti_paga(), verbose = TRUE)
SLINGmodelSub2 <- infer_trajectory(TIsub2.dyn, ti_slingshot(), verbose=TRUE)

png("./sub2TImethodsPAGA.png", units = "in", width = 10, height = 8, res = 500)
plot_dimred(PAGAmodelSub2, dimred =  TIsub2@reductions$umap@cell.embeddings, color_cells = "grouping", size_cells = 1, grouping = TIsub2$orig_ident, label_milestones = TRUE) + ggtitle("PAGA TI Method")
dev.off()
png("./sub2TImethodsSLING.png", units = "in", width = 10, height = 8, res = 500)
plot_dimred(SLINGmodelSub2, dimred =  TIsub2@reductions$umap@cell.embeddings, color_cells = "grouping", size_cells = 1, grouping = TIsub2$orig_ident) + ggtitle("Slingshot TI Method")
dev.off()


# Meristem + Vascular tissue
TIsub3 <- subset(WTdata, idents = c("Meristem", "Stele", "Protoxylem", "Phloem"))
TIsub3$orig_ident <- TIsub3@active.ident
TIsub3 <- SCTransform(TIsub3, verbose = TRUE)
TIsub3 <- RunPCA(TIsub3, verbose = TRUE)
library(reticulate)
TIsub3 <- RunUMAP(TIsub3, dims = 1:15, verbose = TRUE)
TIsub3 <- FindNeighbors(TIsub3, dims = 1:15, verbose = TRUE)
TIsub3 <- FindClusters(TIsub3, verbose = TRUE, resolution = 0.8)

png("./sub3Recluster.newclust.png", units = "in", width = 10, height = 8, res = 500)
DimPlot(TIsub3, reduction = "umap", group.by= "ident") + ggtitle("Group by New Clusters")
dev.off()
png("./sub3Recluster.png", units = "in", width = 10, height = 8, res = 500)
DimPlot(TIsub3, reduction = "umap", group.by = "orig_ident") + ggtitle("Group by Putative Cell Identity")
dev.off()

            ## TI done only with 3000 HV genes ##
TIsub3 <- FindVariableFeatures(TIsub3, selection.method = "vst", nfeatures = 3000)
varGenes3 <- VariableFeatures(TIsub3)
all_expression3 <-  Matrix::t(TIsub3@assays$SCT@scale.data)
all_counts3 <- Matrix::t(TIsub3$SCT@counts)
vargenes_expression3 <- all_expression3[,colnames(all_expression3) %in% varGenes3]
vargenes_counts3 <-  all_counts3[,colnames(all_counts3) %in% varGenes3]

TIsub3.dyn <- wrap_expression(
  expression = vargenes_expression3,
  counts = vargenes_counts3
)
TIsub3.dyn <- add_grouping(TIsub3.dyn, TIsub3$orig_ident)
TIsub3.dyn <- add_prior_information(TIsub3.dyn, start_id = names(TIsub3$orig_ident[TIsub3$orig_ident == "Meristem"]))
PAGAmodelSub3 <- infer_trajectory(TIsub3.dyn, ti_paga(), verbose = TRUE)
SLINGmodelSub3 <- infer_trajectory(TIsub3.dyn, ti_slingshot(), verbose=TRUE)

png("./sub3TImethodsPAGA.png", units = "in", width = 10, height = 8, res = 500)
plot_dimred(PAGAmodelSub3, dimred =  TIsub3@reductions$umap@cell.embeddings, color_cells = "grouping", size_cells = 1, grouping = TIsub3$orig_ident) + ggtitle("PAGA TI Method")
dev.off()
png("./sub3TImethodsSLING.png", units = "in", width = 10, height = 8, res = 500)
plot_dimred(SLINGmodelSub3, dimred =  TIsub3@reductions$umap@cell.embeddings, color_cells = "grouping", size_cells = 1, grouping = TIsub3$orig_ident) + ggtitle("Slingshot TI Method")
dev.off()


##=======================================================
## Discussion with Tessa 7/24/19
##=======================================================
# 1) slingshot was better than PAGA
# 2) cell state <=> gene expression  --> vector of numbers (20k = number of genes)
# 3) plastid state == the expression of plastid genes (vector of the expression of plastid genes by cells)
      # 3a) plastid state refers to the "placement" of the cell's plastid along its developmental trajectory
      # 3b) infer plastid state based on plastid gene expression levels
# 4) do pseudotime analysis without the plastid genes to examine the progression (TI) of cell state alone
# 5) do pseudotime analysis with only the plastid genes to examine the progression (TI) of plastid state alone
# 6) compare those two pseudotime analyses (diversity of the pseudotime statistics)
      # 6a) if we choose a small "chunk" of the pseudotime and gather cell state and plastid state from the 2 pseudotimes,
          # then we can look at the distance between those two vectors (state = vector)
      # 6b) calculate Euclidian distance between 2 vectors (suggest Seriation to Tessa tomorrow)
            # 6bi) if cell state is defined by vector alpha and plastid state is defined by beta, then we want to find dist(alpha, beta)
            # 6bii) dist() in base R and ser_dist() in Seriation
# 7) Look at a cell type where plastid differentiation is known not to take place as much (i.e. hair cells)
      # 7a) within hair cell cluster, as hair cells develop, plastid "score" should stay more steady than in meristem
      # 7b) plastid state vector should be relatively "constant" throughout psedotime
# 8) To get that vector for cell and plastid state, we will use a function called "scoreGenes"
      # 8a) In that small "chunk" of the pseudotime, what is the expression of genes of interest ("plastid score")
# 9) Using pseudotime graphs --> see genes that are changing in expression levels (i.e. proplastid "stuff")
# 10) The advantage of having plastid clusters identified would be that if we are working with merstem cells, for example, we can know
    # exactly which genes we should be scoring for the plastid vector.
      # 10a) If I can't ID the clusters, what I will do is I will allign the plastid clusters to the cell type clusters and use the markers
           # of the plastid cluster that alligns with the cell types cluster being "pseudotime-ed"




##=======================================================
## Discussion with Tessa 7/25/19
##=======================================================
# Defining cell state by finding average expression in given cluster and using merker genes
avgGeneExp <- AverageExpression(WTdata)

# Define plastid states by matching cell clusters to plastid clusters
# Visualize location of certain cell types in plastid clusters
# First, on re-normalized plastid clusters
# Meristem cells on plastid clusters to hopefully identify proplastid cluster
x <-  names(WTdata@active.ident[WTdata@active.ident == "Meristem"])
png("./meristemPlastPlot.png", units = "in", width = 10, height = 10, res = 500)
DimPlot(WTdata.psub, reduction = "umap", cells.highlight = x)
dev.off()
# CLUSTER 6 = PROPLASTID (MERISTEM)

# Columella cells on plastid clusters to hopefully identify amyloplast cluster
x <-  names(WTdata@active.ident[WTdata@active.ident == "Columella Root Cap"])
png("./columellaPlastPlot.png", units = "in", width = 10, height = 10, res = 500)
DimPlot(WTdata.psub, reduction = "umap", cells.highlight = x)
dev.off()
# CLUSTER 3 = AMYLOPLAST (COLUMELLA ROOT CAP)

# Stele cells on plastid clusters to hopefully identify some plastid type cluster
x <-  names(WTdata@active.ident[WTdata@active.ident == "Stele"])
png("./stelePlastPlot.png", units = "in", width = 10, height = 10, res = 500)
DimPlot(WTdata.psub, reduction = "umap", cells.highlight = x)
dev.off()
# CLUSTER 0 = *whatever is found in stele*

# Hair cells on plastid clusters to hopefully identify some plastid type cluster
x <-  names(WTdata@active.ident[WTdata@active.ident == "Hair Cells"])
png("./hairPlastPlot.png", units = "in", width = 10, height = 10, res = 500)
DimPlot(WTdata.psub, reduction = "umap", cells.highlight = x)
dev.off()
# CLUSTERS 4,9,7 = *whatever is found in hair cells*

# Non Hair cells on plastid clusters to hopefully identify some plastid type cluster
x <-  names(WTdata@active.ident[WTdata@active.ident == "Non Hair Cells"])
png("./nonhairPlastPlot.png", units = "in", width = 10, height = 10, res = 500)
DimPlot(WTdata.psub, reduction = "umap", cells.highlight = x)
dev.off()
# CLUSTER 2 = *whatever is found in non hair cells*

# LRC cells on plastid clusters to hopefully identify some plastid type cluster
x <-  names(WTdata@active.ident[WTdata@active.ident == "Lateral Root Cap"])
png("./lrcPlastPlot.png", units = "in", width = 10, height = 10, res = 500)
DimPlot(WTdata.psub, reduction = "umap", cells.highlight = x)
dev.off()
# CLUSTER 8 = *whatever is found in LRC cells*

# Cortex cells on plastid clusters to hopefully identify some plastid type cluster
x <-  names(WTdata@active.ident[WTdata@active.ident == "Cortex"])
png("./cortexPlastPlot.png", units = "in", width = 10, height = 10, res = 500)
DimPlot(WTdata.psub, reduction = "umap", cells.highlight = x)
dev.off()
# CLUSTER 5,10 = *whatever is found in cortex cells*

# Endodermis cells on plastid clusters to hopefully identify some plastid type cluster
x <-  names(WTdata@active.ident[WTdata@active.ident == "Endodermis"])
png("./endodermisPlastPlot.png", units = "in", width = 10, height = 10, res = 500)
DimPlot(WTdata.psub, reduction = "umap", cells.highlight = x)
dev.off()
# CLUSTER 1 = *whatever is found in endodermis cells*

## NOTE THAT PHLOEM AND PROTOXYLEM CELLS WERE NOT EXAMINED AS THEY ARE TECHNICALLY PART OF THE STELE
x <- WTdata.psub
x@active.ident <- WTdata@active.ident
png("./plastOrigIdent.png", units = "in", width = 10, height = 8, res = 600)
DimPlot(x, reduction = "umap", group.by= "ident", label.size = 6,label = TRUE)
dev.off()

png("./plastClust.png", units = "in", width = 10, height = 8, res = 600)
DimPlot(WTdata.psub, reduction = "umap", label.size = 8,label = TRUE, no.legend = TRUE)
dev.off()


## I will now be working with the amyloplast and proplastid clusters sicne we know that they matched with the columella and meristem clusters, respectively
## First, I will find the marker genes for proplastids by examining marker genes in cluster 6 of WTdata.psub
WTdata.psub.markers <- FindAllMarkers(WTdata.psub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "roc")
proplastMarkers <- WTdata.psub.markers[WTdata.psub.markers$cluster == 6,]
## RESULT --> able to identify palstid-specific marker genes for proplastids nad amyloplasts for sure
write.csv(proplastMarkers, "./proplastidMarkerGenes.csv")
## Also able to find marker genes for amyloplasts by examining marker genes in cluster 3 of WTdata.psub
amyloplastMarkers <- WTdata.psub.markers[WTdata.psub.markers$cluster == 3,]
## RESULT --> able to identify palstid-specific marker genes for proplastids nad amyloplasts for sure
write.csv(amyloplastMarkers, "./amyloplastMarkerGenes.csv")

##=======================================================================================
## Pseudotime Analysis with Amyloplast/Columella cells Columella Root Cap Initials (QC)
##=======================================================================================
nonPlast <- row.names(WTdata)[!is.element(row.names(WTdata), gene.uncur)]
QCcol <-  names(WTdata@active.ident[WTdata@active.ident == "Columella Root Cap" | WTdata@active.ident == "Quiescent Center"])
ptimeNoPlast.columellaQC <- subset(WTdata, features = nonPlast, cells = QCcol)

ptimeNoPlast.columellaQC <- FindVariableFeatures(ptimeNoPlast.columellaQC, nFeatures = 5000)
varGenes.columella <- VariableFeatures(ptimeNoPlast.columellaQC)
allExp.col <-  Matrix::t(ptimeNoPlast.columellaQC@assays$SCT@scale.data)
allCount.col <- Matrix::t(ptimeNoPlast.columellaQC$SCT@counts)
varExp.col <- allExp.col[,colnames(allExp.col) %in% varGenes.columella]
varCount.col <-  varExp.col[,colnames(varExp.col) %in% varGenes.columella]

columellaNoPlast.dyn <- wrap_expression(
  expression = varExp.col,
  counts = varCount.col
)
columellaNoPlast.dyn <- add_grouping(columellaNoPlast.dyn, ptimeNoPlast.columellaQC@active.ident)
initials <- names(WTdata@active.ident[WTdata@active.ident == "Quiescent Center"])
columellaNoPlast.dyn <- add_prior_information(columellaNoPlast.dyn, start_id = initials)

ptimeColumella <- infer_trajectory(columellaNoPlast.dyn, ti_slingshot(), verbose = TRUE)
ptimeColumella.simp <- simplify_trajectory(ptimeColumella)

png("./columellaTIpath.png", units = "in", width = 10, height = 9, res = 600)
plot_dimred(ptimeColumella, label_milestones = TRUE, color_cells = "grouping", size_cells = 1.5, grouping = ptimeNoPlast.columellaQC@active.ident) + ggtitle("Slingshot TI Method used with Columella Root Cap and QC/Root Cap Initial Cells")
dev.off()

png("./columellaTIpathSimp.png", units = "in", width = 10, height = 9, res = 600)
plot_dimred(ptimeColumella.simp, label_milestones = TRUE, color_cells = "grouping", size_cells = 1.5, grouping = ptimeNoPlast.columellaQC@active.ident) + ggtitle("Simplified Slingshot TI Method used with Columella Root Cap and QC/Root Cap Initial Cells")
dev.off()

## Adding no-plastid pseudotime back to WTdata.col
WTdata.col <- subset(WTdata, cells = QCcol)
progressions <- ptimeColumella.simp$progressions
#sorting progressions
matrix(progressions[,1], nrow = 784, ncol = 1) -> x
x <- as.data.frame(x)
x$"2" <- 1:784
y <- 0
for(i in 1:784){
  y <- c(y, which(row.names(WTdata.col@meta.data) == x$V1[i]))
}
y <- y[2:785]
x$"2" <- y
progressions <- progressions[order(x$"2"),]
WTdata.col <- AddMetaData(WTdata.col, progressions, col.name = colnames(progressions))
WTdata.col@meta.data$cell_id <- progressions$cell_id
WTdata.col@meta.data$to <- progressions$to
WTdata.col@meta.data$from <- progressions$from
WTdata.col@meta.data$percentage <- progressions$percentage


## Calculating plastid gene "score" using cellcyclescore
WTdata.col <- CellCycleScoring(WTdata.col, s.features = proplastMarkers$gene, g2m.features = amyloplastMarkers$gene, set.ident = TRUE)
names(WTdata.col@meta.data)[names(WTdata.col@meta.data) == "S.Score"] <- "proplastScore"
names(WTdata.col@meta.data)[names(WTdata.col@meta.data) == "G2M.Score"] <- "amyloplastScore"



col16 <- WTdata.col@meta.data[WTdata.col@meta.data$from == 1 & WTdata.col@meta.data$to == 6,]
col16 <- col16[order(col16$percentage),]

GEptime16 <- function(x){
    GE <- Matrix::t(WTdata.col[["SCT"]][row.names(WTdata.col) == x])
    GE.16 <- GE[row.names(GE) %in% col16$cell_id,]
    GE.16 <- t(t(GE.16))
    percentages <- as.matrix(col16$percentage)
    row.names(percentages) <- col16$cell_id
    ptimeGE <- merge(percentages, GE.16, by = 0)
    # pTime plot for gene
    plot(ptimeGE[,2], ptimeGE[,3])
}

png("./col16.png", units = "in", width = 10, height = 8, res = 600)
par(mar = c(5.1, 4.1, 4.1, 13))
plot(col16$percentage, col16$proplastScore, col=ifelse(col16$old.ident == "Columella Root Cap", "lightgreen", "lightblue"), pch = 16, cex = 0.9, ylab = "Plastid Gene Expression Score", xlab = "Pseudotime", main = "Plastid Gene Expression Along Pseudotime form Milestone 1 to 6")
points(col16$percentage, col16$amyloplastScore, col=ifelse(col16$old.ident == "Columella Root Cap", "lightgreen", "lightblue"), pch = 17, cex = 0.9)
lines(smooth.spline(col16$proplastScore ~ col16$percentage, spar = 0.7), lwd=3, col = 'navyblue')
lines(smooth.spline(col16$amyloplastScore ~ col16$percentage, spar = 0.7), lwd=3, col = 'red')
abline(lm(col16$proplastScore ~ col16$percentage), lwd = 3, col = "navyblue")    ## Very wierd that proplastid expression goes up and amyloplast is opposite
abline(lm(col16$amyloplastScore ~ col16$percentage), col = "red", lwd = 3)
legend("right",  inset=c(-0.33,0), xpd = TRUE,  legend = c("Proplastid Exp in QC", "Proplastid Exp in CRC", "Amyloplast Exp in QC", "Amyloplast Exp in CRC", "Proplastid Exp Regression", "Amyloplast Exp Regression"), pch = c(16,16,17,17, NA, NA), lwd = c(NA,NA,NA,NA, 3, 3), col = c("lightblue", "lightgreen","lightblue", "lightgreen", "navyblue", "red"), cex = 0.8, bg = "white")
dev.off()


## Looking at Pseudotime from M6 to M4
col64 <- WTdata.col@meta.data[WTdata.col@meta.data$from == 6 & WTdata.col@meta.data$to == 4,]
col64 <- col64[order(col64$percentage),]
png("./col64.png", units = "in", width = 10, height = 8, res = 600)
par(mar = c(5.1, 4.1, 4.1, 13))
plot(col64$percentage, col64$proplastScore, col=ifelse(col64$old.ident == "Columella Root Cap", "lightgreen", "lightblue"), pch = 16, cex = 0.9, ylab = "Plastid Gene Expression Score", xlab = "Pseudotime", main = "Plastid Gene Expression Along Pseudotime form Milestone 6 to 4")
points(col64$percentage, col64$amyloplastScore, col=ifelse(col64$old.ident == "Columella Root Cap", "lightgreen", "lightblue"), pch = 17, cex = 0.9)
abline(lm(col64$proplastScore ~ col64$percentage), lwd = 3, col = "navyblue")
abline(lm(col64$amyloplastScore ~ col64$percentage), col = "red", lwd = 3)
legend("right",  inset=c(-0.33,0), xpd = TRUE,  legend = c("Proplastid Exp in QC", "Proplastid Exp in CRC", "Amyloplast Exp in QC", "Amyloplast Exp in CRC", "Proplastid Exp Regression", "Amyloplast Exp Regression"), pch = c(16,16,17,17, NA, NA), lwd = c(NA,NA,NA,NA, 3, 3), col = c("lightblue", "lightgreen","lightblue", "lightgreen", "navyblue", "red"), cex = 0.8, bg = "white")
dev.off()



## Looking at Pseudotime from M6 to M3
col63 <- WTdata.col@meta.data[WTdata.col@meta.data$from == 6 & WTdata.col@meta.data$to == 3,]
col63 <- col63[order(col63$percentage),]
png("./col63.png", units = "in", width = 10, height = 8, res = 600)
par(mar = c(5.1, 4.1, 4.1, 13))
plot(col63$percentage, col63$proplastScore, col=ifelse(col63$old.ident == "Columella Root Cap", "lightgreen", "lightblue"), pch = 16, cex = 0.9, ylab = "Plastid Gene Expression Score", xlab = "Pseudotime", main = "Plastid Gene Expression Along Pseudotime form Milestone 6 to 3")
points(col63$percentage, col63$amyloplastScore, col=ifelse(col63$old.ident == "Columella Root Cap", "lightblue", "lightgreen"), pch = 17, cex = 0.9)
abline(lm(col63$proplastScore ~ col63$percentage), lwd = 3, col = "navyblue")
abline(lm(col63$amyloplastScore ~ col63$percentage), col = "red", lwd = 3)
legend("right",  inset=c(-0.33,0), xpd = TRUE,  legend = c("Proplastid Exp in QC", "Proplastid Exp in CRC", "Amyloplast Exp in QC", "Amyloplast Exp in CRC", "Proplastid Exp Regression", "Amyloplast Exp Regression"), pch = c(16,16,17,17, NA, NA), lwd = c(NA,NA,NA,NA, 3, 3), col = c("lightblue", "lightgreen","lightblue", "lightgreen", "navyblue", "red"), cex = 0.8, bg = "white")
dev.off()









##=======================================================================================
## Pseudotime Analysis with Amyloplast/Columella cells and Meristem/Proplastid
##=======================================================================================
nonPlast <- row.names(WTdata)[!is.element(row.names(WTdata), gene.uncur)]
meriCol <-  names(WTdata@active.ident[WTdata@active.ident == "Columella Root Cap" | WTdata@active.ident == "Meristem"])
ptimeNoPlast.meriCol <- subset(WTdata, features = nonPlast, cells = meriCol)

ptimeNoPlast.meriCol <- FindVariableFeatures(ptimeNoPlast.meriCol, nFeatures = 5000)
varGenes.meriCol <- VariableFeatures(ptimeNoPlast.meriCol)
allExp.meriCol <-  Matrix::t(ptimeNoPlast.meriCol@assays$SCT@scale.data)
allCount.meriCol <- Matrix::t(ptimeNoPlast.meriCol$SCT@counts)
varExp.meriCol <- allExp.meriCol[,colnames(allExp.meriCol) %in% varGenes.meriCol]
varCount.meriCol <-  varExp.meriCol[,colnames(varExp.meriCol) %in% varGenes.meriCol]

meriColNoPlast.dyn <- wrap_expression(
  expression = varExp.meriCol,
  counts = varCount.meriCol
)
meriColNoPlast.dyn <- add_grouping(meriColNoPlast.dyn, ptimeNoPlast.meriCol@active.ident)
initials.meriCol <- names(ptimeNoPlast.meriCol@active.ident[ptimeNoPlast.meriCol@active.ident == "Meristem"])
meriColNoPlast.dyn <- add_prior_information(meriColNoPlast.dyn, start_id = initials.meriCol)

ptimeMeriCol <- infer_trajectory(meriColNoPlast.dyn, ti_slingshot(), verbose = TRUE)
ptimeMeriCol.simp <- simplify_trajectory(ptimeMeriCol)

png("./meriColTIpath.png", units = "in", width = 10, height = 9, res = 600)
plot_dimred(ptimeMeriCol, label_milestones = TRUE, color_cells = "grouping", size_cells = 1.5, grouping = ptimeNoPlast.meriCol@active.ident) + ggtitle("Slingshot TI Method used with Mersitem and Columella Root Cap Cells")
dev.off()

png("./meriColTIpathSimp.png", units = "in", width = 10, height = 9, res = 600)
plot_dimred(ptimeMeriCol.simp, label_milestones = TRUE, color_cells = "grouping", size_cells = 1.9, grouping = ptimeNoPlast.meriCol@active.ident) + ggtitle("Simplified Slingshot TI Method used with Meristem and Columella Root Cap Cells")
dev.off()

## Adding no-plastid pseudotime back to WTdata.col
WTdata.meriCol <- subset(WTdata, cells = meriCol)
progressions <- ptimeMeriCol.simp$progressions
#sorting progressions
matrix(progressions[,1], nrow = 1057, ncol = 1) -> x
x <- as.data.frame(x)
x$"2" <- 1:1057
y <- 0
for(i in 1:1057){
  y <- c(y, which(row.names(WTdata.meriCol@meta.data) == x$V1[i]))
}
y <- y[2:1058]
x$"2" <- y
progressions <- progressions[order(x$"2"),]
WTdata.meriCol <- AddMetaData(WTdata.meriCol, progressions, col.name = colnames(progressions))
WTdata.meriCol@meta.data$cell_id <- progressions$cell_id
WTdata.meriCol@meta.data$to <- progressions$to
WTdata.meriCol@meta.data$from <- progressions$from
WTdata.meriCol@meta.data$percentage <- progressions$percentage


## Calculating plastid gene "score" using cellcyclescore
WTdata.meriCol <- CellCycleScoring(WTdata.meriCol, s.features = proplastMarkers$gene, g2m.features = amyloplastMarkers$gene, set.ident = TRUE)
names(WTdata.meriCol@meta.data)[names(WTdata.meriCol@meta.data) == "S.Score"] <- "proplastScore"
names(WTdata.meriCol@meta.data)[names(WTdata.meriCol@meta.data) == "G2M.Score"] <- "amyloplastScore"

meriCol12 <- WTdata.meriCol@meta.data[WTdata.meriCol@meta.data$from == 1 & WTdata.meriCol@meta.data$to == 2,]
meriCol12 <- meriCol12[order(meriCol12$percentage),]

png("./meriCol12.png", units = "in", width = 10, height = 8, res = 600)
plot(meriCol12$percentage, meriCol12$proplastScore, col=ifelse(meriCol12$old.ident == "Columella Root Cap", "darkcyan", "goldenrod3"), pch = 16, cex = 1.6, ylab = "Plastid Gene Expression Score", xlab = "Pseudotime", main = "Plastid Gene Expression Along Pseudotime")
points(meriCol12$percentage, meriCol12$amyloplastScore, col=ifelse(meriCol12$old.ident == "Columella Root Cap", "darkcyan", "goldenrod3"), pch = 24, cex = 1.6)
abline(lm(meriCol12$proplastScore ~ meriCol12$percentage), lwd = 3, col = "navyblue", lty = 2)
abline(lm(meriCol12$amyloplastScore ~ meriCol12$percentage), col = "orange2", lwd = 3, lty = 2)
dev.off()

png("./meriCol12prop.png", units = "in", width = 10, height = 8, res = 600)
plot(meriCol12$percentage, meriCol12$proplastScore, col=ifelse(meriCol12$old.ident == "Columella Root Cap", "darkcyan", "goldenrod3"), pch = 16, cex = 1.6, ylab = "Plastid Gene Expression Score", xlab = "Pseudotime", main = "Plastid Gene Expression Along Pseudotime")
abline(lm(meriCol12$proplastScore ~ meriCol12$percentage), lwd = 3, col = "navyblue", lty = 2)
dev.off()



## Making sure graph is not just going from A to B --> dividing graph into 3 sections and analyzing slope
x <- meriCol12$percentage; y <- meriCol12$proplastScore
plot(x,y)
x1 <- x[x <= 0.20] # first 419 of x
x2 <- x[x > 0.20 & x <= 0.70] # 420 to 478
x3 <- x[x > 0.70] # 479 to end
y1 <- y[1:419]
y2 <- y[420:478]
y3 <- y[479:1057]

png("./divided3ProplastMeriCol.png", units = "in", width = 15, height = 6, res = 600)
par(mfrow=c(1,3))
plot(x1,y1, ylim = c(-0.8, 0.8), main = "Proplastid Gene Expression from 0% to 20% Pseudotime", pch = 16, cex = 1.1, xlab = "Pseudotime (0% to 20%)", ylab = "Proplastid Gene Expression", col = "darkgray"); abline(lm(y1 ~ x1), lwd = 2, col = "navyblue")
plot(x2,y2, ylim = c(-0.8, 0.8), main = "Proplastid Gene Expression from 20% to 70% Pseudotime", pch = 16, cex = 1.1, xlab = "Pseudotime (20% to 70%)", ylab = "Proplastid Gene Expression", col = "darkgray"); abline(lm(y2 ~ x2), lwd = 2, col = "navyblue")
plot(x3,y3, ylim = c(-0.8, 0.8), main = "Proplastid Gene Expression from 70% to 100% Pseudotime", pch = 16, cex = 1.1, xlab = "Pseudotime (70% to 100%)", ylab = "Proplastid Gene Expression", col = "darkgray"); abline(lm(y3 ~ x3), lwd = 2, col = "navyblue")
dev.off()


a <- meriCol12$percentage; b <- meriCol12$amyloplastScore
plot(a,b)
a1 <- a[a <= 0.20] # first 419 of a
a2 <- a[a > 0.20 & a <= 0.70] # 420 to 478
a3 <- a[a > 0.70] # 479 to end
b1 <- b[1:419]
b2 <- b[420:478]
b3 <- b[479:1057]

png("./divided3AmyloplastMeriCol.png", units = "in", width = 15, height = 6, res = 600)
par(mfrow=c(1,3))
plot(a1,b1, ylim = c(-0.8, 0.8), main = "Amyloplast Gene Expression from 0% to 20% Pseudotime", pch = 16, cex = 1.1, xlab = "Pseudotime (0% to 20%)", ylab = "Amyloplast Gene Expression", col = "darkgray"); abline(lm(b1 ~ a1), lwd = 2, col = "navyblue")
plot(a2,b2, ylim = c(-0.8, 0.8), main = "Amyloplast Gene Expression from 20% to 70% Pseudotime", pch = 16, cex = 1.1, xlab = "Pseudotime (20% to 70%)", ylab = "Amyloplast Gene Expression", col = "darkgray"); abline(lm(b2 ~ a2), lwd = 2, col = "navyblue")
plot(a3,b3, ylim = c(-0.8, 0.8), main = "Amyloplast Gene Expression from 70% to 100% Pseudotime", pch = 16, cex = 1.1, xlab = "Pseudotime (70% to 100%)", ylab = "Amyloplast Gene Expression", col = "darkgray"); abline(lm(b3 ~ a3), lwd = 2, col = "navyblue")
dev.off()







##=======================================================================================
## Pseudotime Analysis with Hair Cells
##=======================================================================================
nonPlast <- row.names(WTdata)[!is.element(row.names(WTdata), gene.uncur)]
hair <-  names(WTdata@active.ident[WTdata@active.ident == "Hair Cells"])
ptimeNoPlast.hair <- subset(WTdata, features = nonPlast, cells = hair)

ptimeNoPlast.hair <- FindVariableFeatures(ptimeNoPlast.hair, nFeatures = 5000)
varGenes.hair <- VariableFeatures(ptimeNoPlast.hair)
allExp.hair <-  Matrix::t(ptimeNoPlast.hair@assays$SCT@scale.data)
allCount.hair <- Matrix::t(ptimeNoPlast.hair$SCT@counts)
varExp.hair <- allExp.hair[,colnames(allExp.hair) %in% varGenes.hair]
varCount.hair <-  varExp.hair[,colnames(varExp.hair) %in% varGenes.hair]

hairNoPlast.dyn <- wrap_expression(
  expression = varExp.hair,
  counts = varCount.hair
)
hairNoPlast.dyn <- add_grouping(hairNoPlast.dyn, ptimeNoPlast.hair@active.ident)
initials.hair <- names(ptimeNoPlast.hair@active.ident[ptimeNoPlast.hair@active.ident == "Hair Cells"])
hairNoPlast.dyn <- add_prior_information(hairNoPlast.dyn, start_id = initials.hair)

ptimehair <- infer_trajectory(hairNoPlast.dyn, ti_slingshot(), verbose = TRUE)
ptimehair.simp <- simplify_trajectory(ptimehair)

png("./hairTIpath.png", units = "in", width = 10, height = 9, res = 600)
plot_dimred(ptimehair, label_milestones = TRUE, color_cells = "pseudotime", size_cells = 1.5) + ggtitle("Slingshot TI Method used with Hair Cells")
dev.off()

png("./hairTIpathSimp.png", units = "in", width = 10, height = 9, res = 600)
plot_dimred(ptimehair.simp, label_milestones = TRUE, color_cells = "pseudotime", size_cells = 1.5) + ggtitle("Simplified Slingshot TI Method used with Hair Cells")
dev.off()

## Adding no-plastid pseudotime back to WTdata.col
WTdata.hair <- subset(WTdata, cells = hair)
progressions <- ptimehair.simp$progressions
#sorting progressions
matrix(progressions[,1], nrow = 1447, ncol = 1) -> x
x <- as.data.frame(x)
x$"2" <- 1:1447
y <- 0
for(i in 1:1447){
  y <- c(y, which(row.names(WTdata.hair@meta.data) == x$V1[i]))
}
y <- y[2:1448]
x$"2" <- y
progressions <- progressions[order(x$"2"),]
WTdata.hair <- AddMetaData(WTdata.hair, progressions, col.name = colnames(progressions))
WTdata.hair@meta.data$cell_id <- progressions$cell_id
WTdata.hair@meta.data$to <- progressions$to
WTdata.hair@meta.data$from <- progressions$from
WTdata.hair@meta.data$percentage <- progressions$percentage


## Calculating plastid gene "score" using cellcyclescore
WTdata.hair <- CellCycleScoring(WTdata.hair, s.features = proplastMarkers$gene, g2m.features = amyloplastMarkers$gene, set.ident = TRUE)
names(WTdata.hair@meta.data)[names(WTdata.hair@meta.data) == "S.Score"] <- "proplastScore"
names(WTdata.hair@meta.data)[names(WTdata.hair@meta.data) == "G2M.Score"] <- "amyloplastScore"

hair34 <- WTdata.hair@meta.data[WTdata.hair@meta.data$from == 3 & WTdata.hair@meta.data$to == 4,]
hair34 <- hair34[order(hair34$percentage),]

png("./hair34.png", units = "in", width = 10, height = 8, res = 600)
par(mar = c(5.1, 4.1, 4.1, 13))
plot(hair34$percentage, hair34$proplastScore, col = "lightblue", pch = 16, cex = 0.9, ylab = "Plastid Gene Expression Score", xlab = "Pseudotime", main = "Plastid Gene Expression Along Pseudotime")
abline(lm(hair34$proplastScore ~ hair34$percentage), lwd = 3, col = "navyblue")    ## Very wierd that proplastid expression goes up and amyloplast is opposite
legend("right",  inset=c(-0.33,0), xpd = TRUE,  legend = c("Proplastid Exp", "Amyloplast Exp", "Proplastid Exp Regression", "Amyloplast Exp Regression"), pch = c(16,16, NA, NA), lwd = c(NA,NA, 3, 3), col = c("lightblue", "lightgreen", "navyblue", "red"), cex = 0.8, bg = "white")
dev.off()
