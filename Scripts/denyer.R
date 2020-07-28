##==================================================================
## Working with Denyer data and marker genes
##==================================================================
denyer.counts <- read.csv("./DenyerMarkerGeneData/GSE123818_Root_single_cell_wt_datamatrix.csv")
rownames(denyer.counts) <- toupper(denyer.counts$X)
denyer.counts <- denyer.counts[,-1]
denyer <- CreateSeuratObject(counts = denyer.counts)

library(sctransform)
denyer <- SCTransform(denyer, verbose = TRUE)
VlnPlot(denyer, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
denyer <- subset(denyer, subset = nFeature_RNA > 1000 & nFeature_RNA < 12500 & nCount_RNA < 200000)

denyer <- FindVariableFeatures(denyer, selection.method = "vst", nfeatures = 4000)
top10.denyer <- head(VariableFeatures(denyer), 10)
LabelPoints(VariableFeaturePlot(denyer), points = top10.denyer, repel = TRUE)

denyer <- RunPCA(denyer, verbose = TRUE)
denyer <- RunTSNE(denyer, verbose = TRUE)
library(reticulate)
denyer <- RunUMAP(denyer, dims = 1:50, verbose = TRUE)
denyer <- FindNeighbors(denyer, dims = 1:50, verbose = TRUE)
denyer <- FindClusters(denyer, verbose = TRUE, resolution = 0.8)   #clustering with high resolution and will regroup manually if necessary

denyer.markers <- FindAllMarkers(denyer, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
denyer.markers <- denyer.markers[denyer.markers$pct.1 >= 0.1 & denyer.markers$pct.2 <= 0.1,]

DimPlot(denyer, reduction = "umap")

# loading marker data
c0 <- read.csv("./DenyerMarkerGeneData/c0markers.csv", header = T)  #mature
c1 <- read.csv("./DenyerMarkerGeneData/c1markers.csv", header = T)  #mature
c2 <- read.csv("./DenyerMarkerGeneData/c2markers.csv", header = T)  #meristem
c3 <- read.csv("./DenyerMarkerGeneData/c3markers.csv", header = T)  #non hair cell
c4 <- read.csv("./DenyerMarkerGeneData/c4markers.csv", header = T)  #stele
c5 <- read.csv("./DenyerMarkerGeneData/c5markers.csv", header = T)  #hair cell
c6 <- read.csv("./DenyerMarkerGeneData/c6markers.csv", header = T)  #meristem
c7 <- read.csv("./DenyerMarkerGeneData/c7markers.csv", header = T)  #meristem
c8 <- read.csv("./DenyerMarkerGeneData/c8markers.csv", header = T)  #meristem
c9 <- read.csv("./DenyerMarkerGeneData/c9markers.csv", header = T)  #cortex
c10 <- read.csv("./DenyerMarkerGeneData/c10markers.csv", header = T)  #hair cell
c11 <- read.csv("./DenyerMarkerGeneData/c11markers.csv", header = T)  #QC/columella
c12 <- read.csv("./DenyerMarkerGeneData/c12markers.csv", header = T)  #xylem
c13 <- read.csv("./DenyerMarkerGeneData/c13markers.csv", header = T)  #endodermis
c14 <- read.csv("./DenyerMarkerGeneData/c14markers.csv", header = T)  #mature
c2sub <- read.csv("./DenyerMarkerGeneData/c2subMarkers.csv", header = T)
c4sub <- read.csv("./DenyerMarkerGeneData/c4subMarkers.csv", header = T)
c6sub <- read.csv("./DenyerMarkerGeneData/c6subMarkers.csv", header = T)
c8sub <- read.csv("./DenyerMarkerGeneData/c8subMarkers.csv", header = T)
c11sub <- read.csv("./DenyerMarkerGeneData/c11subMarkers.csv", header = T)

dMark <- rbind(c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14)


table(dMark$Cluster[dMark$Gene %in% denyer.markers$gene[denyer.markers$cluster == 0]])
table(dMark$Cluster[dMark$Gene %in% denyer.markers$gene[denyer.markers$cluster == 1]])
table(dMark$Cluster[dMark$Gene %in% denyer.markers$gene[denyer.markers$cluster == 2]])
table(dMark$Cluster[dMark$Gene %in% denyer.markers$gene[denyer.markers$cluster == 3]])
table(dMark$Cluster[dMark$Gene %in% denyer.markers$gene[denyer.markers$cluster == 4]])
table(dMark$Cluster[dMark$Gene %in% denyer.markers$gene[denyer.markers$cluster == 5]])
table(dMark$Cluster[dMark$Gene %in% denyer.markers$gene[denyer.markers$cluster == 6]])
table(dMark$Cluster[dMark$Gene %in% denyer.markers$gene[denyer.markers$cluster == 7]])
table(dMark$Cluster[dMark$Gene %in% denyer.markers$gene[denyer.markers$cluster == 8]])
table(dMark$Cluster[dMark$Gene %in% denyer.markers$gene[denyer.markers$cluster == 9]])
table(dMark$Cluster[dMark$Gene %in% denyer.markers$gene[denyer.markers$cluster == 10]])
table(dMark$Cluster[dMark$Gene %in% denyer.markers$gene[denyer.markers$cluster == 11]])
table(dMark$Cluster[dMark$Gene %in% denyer.markers$gene[denyer.markers$cluster == 12]])
table(dMark$Cluster[dMark$Gene %in% denyer.markers$gene[denyer.markers$cluster == 13]])
table(dMark$Cluster[dMark$Gene %in% denyer.markers$gene[denyer.markers$cluster == 14]])
table(dMark$Cluster[dMark$Gene %in% denyer.markers$gene[denyer.markers$cluster == 15]])
table(dMark$Cluster[dMark$Gene %in% denyer.markers$gene[denyer.markers$cluster == 16]])


## Cluster ID's Found...
#Cluster0 = mature
#Cluster1 = meristem
#Cluster2 = meristem
#Cluster3 = mature
#Cluster4 = meristem
#Cluster5 = hair cell
#Cluster6 = hair cell
#Cluster7 = meristem
#Cluster8 = hair cell
#Cluster9 = xylem
#Cluster10 = non hair cell
#Cluster11 = stele
#Cluster12 = endodermis
#Cluster13 = stele
#Cluster14 = cortex
#Cluster15 = QC/columella
#Cluster16 = mature

# Mersitematic cells from denyer's cluster 2 seem to localize to cortex/endodermis --> can be used in pseudotime analysis
meri2 <- names(denyer$seurat_clusters[denyer$seurat_clusters == 1])

newClusterIDs <- c("Mature", "Meristem", "Meristem", "Mature", "Meristem", "Hair Cells",
"Hair Cells", "Meristem", "Hair Cells", "Xylem", "Non Hair Cells", "Stele", "Endodermis",
"Stele", "Cortex", "QC/Columella", "Mature")
names(newClusterIDs) <- levels(denyer)
denyer <- RenameIdents(denyer, newClusterIDs)


##=================================================================================
## Connecting work with Denyer data to Ryu data
##=================================================================================

## Seeing if cell type marker genes allign between Ryu and denyer clusters
#Cortex
cortexMark.ryu <- FindMarkers(WTdata, ident.1 = "Cortex", only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
cortexMark.denyer <- FindMarkers(denyer, ident.1 = "Cortex", only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
cortexR <- row.names(cortexMark.ryu)
cortexD <- row.names(cortexMark.denyer)
for(i in 1:length(cortexD)){cortexD <- c(cortexD, as.character(genes$V2[genes$V1 == cortexD[i]]))}
cortexD <- cortexD[round(length(cortexD)/2+1):length(cortexD)]
length(which(cortexR %in% cortexD))/length(cortexR)

#Meristem
meristemMark.ryu <- FindMarkers(WTdata, ident.1 = "Meristem", only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
meristemMark.denyer <- FindMarkers(denyer, ident.1 = "Meristem", only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
meristemR <- row.names(meristemMark.ryu)
meristemD <- row.names(meristemMark.denyer)
for(i in 1:length(meristemD)){meristemD <- c(meristemD, as.character(genes$V2[genes$V1 == meristemD[i]]))}
meristemD <- meristemD[round(length(meristemD)/2+1):length(meristemD)]
length(which(meristemR %in% meristemD))/length(meristemR)

#Stele
steleMark.ryu <- FindMarkers(WTdata, ident.1 = "Stele", only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
steleMark.denyer <- FindMarkers(denyer, ident.1 = "Stele", only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
steleR <- row.names(steleMark.ryu)
steleD <- row.names(steleMark.denyer)
for(i in 1:length(steleD)){steleD <- c(steleD, as.character(genes$V2[genes$V1 == steleD[i]]))}
steleD <- steleD[round(length(steleD)/2+1):length(steleD)]
length(which(steleR %in% steleD))/length(steleR)

#Endodermis
endoMark.ryu <- FindMarkers(WTdata, ident.1 = "Endodermis", only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
endoMark.denyer <- FindMarkers(denyer, ident.1 = "Endodermis", only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
endoR <- row.names(endoMark.ryu)
endoD <- row.names(endoMark.denyer)
for(i in 1:length(endoD)){endoD <- c(endoD, as.character(genes$V2[genes$V1 == endoD[i]]))}
endoD <- endoD[round(length(endoD)/2+1):length(endoD)]
length(which(endoR %in% endoD))/length(endoR)

#QC/Columella
qcColMark.ryu <- FindMarkers(WTdata, ident.1 = "QC/Columella", only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
qcColMark.denyer <- FindMarkers(denyer, ident.1 = "QC/Columella", only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
qcColR <- row.names(qcColMark.ryu)
qcColD <- row.names(qcColMark.denyer)
for(i in 1:length(qcColD)){qcColD <- c(qcColD, as.character(genes$V2[genes$V1 == qcColD[i]]))}
qcColD <- qcColD[round(length(qcColD)/2+1):length(qcColD)]
length(which(qcColR %in% qcColD))/length(qcColR)

#Hair cells
hairMark.ryu <- FindMarkers(WTdata, ident.1 = "Hair Cells", only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
hairMark.denyer <- FindMarkers(denyer, ident.1 = "Hair Cells", only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
hairR <- row.names(hairMark.ryu)
hairD <- row.names(hairMark.denyer)
for(i in 1:length(hairD)){hairD <- c(hairD, as.character(genes$V2[genes$V1 == hairD[i]]))}
hairD <- hairD[round(length(hairD)/2+1):length(hairD)]
length(which(hairR %in% hairD))/length(hairR)

#Non hair cells
nhairMark.ryu <- FindMarkers(WTdata, ident.1 = "Non Hair Cells", only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
nhairMark.denyer <- FindMarkers(denyer, ident.1 = "Non Hair Cells", only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
nhairR <- row.names(nhairMark.ryu)
nhairD <- row.names(nhairMark.denyer)
for(i in 1:length(nhairD)){nhairD <- c(nhairD, as.character(genes$V2[genes$V1 == nhairD[i]]))}
nhairD <- nhairD[round(length(nhairD)/2+1):length(nhairD)]
length(which(nhairR %in% nhairD))/length(nhairR)

## Finding LRC in denyer data
