##==================================================================
## Working with Jean-Baptiste data and marker genes
##==================================================================
list.files("./Jean-Baptiste/control1/")
list.files("./Jean-Baptiste/control2/")   # [1] "barcodes.tsv" "genes.tsv"    "matrix.mtx"

JB1.counts <-  Read10X(data.dir = "./Jean-Baptiste/control1/")
JB2.counts <-  Read10X(data.dir = "./Jean-Baptiste/control2/")

JB1.obj <- CreateSeuratObject(counts = JB1.counts)
JB2.obj <- CreateSeuratObject(counts = JB2.counts)

JBdata <- merge(JB1.obj, y = JB2.obj, add.cell.ids = c("JB1", "JB2"), project = "plantRNA")

#VlnPlot(JBdata, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
JBdata <- subset(JBdata, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & nCount_RNA > 200 & nCount_RNA < 200000)
JBdata[["percent.mt"]] <- PercentageFeatureSet(JBdata, pattern = "^ATM")
#VlnPlot(JBdata, features = "percent.mt")
JBdata <- subset(JBdata, subset = percent.mt < 10)

library(sctransform)
options(future.globals.maxSize = 6000 * 1024^2)
JBdata <- SCTransform(JBdata, verbose = TRUE)

JBdata <- RunPCA(JBdata, verbose = TRUE)
JBdata <- RunTSNE(JBdata, verbose = TRUE)
library(reticulate)

JBdata <- JackStraw(JBdata, num.replicate = 100)
JBdata <- ScoreJackStraw(JBdata, dims = 1:20)
JackStrawPlot(JBdata, dims = 1:20)

JBdata <- RunUMAP(JBdata, dims = 1:15, verbose = TRUE)
JBdata <- FindNeighbors(JBdata, dims = 1:15, verbose = TRUE)
JBdata <- FindClusters(JBdata, verbose = TRUE, resolution = 0.8)   #clustering with high resolution and will regroup manually if necessary

JBdata.markers <- FindAllMarkers(JBdata, only.pos = TRUE, logfc.threshold = 0.25, test.use = "bimod")
JBdata.markers <- JBdata.markers[JBdata.markers$pct.1 >= 0.1 & JBdata.markers$pct.2 <= 0.1,]

DimPlot(JBdata, reduction = "umap")


readRDS(file= "/home/ahmad/scRNA-seqAnalysis/Jean-Baptiste/GSM3440205_Control_1_cds.rds")
read.table(file = "/home/ahmad/scRNA-seqAnalysis/Jean-Baptiste/GSM3440205_Control_1_pData.tsv", sep = '\t')
