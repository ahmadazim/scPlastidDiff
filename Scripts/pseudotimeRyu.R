##=======================================================================================
## Running pseudotime on Ryu dataset
## Ideas:
##    Will try between meristem and other cell types to see if we can remove gap
##    Try to find genes that are indicative of plastid and cell state
##        Track the expression of such genes to see if it can suggest anything about cell state and plastid state
##    Try to address question of synchronicity between cell and plastid dofferentiation
##=======================================================================================
# First, I installed docker on my VM by typing this into debian
#       curl -fsSL https://get.docker.com | sh;
#       sudo service docker start
#       sudo usermod -a -G docker $USER
# install the following packages if needed...
library(Seurat)
library(dyno)
library(dplyr)
library(tidyverse)


#=======================================================
# First try between Meristem and Root Cap Cells
# Ended up not doing this because of gap -- ignore
#=======================================================
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
varCount.meriRC <-  allCount.meriRC[,colnames(allCount.meriRC) %in% varGenes.meriRC]

meriRCNoPlast.dyn <- wrap_expression(
  expression = varExp.meriRC,
  counts = varCount.meriRC
)
meriRCNoPlast.dyn <- add_grouping(meriRCNoPlast.dyn, ptimeNoPlast.meriRC@active.ident)
initials.meriRC <- names(ptimeNoPlast.meriRC@active.ident[ptimeNoPlast.meriRC@active.ident == "Meristem"])
meriRCNoPlast.dyn <- add_prior_information(meriRCNoPlast.dyn, start_id = initials.meriRC)

ptimemeriRC <- infer_trajectory(meriRCNoPlast.dyn, ti_scorpius(), verbose = TRUE)
ptimemeriRC.simp <- simplify_trajectory(ptimemeriRC)

png("./NDmeriRCTI.png", units = "in", width = 10, height = 9, res = 600)
plot_dimred(ptimemeriRC, label_milestones = TRUE, color_cells = "grouping", size_cells = 1.5, grouping = ptimeNoPlast.meriRC@active.ident) + ggtitle("Slingshot TI Method used with Mersitem and Root Cap Cells")
dev.off()

png("./NDmeriRCTISimp.png", units = "in", width = 10, height = 9, res = 600)
plot_dimred(ptimemeriRC.simp, label_milestones = TRUE, color_cells = "grouping", size_cells = 1.9, grouping = ptimeNoPlast.meriRC@active.ident) + ggtitle("Simplified Slingshot TI Method used with Meristem and Root Cap Cells")
dev.off()



#==================================================================
# Trying between Meristem and Endodermis
# The difference between this and the next section is that this
#   pseudotime is run on the meristem and endodermis cells without
#   first re-normalizing those cell separately
#===================================================================
# Assumed that you have the Seurat object named "NDpWTdata", which is all the wild-type (WTdata) cells from the Ryu dataset, worked with post-RSI (pWTdata), with No Doublets included in the data (NDpWTdata)
meriEndo <-  names(NDpWTdata@active.ident[NDpWTdata@active.ident == "Endodermis" | NDpWTdata@active.ident == "Meristem"])
ptimeNoPlast.meriEndo <- subset(NDpWTdata, features = nonPlast, cells = meriEndo)  # only including meristematic/endodermal cells "nonPlast" genes in pseudotime object (independent of plastid development)

ptimeNoPlast.meriEndo <- FindVariableFeatures(ptimeNoPlast.meriEndo)
varGenes.meriEndo <- VariableFeatures(ptimeNoPlast.meriEndo)
allExp.meriEndo <-  Matrix::t(ptimeNoPlast.meriEndo@assays$SCT@scale.data)   # grabbing scaled EXPRESSION matrix (genes as columns and cells as rows)
allCount.meriEndo <- Matrix::t(ptimeNoPlast.meriEndo$SCT@counts)  # grabbing scaled COUNTS matrix (genes as columns and cells as rows)
varExp.meriEndo <- allExp.meriEndo[,colnames(allExp.meriEndo) %in% varGenes.meriEndo]   # scaled expression with only 3000 variable features
varCount.meriEndo <-  allCount.meriEndo[,colnames(allCount.meriEndo) %in% varGenes.meriEndo]   # scaled counts with only 3000 variable features

# Creating (called "wrapping") the dynverse object
meriEndoNoPlast.dyn <- wrap_expression(
  expression = varExp.meriEndo,
  counts = varCount.meriEndo
)
meriEndoNoPlast.dyn <- add_grouping(meriEndoNoPlast.dyn, ptimeNoPlast.meriEndo@active.ident)    # adding cell identity (endodermis or meristem) to dynverse object as a grouping
initials.meriEndo <- names(ptimeNoPlast.meriEndo@active.ident[ptimeNoPlast.meriEndo@active.ident == "Meristem"])    # list of cell names which should be starting point of trajectory
meriEndoNoPlast.dyn <- add_prior_information(meriEndoNoPlast.dyn, start_id = initials.meriEndo)    # adding starting cell names to dynverse object

ptimemeriEndo <- infer_trajectory(meriEndoNoPlast.dyn, ti_scorpius(), verbose = TRUE)   # actually inferring the trajectory with the SCORPIUS TI method
ptimemeriEndo.simp <- simplify_trajectory(ptimemeriEndo)    # simplifying the trajectory (so we only have a beginning point and end point)

# Vizualizing the trjectories
png("./NDmeriEndoTI.png", units = "in", width = 10, height = 9, res = 600)
plot_dimred(ptimemeriEndo, label_milestones = TRUE, color_cells = "grouping", size_cells = 1.5, grouping = ptimeNoPlast.meriEndo@active.ident) + ggtitle("SCORPIUS TI Method used with Mersitem and Endodermis (Original Embedding)")
dev.off()



##===========================================================================================
## Normalizing separate dataset of only meristematic and endodermal cells before pseudotime
##===========================================================================================
meriEndo.names <-  names(NDpWTdata@active.ident[NDpWTdata@active.ident == "Endodermis" | NDpWTdata@active.ident == "Meristem"])   # endodermal and meristematic cells names
meriEndo.obj <- merge(WT1.obj, y = c(WT2.obj, WT3.obj), add.cell.ids = c("R1", "R2", "R3"), project = "plantRNA")   # creating Seurat object with all Ryu data first
nonPgene <- row.names(meriEndo.obj)[which(row.names(meriEndo.obj) %in% gene.uncur == F)]   # nonPgene are all genes not included in "gene.uncur" (list of uncurated plastid genes list)
nonPgene <- nonPgene[!(nonPgene %in% as.character(genes$V2[27207:27416]))]   # also taking out all genes that start with MT- or CT- which are at the end of the gene list ("genes")
meriEndo.obj <- subset(meriEndo.obj, features = nonPgene, cells = meriEndo.names)   # subsetting Seurat object with only meri/endo genes and non plastid genes

# Basic quality control and clustering workflow
meriEndo.obj[["percent.mt"]] <- PercentageFeatureSet(meriEndo.obj, pattern = "^MT-")
VlnPlot(meriEndo.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
meriEndo.obj <- subset(meriEndo.obj, subset = nFeature_RNA > 2000 & nFeature_RNA < 8000 & percent.mt < 0.0075 & nCount_RNA < 75000)
meriEndo.obj <- SCTransform(meriEndo.obj, verbose = TRUE)
meriEndo.obj <- RunPCA(meriEndo.obj, verbose = TRUE)
meriEndo.obj <- RunTSNE(meriEndo.obj, verbose = TRUE)
meriEndo.obj <- RunUMAP(meriEndo.obj, dims = 1:30, verbose = TRUE)
meriEndo.obj <- FindNeighbors(meriEndo.obj, dims = 1:30, verbose = TRUE)
meriEndo.obj <- FindClusters(meriEndo.obj, verbose = TRUE, resolution = 0.1) # note the low resolution of 0.1 to only capture 2 clusters of meristem and endodermis

# Labeling meristem and endodermis clusters
newClusterIDs <- c("Meristem", "Endodermis")
names(newClusterIDs) <- levels(meriEndo.obj)
meriEndo.obj <- RenameIdents(meriEndo.obj, newClusterIDs)

# Vizualizing UMAP plot of meristem and endodermis cells
png("meriEndoUMAP.png", units = "in", width = 10, height = 8, res = 600)
DimPlot(meriEndo.obj, reduction = "umap", label = TRUE, label.size = 7, pt.size = 0.5)
dev.off()




##==============================================================================
## Actually running pseudotime between mersitem and endodermis
## For annotations, refer to lines 64-94
##==============================================================================
varGenes.x <- VariableFeatures(meriEndo.obj)
allExp.x <-  Matrix::t(meriEndo.obj@assays$SCT@scale.data)
allCount.x <- Matrix::t(meriEndo.obj$SCT@counts)
varExp.x <- allExp.x[,colnames(allExp.x) %in% varGenes.x]
varCount.x <-  all.Count.x[,colnames(all.Count.x) %in% varGenes.x]

xNoPlast.dyn <- wrap_expression(
    expression = varExp.x,
    counts = varCount.x
)

xNoPlast.dyn <- add_grouping(xNoPlast.dyn, meriEndo.obj@active.ident)
initials.x <- names(meriEndo.obj@active.ident[meriEndo.obj@active.ident == "Meristem"])
xNoPlast.dyn <- add_prior_information(xNoPlast.dyn, start_id = initials.x)
ptimemeriEndo.norm <- infer_trajectory(xNoPlast.dyn, ti_slingshot(), verbose = TRUE)   # Using Slingshot TI method
ptimemeriEndo.norm.simp <- simplify_trajectory(ptimemeriEndo.norm)

# Getting pseudotime progression percentages from dyno object
slingProg <- ptimemeriEndo.norm.simp$progressions
slingProg <- slingProg[order(slingProg$percentage),]

# Vizulaizing simplified trajectory (slingshot TI method)
png("meriEndoNormSimp.png", units = "in", width = 10, height = 8, res = 600)
plot_dimred(ptimemeriEndo.norm.simp, label_milestones = TRUE, color_cells = "grouping", size_cells = 1.5, grouping = meriEndo.obj@active.ident) + ggtitle("Slingshot TI Method used with Re-Normalized Meristem and Endodermis")
dev.off()

# Trying scorpius ti method
meriEndoPtime.scorpius <- infer_trajectory(xNoPlast.dyn, ti_scorpius(), verbose = TRUE)
png("meriEndoNormSCORPIUS.png", units = "in", width = 10, height = 8, res = 600)    # Vizualizing trajectory
plot_dimred(meriEndoPtime.scorpius, label_milestones = TRUE, color_cells = "grouping", size_cells = 1.5, grouping = meriEndo.obj@active.ident) + ggtitle("SCORPIUS TI Method used with Re-Normalized Meristem and Endodermis")
dev.off()
scorpProg <- meriEndoPtime.scorpius$progressions
scorpProg <- scorpProg[order(scorpProg$percentage),]   # getting progression from scorpius TI method


# Trying embeddr ti method
meriEndoPtime.embeddr <- infer_trajectory(xNoPlast.dyn, ti_embeddr(), verbose = TRUE)
png("meriEndoNormEmbeddr.png", units = "in", width = 10, height = 8, res = 600)     # Vizualizing trajectory
plot_dimred(meriEndoPtime.embeddr, label_milestones = TRUE, color_cells = "grouping", size_cells = 1.5, grouping = meriEndo.obj@active.ident) + ggtitle("Embeddr TI Method used with Re-Normalized Meristem and Endodermis")
dev.off()
embProg <- meriEndoPtime.embeddr$progressions
embProg <- embProg[order(embProg$percentage),]    # getting progression from Embeddr TI method


# Comparing scorpius to slingshot to embeddr by plotting cell ordering (x-axis) against its progression percentage (y-axis)
png("scorpVSslingVSembeddr.png", units = "in", width = 10, height = 8, res = 600)
plot(slingProg$percentage, col = "red", ylab = "Percentage Along Pseudotime Axis", xlab = "Cell Index", main = "Comparison of Pseudotime Progression with SCORPIUS versus Slingshot")
points(embProg$percentage, col = "navyblue")
points(scorpProg$percentage, col = "black")
legend(1050, 0.20, col = c("red", "navyblue", "black"), legend = c("Slingshot", "Embeddr", "SCORPIUS"), pch = c(16, 16, 16))
dev.off()
# CHOSE SCORPIUS AS TI INFERENCE METHOD --> most stable increasing along ptime



##================================================================================================
## Creating meriEndo.plast to...
##    First, see if I can score cells based on previously identified  proplastid genes
##    Also, will serve as a re-normalized plastid clustering (did same with cell)
##    Allowing me to analyze plastid gene expression independently of non-plastid gene expresion
##================================================================================================
meriEndoPlast <- merge(WT1.obj, y = c(WT2.obj, WT3.obj), add.cell.ids = c("R1", "R2", "R3"), project = "plantRNA")   # first creating Seurat object with all Ryu data
meriEndoPlast <- subset(meriEndoPlast, features = gene.uncur, cells = meriEndo.names)   # subset to only include meristematic/endodermal cells and only plastid genes ("gene.uncur")
meriEndoPlast[["percent.mt"]] <- PercentageFeatureSet(meriEndoPlast, pattern = "^MT-")    # basic quality control, clustering workflow...
meriEndoPlast <- SCTransform(meriEndoPlast, verbose = TRUE)
meriEndoPlast <- RunPCA(meriEndoPlast, verbose = TRUE)
meriEndoPlast <- RunTSNE(meriEndoPlast, verbose = TRUE)
meriEndoPlast <- RunUMAP(meriEndoPlast, dims = 1:30, verbose = TRUE)
meriEndoPlast <- FindNeighbors(meriEndoPlast, dims = 1:30, verbose = TRUE)
meriEndoPlast <- FindClusters(meriEndoPlast, verbose = TRUE, resolution = 1)

# Vizualizing UMAP plot of meriEndoPlast (plastid classes among meri/endo cells)
png("meriEndoPlastUMAP.png", units = "in", width = 10, height = 8, res = 600)
DimPlot(meriEndoPlast, reduction = "umap", label = TRUE, label.size = 7, pt.size = 0.5)
dev.off()

# First plotting expression of 2 known proplastid genes against psuedotime
# THIS GOT ME NOWHERE SO YOU CAN SKIP TO LINE 258 IF YOU WANT
proplastMarkers <- c("LPA3","PAC")  # only 2 known proplastid marker genes

meriEndoPlast <- CellCycleScoring(meriEndoPlast, s.features = proplastMarkers, g2m.features = amyloplastMarkers$gene, set.ident = TRUE)    # Note that I am using amyloplast markers, but that odesnt really matter --> I just need to give another set of genes (again you can skip this section)
names(meriEndoPlast@meta.data)[names(meriEndoPlast@meta.data) == "S.Score"] <- "proplastScore"
names(meriEndoPlast@meta.data)[names(meriEndoPlast@meta.data) == "G2M.Score"] <- "amyloplastScore"     # refer to comment on line 216

proplastScore <- as.data.frame(meriEndoPlast@meta.data[,"proplastScore"])    # scores for proplastid expression from CellCycleScoring()
row.names(proplastScore) <- row.names(meriEndoPlast@meta.data)   # cell names as row.names
proplastScore$extra <- 0   # just to keep as data.frame and not list

scorpProg <- meriEndoPtime.scorpius$progressions
scorpProg <- scorpProg[order(scorpProg$percentage),]     # just getting progression percentages
row.names(scorpProg) <- scorpProg$cell_id

proplastScore <- proplastScore[row.names(proplastScore) %in% row.names(scorpProg),]    # only using cells common between proplastScore and scorpProg data frames
proplastScore$extra <- NULL
scorpProg <- scorpProg[row.names(scorpProg) %in% row.names(proplastScore),]   # only using cells common between proplastScore and scorpProg data frames
scorpProg$cell_id <- NULL; scorpProg$from <- NULL; scorpProg$to <- NULL    # removing columns that will not be used later

proplastScore$id <- row.names(proplastScore); row.names(proplastScore) <- NULL    # making row names into a column ("id")
scorpProg$id <- row.names(scorpProg); row.names(scorpProg) <- NULL
linkage <- merge(scorpProg, proplastScore, by = "id")   # merging progressions with expression of proplastid genes
colnames(linkage) <- c("id", "Pseudotime", "ProplastidScore")

plot(linkage$Pseudotime, linkage$ProplastidScore)    # plotting ptime vs expression of the 2 proplastid markers
lines(smooth.spline(linkage$Pseudotime, linkage$ProplastidScore, spar = 0.95), lwd=4, col = 'red')

# For some reason, proplastid expression is not decreasing along pseudotime....
# Ideas: maybe check where paseudotime is actually starting and make sure it
#        is at meristem. Also, find better marker genes as these could just be
#        trashy genes for proplastid expression. Maybe look at plastid active
#        differentiation genes. Could also be because I used amyloplast genes,
#        which down-regulated some proplastid genes.
# Also, these genes really have nothing to do with proplastids; the were just differentially expressed --> find better genes
# CONSENSUS: bad work done. Try something else
# let's try something else to get better genes...





##=======================================================================================
## Using dynfeature to see cell-type features that change anywhere along the trajectory
##=======================================================================================
# THIS IS ONLY FOR CELL-TYPE GENES (no plastid genes yet -- just trying this method out and validating that it works)
overall_feature_importances <- dynfeature::calculate_overall_feature_importance(meriEndoPtime.scorpius, expression_source= t(meriEndo.obj@assays$SCT@scale.data))   # basically calculating how important each gene contributes to shift from meristem to endodermis along pseudotime
features <- overall_feature_importances %>% top_n(40, importance) %>% pull(feature_id)   # getting top 40 genes that contribute to development from meristem to endodermis (most dramatic change in expression along ptime)
  # plot_dimred(meriEndoPtime.scorpius) for trajectory again
png("TIheatmap.png", units = "in", width = 16, height = 10, res = 700)   # plotting heatmap to show top 20 genes that changed dramatically over ptime from meristem to endodermis
dynplot::plot_heatmap(meriEndoPtime.scorpius, expression_source= t(meriEndo.obj@assays$SCT@scale.data))  #red represents high expression of gene while blue represents low expression
dev.off()

# After looking at top 20 important genes...
#   meristem (beiginning) shows high expression of genes involved in differentiation and formation of cell parts
#   endodermis (end) shows high expression of genes chracteristic to endodermis, especially casparian strip
#   Next, see if I can do this same thing but based on plastid gene expression profiles


# Validating use of pseudotime to evaluate development by showing with cell type genes
# Will score genes from heatmap (red to blue as proplastid, blue to red as endodermis) and plot against pseudotime just for validation
endoCell <- c("DIR24","CASP1","PER72", "DIR18", "AT1G61590")  # genes that go from blue to red along pseudotime
meriCell <- c("SMT3", "ATHB-20", "ACT7","AT3G44870","AGP4","PIP2-7","PIP2-2","AGP9","TUBA4","CYSD2","PER3","PER39","PME18", "PAP4", "PER45")  # genes that go from red to blue along pseudotime

meriEndo.obj <- CellCycleScoring(meriEndo.obj, s.features = endoCell, g2m.features = meriCell, set.ident = TRUE)   # scoring cells based on endo/meri genes
names(meriEndo.obj@meta.data)[names(meriEndo.obj@meta.data) == "S.Score"] <- "EndoCellGE"     # renaming meta data columns
names(meriEndo.obj@meta.data)[names(meriEndo.obj@meta.data) == "G2M.Score"] <- "MeriCellGE"

endoScore <- as.data.frame(meriEndo.obj@meta.data[,"EndoCellGE"])
meriScore <- as.data.frame(meriEndo.obj@meta.data[,"MeriCellGE"])
cellScores <- cbind(endoScore, meriScore)     # merging meri and endo scores
colnames(cellScores) <- c("endoScore", "meriScore")
cellScores$id <- row.names(meriEndo.obj@meta.data)
# scorpProg <- meriEndoPtime.scorpius$progressions
# scorpProg <- scorpProg[order(scorpProg$percentage),]
cellScorePtime <- merge(scorpProg, cellScores, by = "id")   # merging scorpius prgressions wth endo/meri scores

png("./validatePtimeCell.png", units = "in", width = 10, height = 8, res = 600)     # plotting scores versus pseudotime (again, just validation since I used cell-type genes)
plot(cellScorePtime$percentage, cellScorePtime$endoScore, main = "Cell-Type Gene Expression Along Pseudotime", ylab = "Cell-Type Gene Expression", xlab = "Progression Along Pseudotime", col = "navy blue", pch = 16, ylim = c(-4,4))
abline(lm(cellScorePtime$endoScore ~ cellScorePtime$percentage), lwd = 3, col = "blue", lty = 2)
points(cellScorePtime$percentage, cellScorePtime$meriScore, col = "dark red", pch = 16)
abline(lm(cellScorePtime$meriScore ~ cellScorePtime$percentage), lwd = 3, col = "red", lty = 2)
dev.off()


##===========================================================================================
## Using dynfeature to see plastid-type features that change anywhere along the trajectory
## Trying to find plastid genes that dramatically change in expresion
##    somehere along the trajectory by feeding in plastid gene expression
##    instead of cell type gene expression
##===========================================================================================
pExp <- t(meriEndoPlast@assays$SCT@scale.data)    # sclaed expression of PLASTID genes from meriEndoPlast object
pExp <- pExp[row.names(pExp) %in% meriEndoPtime.scorpius$cell_ids,]     # only using cell names that are in the pseudotime ordering
overall_pfeature_importances <- dynfeature::calculate_overall_feature_importance(meriEndoPtime.scorpius, expression_source= pExp)   # see line 259
pfeatures <- overall_feature_importances %>% top_n(40, importance) %>% pull(feature_id)   # see line 260
  # plot_dimred(meriEndoPtime.scorpius)
png("TIheatmapPlastid.png", units = "in", width = 16, height = 10, res = 700)
dynplot::plot_heatmap(meriEndoPtime.scorpius, expression_source= pExp)  #red represents high expression of gene while blue represents low expression
dev.off()

# Looking at proplastid and "endodermis-plastid" gene expression versus cell development (progression along pseudotime)
endoPlast <- c("NCED3", "NUDT17", "GSTF8", "E1-BETA-2")  # genes that go from blue to red along pseudotime (endodermal plastid genes)
meriPlast <- c("PSP","GLT1","KAS1","CYP74A","BCCP2","CAC2","MOD1","PDH-E1 BETA","LOX4","PYD1","IPP1","BCCP1","DAD1","PPA6", "DRP1C", "AT2G31670")  # genes that go from red to blue along pseudotime (proplastid genes)

meriEndoPlast <- CellCycleScoring(meriEndoPlast, s.features = endoPlast, g2m.features = meriPlast, set.ident = TRUE)   # givng each cell a proplastid and endodermal score
names(meriEndoPlast@meta.data)[names(meriEndoPlast@meta.data) == "S.Score"] <- "EndoPlastGE"
names(meriEndoPlast@meta.data)[names(meriEndoPlast@meta.data) == "G2M.Score"] <- "MeriPlastGE"

endoScorePlast <- as.data.frame(meriEndoPlast@meta.data[,"EndoPlastGE"])
meriScorePlast <- as.data.frame(meriEndoPlast@meta.data[,"MeriPlastGE"])
plastScores <- cbind(endoScorePlast, meriScorePlast)
colnames(plastScores) <- c("endoScore", "meriScore")
plastScores$id <- row.names(meriEndoPlast@meta.data)
# scorpProg <- meriEndoPtime.scorpius$progressions
# scorpProg <- scorpProg[order(scorpProg$percentage),]
plastScorePtime <- merge(scorpProg, plastScores, by = "id")   # merging scores and progressions

png("./ptimePlastid.png", units = "in", width = 10, height = 8, res = 600)    # plotting plastid dev (pExp) vs cell dev (ptime)
plot(plastScorePtime$percentage, plastScorePtime$endoScore, main = "Plastid-Type Gene Expression Along Pseudotime", ylab = "Plastid-Type Gene Expression", xlab = "Progression Along Pseudotime", col = "navy blue", pch = 16)
abline(lm(plastScorePtime$endoScore ~ plastScorePtime$percentage), lwd = 3, col = "blue", lty = 2)
points(plastScorePtime$percentage, plastScorePtime$meriScore, col = "dark red", pch = 16)
abline(lm(plastScorePtime$meriScore ~ plastScorePtime$percentage), lwd = 3, col = "red", lty = 2)
dev.off()

summary(lm(plastScorePtime$endoScore ~ plastScorePtime$percentage))    # seeing if slopes are significant
summary(lm(plastScorePtime$meriScore ~ plastScorePtime$percentage))




##==============================================================================
## Plotting expression of single plastid genes along pseudotime
##==============================================================================
# NCED3
which(row.names(meriEndoPlast@assays$SCT@scale.data) == "NCED3")   # finding which row of expression matrix is this plastid gene
NCED3 <- as.data.frame(meriEndoPlast@assays$SCT@scale.data[571,])   # expression data for each cell
NCED3$ids <- row.names(NCED3)    # making cell names as a column not row name
colnames(NCED3) <- c("GE", "id")
row.names(NCED3) <- NULL
ptNCED3 <- merge(NCED3, scorpProg, by = "id")   # merging expression data with ptime progression
png("./NCED3.png", units = "in", width = 10, height = 6, res = 400)    # plotting expression of genes versus ptime progression
plot(ptNCED3$percentage, ptNCED3$GE, pch = 16, col= "gray", main = "NCED3", xlab = "Pseudotime", ylab = "Gene Expression", cex = 1.1)
lines(smooth.spline(ptNCED3$percentage, ptNCED3$GE, spar = 0.95), lwd=4, col = 'red')
dev.off()

# IPP1  (annotations above applied below)
which(row.names(meriEndoPlast@assays$SCT@scale.data) == "IPP1")
IPP1 <- as.data.frame(meriEndoPlast@assays$SCT@scale.data[822,])
IPP1$ids <- row.names(IPP1)
colnames(IPP1) <- c("GE", "id")
row.names(IPP1) <- NULL
ptIPP1 <- merge(IPP1, scorpProg, by = "id")
png("./IPP1.png", units = "in", width = 10, height = 6, res = 400)
plot(ptIPP1$percentage, ptIPP1$GE, pch = 16, col= "gray", main = "IPP1", xlab = "Pseudotime", ylab = "Gene Expression", cex = 1.1)
lines(smooth.spline(ptIPP1$percentage, ptIPP1$GE, spar = 0.95), lwd=4, col = 'red')
dev.off()

# LOX4
which(row.names(meriEndoPlast@assays$SCT@scale.data) == "LOX4")
LOX4 <- as.data.frame(meriEndoPlast@assays$SCT@scale.data[251,])
LOX4$ids <- row.names(LOX4)
colnames(LOX4) <- c("GE", "id")
row.names(LOX4) <- NULL
ptLOX4 <- merge(LOX4, scorpProg, by = "id")
png("./LOX4.png", units = "in", width = 10, height = 6, res = 400)
plot(ptLOX4$percentage, ptLOX4$GE, pch = 16, col= "gray", main = "LOX4", xlab = "Pseudotime", ylab = "Gene Expression", cex = 1.1)
lines(smooth.spline(ptLOX4$percentage, ptLOX4$GE, spar = 0.95), lwd=4, col = 'red')
dev.off()

# KAS1
which(row.names(meriEndoPlast@assays$SCT@scale.data) == "KAS1")
KAS1 <- as.data.frame(meriEndoPlast@assays$SCT@scale.data[817,])
KAS1$ids <- row.names(KAS1)
colnames(KAS1) <- c("GE", "id")
row.names(KAS1) <- NULL
ptKAS1 <- merge(KAS1, scorpProg, by = "id")
png("./KAS1.png", units = "in", width = 10, height = 6, res = 400)
plot(ptKAS1$percentage, ptKAS1$GE, pch = 16, col= "gray", main = "KAS1", xlab = "Pseudotime", ylab = "Gene Expression", cex = 1.1)
lines(smooth.spline(ptKAS1$percentage, ptKAS1$GE, spar = 0.95), lwd=4, col = 'red')
dev.off()

# NUDT17
which(row.names(meriEndoPlast@assays$SCT@scale.data) == "NUDT17")
NUDT17 <- as.data.frame(meriEndoPlast@assays$SCT@scale.data[316,])
NUDT17$ids <- row.names(NUDT17)
colnames(NUDT17) <- c("GE", "id")
row.names(NUDT17) <- NULL
ptNUDT17 <- merge(NUDT17, scorpProg, by = "id")
png("./NUDT17.png", units = "in", width = 10, height = 6, res = 400)
plot(ptNUDT17$percentage, ptNUDT17$GE, pch = 16, col= "gray", main = "NUDT17", xlab = "Pseudotime", ylab = "Gene Expression", cex = 1.1)
lines(smooth.spline(ptNUDT17$percentage, ptNUDT17$GE, spar = 0.95), lwd=4, col = 'red')
dev.off()

# GSTF8
which(row.names(meriEndoPlast@assays$SCT@scale.data) == "GSTF8")
GSTF8 <- as.data.frame(meriEndoPlast@assays$SCT@scale.data[514,])
GSTF8$ids <- row.names(GSTF8)
colnames(GSTF8) <- c("GE", "id")
row.names(GSTF8) <- NULL
ptGSTF8 <- merge(GSTF8, scorpProg, by = "id")
png("./GSTF8.png", units = "in", width = 10, height = 6, res = 400)
plot(ptGSTF8$percentage, ptGSTF8$GE, pch = 16, col= "gray", main = "GSTF8", xlab = "Pseudotime", ylab = "Gene Expression", cex = 1.1)
lines(smooth.spline(ptGSTF8$percentage, ptGSTF8$GE, spar = 0.95), lwd=4, col = 'red')
dev.off()


# geneName  --> use for any gene to do same thing as above
plotGEptime <- function(gene){
  geneName <- as.data.frame(meriEndoPlast@assays$SCT@scale.data[(which(row.names(meriEndoPlast@assays$SCT@scale.data) == gene)),])
  geneName$ids <- row.names(geneName)
  colnames(geneName) <- c("GE", "id")
  row.names(geneName) <- NULL
  ptgeneName <- merge(geneName, scorpProg, by = "id")
  plot(ptgeneName$percentage, ptgeneName$GE, pch = 16, col= "gray", main = gene, xlab = "Pseudotime", ylab = "Gene Expression", cex = 1.1)
  lines(smooth.spline(ptgeneName$percentage, ptgeneName$GE, spar = 1), lwd=4, col = 'red')
}




##========================================================================
## Checking if plastid genes are just "housekeeping" genes
## Are those genes changing dramatically just because they are being
##      downregulated by the high expression of endodermal genes toward
##      the end of the trajectory??
## Try to see how much other clusters express those genes....
##========================================================================
png("./propGenesClust.png", units = "in", width = 20, height = 20, res = 1000)
FeaturePlot(NDpWTdata.psub, features = meriPlast)
dev.off()

png("./endoGenesClust.png", units = "in", width = 8, height = 8, res = 1000)
FeaturePlot(NDpWTdata.psub, features = endoPlast)
dev.off()


## Making a function to get expression of each a gene within a cluster
pGeneExp <- NDpWTdata.psub@assays$SCT@scale.data

meanExp.byClust <- function(gene){
  gene.exp <- as.data.frame(pGeneExp[row.names(pGeneExp) == gene,])
  gene.exp$ids <- row.names(gene.exp); row.names(gene.exp) <- NULL; colnames(gene.exp) <- c("Exp", "ids")
  pIdent <- as.data.frame(NDpWTdata.psub@active.ident)
  pIdent$ids <- row.names(pIdent); row.names(pIdent) <- NULL; colnames(pIdent) <- c("ident", "ids")
  gene.expId <- merge(pIdent, gene.exp, by = "ids")
  meansGene <- gene.expId %>%
    group_by(ident) %>%
    summarize(average = mean(Exp))
  as.data.frame(meansGene)
}

propGE <- data.frame(ident = 1:13, Exp = 0)
for(i in meriPlast){
  x <- meanExp.byClust(i)
  colnames(x) <- c("ident", i)
  propGE <- merge(propGE, x, by = "ident")
}
propGE <- propGE[,-(1:2)]


endoGE <- data.frame(ident = 1:13, Exp = 0)
for(i in endoPlast){
  x <- meanExp.byClust(i)
  colnames(x) <- c("ident", i)
  endoGE <- merge(endoGE, x, by = "ident")
}
endoGE <- endoGE[,-(1:2)]



##============================================
## Comparing embeddings (original vs Re-...)
##============================================
scorpProg.orig <- ptimemeriEndo$progressions
scorpProg.orig <- scorpProg.orig[order(scorpProg.orig$percentage),]
row.names(scorpProg.orig) <- scorpProg.orig$cell_id
scorpProg.orig$cell_id <- NULL; scorpProg.orig$from <- NULL; scorpProg.orig$to <- NULL
scorpProg.orig$id <- row.names(scorpProg.orig); row.names(scorpProg.orig) <- NULL
plastScorePtime.orig <- merge(scorpProg.orig, plastScores, by = "id")

png("./compareEmbeddings.png", units = "in", width = 21, height = 8, res = 600)
par(mfrow = c(1,2))
plot(plastScorePtime.orig$percentage, plastScorePtime.orig$endoScore, main = "Original Embedding", ylab = "Plastid-Type Gene Expression", xlab = "Progression Along Pseudotime", col = "navy blue", pch = 16, ylim = c(-2,3), cex = 1)
points(plastScorePtime.orig$percentage, plastScorePtime.orig$meriScore, col = "dark red", pch = 16, cex = 1)
abline(lm(plastScorePtime.orig$meriScore ~ plastScorePtime.orig$percentage), lwd = 3, col = "red", lty = 2)
abline(lm(plastScorePtime.orig$endoScore ~ plastScorePtime.orig$percentage), lwd = 3, col = "blue", lty = 2)
plot(plastScorePtime$percentage, plastScorePtime$endoScore, main = "Re-Normalized, Re-Scaled, Re-Clustered Embedding", ylab = "Plastid-Type Gene Expression", xlab = "Progression Along Pseudotime", col = "navy blue", pch = 16)
abline(lm(plastScorePtime$endoScore ~ plastScorePtime$percentage), lwd = 3, col = "blue", lty = 2)
points(plastScorePtime$percentage, plastScorePtime$meriScore, col = "dark red", pch = 16)
abline(lm(plastScorePtime$meriScore ~ plastScorePtime$percentage), lwd = 3, col = "red", lty = 2)
dev.off()


png("./compareEmbeddings_progression.png", units = "in", width = 10, height = 10, res = 600)
plot(sort(ptimemeriEndo$progressions$percentage), col = "navy blue", main = "\"Gradualness\" of Pseudotime Progression", xlab = "Cell Ordering", ylab = "Progression Along Pseudotime")
points(sort(plastScorePtime$percentage), col = "dark red")
legend("bottomright", legend = c("Original Embedding", "Re-Normalized, Re-Scaled, Re-Clustered Embedding"), pch = 16, col = c("blue", "red"))
dev.off()

x <- t(NDpWTdata.psub@assays$SCT@scale.data)
x <- x[row.names(x) %in% ptimemeriEndo$cell_ids,]
dynplot::plot_heatmap(ptimemeriEndo, expression_source= x)









## Running pseudotime with plastid genes only (2/15/20) RE-EMBEDDING
meriEndo.names <-  names(NDpWTdata@active.ident[NDpWTdata@active.ident == "Endodermis" | NDpWTdata@active.ident == "Meristem"])   # endodermal and meristematic cells names
meriEndoPlast.obj <- merge(WT1.obj, y = c(WT2.obj, WT3.obj), add.cell.ids = c("R1", "R2", "R3"), project = "plantRNA")   # creating Seurat object with all Ryu data first
meriEndoPlast.obj <- subset(meriEndoPlast.obj, features = gene.uncur, cells = meriEndo.names)   # subsetting Seurat object with only meri/endo genes and non plastid genes

# Basic quality control and clustering workflow
meriEndoPlast.obj[["percent.mt"]] <- PercentageFeatureSet(meriEndoPlast.obj, pattern = "^MT-")
VlnPlot(meriEndoPlast.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
meriEndoPlast.obj <- subset(meriEndoPlast.obj, subset = nFeature_RNA < 480 & nCount_RNA < 2500)
meriEndoPlast.obj <- SCTransform(meriEndoPlast.obj, verbose = TRUE)
meriEndoPlast.obj <- RunPCA(meriEndoPlast.obj, verbose = TRUE)
meriEndoPlast.obj <- RunTSNE(meriEndoPlast.obj, verbose = TRUE)
meriEndoPlast.obj <- RunUMAP(meriEndoPlast.obj, dims = 1:30, verbose = TRUE)
meriEndoPlast.obj <- FindNeighbors(meriEndoPlast.obj, dims = 1:30, verbose = TRUE)
meriEndoPlast.obj <- FindClusters(meriEndoPlast.obj, verbose = TRUE, resolution = 0.5)

# Labeling meristem and endodermis clusters
newClusterIDs <- c("Meristem", "Endodermis")
names(newClusterIDs) <- levels(meriEndoPlast.obj)
meriEndoPlast.obj <- RenameIdents(meriEndoPlast.obj, newClusterIDs)

# Vizualizing UMAP plot of meristem and endodermis cells
DimPlot(meriEndoPlast.obj, reduction = "umap", label = TRUE, label.size = 7, pt.size = 0.5)


ptimePlast.meriEndo <- subset(meriEndoPlast.obj, features = gene.uncur)

ptimePlast.meriEndo <- FindVariableFeatures(ptimePlast.meriEndo)
varPlastGenes.meriEndo <- VariableFeatures(ptimePlast.meriEndo)
allPlastExp.meriEndo <-  Matrix::t(ptimePlast.meriEndo@assays$SCT@scale.data)
allPlastCount.meriEndo <- Matrix::t(ptimePlast.meriEndo$SCT@counts)
varPlastExp.meriEndo <- allPlastExp.meriEndo[,colnames(allPlastExp.meriEndo) %in% varPlastGenes.meriEndo]
varPlastCount.meriEndo <-  allPlastCount.meriEndo[,colnames(allPlastCount.meriEndo) %in% varPlastGenes.meriEndo]

# Creating (called "wrapping") the dynverse object
meriEndoPlast.dyn <- wrap_expression(
  expression = varPlastExp.meriEndo,
  counts = varPlastCount.meriEndo
)
meriEndoPlast.dyn <- add_grouping(meriEndoPlast.dyn, ptimePlast.meriEndo@active.ident)
meriEndoPlast.dyn <- add_prior_information(meriEndoPlast.dyn, start_id = meriNames)

ptimePlastMeriEndo <- infer_trajectory(meriEndoPlast.dyn, ti_scorpius(), verbose = TRUE)
ptimePlastMeriEndo.simp <- simplify_trajectory(ptimePlastMeriEndo)

# Vizualizing the trjectories
plot_dimred(ptimePlastMeriEndo, label_milestones = TRUE, color_cells = "grouping", size_cells = 1.5, grouping = ptimeNoPlast.meriEndo@active.ident) + ggtitle("SCORPIUS TI Method used with Mersitem and Endodermis PLASTID GENES ONLY (Re-embedded)")
