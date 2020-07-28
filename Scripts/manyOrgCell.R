library(Seurat)
library(dyno)
library(dplyr)
library(tidyverse)

cellState <- ptimemeriEndo$progressions$percentage  ## CELL STATE
meriEndoPtime.scorpius  ## Ptime trajectory for cell state

geneList <- function(filePath){
  a <- read.csv(filePath, header=TRUE)
  a$Accession <- substr(a$Accession, 1, 9)
  b <- unique(a$Accession)
  # Convert locus names to gene names
  geneNames <- read.table(file = "/home/ahmad/Ryu/WT1/genes.tsv", sep = '\t', header = F)
  a <- c("remove")
  for(i in 1:length(b)){
    index <- match(b[i], geneNames$V1)
    a <- c(b, as.character(geneNames$V2[index]))
  }
  a <- as.character(na.omit(a[-1]))
  a
}

##=================================================
## Plastid development
##=================================================

# First, prepare list of plastid genes
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


# Next, create meriEndoPlast object
meriEndoPlast <- merge(WT1.obj, y = c(WT2.obj, WT3.obj), add.cell.ids = c("R1", "R2", "R3"), project = "plantRNA")   # first creating Seurat object with all Ryu data
meriEndoPlast <- subset(meriEndoPlast, features = gene.uncur, cells = meriEndo.names)   # subset to only include meristematic/endodermal cells and only plastid genes ("gene.uncur")
meriEndoPlast[["percent.mt"]] <- PercentageFeatureSet(meriEndoPlast, pattern = "^MT-")    # basic quality control, clustering workflow...
meriEndoPlast <- SCTransform(meriEndoPlast, verbose = TRUE)
meriEndoPlast <- RunPCA(meriEndoPlast, verbose = TRUE)
meriEndoPlast <- RunTSNE(meriEndoPlast, verbose = TRUE)
meriEndoPlast <- RunUMAP(meriEndoPlast, dims = 1:30, verbose = TRUE)
meriEndoPlast <- FindNeighbors(meriEndoPlast, dims = 1:30, verbose = TRUE)
meriEndoPlast <- FindClusters(meriEndoPlast, verbose = TRUE, resolution = 1)

#DimPlot(meriEndoPlast, reduction = "umap", label = TRUE, label.size = 7, pt.size = 0.5) # Vizualizing UMAP plot of meriEndoPlast (plastid classes among meri/endo cells)

# Then, get list of dramtically-changing plastid genes that will be used to assess plastid state
pExp <- t(meriEndoPlast@assays$SCT@scale.data)    # sclaed expression of PLASTID genes from meriEndoPlast object
pExp <- pExp[row.names(pExp) %in% meriEndoPtime.scorpius$cell_ids,]     # only using cell names that are in the pseudotime ordering
#overall_pfeature_importances <- dynfeature::calculate_overall_feature_importance(meriEndoPtime.scorpius, expression_source= pExp) # basically calculating how important each gene contributes to shift from meristem to endodermis along pseudotime
plastHeatmap <- dynplot::plot_heatmap(meriEndoPtime.scorpius, expression_source= pExp, features_oi = 50)
pfeatures.table  <- plastHeatmap$assemble$plots[[4]]$data

a <- matrix(0, nrow = 30, ncol = 3)
for(i in 1:30){
  x <- pfeatures.table[pfeatures.table$feature_id == unique(pfeatures.table$feature_id)[i],]
  lm.x <- lm(expression  ~ percentage, data = x)
  a[i,1] <- unique(x$feature_id)
  a[i,2] <- coef(lm.x)[2]
  a[i,3] <- summary(lm.x)$coefficients[2,4]
}
a[order(a[,2]),]
goi <-  a[a[,2] > 0,1] # goi = genes of interest

# Finally, looking at endodermis-plastid gene expression versus cell development (progression along pseudotime)
      # Average GE of goi
plastScores <- AddModuleScore(object = meriEndoPlast, features = goi, ctrl = 16, name = "Feature")
plastScores <- plastScores@meta.data[, grep("^Feature", colnames(plastScores@meta.data))]
plastScores$Avg <- rowMeans(plastScores)
plastScores$id <- row.names(plastScores)
plastScorePtime <- merge(scorpProg, plastScores, by = "id")   # merging scores and progressions

# Plotting cell state vs plastid state
plot(plastScorePtime$percentage, plastScorePtime$Avg, main = "Plastid-Type Gene Expression Along Pseudotime", ylab = "Plastid-Type Gene Expression", xlab = "Progression Along Pseudotime", col = "navy blue", pch = 16)
abline(lm(plastScorePtime$Avg ~ plastScorePtime$percentage), lwd = 3, col = "blue", lty = 2)














##=================================================
## Vacuole development
##=================================================

# First, prepare list of vacuole genes
vacGenes <- geneList("./vacuoleGenes.csv")


# Next, create meriEndoVac object
meriEndoVac <- merge(WT1.obj, y = c(WT2.obj, WT3.obj), add.cell.ids = c("R1", "R2", "R3"), project = "plantRNA")   # first creating Seurat object with all Ryu data
meriEndoVac <- subset(meriEndoVac, features = vacGenes, cells = meriEndo.names)   # subset to only include meristematic/endodermal cells and only vacuole genes ("VacGenes")
meriEndoVac[["percent.mt"]] <- PercentageFeatureSet(meriEndoVac, pattern = "^MT-")    # basic quality control, clustering workflow...
meriEndoVac <- SCTransform(meriEndoVac, verbose = TRUE)
meriEndoVac <- RunPCA(meriEndoVac, verbose = TRUE)
meriEndoVac <- RunUMAP(meriEndoVac, dims = 1:10, verbose = TRUE)
meriEndoVac <- FindNeighbors(meriEndoVac, dims = 1:10, verbose = TRUE)
meriEndoVac <- FindClusters(meriEndoVac, verbose = TRUE, resolution = 1)

# DimPlot(meriEndoVac, reduction = "umap", label = TRUE, label.size = 7, pt.size = 0.5) # Vizualizing UMAP plot of meriEndoVac (vacuole classes among meri/endo cells)

# Then, get list of dramtically-changing vacuole genes that will be used to assess vacuole state
vacExp <- t(meriEndoVac@assays$SCT@scale.data)    # sclaed expression of VACUOLE genes from meriEndoVac object
vacExp <- vacExp[row.names(vacExp) %in% meriEndoPtime.scorpius$cell_ids,]     # only using cell names that are in the pseudotime ordering
# overall_pfeature_importances <- dynfeature::calculate_overall_feature_importance(meriEndoPtime.scorpius, expression_source= vacExp) # basically calculating how important each gene contributes to shift from meristem to endodermis along pseudotime
vacHeatmap <- dynplot::plot_heatmap(meriEndoPtime.scorpius, expression_source= vacExp, features_oi = 50)
vacfeatures.table  <- vacHeatmap$assemble$plots[[4]]$data

a <- matrix(0, nrow = length(unique(vacfeatures.table$feature_id)), ncol = 3)
for(i in 1:length(unique(vacfeatures.table$feature_id))){
  x <- vacfeatures.table[vacfeatures.table$feature_id == unique(vacfeatures.table$feature_id)[i],]
  lm.x <- lm(expression  ~ percentage, data = x)
  a[i,1] <- unique(x$feature_id)
  a[i,2] <- coef(lm.x)[2]
  a[i,3] <- summary(lm.x)$coefficients[2,4]
}
a[order(a[,2]),]
goi.vac <-  a[,1] # goi = genes of interest

# Finally, looking at endodermis-vacuole gene expression versus cell development (progression along pseudotime)
      # Average GE of goi.vac
vacScores <- AddModuleScore(object = meriEndoVac, features = goi.vac, ctrl = 2, name = "Feature", nbin=2)
vacScores <- vacScores@meta.data[, grep("^Feature", colnames(vacScores@meta.data))]
vacScores$Avg <- rowMeans(vacScores)
vacScores$id <- row.names(vacScores)
vacScorePtime <- merge(scorpProg, vacScores, by = "id")   # merging scores and progressions

# Plotting cell state vs vacuole state
plot(vacScorePtime$percentage, vacScorePtime$Avg, main = "Vacuole-Type Gene Expression Along Pseudotime", ylab = "Vacuole-Type Gene Expression", xlab = "Progression Along Pseudotime", col = "navy blue", pch = 16)
abline(lm(vacScorePtime$Avg ~ vacScorePtime$percentage), lwd = 3, col = "blue", lty = 2)











##=================================================
## Peroxisome development
##=================================================
# First, prepare list of peroxisome genes
peroxGenes <- geneList("./peroxGenes.csv")

# Next, create meriEndoPerox object
meriEndoPerox <- merge(WT1.obj, y = c(WT2.obj, WT3.obj), add.cell.ids = c("R1", "R2", "R3"), project = "plantRNA")   # first creating Seurat object with all Ryu data
meriEndoPerox <- subset(meriEndoPerox, features = peroxGenes, cells = meriEndo.names)   # subset to only include meristematic/endodermal cells and only peroxisome genes ("gene.uncur")
# meriEndoPerox[["percent.mt"]] <- PercentageFeatureSet(meriEndoPerox, pattern = "^MT-")    # basic quality control, clustering workflow...
meriEndoPerox <- SCTransform(meriEndoPerox, verbose = TRUE)
meriEndoPerox <- RunPCA(meriEndoPerox, verbose = TRUE, approx=FALSE)
meriEndoPerox <- RunTSNE(meriEndoPerox, verbose = TRUE)
meriEndoPerox <- RunUMAP(meriEndoPerox, dims = 1:30, verbose = TRUE)
meriEndoPerox <- FindNeighbors(meriEndoPerox, dims = 1:30, verbose = TRUE)
meriEndoPerox <- FindClusters(meriEndoPerox, verbose = TRUE, resolution = 1)

# DimPlot(meriEndoPerox, reduction = "umap", label = TRUE, label.size = 7, pt.size = 0.5) # Vizualizing UMAP plot of meriEndoPerox (peroxisome classes among meri/endo cells)

# Then, get list of dramtically-changing peroxisome genes that will be used to assess peroxisome state
peroxExp <- t(meriEndoPerox@assays$SCT@scale.data)    # sclaed expression of peroxisome genes from meriEndoPerox object
peroxExp <- peroxExp[row.names(peroxExp) %in% meriEndoPtime.scorpius$cell_ids,]     # only using cell names that are in the pseudotime ordering
# overall_pfeature_importances <- dynfeature::calculate_overall_feature_importance(meriEndoPtime.scorpius, expression_source= peroxExp) # basically calculating how important each gene contributes to shift from meristem to endodermis along pseudotime
PeroxHeatmap <- dynplot::plot_heatmap(meriEndoPtime.scorpius, expression_source= peroxExp, features_oi = 50)
perox_features.table  <- PeroxHeatmap$assemble$plots[[4]]$data

#START HERE
a <- matrix(0, nrow = 30, ncol = 3)
for(i in 1:30){
  x <- perox_features.table[perox_features.table$feature_id == unique(perox_features.table$feature_id)[i],]
  lm.x <- lm(expression  ~ percentage, data = x)
  a[i,1] <- unique(x$feature_id)
  a[i,2] <- coef(lm.x)[2]
  a[i,3] <- summary(lm.x)$coefficients[2,4]
}
a[order(a[,2]),]
a <- a[-a[,3] > 0.05,]
goi.perox <-  a[,1] # goi = genes of interest

# Finally, looking at endodermis-peroxisome gene expression versus cell development (progression along pseudotime)
      # Average GE of goi
peroxScores <- AddModuleScore(object = meriEndoPerox, features = goi.perox, ctrl = 5, name = "Feature", nbin = 3)
peroxScores <- peroxScores@meta.data[, grep("^Feature", colnames(peroxScores@meta.data))]
peroxScores$Avg <- rowMeans(peroxScores)
peroxScores$id <- row.names(peroxScores)
peroxScorePtime <- merge(scorpProg, peroxScores, by = "id")   # merging scores and progressions

# Plotting cell state vs peroxisome state
plot(peroxScorePtime$percentage, peroxScorePtime$Avg, main = "peroxisome-Type Gene Expression Along Pseudotime", ylab = "peroxisome-Type Gene Expression", xlab = "Progression Along Pseudotime", col = "navy blue", pch = 16)
abline(lm(peroxScorePtime$Avg ~ peroxScorePtime$percentage), lwd = 3, col = "blue", lty = 2)












##=================================================
## Plasma Membrane development
##=================================================
# First, prepare list of plasma membrane genes
PmemGenes <- geneList("./pmemGenes.csv")

# Next, create meriEndoPmem object
meriEndoPmem <- merge(WT1.obj, y = c(WT2.obj, WT3.obj), add.cell.ids = c("R1", "R2", "R3"), project = "plantRNA")   # first creating Seurat object with all Ryu data
meriEndoPmem <- subset(meriEndoPmem, features = PmemGenes, cells = meriEndo.names)   # subset to only include meristematic/endodermal cells and only plasma membrane genes ("gene.uncur")
# meriEndoPmem[["percent.mt"]] <- PercentageFeatureSet(meriEndoPmem, pattern = "^MT-")    # basic quality control, clustering workflow...
meriEndoPmem <- SCTransform(meriEndoPmem, verbose = TRUE)
meriEndoPmem <- RunPCA(meriEndoPmem, verbose = TRUE, approx=FALSE)
meriEndoPmem <- RunTSNE(meriEndoPmem, verbose = TRUE)
meriEndoPmem <- RunUMAP(meriEndoPmem, dims = 1:30, verbose = TRUE)
meriEndoPmem <- FindNeighbors(meriEndoPmem, dims = 1:30, verbose = TRUE)
meriEndoPmem <- FindClusters(meriEndoPmem, verbose = TRUE, resolution = 1)

# DimPlot(meriEndoPmem, reduction = "umap", label = TRUE, label.size = 7, pt.size = 0.5) # Vizualizing UMAP plot of meriEndoPmem (plasma membrane classes among meri/endo cells)

# Then, get list of dramtically-changing plasma membrane genes that will be used to assess plasma membrane state
PmemExp <- t(meriEndoPmem@assays$SCT@scale.data)    # sclaed expression of plasma membrane genes from meriEndoPmem object
PmemExp <- PmemExp[row.names(PmemExp) %in% meriEndoPtime.scorpius$cell_ids,]     # only using cell names that are in the pseudotime ordering
# overall_pfeature_importances <- dynfeature::calculate_overall_feature_importance(meriEndoPtime.scorpius, expression_source= PmemExp) # basically calculating how important each gene contributes to shift from meristem to endodermis along pseudotime
PmemHeatmap <- dynplot::plot_heatmap(meriEndoPtime.scorpius, expression_source= PmemExp, features_oi = 50)
Pmem_features.table  <- PmemHeatmap$assemble$plots[[4]]$data

#START HERE
a <- matrix(0, nrow = 30, ncol = 3)
for(i in 1:30){
  x <- Pmem_features.table[Pmem_features.table$feature_id == unique(Pmem_features.table$feature_id)[i],]
  lm.x <- lm(expression  ~ percentage, data = x)
  a[i,1] <- unique(x$feature_id)
  a[i,2] <- coef(lm.x)[2]
  a[i,3] <- summary(lm.x)$coefficients[2,4]
}
a[order(a[,2]),]
a <- a[-a[,3] > 0.05,]
goi.Pmem <-  a[,1] # goi = genes of interest

# Finally, looking at endodermis-plasma membrane gene expression versus cell development (progression along pseudotime)
      # Average GE of goi
PmemScores <- AddModuleScore(object = meriEndoPmem, features = goi.Pmem, ctrl = 16, name = "Feature")
PmemScores <- PmemScores@meta.data[, grep("^Feature", colnames(PmemScores@meta.data))]
PmemScores$Avg <- rowMeans(PmemScores)
PmemScores$id <- row.names(PmemScores)
PmemScorePtime <- merge(scorpProg, PmemScores, by = "id")   # merging scores and progressions

# Plotting cell state vs plasma membrane state
plot(PmemScorePtime$percentage, PmemScorePtime$Avg, main = "plasma membrane-Type Gene Expression Along Pseudotime", ylab = "plasma membrane-Type Gene Expression", xlab = "Progression Along Pseudotime", col = "navy blue", pch = 16)
abline(lm(PmemScorePtime$Avg ~ PmemScorePtime$percentage), lwd = 3, col = "blue", lty = 2)









##=================================================
## Mitochondria development
##=================================================
# First, prepare list of mitochondria genes
MitoGenes <- geneList("./mitoGenes.csv")

# Next, create meriEndoMito object
meriEndoMito <- merge(WT1.obj, y = c(WT2.obj, WT3.obj), add.cell.ids = c("R1", "R2", "R3"), project = "plantRNA")   # first creating Seurat object with all Ryu data
meriEndoMito <- subset(meriEndoMito, features = MitoGenes, cells = meriEndo.names)   # subset to only include meristematic/endodermal cells and only mitochondria genes ("gene.uncur")
# meriEndoMito[["percent.mt"]] <- PercentageFeatureSet(meriEndoMito, pattern = "^MT-")    # basic quality control, clustering workflow...
meriEndoMito <- SCTransform(meriEndoMito, verbose = TRUE)
meriEndoMito <- RunPCA(meriEndoMito, verbose = TRUE, approx=FALSE)
meriEndoMito <- RunTSNE(meriEndoMito, verbose = TRUE)
meriEndoMito <- RunUMAP(meriEndoMito, dims = 1:30, verbose = TRUE)
meriEndoMito <- FindNeighbors(meriEndoMito, dims = 1:30, verbose = TRUE)
meriEndoMito <- FindClusters(meriEndoMito, verbose = TRUE, resolution = 1)

# DimPlot(meriEndoMito, reduction = "umap", label = TRUE, label.size = 7, pt.size = 0.5) # Vizualizing UMAP plot of meriEndoMito (mitochondria classes among meri/endo cells)

# Then, get list of dramtically-changing mitochondria genes that will be used to assess mitochondria state
MitoExp <- t(meriEndoMito@assays$SCT@scale.data)    # sclaed expression of mitochondria genes from meriEndoMito object
MitoExp <- MitoExp[row.names(MitoExp) %in% meriEndoPtime.scorpius$cell_ids,]     # only using cell names that are in the pseudotime ordering
# overall_pfeature_importances <- dynfeature::calculate_overall_feature_importance(meriEndoPtime.scorpius, expression_source= MitoExp) # basically calculating how important each gene contributes to shift from meristem to endodermis along pseudotime
MitoHeatmap <- dynplot::plot_heatmap(meriEndoPtime.scorpius, expression_source= MitoExp, features_oi = 50)
Mito_features.table  <- MitoHeatmap$assemble$plots[[4]]$data

#START HERE
a <- matrix(0, nrow = 30, ncol = 3)
for(i in 1:30){
  x <- Mito_features.table[Mito_features.table$feature_id == unique(Mito_features.table$feature_id)[i],]
  lm.x <- lm(expression  ~ percentage, data = x)
  a[i,1] <- unique(x$feature_id)
  a[i,2] <- coef(lm.x)[2]
  a[i,3] <- summary(lm.x)$coefficients[2,4]
}
a[order(a[,2]),]
a <- a[-a[,3] > 0.05,]
goi.Mito <-  a[,1] # goi = genes of interest

# Finally, looking at endodermis-mitochondria gene expression versus cell development (progression along pseudotime)
      # Average GE of goi
MitoScores <- AddModuleScore(object = meriEndoMito, features = goi.Mito, ctrl = 16, name = "Feature")
MitoScores <- MitoScores@meta.data[, grep("^Feature", colnames(MitoScores@meta.data))]
MitoScores$Avg <- rowMeans(MitoScores)
MitoScores$id <- row.names(MitoScores)
MitoScorePtime <- merge(scorpProg, MitoScores, by = "id")   # merging scores and progressions

# Plotting cell state vs mitochondria state
plot(MitoScorePtime$percentage, MitoScorePtime$Avg, main = "mitochondria-Type Gene Expression Along Pseudotime", ylab = "mitochondria-Type Gene Expression", xlab = "Progression Along Pseudotime", col = "navy blue", pch = 16)
abline(lm(MitoScorePtime$Avg ~ MitoScorePtime$percentage), lwd = 3, col = "blue", lty = 2)
