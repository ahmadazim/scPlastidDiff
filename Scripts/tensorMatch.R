##===========================================================================================================================
## Thoughts after seeing cluster matching...
#    - it's clear that NEPG can in fact be used to identify cell types because plastid types are "localizing" to cell types
#    - if cell types are reclustered with NEPG, do we get a blob on UMAP (all same type but just high resolution) or do we
#      we get n number of far away communities on plot
##===========================================================================================================================

# Preparing data sets (doublets included and normalized)
t(pWTdata.psub@assays$SCT@scale.data) -> plastGE
x <- data.frame(ids = names(Idents(pWTdata)), ct = as.character(Idents(pWTdata)))
x[x$ct == "Protoxylem", "ct"] <- "Xylem"
x[x$ct == "Protophloem", "ct"] <- "Phloem"
plastGE[row.names(plastGE) %in% x[,1],] -> plastGE
x[x$ids %in% row.names(plastGE),] -> x
plastGE <- as.data.frame(plastGE)
plastGE$ids <- row.names(plastGE)

plastML <- merge(x, plastGE)
plastML$ids <- NULL
write.csv(plastML, file = "./plastML.csv")


# Preparing data sets (doublets included and unnormalized)
t(pWTdata.psub@assays$RNA@data) -> unnorm
ident <- data.frame(ids = names(Idents(pWTdata)), ct = as.character(Idents(pWTdata)))
ident[ident$ct == "Protoxylem", "ct"] <- "Xylem"
ident[ident$ct == "Protophloem", "ct"] <- "Phloem"
unnorm[row.names(unnorm) %in% ident[,1],] -> unnorm
ident[ident$ids %in% row.names(unnorm),] -> ident
unnorm <- as.data.frame(unnorm)
unnorm$ids <- row.names(unnorm)

unnorm <- merge(ident, unnorm)
unnorm$ids <- NULL
write.csv(unnorm, file = "./unnormML.csv")


# Preparing data sets (doublets REMOVED and normalized)
t(NDpWTdata.psub@assays$SCT@scale.data) -> NDplastGE
x <- data.frame(ids = names(Idents(NDpWTdata)), ct = as.character(Idents(NDpWTdata)))
x[x$ct == "Protophloem", "ct"] <- "Phloem"
NDplastGE[row.names(NDplastGE) %in% x[,1],] -> NDplastGE
x[x$ids %in% row.names(NDplastGE),] -> x
NDplastGE <- as.data.frame(NDplastGE)
NDplastGE$ids <- row.names(NDplastGE)

NDplastML <- merge(x, NDplastGE)
NDplastML$ids <- NULL
write.csv(NDplastML, file = "./NDplastML.csv")




##===============================================================================================================
## Use Regression to classify plastid gene expression profiles into cell type classes with DOUBLETS INCLUDED!!!
## Read/Visualize/Section data into training, validation and test sets
##===============================================================================================================
plastGE$ids <- NULL
pca <- prcomp(plastGE, scale = TRUE)
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.per, main = "Scree Plot", xlab = "PCs", ylab = "Percent Variation")

lmData <- pca$x[,1:100]
lmData <- data.frame(lmData)
lmData$ids <- row.names(lmData)
lmData <- merge(lmData, x)
levels(lmData$ct)
# [1] "Cortex"         "Endodermis"     "Hair Cells"     "Meristem"
# [5] "Non Hair Cells" "Phloem"         "Protophloem"    "Protoxylem"
# [9] "Root Cap Cells" "Stele"          "Xylem"
lmData$ct <- as.numeric(lmData$ct)
trainPct <- 0.7
N <- nrow(lmData)
trainRows <- sample(1:N, size = round(trainPct*N), replace = F)
testRows <- (1:N)[!is.element(1:N, trainRows)]
train <- lmData[trainRows,]
test <- lmData[testRows,]

train.lm <- lm(ct ~ ., data = train[,-1])
summary(train.lm)

test$yhat <- predict(train.lm, newdata = test[,-1])
cor(test[,c("ct", "yhat")])
## 86% accurate!!!



##=================================================================================================================================
## Use TensorFlow to classify plastid gene expression profiles into cell type classes with DOUBLETS INCLUDED!!! and NORMALIZED!!!
## Read/Visualize/Section data into training, validation and test sets
##=================================================================================================================================
library(data.table)
plastData <- fread('./plastML.csv', header = T, sep = ',')
plastData[sample(1:nrow(plastData), 10), 1:5]  # look at 10 random rows
plastData$V1 <- NULL

## data set ratios
trnRatio <- 0.7; valRatio <- 0.2; tstRatio <- 0.1; N <- nrow(plastData)
trnN <- round(trnRatio * N); valN <- round(valRatio * N); tstN = N - trnN - valN;

rows <- 1:N
trnRows <- sample(rows, size=trnN); rows <- rows[-trnRows]  # values and vector indexes are the same
valRows <- sample(rows, size=valN); rows <- rows[!is.element(rows, valRows)]  # values are not the indexes
tstRows <- rows  # whatever remains is used as the test dataset --

write.csv(plastData[trnRows,], file = "./trnData.csv", row.names=F)
write.csv(plastData[valRows,], file = "./valData.csv", row.names=F)
write.csv(plastData[tstRows,], file = "./tstData.csv", row.names=F)




##============================================================================================================================
## Use TensorFlow to classify plastid gene expression profiles into cell type classes DOUBLETS INCLUDED!! and UNNORMALIZED!!!
## Read/Visualize/Section data into training, validation and test sets
##============================================================================================================================
library(data.table)
plastDataUN <- fread('./unnormML.csv', header = T, sep = ',')
plastDataUN[sample(1:nrow(plastDataUN), 10), 1:5]  # look at 10 random rows
plastDataUN$V1 <- NULL

## data set ratios
trnRatio <- 0.7; valRatio <- 0.2; tstRatio <- 0.1; Nun <- nrow(plastDataUN)
UNtrnN <- round(trnRatio * Nun); UNvalN <- round(valRatio * Nun); UNtstN = Nun - UNtrnN - UNvalN;

rowsUN <- 1:Nun
UNtrnRows <- sample(rowsUN, size=UNtrnN); rowsUN <- rowsUN[-UNtrnRows] # values and vector indexes are the same
UNvalRows <- sample(rowsUN, size=UNvalN); rowsUN <- rowsUN[!is.element(rowsUN, UNvalRows)] # values are not the indexes
UNtstRows <- rowsUN # whatever remains is used as the test dataset --

write.csv(plastDataUN[UNtrnRows,], file = "./trnDataUN.csv", row.names=F)
write.csv(plastDataUN[UNvalRows,], file = "./valDataUN.csv", row.names=F)
write.csv(plastDataUN[UNtstRows,], file = "./tstDataUN.csv", row.names=F)




##================================================================================================================================
## Use TensorFlow to classify plastid gene expression profiles into cell type classes with DOUBLETS REMOVED!!! and NORMALIZED!!!
## Read/Visualize/Section data into training, validation and test sets
##================================================================================================================================
library(data.table)
NDplastData <- fread('./NDplastML.csv', header = T, sep = ',')
NDplastData[sample(1:nrow(NDplastData), 10), 1:5]  # look at 10 random rows
NDplastData$V1 <- NULL

## data set ratios
trnRatio <- 0.7; valRatio <- 0.2; tstRatio <- 0.1; Nnd <- nrow(NDplastData)
NDtrnN <- round(trnRatio * Nnd); NDvalN <- round(valRatio * Nnd); NDtstN = Nnd - NDtrnN - NDvalN;

rowsND <- 1:Nnd
NDtrnRows <- sample(rowsND, size=NDtrnN); rowsND <- rowsND[-NDtrnRows] # values and vector indexes are the same
NDvalRows <- sample(rowsND, size=NDvalN); rowsND <- rowsND[!is.element(rowsND, NDvalRows)] # values are not the indexes
NDtstRows <- rowsND # whatever remains is used as the test dataset --

write.csv(NDplastData[NDtrnRows,], file = "/home/ahmad/TensorFlow/trnDataND.csv", row.names=F)
write.csv(NDplastData[NDvalRows,], file = "/home/ahmad/TensorFlow/valDataND.csv", row.names=F)
write.csv(NDplastData[NDtstRows,], file = "/home/ahmad/TensorFlow/tstDataND.csv", row.names=F)




##===========================================================================================================================================
## Use TensorFlow to classify cell gene expression profiles into cell type classes as validation (DOUBLETS REMOVED!!! and NORMALIZED!!!)
## Read/Visualize/Section data into training, validation and test sets
##===========================================================================================================================================
# Preparing data sets -- doublets REMOVED, normalized, cell type gene expression profiles (to be used for validation)
t(NDpWTdata@assays$SCT@scale.data) -> NDcellGE
ctype <- data.frame(ids = names(Idents(NDpWTdata)), ct = as.character(Idents(NDpWTdata)))
ctype[ctype$ct == "Protophloem", "ct"] <- "Phloem"
NDcellGE[row.names(NDcellGE) %in% ctype[,1],] -> NDcellGE
ctype[ctype$ids %in% row.names(NDcellGE),] -> ctype
NDcellGE <- as.data.frame(NDcellGE)
NDcellGE$ids <- row.names(NDcellGE)

NDcellML <- merge(ctype, NDcellGE, by = "ids")
NDcellML$ids <- NULL
write.csv(NDcellML, file = "./NDcellML.csv")

library(data.table)
NDcellData <- fread('./NDcellML.csv', header = T, sep = ',')
NDcellData[sample(1:nrow(NDcellData), 10), 1:5]  # look at 10 random rows
NDcellData$V1 <- NULL

## data set ratios
trnRatio <- 0.7; valRatio <- 0.2; tstRatio <- 0.1; cNnd <- nrow(NDcellData)
cNDtrnN <- round(trnRatio * cNnd); cNDvalN <- round(valRatio * cNnd); cNDtstN = cNnd - cNDtrnN - cNDvalN;

CrowsND <- 1:cNnd
cNDtrnRows <- sample(CrowsND, size=cNDtrnN); CrowsND <- CrowsND[-cNDtrnRows] # values and vector indexes are the same
cNDvalRows <- sample(CrowsND, size=cNDvalN); CrowsND <- CrowsND[!is.element(CrowsND, cNDvalRows)] # values are not the indexes
cNDtstRows <- CrowsND # whatever remains is used as the test dataset --

write.csv(NDcellData[cNDtrnRows,], file = "./trnDataNDcell.csv", row.names=F)
write.csv(NDcellData[cNDvalRows,], file = "./valDataNDcell.csv", row.names=F)
write.csv(NDcellData[cNDtstRows,], file = "./tstDataNDcell.csv", row.names=F)



##==============================================================================
## Using variable features to classify
NDpWTdata.lf <- FindVariableFeatures(NDpWTdata, nfeatures = 1000)
NDcellML.var <- NDcellML[,c(which(colnames(NDcellML) %in% NDpWTdata.lf@assays$SCT@var.features), 1)]
write.csv(NDcellML.var, file = "./NDcellML_var.csv")

## Using diffExp features to classify
NDcellML.diff <- NDcellGE[,c(colnames(NDcellGE) %in% NDpWTdata.markers$gene, 3001)]
row.names(NDcellML.diff) <- NULL
NDcellML.diff <- merge(NDcellML.diff, ctype, by = "ids")
NDcellML.diff$ids <- NULL
write.csv(NDcellML.diff, file = "./NDcellML_diff.csv")


library(data.table)
NDcellData <- fread('./NDcellML_diff.csv', header = T, sep = ',')
NDcellData[sample(1:nrow(NDcellData), 10), 1:5]  # look at 10 random rows
NDcellData$V1 <- NULL

## Data set ratios
trnRatio <- 0.7; valRatio <- 0.2; tstRatio <- 0.1; cNnd <- nrow(NDcellData)
cNDtrnN <- round(trnRatio * cNnd); cNDvalN <- round(valRatio * cNnd); cNDtstN = cNnd - cNDtrnN - cNDvalN;

CrowsND <- 1:cNnd
cNDtrnRows <- sample(CrowsND, size=cNDtrnN); CrowsND <- CrowsND[-cNDtrnRows] # values and vector indexes are the same
cNDvalRows <- sample(CrowsND, size=cNDvalN); CrowsND <- CrowsND[!is.element(CrowsND, cNDvalRows)] # values are not the indexes
cNDtstRows <- CrowsND # whatever remains is used as the test dataset --

write.csv(NDcellData[cNDtrnRows,], file = "./trnDataNDcell_diff.csv", row.names=F)
write.csv(NDcellData[cNDvalRows,], file = "./valDataNDcell_diff.csv", row.names=F)
write.csv(NDcellData[cNDtstRows,], file = "./tstDataNDcell_diff.csv", row.names=F)


## Trying to use regression to classify
NDcellGE$ids <- NULL
pca <- prcomp(NDcellGE, scale = TRUE)
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.per, main = "Scree Plot", xlab = "PCs", ylab = "Percent Variation")

lmData <- pca$x[,1:25]
lmData <- data.frame(lmData)
lmData$ids <- row.names(lmData)
lmData <- merge(lmData, ctype)
levels(lmData$ct)
# [1] "Cortex"         "Endodermis"     "Hair Cells"     "Meristem"
# [5] "Non Hair Cells" "Phloem"         "Protophloem"    "Root Cap Cells"
# [9] "Stele"          "Xylem"
lmData$ct <- as.numeric(lmData$ct)
trainPct <- 0.7
N <- nrow(lmData)
trainRows <- sample(1:N, size = round(trainPct*N), replace = F)
testRows <- (1:N)[!is.element(1:N, trainRows)]
train <- lmData[trainRows,]
test <- lmData[testRows,]

train.lm <- lm(ct ~ ., data = train[,-1])
summary(train.lm)

test$yhat <- predict(train.lm, newdata = test[,-1])
cor(test[,c("ct", "yhat")])
## 93% accurate based on first 25 PCs

## Using PCs in TensorFlow
lmData <- pca$x[,1:200]
lmData <- data.frame(lmData)
lmData$ids <- row.names(lmData)
lmData <- merge(lmData, ctype, by = "ids")
lmData$ids <- NULL
levels(lmData$ct)
# [1] "Cortex"         "Endodermis"     "Hair Cells"     "Meristem"
# [5] "Non Hair Cells" "Phloem"         "Protophloem"    "Root Cap Cells"
# [9] "Stele"          "Xylem"

trnRatio <- 0.7; valRatio <- 0.2; tstRatio <- 0.1; cNnd <- nrow(lmData)
cNDtrnN <- round(trnRatio * cNnd); cNDvalN <- round(valRatio * cNnd); cNDtstN = cNnd - cNDtrnN - cNDvalN;

CrowsND <- 1:cNnd
cNDtrnRows <- sample(CrowsND, size=cNDtrnN); CrowsND <- CrowsND[-cNDtrnRows] # values and vector indexes are the same
cNDvalRows <- sample(CrowsND, size=cNDvalN); CrowsND <- CrowsND[!is.element(CrowsND, cNDvalRows)] # values are not the indexes
cNDtstRows <- CrowsND # whatever remains is used as the test dataset --

write.csv(lmData[cNDtrnRows,], file = "./trnDataNDcell_PCs.csv", row.names=F)
write.csv(lmData[cNDvalRows,], file = "./valDataNDcell_PCs.csv", row.names=F)
write.csv(lmData[cNDtstRows,], file = "./tstDataNDcell_PCs.csv", row.names=F)



## Creating ND plastid data set with only one cell type: meristem
# NDplastData.meri <- NDplastData[NDplastData$ct == "Meristem",]
# write.csv(NDplastData.meri, file = "/home/ahmad/TensorFlow/meriPlastData.csv", row.names=F)
## Done in TensorFlow instead


## Running classifier on aberant cells types (doublets)
pWTdata.doublets <- subset(pWTdata.psub, cells = doublets, features = colnames(NDplastData)[2:1177])
t(pWTdata.doublets@assays$SCT@scale.data) -> doubletsPlast
x <- data.frame(ids = names(Idents(pWTdata)), ct = as.character(Idents(pWTdata)))
x[x$ct == "Protophloem", "ct"] <- "Phloem"
doubletsPlast[row.names(doubletsPlast) %in% x[,1],] -> doubletsPlast
x[x$ids %in% row.names(doubletsPlast),] -> x
doubletsPlast <- as.data.frame(doubletsPlast)
doubletsPlast$ids <- row.names(doubletsPlast)

plastML.doublets <- merge(x, doubletsPlast)
plastML.doublets$ids <- NULL
write.csv(plastML.doublets, file = "/home/ahmad/TensorFlow/doublets.csv")
