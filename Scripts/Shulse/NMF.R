# Let’s say that we have 2 plastid clusters and 3 cell clusters. Let’s also assume 
# that plastid cluster 1 is more common with cell type 3. The 2nd plastid cluster 
# occurs equally in cell types 1 and 2. 

# The following matrix with 9 cells and two columns for plastid and cell clustering 
# are assumed. Note that plastid cluster 1 goes with cell type 3 as assumed. 

a <- matrix(
c(2,	1,
  2,	1,
  2,	1,
  2,	2,
  2,	2,
  2,	2,
  1,	3,
  1,	3,
  1,	3),
byrow = TRUE,
nrow = 9, ncol = 2)

# Now, let’s pretend that we don’t know which cell types go with which plastid types. 
# So we will just shuffle the cells so that we don’t see a pattern! Note that with 
# thousands of cells, it is not possible to just realize a pattern as in our 9 cells here!
a <- a[sample(1:nrow(a)),]  # randomize the matrix to lose the pattern 

# We then use NMF to extract pattern information:
install.packages("NMF")
library(NMF)

fct <- nmf(a, rank=2)
# get factor matrics w and h where a = w %*% h 
w <-  basis(fct)
h <- coef(fct)

# Reconstruct based on plastid types in column 1 of w and row 1 of h
x <- w[order(w[,1]),] %*% h[,order(h[,1])]
# This should be an ‘ordered’ reconstructed 'a' matrix that aligns plastid clusters 
# nicely against cell types.

# custom palette
my_palette <- colorRampPalette(c("white", "yellow", "red"))(n = 9)
# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(1,1.99,length=3),  # for red
               seq(2,2.99,length=3),              # for yellow
               seq(3,3.99,length=3),
               seq(4, 4.99))              # for green

# Run heat map on randomized matrix 'a'
heatmap.2(a,
          density.info="none",  
          trace="none",         
          margins =c(12,9),    
          #col=my_palette,       
          #breaks=col_breaks,   
          dendrogram='none',     
          Rowv=FALSE,
          Colv=FALSE) 

# Run heat map on reconstructed matrix 'x' 
# -- results are diff from(heatmap(a)
heatmap.2(x,
          density.info="none",  
          trace="none",         
          margins =c(12,9),    
          #col=my_palette,       
          #breaks=col_breaks,   
          dendrogram='none',     
          Rowv=FALSE,
          Colv=FALSE) 
# The heatmap shows that 
#  -- plastid type 1 was highly emphasized in cell type 3
#  -- plastid type 2 was vague between cell types 1 & 2
#  -- plastid type 1 did not exist in cell type 1 
#  -- plastid type 2 did not exist in cell type 3 
#  -- plastid type 1 was not supposed to be completely absent from cell type 2 



## Shulse data
## ============
x <- as.data.frame(Idents(d.obj)); names(x) <- "cells";
y <- as.data.frame(Idents(dp.obj)); names(y) <- "plastids";
x$genes <- row.names(x); y$genes <- row.names(y);
z <- merge(x, y); 
for(i in 2:3) z[,i] <- as.numeric(as.character(z[,i]))
a.mat <- as.matrix(z[,2:3])

# We then use NMF to extract pattern information:
# install.packages("NMF")
library(NMF)

fct <- nmf(a.mat+1, rank=2, .options='v')
# get factor matrics w and h where a = w %*% h 
w <-  basis(fct)
h <- coef(fct)

# Reconstruct based on plastid types in column 1 of w and row 1 of h
x <- round(w[order(w[,1]),] %*% h[,order(h[,1])])
# This should be an ‘ordered’ reconstructed 'a' matrix that aligns plastid clusters 
# nicely against cell types.

library(gplots)

# Run heat map on matched clusters without ordering
heatmap.2(a.mat,
          density.info="none",  
          trace="none",         
          margins =c(12,9),    
          #col=my_palette,       
          #breaks=col_breaks,   
          dendrogram='none',     
          Rowv=FALSE,
          Colv=FALSE) 

# Run heat map on reconstructed matrix 'x' 
# -- results are diff from(heatmap(a)
heatmap.2(x,
          density.info="none",  
          trace="none",         
          margins =c(12,9),    
          #col=my_palette,       
          #breaks=col_breaks,   
          dendrogram='none',     
          Rowv=FALSE,
          Colv=FALSE) 

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
  'Endodermis I-II',
  'Non-hair cells III',
  'Columella',
  'Cortex I',
  'Hair cells I',
  'Cortex IV',
  'Cortex II',
  'Non-hair cells II',
  'Hair cells II',
  'Cortex II',
  'Stele',
  'Endodermis I-II',
  'Non-hair cells I')
  
matchFreq$cellType <- cellType[matchFreq$Cell+1]
