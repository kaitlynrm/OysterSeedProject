source("Documents/robertslab/labnotebook/analysis/clustering/biostats.R")

# https://www.bioconductor.org/install/
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install("limma")

library(limma)

#example
A <- runif(1000,4,16)
y <- A + matrix(rnorm(1000*3,sd=0.2),1000,3)
status <- rep(c(0,-1,1),c(950,40,10))
y[,1] <- y[,1] + status
plotMA(y, array=1, values=c(-1,1), hl.col=c("blue","red"))

#silo 3 vs silo 9 on day 13
ma.log <- NSAF.log[,c(1:12)]
ma.log$d3.mean <- rowMeans(ma.log[,c('avg4', 'avg8')], na.rm=TRUE)
ma.log$d3.fc <- foldchange(ma.log$avg4, ma.log$avg8)

ma.plot <- ma.log[,c(13:14)]



plotMA(ma.plot, array=1, hl.col=c("blue", "red"))
