
setwd("Documents/robertslab/labnotebook/analysis/clustering/silo3_9/")

#load in tab delimited file
bp.direct <- read.delim("all-proteins/enrichment/bp-direct.txt", header = TRUE, sep = "\t", dec = ".")

#split columns
library(tidyr)
split <- separate(bp.direct, 'Term', c("Accession", "Term"), sep="~")
split2 <- separate(split, 'Genes', paste("Genes", 1:3, sep="_"), sep=",", extra="drop")

#just terms and accession codes
terms <- split2[,c(3,7:9)]

#gather genes
gathered <- gather(terms, key="number", value="genes", c(2:4))
gathered <- gathered[,c(1,3)]

#load in unique proteins
unq.prot <- read.csv("all-proteins/unique-clus-prot-silo3_9.csv")
unq.acces <- unq.prot[,c(1,3:11,38)]
colnames(unq.acces)[11] <- "genes"

#get terms for unique proteins
merge <- merge(unq.acces, gathered, by = "genes") 
sort <- merge[order(merge$Term),]

n <- length(sort)

#silo9
silo9 <- sort[seq(1, n, 2),]
rownames(silo9) <- silo9[,12]
silo9 <- silo9[,c(4:11)]
colnames(silo9) <- c("0", "3", "5", "7", "9", "11", "13","15")

#silo3
silo3 <- sort[seq(2,n,2),]
rownames(silo3) <- silo3[,12]
silo3 <- silo3[,c(4:11)]
colnames(silo3) <- c("0", "3", "5", "7", "9", "11", "13","15")

#heatmap
setwd("all-proteins/enrichment/")

#silo9
jpeg(filename = "heatmap-silo9.jpeg", width = 1400, height = 550)
heatmap(as.matrix(silo9), scale = "column", col=heat.colors(256), main="Silo 9 ", Rowv=NA, Colv=NA)
dev.off()

jpeg(filename = "heatmap-silo9-clus.jpeg", width = 1400, height = 550)
heatmap(as.matrix(silo9), scale = "column", col=heat.colors(256), main="Silo 9 ", Rowv=NA)
dev.off()


#silo3
jpeg(filename = "heatmap-silo3.jpeg", width = 1400, height = 550)
heatmap(as.matrix(silo3), scale = "column", col=heat.colors(256), main="Silo 3 ", Rowv=NA, Colv=NA)
dev.off()

jpeg(filename = "heatmap-silo3-clus.jpeg", width = 1400, height = 550)
heatmap(as.matrix(silo3), scale = "column", col=heat.colors(256), main="Silo 3 ", Rowv=NA)
dev.off()






