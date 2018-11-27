#Heatmap

#organize data
unq.prot <- read.csv("unique-clus-prot-silo3_9.csv")
unq.prot <- unq.prot[,c(1:11)]
unq.prot$Protein <-  paste(unq.prot$Silo, "_", unq.prot$Protein, sep = "")
colnames(unq.prot) <- c("Protein", "Cluster","Silo", "0", "3", "5", "7", "9", "11", "13", "15")
row.names(unq.prot) <- unq.prot$Protein
unq.prot <- unq.prot[,-c(1:3)]

#Note that all dendrograms are done with eucledian distance matricies

jpeg(filename = "heatmap-protclus.jpeg", width = 1000, height = 1000)
heatmap(as.matrix(unq.prot), scale = "column", col=heat.colors(256), main="Unique Proteins", Colv=NA)
dev.off()

jpeg(filename = "heatmap-dayclus.jpeg", width = 1000, height = 1000)
heatmap(as.matrix(unq.prot), scale = "column", col=heat.colors(256), main="Unique Proteins", Rowv=NA)
dev.off()

jpeg(filename = "heatmap-allclus.jpeg", width = 1000, height = 1000)
heatmap(as.matrix(unq.prot), scale = "column", col=heat.colors(256), main="Unique Proteins")
dev.off()

#Seperate silos
unq.prot2 <- read.csv("unique-clus-prot-silo3_9.csv")
unq.prot2 <- unq.prot2[,c(1:11)]
unq.prot2 <- unq.prot2[,-c(2)]
colnames(unq.prot2) <- c("Protein", "Silo", "0", "3", "5", "7", "9", "11", "13", "15")

#seperate silos
silo3 <- unq.prot2[ which(unq.prot2$Silo == '3'),]
silo9 <- unq.prot2[ which(unq.prot2$Silo == '9'),]

row.names(silo3) <- silo3$Protein
row.names(silo9) <- silo9$Protein

silo3 <- silo3[,-c(1:2)]
silo9 <- silo9[,-c(1:2)]

jpeg(filename = "heatmap-silo3.jpeg", width = 1000, height = 1000)
heatmap(as.matrix(silo3), scale = "column", col=heat.colors(256), main="Silo 3", Colv=NA)
dev.off()

jpeg(filename = "heatmap-silo9.jpeg", width = 1000, height = 1000)
heatmap(as.matrix(silo9), scale = "column", col=heat.colors(256), main="Silo 9", Colv=NA)
dev.off()

#Try to seperate silos in same heatmap- orginally I thought htis was a good idea but because we know that these proteins cluster differently,
#it seems quite redundant to put them together and cluster them again

#Reformat like for metboanalyst
metboanalyst <- read.csv("../../../MetboAnalyst/silo3_9-Metbo_format.csv")
library(tidyr)

test <- gather(unq.prot2, key = "Time", value = "Abundance", c(3:10))
test <- spread(test, "Protein", "Abundance")

class(test$Time)
test$Time <- as.numeric(test$Time)
test2 <- test[order(test$Time),]
test2$Sample <- paste("S", test2$Silo,"T", test2$Time, sep = "")
test3 <- test2[,c(ncol(test2),1:ncol(test2)-1)]

write.csv(test3, "../../../../../metbo-unqprot-test.csv", row.names = FALSE)

#enter data into metboanalyst

#copied code from metboanalyst

mSet<-InitDataObjects("pktable", "ts", FALSE)

mSet<-SetDesignType(mSet, "time")

mSet<-Read.TextData(mSet, "Replacing_with_your_file_path", "rowts", "disc");

mSet<-SanityCheckData(mSet)

mSet<-ReplaceMin(mSet);

mSet<-Normalization(mSet, "NULL", "NULL", "MeanCenter", ratio=FALSE, ratioNum=20)

mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)

mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)

mSet<-PlotHeatMap2(mSet, "heatmap2_0_", "png", 72, width=NA, "euclidean","ward.D","bwm","overview", F, 1, F, F)












