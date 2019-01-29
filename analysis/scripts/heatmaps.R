#Heatmap

#organize data
unq.prot <- read.csv("Documents/robertslab/labnotebook/analysis/clustering/NSAF-s3s9/travg-uniqueproteins.csv", row.names = 1)

#seperate silos
silo3 <- unq.prot[ which(unq.prot$Silo == '3'),]
silo9 <- unq.prot[ which(unq.prot$Silo == '9'),]

row.names(silo3) <- silo3$ID
row.names(silo9) <- silo9$ID

silo3 <- silo3[,-c(1:3)]
silo9 <- silo9[,-c(1:3)]

s3s9 <- cbind(silo3, silo9)

colnames(s3s9) <- c(paste("23", colnames(s3s9[,c(1:6)]), sep=""),  paste("29", colnames(s3s9[,c(7:12)]), sep=""))

#add back in day 0 for the heatmap
#load in and format original data
NSAF.filtered <- read.csv("Documents/robertslab/labnotebook/data/NSAF/silo3and9_nozerovals_AVGs.csv")
rownames(NSAF.filtered) <- paste(NSAF.filtered$day, NSAF.filtered$temp, sep="_")
NSAF.filtered <- NSAF.filtered[,-c(1:3)]

NSAF.trans <- data.frame(t(NSAF.filtered))

#remove contaminant proteins
NSAF.trans <- subset(NSAF.trans, grepl(paste('CHOYP', collapse="|"), rownames(NSAF.trans)))
which(grepl('ALBU_BOVIN', rownames(NSAF.trans))) #ensure contaminants are gone

#set gene names as row names
annotations <- read.csv("Documents/robertslab/labnotebook/data/allsilos-tag_and_annot.csv")
annotations <- annotations[,c(1,67)]
annotations$Protein.ID <- sub("\\|", ".", annotations$Protein.ID)

merge <- merge(NSAF.trans, annotations, by.x = "row.names", by.y="Protein.ID", all.x=TRUE)
merge$Gene.names <- as.character(unlist(merge$Gene.names))
class(merge$Gene.names)
merge$Names <- ifelse(merge$Gene.names == "", yes=merge$Row.names, no= merge$Gene.names)
merge$Names <- ifelse(merge$Names == "None", yes=merge$Row.names, no= merge$Names)
merge$Names <- ifelse(is.na(merge$Names) == TRUE, yes=merge$Row.names, no=merge$Names)
any(is.na(merge$Names))

day0 <-merge[,c(16,2)]

#merge day 0
s3s9$Names <- rownames(s3s9)

library(dplyr)
s3s9 <- merge(day0, s3s9, by.x="Names", by.y="Names") 
rownames(s3s9) <- s3s9$Names
s3s9 <- s3s9[,-c(1)]
colnames(s3s9)[1] <- '16D0'

library(RColorBrewer)
hmcol<-brewer.pal(11,"RdBu")

library(heatmap3)
jpeg(filename = "Documents/robertslab/labnotebook/analysis/clustering/NSAF-s3s9/travg-day0-silo9-unqprot-heatmap.jpeg", width = 1500, height = 1000)
heatmap3(as.matrix(s3s9), method = "average", scale = "row", col=heat.colors(50), main="Proteins from both temperatures that differentially cluster", Colv=NA)
dev.off()

##############################################################################

#Note that all dendrograms are done with eucledian distance matricies
library(heatmap3)
heatmap3(as.matrix(unq.prot), method = "average", scale = "row", col=heat.colors(256), main="Proteins in Silo 3 and 9 that differentially hierarchically cluster", Colv=NA)


### below are indvidual silo heatmaps and metbo reformatting  ##

#proteins clustered
jpeg(filename = "Documents/robertslab/labnotebook/analysis/clustering/NSAF-s3s9/50-silo9-unqprot-heatmap.jpeg", width = 1000, height = 1000)
heatmap(as.matrix(silo9), scale = "row", col=heat.colors(256), main="Unique Proteins", Colv=NA)
dev.off()

#days clsutered
jpeg(filename = "heatmap-dayclus.jpeg", width = 1000, height = 1000)
heatmap(as.matrix(unq.prot), scale = "column", col=heat.colors(256), main="Unique Proteins", Rowv=NA)
dev.off()

#no clustering
jpeg(filename = "heatmap-allclus.jpeg", width = 1000, height = 1000)
heatmap(as.matrix(unq.prot), scale = "column", col=heat.colors(256), main="Unique Proteins")
dev.off()

#Seperate silos
unq.prot2 <- read.csv("no0-unique-clus-prot-silo3_9.csv")
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

#silo 3/23C
jpeg(filename = "heatmap-silo3.jpeg", width = 1500, height = 1000)
heatmap(as.matrix(silo3), scale = "column", col=heat.colors(256), main="Silo 3", Colv=NA, Rowv=NA)
dev.off()

#silo 9/29C
jpeg(filename = "heatmap-silo9.jpeg", width = 1500, height = 1000)
heatmap(as.matrix(silo9), scale = "column", col=heat.colors(256), main="Silo 9", Colv=NA, Rowv=NA)
dev.off()






#Try to seperate silos in same heatmap- orginally I thought this was a good idea but because we know that these proteins cluster differently,
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





