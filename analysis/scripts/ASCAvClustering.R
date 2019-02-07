
#ASCA proteins
asca.prot <- read.csv("Downloads/ASCA_TempAffectedProteins.csv", stringsAsFactors = FALSE)
asca.prot <- t(asca.prot)
colnames(asca.prot) <- asca.prot[1,]
asca.prot <- asca.prot[-1,]
asca.prot <- data.frame(asca.prot)

asca.prot$ID <- rownames(asca.prot)

#Clustering proteins
hc.prot <- read.csv("Documents/robertslab/labnotebook/analysis/clustering/silo3_9-NSAF/techreps-avgs/silo3and9_nozerovals_AVGs/unq-prot-anno.csv", stringsAsFactors = FALSE)

silo3 <- hc.prot[which(hc.prot$Silo == 3),]
colnames(silo3) <- paste("T23", colnames(silo3), sep = "")

silo9 <- hc.prot[which(hc.prot$Silo == 9),]
colnames(silo9) <- paste("T29", colnames(silo9), sep = "")

s3s9 <- cbind(silo3,silo9)

rownames(s3s9) <- silo3$T23ID

s3s9 <- s3s9[,-c(1:3,11:14,22)]

colnames(s3s9)[1] <- "T16D0"

s3s9$ID <- rownames(s3s9)

#proteins that differ 
library(dplyr)
clus.difprot <- anti_join(s3s9, asca.prot, by="ID") #Clustering had 8 proteins not IDed by ASCA
asca.difprot <- anti_join(asca.prot, s3s9, by="ID") #ASCA had 113 proteins not IDed by cluster
                                                    #121 proteins do not match

#proteins identified by both ASCA and clustering
same.prot <- merge(s3s9, asca.prot, by="ID") #25 proteins
rownames(same.prot) <- same.prot[,1]
same.prot <- same.prot[,-c(1)]

same.prot <- same.prot[,-c(16:28,9)]

matrix <- as.matrix(same.prot)
class(matrix)

library(RColorBrewer)
hmcol<-brewer.pal(11,"RdBu")

library(heatmap3)
heatmap3(as.matrix(same.prot), method = "average", scale = "row", col=heat.colors(50), main="Proteins identified by ASCA and Clustering", Colv=NA)
