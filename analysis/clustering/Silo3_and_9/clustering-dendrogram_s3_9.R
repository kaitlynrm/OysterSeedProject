###Kmeans clusering for multiple silos###

#Create new dataframe
#deliminate silo in protein name and combine silos
system("awk '$1=\"3_\"$1' /home/srlab/Documents/Kaitlyn/Github/OysterSeedProject/analysis/kmeans/silo3/silo3.csv >> /home/srlab/Documents/Kaitlyn/Github/OysterSeedProject/analysis/kmeans/Silo3_and_9/silo3_9.csv")
system("awk '$1=\"9_\"$1' /home/srlab/Documents/Kaitlyn/Github/OysterSeedProject/analysis/kmeans/silo9/silo9.csv >> /home/srlab/Documents/Kaitlyn/Github/OysterSeedProject/analysis/kmeans/Silo3_and_9/silo3_9.csv")

setwd("Documents/Kaitlyn/Github/OysterSeedProject/analysis/kmeans/Silo3_and_9/")

#Load in NSAF data
silo3.9 <- read.csv(file = "silo3_9.csv")
colnames(silo3.9)[colnames(silo3.9)=="X3_S3.Protein"] <- "Protein"

#remove silo header in datafame that was inserted into the middle of the dataframe from combining silos
silo3.9 <- silo3.9[-(which(grepl('9_S9-Protein', silo3.9$Protein))),]

#remove contaminant proteins (22 total)
proteins3_9<-subset(silo3.9, grepl(paste('CHOYP', collapse="|"), silo3.9$Protein))
which(grepl('ALBU_BOVIN', proteins3_9$Protein)) #ensure contaminants are gone

#Organize and save edited dataframe
rownames(proteins3_9)<-proteins3_9$Protein
silo3_9<-proteins3_9[,2:9]
names(silo3_9) <- c("0", "3", "5", "7", "9", "11", "13", "15")

#save silo3_9 without contaminant proteins
write.csv(silo3_9, file = "silo3_9-edited.csv")

##############################################################################################################################
###begin clustering###
silo3_9 <- read.csv("silo3_9-edited.csv", row.names = 1)
source("../biostats.R")

#use euclidean dissimilarity for clustering
library(vegan)
nsaf.euc<-vegdist(silo3_9, method='euclidean')

#average clustering method to cluster the data
library(cluster)
clust.avg<-hclust(nsaf.euc, method='average')

coef.hclust(clust.avg)
#agglomerate coefficent = 0.9959262 which means my proteins are more likely to be added into a new cluster

#cophenetic correlation
#how well cluster hierarchy represents original object-by-object dissimilarity space
cor(nsaf.euc, cophenetic(clust.avg))
#Above 0.75 is good for a cophenetic correlation; 0.6299225 for ward.D2, 0.9433488 for average

#Scree plot
hclus.scree(clust.avg)

jpeg(filename = "s3_9_scree.jpeg", width = 1000, height = 1000)
hclus.scree(clust.avg)
dev.off()

#Look for the elbow/inflection point on the scree plot and you can estimate number of clusters. But  it seems that this information cannot be pulled from the scree plot. (less than 500, maybe around 300?)

#cut dendrogram at selected height (example is given for 0.5) based on what looks reasonable because SCIENCE
plot(clust.avg)
rect.hclust(clust.avg, h=100)

jpeg(filename = "s3_9_dendrogram.jpeg", width = 1000, height = 1000)
plot(clust.avg)
rect.hclust(clust.avg, h=100)
dev.off()

#this looks reasonable
clust.class<-cutree(clust.avg, h=100)
max(clust.class) #41 clusters

#Cluster Freq table
silo3_9.freq <- data.frame(table(clust.class))

write.csv(silo3_9.freq, file = "silo3_9-freq.csv", row.names = FALSE)

#Make df
silo3_9.clus <- data.frame(clust.class)
names <- rownames(silo3_9.clus)
silo3_9.clus <- cbind(names, silo3_9.clus)
rownames(silo3_9.clus) <- NULL
colnames(silo3_9.clus)[1] <- "S3_9.Protein"
colnames(silo3_9.clus)[2] <- "Cluster"
silo3_9.norownames <- read.csv("silo3_9-edited.csv")

silo3_9.all <- merge(silo3_9.clus, silo3_9.norownames, by.x = "S3_9.Protein", by.y = "X")

#this gives matrix of 2 columns, first with proteins second with cluster assignment
#Line plots for each cluster
library(ggthemes)
library(reshape)
library(ggplot2)

melted_all_s3_9<-melt(silo3_9.all, id.vars=c('S3_9.Protein', 'Cluster'))

ggplot(melted_all_s3_9, aes(x=variable, y=value, group=S3_9.Protein)) +geom_line(alpha=0.1) + theme_bw() +
  facet_wrap(~Cluster, scales='free_y') + labs(x='Time Point', y='Normalized Spectral Abundance Factor')

jpeg(filename = "silo3_9clus_lineplots.jpeg", width = 1000, height = 1000)
ggplot(melted_all_s3_9, aes(x=variable, y=value, group=S3_9.Protein)) +geom_line(alpha=0.1) + theme_bw() +
  facet_wrap(~Cluster, scales='free_y') + labs(x='Time Point', y='Normalized Spectral Abundance Factor')
dev.off()

#Merge Silo clusters with Silo annotated and tagged datasheet
library(dplyr)

raw <- read.csv("../../../raw_data/os-allsilos-stats_and_annot.csv")
silo3_9.annotated <- raw[,-c(3, 6, 9, 12, 15, 18, 21)]
silo3_9.clus$Silo <- substr(silo3_9.clus$S3_9.Protein,1,1)
class(silo3_9.clus$S3_9.Protein)
silo3_9.clus$S3_9.Protein <- as.character(unlist(silo3_9.clus$S3_9.Protein))
silo3_9.clus$S3_9.Protein <- substr(silo3_9.clus$S3_9.Protein, 3, nchar(silo3_9.clus$S3_9.Protein))
merge <- merge(silo3_9.annotated, silo3_9.clus, by.x = "Protein.ID", by.y = "S3_9.Protein")

#reorganize
merge2 <- merge %>% select(Silo, everything())
merge3 <- merge2 %>% select(Cluster, everything())

write.csv(merge3, file = "silo3_9-anno_clus.csv", row.names = FALSE)

###Now I need to remove or subset proteins that are in the same cluster for both silos

#test with anti_join - didn't work
#s3.test <- filter(merge3, Silo == 3)
#s9.test <- filter(merge3, Silo == 9)

#s9.prot <- anti_join(s9.test, s3.test, by="Cluster")
#s3.prot <- anti_join(s3.test, s9.test, by="Cluster")

#test with duplicated
library(data.table)
unique.prot <- merge3[!(duplicated(merge3[c("Protein.ID", "Cluster")]) | duplicated(merge3[c("Protein.ID", "Cluster")], fromLast = TRUE)), ]

#determine if duplicates were removed
anyDuplicated(merge3[,c("Protein.ID","Cluster")]) #returns first duplicated rows, in this case columns 1 and 2 are duplicated so 2 is returned
anyDuplicated(test[,c("Protein.ID","Cluster")]) #returns 0 because there are no duplicates

write.csv(unique.prot, file = "unique-clus-prot-silo3_9.csv", row.names = FALSE)


