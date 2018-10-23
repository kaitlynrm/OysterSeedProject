### Hierarchical clusering for multiple silos###

#create dataframe from raw file instead of merging two seperated files
#






###############################################################################################################################
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


### Start here if dataframe combining silos has already been created! ###



silo3_9 <- read.csv("silo3_9-edited.csv", row.names = 1)
colnames(silo3_9) <- c(0, 3, 5, 7, 9, 11, 13, 15)
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
colnames(silo3_9.all) <- c("Protein", "Cluster", "0", "3", "5", "7", "9", "11", "13", "15")

#this gives matrix of 2 columns, first with proteins second with cluster assignment
#Line plots for each cluster
library(ggthemes)
library(reshape)
library(ggplot2)

melted_all_s3_9<-melt(silo3_9.all, id.vars=c('Protein', 'Cluster'))

jpeg(filename = "silo3_9clus_lineplots.jpeg", width = 1000, height = 1000)
ggplot(melted_all_s3_9, aes(x=variable, y=value, group=Protein)) +geom_line(alpha=0.8) + theme_bw() +
  facet_wrap(~Cluster, scales='free_y') + labs(x='Time Point', y='Normalized Spectral Abundance Factor')
dev.off()

#Line plot with color for each silo

jpeg(filename = "bycolour-silo3_9clus_lineplots.jpeg", width = 1000, height = 1000)

ggplot(melted_all_s3_9, aes(x=variable, y=value, group=Protein, color = substr(melted_all_s3_9$Protein, 0, 1))) +geom_line(alpha=0.8) + theme_bw() +
         facet_wrap(~Cluster, scales='free_y') + labs(x='Time Point', y='Normalized Spectral Abundance Factor')

dev.off()

#Seperate Silo from protein name
silo3_9.edit <- silo3_9.all
silo3_9.edit$Silo <- substr(silo3_9.edit$Protein,1,1)

#Remove protein silo notation and organize
class(silo3_9.edit$Protein)
silo3_9.edit$Protein <- as.character(unlist(silo3_9.edit$Protein))
silo3_9.edit$Protein <- substr(silo3_9.edit$Protein, 3, nchar(silo3_9.edit$Protein))
silo3_9.edit <- silo3_9.edit %>% select(Silo, everything())
silo3_9.edit <- silo3_9.edit %>% select(Cluster, everything())

#Obtain annotated Silo3_9 sheet with clusters
all.annotated <- read.csv("../../../raw_data/os-allsilos-stats_and_annot.csv")
annotations <- all.annotated[,-c(2:23)]

silo3_9.annot <- merge(silo3_9.edit, annotations, by.x = "Protein", by.y = "Protein.ID")
write.csv(final.unique.prot, file = "silo3_9-clus-annot.csv", row.names = FALSE)

#Parse out unique proteins
library(data.table)
unique.prot <- silo3_9.edit[!(duplicated(silo3_9.edit[c("Protein", "Cluster")]) | duplicated(silo3_9.edit[c("Protein", "Cluster")], fromLast = TRUE)), ]

#determine if duplicates were removed
anyDuplicated(silo3_9.edit[,c("Protein","Cluster")]) #returns first duplicated rows
anyDuplicated(unique.prot[,c("Protein","Cluster")]) #returns 0 because there are no duplicates

#Annotate unique proteins
final.unique.prot <- merge(unique.prot, annotations, by.x = "Protein", by.y = "Protein.ID")

write.csv(final.unique.prot, file = "unique-clus-prot-silo3_9.csv", row.names = FALSE)

#Remove annotations
s39.unq.abudance <- final.unique.prot[,-c(12:64)]

#Now I need to have only one column for each day rather than one column per silo per day
protein.names <-  data.frame(paste(s39.unq.abudance$Silo, "_", s39.unq.abudance$Protein, sep = ""))
colnames(protein.names) <- "Protein"
plot.unq.prot <-  silo3_9.all[which(silo3_9.all$Protein %in% protein.names$Protein),]

#Plot abudances of unique proteins
unq_melted_all_s3_9<-melt(plot.unq.prot, id.vars=c('Protein', 'Cluster'))

jpeg(filename = "bycolour-silo3_9-unq_lineplots.jpeg", width = 1000, height = 1000)

ggplot(unq_melted_all_s3_9, aes(x=variable, y=value, group=Protein, color = substr(unq_melted_all_s3_9$Protein, 0, 1))) +geom_line(alpha=0.8) + theme_bw() +
  facet_wrap(~Cluster, scales='free_y') + labs(x='Time Point', y='Normalized Spectral Abundance Factor')

dev.off()
