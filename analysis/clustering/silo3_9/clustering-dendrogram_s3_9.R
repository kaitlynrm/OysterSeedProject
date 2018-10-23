
### Hierarchical clusering
#REVIEW the DOCUMENT to find your starting place first!!!

#Option 1: create a dataframe to include only detected proteins in analysis
#Option 2: a dataframe with only detected proteins was already created so it just needs to be loaded in first
#Option 3: create a dataframe to include all proteins (even undetected) in analysis
#Option 4: a dataframe with all proteins was already created and just needs to be loaded in 


#################################################################   OPTION 1    #############################################################

#Create new from seperated silos 
#this approach will not include any proteins that were never detected
#deliminate silo in protein name and combine silos with bash
system("awk '$1=\"3_\"$1' /home/srlab/Documents/Kaitlyn/Github/OysterSeedProject/analysis/clustering/silo3/silo3.csv >> /home/srlab/Documents/Kaitlyn/Github/OysterSeedProject/analysis/clustering/silo3_9/silo3_9.csv")
system("awk '$1=\"9_\"$1' /home/srlab/Documents/Kaitlyn/Github/OysterSeedProject/analysis/clustering/silo9/silo9.csv >> /home/srlab/Documents/Kaitlyn/Github/OysterSeedProject/analysis/clustering/silo3_9/silo3_9.csv")

setwd("Documents/Kaitlyn/Github/OysterSeedProject/analysis/clustering/silo3_9/")

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

source("../biostats.R")

##############################################################################################################################


#################################################################   OPTION 2    #############################################################
#Load in data made in bash
setwd("Documents/Kaitlyn/Github/OysterSeedProject/analysis/clustering/silo3_9/")
silo3_9 <- read.csv("silo3_9-edited.csv", row.names = 1)
colnames(silo3_9) <- c(0, 3, 5, 7, 9, 11, 13, 15)

source("../biostats.R")

##############################################################################################################################


#################################################################   OPTION 3    #############################################################

#Create dataframe from raw file instead of merging two seperated files 
#this approach will include all proteins detected in either silo even if it was only in 1 silo
setwd("Documents/robertslab/labnotebook/analysis/clustering/silo3_9/all_proteins")
all.silos <- read.csv("../../../../data/ABACUSdata_only.csv")
colnames(all.silos)[1] <- "Protein"
head(all.silos)

#seperate silos
silo3 <- all.silos[,c(1,2,4,7,10,13,16,19,22)]

silo3$Protein <- sub("^", "3_", silo3$Protein)
colnames(silo3) <- c("Protein", "0", "3", "5", "7", "9", "11", "13", "15")

#load in all dected proteins in silo 3 without nonabundant proteins
silo3.detected <- read.csv("../../silo3/silo3.csv")
str(silo3)
colnames(silo3.detected) <- c("Protein", "0", "3", "5", "7", "9", "11", "13", "15")
silo3.detected$Protein <- sub("^", "3_", silo3.detected$Protein)

#confirm different row lengths
nrow(silo3) #8443
nrow(silo3.detected) #7345
8443-7345 #1098 undetected proteins

#look at list of undected proteins
library(dplyr)
silo3.undetected <- anti_join(silo3, silo3.detected, by = "Protein")
nrow(silo3.undetected) #confirming it is still 1098
View(silo3.undetected)

#so all proteins are included in silo3 (1098 proteins were not detected but will be added to the cluster analysis)

#do the same for silo9
silo9 <- all.silos[,c(1,2,5,8,11,14,17,20,23)]
silo9$Protein <- sub("^", "9_", silo9$Protein)
colnames(silo9) <- c("Protein", "0", "3", "5", "7", "9", "11", "13", "15")

#combine both dataframes
silo3.9 <- rbind(silo3, silo9)

#remove contaminant proteins (22 total)
silo3_9<-subset(silo3.9, grepl(paste('CHOYP', collapse="|"), silo3.9$Protein))
which(grepl('ALBU_BOVIN', silo3_9$Protein)) #ensure contaminants are gone

write.csv(silo3_9, "silo3_9-allproteins.csv")
silo3_9 <- read.csv("silo3_9-all_proteins.csv", row.names = 1)

colnames(silo3_9) <- c("0", "3", "5", "7", "9", "11", "13", "15")
source("../../biostats.R")

##############################################################################################################################


#################################################################   OPTION 4    #############################################################

setwd("Documents/robertslab/labnotebook/analysis/clustering/silo3_9/all_proteins")
silo3_9 <- read.csv("silo3_9-all_proteins.csv", row.names = 1)
colnames(silo3_9) <- c("0", "3", "5", "7", "9", "11", "13", "15")

source("../../biostats.R")

##############################################################################################################################


#Data is loaded in as silo3_9 with proteins as row names and biostats is loaded
#Now start cluster analysis!

#use euclidean dissimilarity for clustering
library(vegan)
nsaf.euc<-vegdist(silo3_9, method='euclidean')

#average clustering method to cluster the data
library(cluster)
clust.avg<-hclust(nsaf.euc, method='average')

coef.hclust(clust.avg)
#agglomerate coefficent = (0.9964979 all) (0.9959262 detected) which means my proteins are more likely to be added into a new cluster

#cophenetic correlation
#how well cluster hierarchy represents original object-by-object dissimilarity space
cor(nsaf.euc, cophenetic(clust.avg))
#Above 0.75 is good for a cophenetic correlation; (0.6299225 for ward.D2, 0.9433488 for average for detected) (0.)

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
silo <- substr(melted_all_s3_9$Protein, 0, 1)

jpeg(filename = "bycolour-silo3_9clus_lineplots.jpeg", width = 1000, height = 1000)

ggplot(melted_all_s3_9, aes(x=variable, y=value, group=Protein, color = silo)) +geom_line(alpha=0.8) + theme_bw() +
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
write.csv(silo3_9.annot, file = "silo3_9-clus-annot.csv", row.names = FALSE)

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
silo <- substr(unq_melted_all_s3_9$Protein, 0, 1)

jpeg(filename = "bycolour-silo3_9-unq_lineplots.jpeg", width = 1000, height = 1000)

ggplot(unq_melted_all_s3_9, aes(x=variable, y=value, group=Protein, color = silo)) +geom_line(alpha=0.8) + theme_bw() +
  facet_wrap(~Cluster, scales='free_y') + labs(x='Time Point', y='Normalized Spectral Abundance Factor')

dev.off()
