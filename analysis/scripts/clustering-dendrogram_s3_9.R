
### Hierarchical clustering

#Option 1: All NSAF values and detected proteins
#Option 2: Filtered NSAF values
#Optional: remove day 0 and day 15

#####################################  OPTION 1: All NSAF values and detected proteins  #################################################################################

NSAF.avg <- read.csv("Documents/robertslab/labnotebook/data/NSAF/")

##############################################################################################################################

#####################################  OPTION 2: Filtered NSAF values  #################################################################################

NSAF.filtered <- read.csv("Documents/robertslab/labnotebook/data/NSAF/silo3and9_nozerovals_AVGs.csv")
rownames(NSAF.filtered) <- paste(NSAF.filtered$day, NSAF.filtered$temp, sep="_")
NSAF.filtered <- NSAF.filtered[,-c(1:3)]

NSAF.trans <- data.frame(t(NSAF.filtered))

#remove contaminant proteins
NSAF.trans <- subset(NSAF.trans, grepl(paste('CHOYP', collapse="|"), rownames(NSAF.trans)))
which(grepl('ALBU_BOVIN', rownames(NSAF.trans))) #ensure contaminants are gone
NSAF.trans <- NSAF.trans[,-c(1)]

#####################################  For genes as names  ####################################################################

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

#seperate silos (for merge file)
s3 <- merge[,c(15,2,4,6,8,10,12)]
s9 <- merge[,c(15,3,5,7,9,11,13)]
#########################################  For Proteins as names  ################################################################

#seperate silos (for NSAF.trans file)
NSAF.trans$Names <- rownames(NSAF.trans)
silo3 <- NSAF.trans[,c(13,1,3,5,7,9,11)]
silo9 <- NSAF.trans[,c(13,2,4,6,8,10,12)]

####################################################################################
#add silo/temp prefix into protein/gene names
silo3$Names <- sub("^", "3_", silo3$Names)
silo9$Names <- sub("^", "9_", silo9$Names)

colnames(silo3) <- c("ID", "D3", "D5", "D7", "D9", "D11", "D13")
colnames(silo9) <- c("ID", "D3", "D5", "D7", "D9", "D11", "D13")

#combine dataframes
silo3.9 <- rbind(silo3,silo9)

# do below lines if genes are names
rownames(silo3.9) <-  make.names(silo3.9$ID, unique=TRUE)
rownames(silo3.9) <- sub("X", "", rownames(silo3.9))

#do below if proteins are names
rownames(silo3.9) <- silo3.9$ID

silo3.9 <- silo3.9[,-c(1)]

#biostats package
source("Documents/robertslab/labnotebook/analysis/clustering/biostats.R")

##############################################################################################################################

#Data is loaded in as silo3_9 with proteins as row names and biostats is loaded
#Now start cluster analysis!

#dissimilarity matrix and clustering
library(vegan)
nsaf.euc<-vegdist(silo3.9, method='euclidean')

library(cluster)
clust.avg<-hclust(nsaf.euc, method='average')

#agglomerate coefficent
coef.hclust(clust.avg) # euclidean = 0.9991113; bray = 0.9542416; euclidean(avg tech reps) = 0.9975317

#cophenetic correlation (want at least 0.75)
cor(nsaf.euc, cophenetic(clust.avg)) # euclidean = 0.9635755; bray = 0.7826124; euclidean(avg tech reps) = 0.9635755

#Scree plot
jpeg(filename = "Documents/robertslab/labnotebook/analysis/clustering/NSAF-s3s9/travg-scree-plot.jpeg", width = 1000, height = 1000)
hclus.scree(clust.avg)
dev.off()

#Look for the elbow/inflection point on the scree plot and you can estimate number of clusters. 
#But  it seems that this information cannot be pulled from the scree plot.

#cut dendrogram at selected height
plot(clust.avg)
rect.hclust(clust.avg, h=250)


jpeg(filename = "Documents/robertslab/labnotebook/analysis/clustering/NSAF-s3s9/travg-dendrogram.jpeg", width = 1000, height = 1000)
plot(clust.avg)
rect.hclust(clust.avg, h=250)
dev.off()

#this looks reasonable
clust.class<-cutree(clust.avg, h=250)
max(clust.class) #33 clusters at 50 height; bray = 22 at 0.5; 28 cluster travg

#Cluster Freq table
silo3_9.freq <- data.frame(table(clust.class))

write.csv(silo3_9.freq, file = "Documents/robertslab/labnotebook/analysis/clustering/NSAF-s3s9/travg-s3s9-freq.csv", row.names = FALSE)

#Make df
silo3_9.clus <- data.frame(clust.class)
names <- rownames(silo3_9.clus)
silo3_9.clus <- cbind(names, silo3_9.clus)
rownames(silo3_9.clus) <- NULL
colnames(silo3_9.clus)[1] <- "ID"
colnames(silo3_9.clus)[2] <- "Cluster"

#merge for abundance values
silo3_9.all <- merge(silo3_9.clus, silo3.9, by.x = "ID", by.y = "row.names")

#this gives matrix of 2 columns, first with proteins second with cluster assignment
#Line plots for each cluster
library(ggthemes)
library(reshape)
library(ggplot2)

melted_all_s3_9<-melt(silo3_9.all, id.vars=c('ID', 'Cluster'))

silo <- substr(melted_all_s3_9$ID, 0, 1)

jpeg(filename = "Documents/robertslab/labnotebook/analysis/clustering/NSAF-s3s9/travg-250-unqabundanceplot.jpeg", width = 1000, height = 1000)

ggplot(melted_all_s3_9, aes(x=variable, y=value, group=ID, color=silo)) +geom_line(alpha=0.8) + theme_bw() +
  facet_wrap(~Cluster, scales='free_y') + labs(x='Time Point', y='Normalized Spectral Abundance Factor')

dev.off()

#Seperate Silo from protein name
silo3_9.edit <- silo3_9.all
silo3_9.edit$Silo <- substr(silo3_9.edit$ID,1,1)

#Remove protein silo notation and organize
class(silo3_9.edit$ID)
silo3_9.edit$ID <- as.character(unlist(silo3_9.edit$ID))
silo3_9.edit$ID <- substr(silo3_9.edit$ID, 3, nchar(silo3_9.edit$ID))

library(dplyr)
silo3_9.edit <- silo3_9.edit %>% select(Silo, everything())
silo3_9.edit <- silo3_9.edit %>% select(Cluster, everything())

#save datasheet with abundance values, clusters, and annotations
#silo3_9.annot <- merge(silo3_9.edit, annotations, by.x = "Protein", by.y = "Protein.ID")
#write.csv(silo3_9.annot, file = "no0-silo3_9-clus-annot.csv", row.names = FALSE)

#Parse out unique proteins
library(data.table)
unique.prot <- silo3_9.edit[!(duplicated(silo3_9.edit[c("ID", "Cluster")]) | duplicated(silo3_9.edit[c("ID", "Cluster")], fromLast = TRUE)), ]

#determine if duplicates were removed
anyDuplicated(silo3_9.edit[,c("ID","Cluster")]) #returns first duplicated rows
anyDuplicated(unique.prot[,c("ID","Cluster")]) #returns 0 because there are no duplicates

#Save datasheet with unique proteins and annotations
#final.unique.prot <- merge(unique.prot, annotations, by.x = "ID", by.y = "Protein.ID")

#sum(final.unique.prot$Silo == "3")
#sum(final.unique.prot$Silo == "9")

#write.csv(final.unique.prot, file = "no0-unique-clus-prot-silo3_9.csv", row.names = FALSE)

write.csv(unique.prot, "Documents/robertslab/labnotebook/analysis/clustering/NSAF-s3s9/travg-protid--uniqueproteins.csv")

#Remove annotations
#s39.unq.abudance <- final.unique.prot[,-c(12:64)]

#Now I need to have only one column for each day rather than one column per silo per day
#protein.names <-  data.frame(paste(s39.unq.abudance$Silo, "_", s39.unq.abudance$Protein, sep = ""))
#colnames(protein.names) <- "Protein"
#plot.unq.prot <-  silo3_9.all[which(silo3_9.all$Protein %in% protein.names$Protein),]

#Plot abudances of unique proteins
plot.unq <- unique.prot #5632 compared to 9180 = 3548 proteins eliminated
plot.unq$S.ID <- paste(plot.unq$Silo, plot.unq$ID, sep="_")
plot.unq <- plot.unq[,-c(2:3)]

unq_melted_all_s3_9<-melt(plot.unq, id.vars=c('S.ID', 'Cluster'))
silo <- substr(unq_melted_all_s3_9$S.ID, 0, 1)

jpeg(filename = "Documents/robertslab/labnotebook/analysis/clustering/NSAF-s3s9/travg250-unique-abundanceplots.jpeg", width = 1000, height = 1000)

ggplot(unq_melted_all_s3_9, aes(x=variable, y=value, group=S.ID, color = silo)) +geom_line(alpha=0.8) + theme_bw() +
  facet_wrap(~Cluster, scales='free_y') + labs(x='Time Point', y='Normalized Spectral Abundance Factor')

dev.off()

###########################

annotations <- read.csv("Documents/robertslab/labnotebook/data/allsilos-tag_and_annot.csv")
annotations <- annotations[,c(1,63)]
annotations$Protein.ID <- sub("\\|", ".", annotations$Protein.ID)

unq.anno <- merge(unique.prot, annotations, by.x="ID", by.y="Protein.ID")

write.csv(unq.anno, "Documents/robertslab/labnotebook/analysis/clustering/NSAF-s3s9/unqprot-travg-anno.csv")

library(dplyr)
count(unq.anno, Entry)

