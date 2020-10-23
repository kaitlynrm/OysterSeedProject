#Heatmaps
#https://www.nceas.ucsb.edu/~frazier/RSpatialGuides/colorPaletteCheatsheet.pdf <- colorpalettes

########################## ASCA ######################################
#Reformatting ASCA proteins for cluster extraction from heatmap
asca.prot <- read.csv("Downloads/ASCA_TempAffectedProteins.csv", stringsAsFactors = FALSE)
rownames(asca.prot) <- asca.prot[,1]
asca.prot <- asca.prot[,-c(1)]
asca.prot <- t(asca.prot)
asca.prot <- data.frame(asca.prot)

names(asca.prot) <- gsub(x = colnames(asca.prot), pattern = "X", replacement = "D")
names(asca.prot) <- gsub(x = colnames(asca.prot), pattern = "\\_", replacement = "T")

#seperate silos
silo3 <- unq.prot[,grepl("23", colnames(unq.prot))]
silo9 <- unq.prot[,grepl("29", colnames(unq.prot))]
s3s9 <- cbind(silo3, silo9)
class(s3s9$D3T23) #columns are factor but should be numeric

################################################################

############################# HClust ###################################
#organize data
unq.prot <- read.csv("Documents/robertslab/labnotebook/analysis/clustering/silo3_9-NSAF/techreps-avgs/silo3and9_nozerovals_AVGs/unq-prot-anno.csv")

#seperate silos
silo3 <- unq.prot[ which(unq.prot$Silo == '3'),]
silo9 <- unq.prot[ which(unq.prot$Silo == '9'),]

row.names(silo3) <- silo3$ID
row.names(silo9) <- silo9$ID

silo3 <- silo3[,-c(1:3,11)]
silo9 <- silo9[,-c(1:4,11)]

s3s9 <- cbind(silo3, silo9)

colnames(s3s9) <- c("16D0", paste("23", colnames(s3s9[,c(2:7)]), sep=""),  paste("29", colnames(s3s9[,c(8:13)]), sep=""))
################################################################

###########################  Set gene names as row names  #####################################
#load in and format original data
NSAF.filtered <- read.csv("Documents/robertslab/labnotebook/data/NSAF/silo3and9_nozerovals_AVGs.csv")
rownames(NSAF.filtered) <- paste(NSAF.filtered$day, NSAF.filtered$temp, sep="_")
NSAF.filtered <- NSAF.filtered[,-c(1:3)]

NSAF.trans <- data.frame(t(NSAF.filtered))

#remove contaminant proteins
NSAF.trans <- subset(NSAF.trans, grepl(paste('CHOYP', collapse="|"), rownames(NSAF.trans)))
which(grepl('ALBU_BOVIN', rownames(NSAF.trans))) #ensure contaminants are gone

#merge tables
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
###############################

#heatmap3
library(RColorBrewer)
hmcol<-brewer.pal(11,"RdBu")

library(heatmap3)
jpeg(filename = "Documents/robertslab/labnotebook/analysis/clustering/silo3_9-NSAF/techreps-avgs/silo3and9_nozerovals_AVGs/heatmap.jpeg", width = 1500, height = 1000)
hm <- heatmap3(as.matrix(s3s9), hclustfun=hclust, distfun = function(x) as.dist(1 - cor(t(x), use = "pa")), method = "average", scale = "row", col=heat.colors(50), main="Proteins from both temperatures that differentially cluster", Colv=NA)
dev.off()

### Extract clades from heatmap
#Make dendrogram from heatmap
hm <- as.dist(1-cor(t(as.matrix(s3s9)), use="pa"))
hm2 <- hclust(hm, method = 'average')
plot(hm2)

jpeg(filename= "Documents/robertslab/labnotebook/analysis/clustering/silo3_9-NSAF/techreps-avgs/silo3and9_nozerovals_AVGs/heatmap-dendrogram.jpeg")
plot(hm2)
dev.off()

#Identify clades in heatmap
# define some clusters
mycl <- cutree(hm2, h=max(hm2$height/1.1))

# get a color palette equal to the number of clusters
clusterCols <- colorRamps::primary.colors(length(unique(mycl)))

# create vector of colors for dendrogram
myClusterSideBar <- clusterCols[mycl]

#get colors for temperature treaments
silos <- as.integer(c(1,2,2,2,2,2,2,3,3,3,3,3,3))

#create vector of colors for treatments
siloCols <- rainbow(length(unique(silos)))

mySiloSidebar <- siloCols[silos]

# choose a color palette for the heat map
myheatcol <- rev(redgreen(75))

# draw the heat map
library(colorRamps)
rdtbl <- blue2red(15)

library(gplots)
heatmap.2(as.matrix(s3s9), distfun = function(x) as.dist(1 - cor(t(x))), main="Differentially Clustered Proteins", Rowv=as.dendrogram(hm2), Colv=NA, dendrogram="row", scale="row", col=rdtbl, density.info="none", trace="none", RowSideColors= myClusterSideBar, ColSideColors = mySiloSidebar)
coords <- locator(1) #select on plot where you want treatment legend; (0.8002414, 0.9834021)
coords2 <- locator(1) #select on plot where you want cluster legend; (0.6019484, 0.9807591)

mySiloSidebar # c("#FF0000FF", "#00FF00FF", "#0000FFFF")
myClusterSideBar # c("#000000", #80FFFF")

jpeg(filename = "Documents/robertslab/labnotebook/analysis/clustering/silo3_9-NSAF/techreps-avgs/silo3and9_nozerovals_AVGs/clades-silos-heatmap.jpeg", width = 612, height = 555)
heatmap.2(as.matrix(s3s9), distfun = function(x) as.dist(1 - cor(t(x))), main="Differentially Clustered Proteins", Rowv=as.dendrogram(hm2), Colv=NA, dendrogram="row", scale="row", col=rdtbl, density.info="none", trace="none", RowSideColors= myClusterSideBar, ColSideColors = mySiloSidebar)
legend(yjust=0.8,box.lty = 0, coords, inset=-0.00000001, legend = c("16C", "23C", "29C"), col=c("#FF0000FF", "#00FF00FF", "#0000FFFF"), pch=15)
legend(yjust=0.8,box.lty = 0, coords2, inset=-0.00000001, legend = c("Cluster 1", "Cluster 2"), col=c("#000000", "#80FFFF"), pch=15)

dev.off()

# cutree returns a vector of cluster membership
# in the order of the original data rows
# examine it
mycl

# examine the cluster membership by it's order
# in the heatmap
mycl[hm2$order]

# grab a cluster
cluster1 <- s3s9[mycl == 1,]
cluster2 <- s3s9[mycl == 2,]

# or add the cluster ID to your data
foo <- cbind(s3s9, clusterID=mycl)

# examine the data with cluster ids attached, and ordered like the heat map
clade.prot <- foo[hm2$order,]

#df w/ averages for each temp and cluster
t23c1 <- rowMeans(cluster1[,c(2:7)])
t29c1 <- rowMeans(cluster1[,c(8:13)])
meanscluster1 <- data.frame(cbind(t23c1,t29c1))
meanscluster1 <- rbind(colMeans(meanscluster1), meanscluster1)
rownames(meanscluster1)[1] <- "Mean Protein Abundance"


t23c2 <- rowMeans(cluster2[,c(2:7)])
t29c2 <- rowMeans(cluster2[,c(8:13)])
meanscluster2 <- data.frame(cbind(t23c2,t29c2))
meanscluster2 <- rbind(colMeans(meanscluster2), meanscluster2)
rownames(meanscluster2)[1] <- "Mean Protein Abundance"

meanprotabun <- cbind(meanscluster1[1,c(1:2)],meanscluster2[1,c(1:2)])
View(meanprotabun)

