#Heatmap

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

###################################
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
###############################

library(RColorBrewer)
hmcol<-brewer.pal(11,"RdBu")

library(heatmap3)
jpeg(filename = "Documents/robertslab/labnotebook/analysis/clustering/silo3_9-NSAF/techreps-avgs/silo3and9_nozerovals_AVGs/heatmap.jpeg", width = 1500, height = 1000)
heatmap3(as.matrix(s3s9), method = "average", scale = "row", col=heat.colors(50), main="Proteins from both temperatures that differentially cluster", Colv=NA)
dev.off()

# Extract clades from heatmap








