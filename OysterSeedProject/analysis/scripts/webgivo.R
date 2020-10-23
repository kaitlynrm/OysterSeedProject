
all.proteins <- read.csv("silo3_9-all_proteins.csv")
colnames(all.proteins)[1] <- "Protein"
all.proteins$silo <- substr(all.proteins$Protein,1,1)
View(all.proteins)

#change class of protein column from factor to character
class(all.proteins$Protein)
all.proteins$Protein <- as.character(unlist(all.proteins$Protein))

#remove the silo identifier from the protein names
all.proteins$Protein <- substr(all.proteins$Protein, 3, nchar(all.proteins$Protein))
all.proteins <- all.proteins[,c(ncol(all.proteins),1:ncol(all.proteins)-1)]

annotations <- read.csv("../../../../data/allsilos-tag_and_annot.csv")
colnames(annotations)[1] <- "Protein"
annotations <- annotations[,-c(2:23)]
View(annotations)

merge <- merge(all.proteins, annotations, by = "Protein")
View(merge)
colnames(merge)[37] <- "Uniprot"

#move last column (Sample) to first column
#alternatively `which(colnames(data)=="Sample")` tells you which column it is
merge <- merge[,c(1:10,37)]
merge <- merge[,c(ncol(merge),1:ncol(merge)-1)]

deduped.data <- unique( yourdata[ , 1:3 ] )

merge1 <- merge[duplicated(merge$Protein), ]
final <- merge1[,1:2]

write.csv(final, "all_prot-uniprot.csv", row.names = FALSE)

#Webgivo test

setwd("Documents/robertslab/")
web <- read.csv("webgivo-test_inR.csv")
web2 <- web[,c(3,14:ncol(web))]
names(web2) <- c("Term", "G1","G2","G3","G4","G5","G6","G7", "G8", "G9", "G10")

library(tidyr)
web3 <- gather(web2, key = "Gene", "Term")

