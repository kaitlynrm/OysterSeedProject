
#load in data
s3s9 <- read.csv("Documents/robertslab/labnotebook/data/NSAF/silo3and9_nozerovals_AVGs.csv")
rownames(s3s9) <- paste("D",s3s9$day,"T",s3s9$temp, sep="")
s3s9 <- t(s3s9)
s3s9 <- s3s9[-c(1:2),]
s3s9 <- s3s9[,-c(1)]

#seperate silos and remove undetected proteins from each silo
silo3 <- data.frame(s3s9[, c(grepl("23", colnames(s3s9)))])
silo3 <- silo3[rowSums(silo3 > 0.6),]

silo9 <- data.frame(s3s9[, c(grepl("29", colnames(s3s9)))])
silo9 <- silo3[rowSums(silo9 > 0.6),]

silo3$Names <- rownames(silo3)
silo9$Names <- rownames(silo9)

s3Names <- silo3$Names
s9Names <- silo9$Names

ven <- cbind(silo3$Names, silo9$Names)
colnames(ven) <- c("Silo3", "Silo9")

#venn diagram
library(gplots)
jpeg(filename = "Documents/robertslab/labnotebook/analysis/venn-diagram/s3s9-quickvenn.jpeg", width = 1500, height = 1000)
venn(list("23°C"=s3Names, "29°C"=s9Names))
dev.off()

library(VennDiagram)
venn <- data.frame(ven)
write.csv(venn,"Documents/robertslab/labnotebook/analysis/venn-diagram/temp-prots.csv")

class(venn)
class(venn$Silo3)

t23 <- list(venn$Silo3)
t29 <- list(venn$Silo9)
venn.plot <- draw.pairwise.venn(list("23°C"=t23, "29°C"=t29), fill=c("darkmagenta", "darkblue"), filename="Documents/robertslab/labnotebook/analysis/venn-diagram/venn-diagram.jpeg")




