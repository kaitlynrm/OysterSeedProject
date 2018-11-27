#Update R (do on Rgui)
install.packages("installr"); 
library(installr)
updateR()

#Set WD
getwd()
setwd("C:/Users/Katty/Documents/robertslab/labnotebook/")

#Install and load packages
install.packages("vegan")
library(vegan)
install.packages("raster")
library(raster)
install.packages("BioStatR")
library(BioStatR)
source("analysis/nmds_R/biostats.R")

########Create NMDS plot with silo 3 and 9################

#upload file
ABACUSdata <- read.csv("~/Documents/robertslab/labnotebook/data/ABACUSdata_only.csv", header=TRUE)
View(ABACUSdata)

#Remove NA
ABACUSdata[is.na(ABACUSdata)] <- 0

#Transpose- switch rows and columns
tABACUSdata <- t(ABACUSdata)
View(tABACUSdata)

#Rename Columns and remove row
colnames(tABACUSdata) <- tABACUSdata[1,]
tABACUSdata = tABACUSdata[-1,]

#Remove Silo 2
silo3and9 <- tABACUSdata[-(seq(from = 2, to = 22, by = 3)), ]

#Convert to numeric
is.numeric(silo3and9)
num.silo3and9 <- apply(silo3and9,c(1,2),as.numeric)
is.numeric(num.silo3and9)

#Make MDS dissimilarity matrix
nmds.silo3and9 <- metaMDS(num.silo3and9, distance = 'euclidean', k = 2, trymax = 3000, autotransform = FALSE)

#Save as jpeg
jpeg(filename = "analysis/nmds_R/nmdss3ands9-color.jpeg", width = 670, height = 650)
fig <- ordiplot(nmds.silo3and9, choices=c(1,2),type="none", display="sites", xlab='Axis 1', ylab='Axis 2', cex=0.5)
points(nmds.silo3and9, "sites", col=c(rep('black',1), rep('red',2), rep('orange',2), rep('palegreen',2), rep('springgreen4',2),
                                      rep('lightskyblue1',2), rep('blue',2),rep('purple',2)),
                                pch=c(rep(19,1), rep(17,1), rep(15,1), rep(17,1), rep(15,1), rep(17,1), rep(15,1), 
                                      rep(17,1), rep(15,1), rep(17,1), rep(15,1), rep(17,1), rep(15,1), rep(17,1), rep(15,1)))
legend("topright", legend=c("pool", "23C-Silo3", "29C-Silo9"), pch=c(19,17,15))
legend("topleft", legend=c("Day 0", "Day 3", "Day 5", "Day 7", "Day 9", "Day 11", "Day 13", "Day 15"), 
        col=c('black', 'red', 'orange', 'palegreen','springgreen4','lightskyblue1','blue','purple'), pch=19)
dev.off()



