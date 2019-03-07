library(dplyr)
library(tidyr)
library(MetStaT)

### Reformat Silo3_9 dataframe ###

#read in Oyster Temp. protein data frame
data <- read.csv("kmeans/Silo3_and_9/silo3_9-edited.csv")
names(data) <- c("Protein", "0", "3", "5", "7", "9", "11", "13", "15")

#create silo column
data$Silo <- substr(data$Protein,1,1)

#change class of protein column from factor to character
class(data$Protein)
data$Protein <- as.character(unlist(data$Protein))

#remove the silo identifier from the protein names
data$Protein <- substr(data$Protein, 3, nchar(data$Protein))

#Change table orientation so that abundances are listed in one column and time points are listed down a column

#`gather()` makes wide datasets long- if columns are values of a variable
#Apply to timepoint columns that have abundance values (keep protein and silo unchanged)
#Time values in c(2:9) will become multiple rows along with a new abundance column
data <- gather(data, key = "Time", value = "abundance", c(2:9))

#`spread()` makes long datasets wide- if there are multiple rows for an observation
#list proteins across the top (note that proteins not present in one silo will become NA values)
data <- spread(data, "Protein", "abundance")
any(is.na(data))
data[is.na(data)] <- 0
any(is.na(data))

#change silo column class to numeric
class(data$Silo)
data$Silo <- as.numeric(data$Silo)

#change time column class to numeric; exclude the X with the substr command
data$Time <- as.numeric(substr(data$Time,1,2))

#format and export table for metaboanalysthttps://stats.stackexchange.com/questions/22329/how-does-centering-the-data-get-rid-of-the-intercept-in-regression-and-pca
data$Sample <- paste("S",data$Silo,"T",data$Time, sep = "")

#move last column (Sample) to first column
#alternatively `which(colnames(data)=="Sample")` tells you which column it is
data <- data[,c(ncol(data),1:ncol(data)-1)]
data <- data[order(data$Time),]

write.csv(data, "MetboAnalyst/silo3_9-Metbo_format.csv", row.names = FALSE, quote = FALSE)

#find max abundance value 
max(data[,4:ncol(data)]) #379

#create matrix to pass to ASCA command as data (abundances are observations)
ASCAX <- as.matrix(data[,-c(1:3)])
#create matrix to pass to ASCA command as levels (time and temp are factors)
ASCAF <- as.matrix(data[,2:3])

### Run ASCA command ###
setwd("Documents/Kaitlyn/Github/OysterSeedProject/analysis/MetboAnalyst/")
ASCA <- ASCA.Calculate(ASCAX, ASCAF, equation.elements = "1,2,12")

## plot loadings of the first two principal components of the first factor (temp.)
ASCA.PlotLoadings(ASCA, ee = "1", pcs="1,2")

jpeg(filename = "temp-loadings-plot.jpeg", width = 1000, height = 1000)
ASCA.PlotLoadings(ASCA, ee = "1", pcs="1,2")
dev.off()

## plot loadings of the first two principal components of the second factor (time)
ASCA.PlotLoadings(ASCA, ee = "2", pcs="1,2")

jpeg(filename = "time-loadings-plot.jpeg", width = 1000, height = 1000)
ASCA.PlotLoadings(ASCA, ee = "2", pcs="1,2")
dev.off()

## plot loadings of the first two principal components of the factor interaction
ASCA.PlotLoadings(ASCA, ee = "12", pcs="1,2")

jpeg(filename = "interaction-loadings-plot.jpeg", width = 1000, height = 1000)
ASCA.PlotLoadings(ASCA, ee = "12", pcs="1,2")
dev.off()

#plot single score plot for the first factor on the first two PCs
ASCA.PlotScores(ASCA, ee = "1", PCs = "1,2")

jpeg(filename = "temp-scores-plot.jpeg", width = 1000, height = 1000)
ASCA.PlotScores(ASCA, ee = "1", PCs = "1,2")
dev.off()

#plot single score plot for second factor on first two PCs
ASCA.PlotScores(ASCA, ee = "2", PCs = "1,2")

jpeg(filename = "time-scores-plot.jpeg", width = 1000, height = 1000)
ASCA.PlotScores(ASCA, ee = "2", PCs = "1,2")
dev.off()

#plot single score plot for interaction on first two PCs
ASCA.PlotScores(ASCA, ee = "12", PCs = "1,2")

jpeg(filename = "interaction-scores-plot.jpeg", width = 1000, height = 1000)
ASCA.PlotScores(ASCA, ee = "12", PCs = "1,2")
dev.off()

## plot scores for the first two PCs and the projections of the data for the second factor (time)
ASCA.PlotScoresPerLevel(ASCA, ee = "2", pcs = "1,2")

jpeg(filename = "time-projected_scores-plot.jpeg", width = 1000, height = 1000)
ASCA.PlotScoresPerLevel(ASCA, ee = "2", pcs = "1,2")
dev.off()

## the data for the first factor (i.e. temp.)
ASCA.PlotScoresPerLevel(ASCA, ee = "1", pcs = "1,2")

jpeg(filename = "temp-projected_scores-plot.jpeg", width = 1000, height = 1000)
ASCA.PlotScoresPerLevel(ASCA, ee = "1", pcs = "1,2")
dev.off()

## the data for the interaction of factors (i.e. temp. and time interactive effect)
ASCA.PlotScoresPerLevel(ASCA, ee = "12", pcs = "1,2")

jpeg(filename = "interaction-projected_scores-plot.jpeg", width = 1000, height = 1000)
ASCA.PlotScoresPerLevel(ASCA, ee = "12", pcs = "1,2")
dev.off()

## Do a permutation test to evaluate the significance to the two factors and the interaction.
ASCA.DoPermutationTest(ASCA, perm=1000)

#output
#1, 2, 12
#0.119, 0.015, 1

#Examine the loading values
loadings <-  as.data.frame(ASCA$`1`$svd$v[,1])
