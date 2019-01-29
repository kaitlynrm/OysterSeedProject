#data
allsilos <- read.csv("Documents/robertslab/labnotebook/data/NSAF/ABACUS_output.tsv.csv")

NSAF.reps <- allsilos[, c(1, grep("ADJNSAF", names(allsilos)))]
colnames(NSAF.reps) <- substr(colnames(NSAF.reps), 18, 20)
colnames(NSAF.reps) <- sub("_.*", "", colnames(NSAF.reps))

#replicate averages
NSAF.reps$avg1 <-  rowMeans(NSAF.reps[c('1','1A')], na.rm=TRUE)
NSAF.reps$avg4 <-  rowMeans(NSAF.reps[c('4','4A')], na.rm=TRUE)
NSAF.reps$avg8 <-  rowMeans(NSAF.reps[c('8','8A')], na.rm=TRUE)
NSAF.reps$avg12 <- rowMeans(NSAF.reps[c('12','12A')], na.rm=TRUE)
NSAF.reps$avg16 <- rowMeans(NSAF.reps[c('16','16A')], na.rm=TRUE)
NSAF.reps$avg20 <- rowMeans(NSAF.reps[c('20','20A')], na.rm=TRUE)
NSAF.reps$avg24 <- rowMeans(NSAF.reps[c('24','24A')], na.rm=TRUE)
NSAF.reps$avg28 <- rowMeans(NSAF.reps[c('28','20A')], na.rm=TRUE)
NSAF.reps$avg32 <-  rowMeans(NSAF.reps[c('32','32A')], na.rm=TRUE)
NSAF.reps$avg36 <-  rowMeans(NSAF.reps[c('36','36A')], na.rm=TRUE)
NSAF.reps$avg40 <-  rowMeans(NSAF.reps[c('40','40A')], na.rm=TRUE)
NSAF.reps$avg44 <-  rowMeans(NSAF.reps[c('44','44A')], na.rm=TRUE)
NSAF.reps$avg48 <-  rowMeans(NSAF.reps[c('48','48A')], na.rm=TRUE)

#NSAF averages
NSAF.avg <- NSAF.reps[,c(1,46:ncol(NSAF.reps))]

#set rownames
rownames(NSAF.avg) <- NSAF.avg[,1]
NSAF.avg <- NSAF.avg[,-1]

write.csv(NSAF.avg, "Documents/robertslab/labnotebook/data/NSAF_avgs.csv")

#logtransform
NSAF.log <- log(NSAF.avg[,c(2:ncol(NSAF.avg))]+1)

#foldchanges
library(gtools)

NSAF.log$fc3 <- foldchange(NSAF.log$avg4, NSAF.avg$avg8)

NSAF.log$fc5 <- foldchange(NSAF.log$avg12, NSAF.avg$avg16)

NSAF.log$fc7 <- foldchange(NSAF.log$avg20, NSAF.avg$avg24)

NSAF.log$fc9 <- foldchange(NSAF.log$avg28, NSAF.avg$avg32)

NSAF.log$fc11 <- foldchange(NSAF.log$avg36, NSAF.avg$avg40)

NSAF.log$fc13 <- foldchange(NSAF.log$avg44, NSAF.avg$avg48)

NSAF.fc <- NSAF.log[,c(13:ncol(NSAF.log))]

#remove all rows containing Na, NaN and -/+Inf

NSAF.rem <- NSAF.fc[Reduce('&', lapply(NSAF.fc, function(NSAF.fc) !is.na(NSAF.fc) & is.finite(NSAF.fc))),]

write.csv(NSAF.rem, "Documents/robertslab/labnotebook/analysis/fold_change/NSAF_foldX.csv")

#make a heatmap
library(heatmap3)
heatmap3(as.matrix(NSAF.rem), method = "average", scale = "row", col=heat.colors(256), main="Fold Change Analysis on All Proteins", Colv=NA)



