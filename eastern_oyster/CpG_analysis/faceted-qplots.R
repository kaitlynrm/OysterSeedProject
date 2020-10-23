### CpG O/E Eastern oyster
cpg <- read.csv("Documents/robertslab/work/CpG/eastern-oyster/ID_CpG_labelled_all.csv")

### single q-q plot
data <- cpg$CL_1
filt.data <- data[data >= 0.001 & data <= 1.5] #Cutting off high and low values
set.seed(101)
mixmdl <- normalmixEM(filt.data, k=2)
plot(mixmdl, which = 2, col2 = c("red", "blue"), xlab2 = " ", main2 = "CL_1", font.main = 3)
lines(density(filt.data), lty=2, lwd=2)

### faceted q-q plots
library(mixtools)

colname.list <- colnames(cpg)

jpeg(filename = "Documents/robertslab/work/faceted_cpg.jpeg", width = 1000, height = 1000)

a <- par(mfrow = c(10, 9)) # r x c plot

for(i in 2:length(colname.list)){
data <- cpg[[i]]
filt.data <- data[data >= 0.001 & data <= 1.5] #Cutting off high and low values
set.seed(101)
mixmdl <- normalmixEM(filt.data, k=2)
plot(mixmdl, which = 2, col2 = c("red", "blue"), xlab2 = " ", main2 = colname.list[[1]][i], font.main = 3)
lines(density(filt.data), lty=2, lwd=2)
print(mean(cpg[[i]]))
}

dev.off()

### pull out stress response genes
GOslim <- read.csv("Documents/robertslab/work/CpG/eastern-oyster/Blastquery-GOslim.csv", header=FALSE)
GOslim$V1 <- gsub(".*?NC","NC", GOslim$V1) #remove any artifacts from excel

stress.GOslim <- GOslim[grep("stress", GOslim$V3), ]
write.csv(merge(stress.GOslim, cpg, by.x="V1", by.y="ID"), "Documents/robertslab/work/CpG/eastern-oyster/stress-annotated.csv")

BPspread <- aggregate(V3~V1,GOslim,FUN=paste)
GOspread <- aggregate(V2~V1,GOslim,FUN=paste)
spread.GOslim <- merge(BPspread, GOspread, by="V1")
class(spread.GOslim$V3)

stress <- spread.GOslim[grep("stress", spread.GOslim$V3), ]

stress.cpg <- merge(stress, cpg, by.x = "V1", by.y="ID")

### scatter plots with mean, median and range of values in original cpg table
data.frame(ID=DF[,1], Means=rowMeans(DF[,-1]))

rownames(cpg) <- cpg$ID
  cpg <- cpg[,-c(1)]
  
#cpg$means <- rowMeans(cpg)
library(matrixStats)
#cpg$median <- rowMedians(as.matrix(cpg))
#cpg$range <- rowRanges(as.matrix(cpg))
#library(Rfast)
#cpg$cov <- rowcvs(cpg, ln = FALSE, unbiased = FALSE)
#cpg <- subset(cpg, select=-c(means,median,range))

stats <- data.frame(colMeans(cpg))
colnames(stats)[1] <- "mean"
stats$median <- colMedians(as.matrix(cpg))
#stats$range <- colRanges(as.matrix(cpg)) #values too small
#stats <- subset(stats, select=-c(range))
stats$skew <- colskewness(as.matrix(cpg))
stats$kurtosis <- colkurtosis(as.matrix(cpg))
stats$variance <- colVars(as.matrix(cpg))
stats$stdev <- colVars(as.matrix(cpg), std = TRUE)

stats$ID <- rownames(stats)

a <- ggplot(stats, aes(x=ID, y=mean)) + geom_point()
b <- ggplot(stats, aes(x=ID, y=median)) + geom_point()
c <- ggplot(stats, aes(x=ID, y=kurtosis)) + geom_point()
d <- ggplot(stats, aes(x=ID, y=skew)) + geom_point()
e <- ggplot(stats, aes(x=ID, y=variance)) + geom_point()
f <- ggplot(stats, aes(x=ID, y=stdev)) + geom_point()

library(gridExtra)
jpeg("Documents/robertslab/work/CpG/eastern-oyster/basicstats.jpeg")
grid.arrange(a,b,c,d,e,f, nrow = 3)
dev.off()

### stress response stats
rownames(stress.cpg) <- stress.cpg$V1
  stress.cpg <- stress.cpg[,-c(1:3)]

stress.stats <- data.frame(colMeans(stress.cpg))
colnames(stress.stats)[1] <- "mean"
stress.stats$median <- colMedians(as.matrix(stress.cpg))
stress.stats$skew <- colskewness(as.matrix(stress.cpg))
stress.stats$kurtosis <- colkurtosis(as.matrix(stress.cpg))
stress.stats$variance <- colVars(as.matrix(stress.cpg))
stress.stats$stdev <- colVars(as.matrix(stress.cpg), std = TRUE)

stress.stats$ID <- rownames(stress.stats)

a <- ggplot(stress.stats, aes(x=ID, y=mean)) + geom_point()
b <- ggplot(stress.stats, aes(x=ID, y=median)) + geom_point()
c <- ggplot(stress.stats, aes(x=ID, y=kurtosis)) + geom_point()
d <- ggplot(stress.stats, aes(x=ID, y=skew)) + geom_point()
e <- ggplot(stress.stats, aes(x=ID, y=variance)) + geom_point()
f <- ggplot(stress.stats, aes(x=ID, y=stdev)) + geom_point()

jpeg("Documents/robertslab/work/CpG/eastern-oyster/stress-basicstats.jpeg")
grid.arrange(a,b,c,d,e,f, nrow = 3)
dev.off()

