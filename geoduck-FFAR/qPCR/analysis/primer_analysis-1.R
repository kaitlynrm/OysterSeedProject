# data
pooled <- read.csv("Desktop/20200210qPCR_data.csv", stringsAsFactors = FALSE)
library(plyr)
pooled$Target <- revalue(pooled$Target, c("NFIP"="NFIP1"))

run1 <- read.csv("Desktop/20200218-knownsamples_run1.csv", stringsAsFactors = FALSE)
run2 <- read.csv("Desktop/20200218-knownsamples_run2 (1).csv", stringsAsFactors = FALSE)

#dfs with target, sample, Cq.mean, Melt temp, Melt height
apooled <- pooled[,c(2,5,7,9,10)]
arun1 <- run1[,c(2,5,8,10,11)]
arun2 <- run2[,c(2,5,8,10,11)]

#subset dfs based on targets
  #only prints or saves as list
  #test <- split(apooled, with(apooled, apooled$`Target-pooled`), drop = TRUE)

for(i in unique(apooled$`Target`)) {
  nam <- paste("pooled", i, sep = ".")
  assign(nam, apooled[apooled$`Target`==i,])
}

for(i in unique(arun1$`Target`)) {
  nam <- paste("run1", i, sep = ".")
  assign(nam, arun1[arun1$`Target`==i,])
}

for(i in unique(arun2$`Target`)) {
  nam <- paste("run2", i, sep = ".")
  assign(nam, arun2[arun2$`Target`==i,])
}

#bind same primer dfs
  #pooled = 17 primers (all)
  #runs = 6 primers (APLP, GSK3B, NFIP1, RPL5, SPTN1, TIF3s6B)
APLP <- rbind(pooled.APLP, run1.APLP)
GSK3B <- rbind(pooled.GSK3B, run1.GSK3B, run2.GSK3B)
NFIP1 <- rbind(pooled.NFIP1, run1.NFIP1)
RPL5 <- rbind(pooled.RPL5, run2.RPL5)
SPTN1 <- rbind(pooled.SPTN1, run2.SPTN1)
TIF3s6B <- rbind(pooled.TIF3s6B, run2.TIF3s6B)

#df with all primers
primers <- rbind(APLP,GSK3B,NFIP1,RPL5,SPTN1,TIF3s6B, pooled.ECHD3,pooled.FEN1,pooled.GLYG,pooled.NSF,pooled.TIF3s10,pooled.TIF3s12,pooled.TIF3s4a,pooled.TIF3s7,`pooled.TIF3s8-1`,`pooled.TIF3s8-2`,pooled.TIF3sF)

#plots
  # 17 total (6 have 30 observastions. 11 have 2 obsv.)
library(ggplot2)
library(ggbeeswarm)

ggplot(primers, aes(Target, Cq.Mean, col = Sample)) + geom_point()

ggplot(primers, aes(Target, Melt.Temperature, col = Sample)) + geom_point()

ggplot(primers, aes(Target, Peak.Height, col = Sample)) + geom_point()

#create M/F df and merge with primer/sample df
Sample <- c('19','21','27','28','31','37','39','43','54','55','57','59','61','Pooled')
sex <- c('F','F','M','M','F','F','F','M','M','F','F','M','F','Pooled')
stage <- c('2','2','2','1','4','6','5','4','3','3','7','4','7','Pooled')
sss <- data.frame(Sample,sex,stage)
View(sss)

prsxst <- merge(primers,sss, by = 'Sample', all=FALSE)
View(prsxst)

all.prsxst <- merge(primers,sss, by = 'Sample', all=FALSE)
View(all.prsxst)

#plots
ggplot(all.prsxst, aes(Target, Cq.Mean, col = sex)) + geom_point() 

ggplot(all.prsxst, aes(Target, Cq.Mean, col = stage)) + geom_point()

ggplot(all.prsxst, aes(sex, Cq.Mean, col = Target)) + geom_point()

is.numeric(all.prsxst$Cq.Mean)
all.prsxst$Cq.Mean <- as.numeric(as.character(all.prsxst$Cq.Mean))

ggplot(all.prsxst, aes(x=Target, y=Cq.Mean, color=sex)) + geom_boxplot(aes(x=Target, y=Cq.Mean, color=sex), width = 0.5, alpha = 0.1)


library(ggbeeswarm)
ggplot(all.prsxst, aes(sex, Cq.Mean, col = Target)) + geom_beeswarm(size = 1, cex = 1) + geom_boxplot(aes(sex, Cq.Mean), width = 0.5, alpha = 0.1)

#dataframe w/ primers tested on knowns (APLP, GSK3B, NFIP1, RPL5, SPTN1, TIF3s6B)
A <-  all.prsxst[all.prsxst$Target == "APLP",]
B <-  all.prsxst[all.prsxst$Target == "GSK3B",]
C <-  all.prsxst[all.prsxst$Target == "NFIP1",]
D <-  all.prsxst[all.prsxst$Target == "RPL5",]
E <-  all.prsxst[all.prsxst$Target == "SPTN1",]
F <-  all.prsxst[all.prsxst$Target == "TIF3s6B",]

pp <- rbind(A,B,C,D,E,F)
View(pp)

pp$Cq.Mean <- as.numeric(as.character(pp$Cq.Mean))

ggplot(pp, aes(sex, Cq.Mean, col = Target)) + geom_beeswarm(size = 3, cex = 3)

ggplot(pp, aes(x=Target, y=Cq.Mean, group=sex)) +  geom_point(aes(shape=sex, color=stage))

ggplot(pp, aes(x=Target, y=Cq.Mean, color=sex)) + geom_boxplot(aes(x=Target, y=Cq.Mean, color=sex), width = 0.5, alpha = 0.1)

is.numeric(pp$Melt.Temperature)
pp$Melt.Temperature <- as.numeric(as.character(pp$Melt.Temperature)) #NAs coerced
pp$Melt.Temperature[is.na(pp$Melt.Temperature)] <- 0

ggplot(pp, aes(x=Target, y=Melt.Temperature, color=stage)) + geom_boxplot(aes(x=Target, y=Melt.Temperature, color=stage), width = 0.5, alpha = 0.1)

is.numeric(pp$Peak.Height)
pp$Peak.Height <- as.numeric(as.character(pp$Peak.Height)) #NAs coerced

ggplot(pp, aes(x=Target, y=Peak.Height, color=stage)) + geom_boxplot(aes(x=Target, y=Peak.Height, color=stage), width = 0.5, alpha = 0.1)

pp$Peak.Height[is.na(pp$Peak.Height)] <- 0

ggplot(pp, aes(x=Target, y=Melt.Temperature, color=sex)) + geom_boxplot(aes(x=Target, y=Melt.Temperature, color=stage), width = 0.5, alpha = 0.1)


#boxplot females by stage

female <- subset(pp, sex == 'F')
View(female)
male <- subset(pp, sex == 'M')


ggplot(male, aes(x=Target, y=Cq.Mean, color=stage)) + geom_point()
  
#geom_boxplot(aes(x=Target, y=Cq.Mean, color=stage), width = 0.5, alpha = 0.1)




