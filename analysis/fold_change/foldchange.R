#data
allsilos <- read.csv("Documents/robertslab/labnotebook/data/allsilos-tag_and_annot.csv")

silo3.9 <- allsilos[,-c(3,6,9,12,15,18,21)]

#logtransform
silo3.9 <- cbind(silo3.9[,c(1,43,59)],log(silo3.9[,3:16], 2))

colnames(silo3.9)[2] <- "Uniprot"

#foldchanges
library(gtools)
silo3.9$fc3 <- foldchange(silo3.9$X9_3,silo3.9$X3_3)

silo3.9$fc5 <- foldchange(silo3.9$X9_5,silo3.9$X3_5)

silo3.9$fc7 <- foldchange(silo3.9$X9_7,silo3.9$X3_7)

silo3.9$fc9 <- foldchange(silo3.9$X9_9,silo3.9$X3_9)

silo3.9$fc11 <- foldchange(silo3.9$X9_11,silo3.9$X3_11)

silo3.9$fc13 <- foldchange(silo3.9$X9_13,silo3.9$X3_13)

silo3.9$fc15 <- foldchange(silo3.9$X9_15,silo3.9$X3_15)

fc <- silo3.9[,c(1:3,18:24)]

#make seperate dataframes
day3 <- fc[,c(1:3,4)]
  
day5 <- fc[,c(1:3,5)]

day7 <- fc[,c(1:3,6)]

day9 <- fc[,c(1:3,7)]

day11 <- fc[,c(1:3,8)]

day13 <- fc[,c(1:3,9)]

day15 <- fc[,c(1:3,10)]

#remove NaN and Inf, and set cutoff fold change
day3 <-  day3[which(!is.na(day3$fc3) & !is.infinite(day3$fc3) & day3$fc3 >= 2 | !is.na(day3$fc3) & !is.infinite(day3$fc3) & day3$fc3 <= -2),]

day5 <-  day5[which(!is.na(day5$fc5) & !is.infinite(day5$fc5) & day5$fc5 >= 2 | !is.na(day5$fc5) & !is.infinite(day5$fc5) & day5$fc5 <= -2),]

day7 <-  day7[which(!is.na(day7$fc7) & !is.infinite(day7$fc7) & day7$fc7 >= 2 | !is.na(day7$fc7) & !is.infinite(day7$fc7) & day7$fc7 <= -2),]

day9 <-  day9[which(!is.na(day9$fc9) & !is.infinite(day9$fc9) & day9$fc9 >= 2 | !is.na(day9$fc9) & !is.infinite(day9$fc9) & day9$fc9 <= -2),]

day11 <-  day11[which(!is.na(day11$fc11) & !is.infinite(day11$fc11) & day11$fc11 >= 2 | !is.na(day11$fc11) & !is.infinite(day11$fc11) & day11$fc11 <= -2),]

day13 <-  day13[which(!is.na(day13$fc13) & !is.infinite(day13$fc13) & day13$fc13 >= 2 | !is.na(day13$fc13) & !is.infinite(day13$fc13) & day13$fc13 <= -2),]

day15 <-  day15[which(!is.na(day15$fc15) & !is.infinite(day15$fc15) & day15$fc15 >= 2 | !is.na(day15$fc15) & !is.infinite(day15$fc15) & day15$fc15 <= -2),]

#save dataframes
setwd("Documents/robertslab/labnotebook/analysis/fold_change/")
write.csv(day3, "fc3.csv")
write.csv(day5, "fc5.csv")
write.csv(day7, "fc7.csv")
write.csv(day9, "fc9.csv")
write.csv(day11, "fc11.csv")
write.csv(day13, "fc13.csv")
write.csv(day15, "fc15.csv")


