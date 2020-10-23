
length <- read.csv("Documents/robertslab/work/FFAR-geoduck/heath-stack-juveniles/20190426-juvenile_measurements.csv")

colnames(length) <- c("Order", "Tank", "Treatment", "Alive", "Dead", "n", "Length", "Width", "Notes")
length <- length[,-c(ncol(length))]

length <- length[-c(8:9, 19:20, 36:37,54:55,70:71,83:84,99:100,114:115),]

H0_T <- length[c(1:7),]
H0_B <- length[c(8:16),]
H1_T <- length[c(17:31),]
H1_B <- length[c(32:47),]
H2_T <- length[c(48:61),]
H2_B <- length[c(62:72),]
H3_T <- length[c(73:86),]
H3_B <- length[c(87:99),]

H0_T$Treatment <- rep("EA", length(H0_T$Treatment))
H0_B$Treatment <- rep("EE", length(H0_B$Treatment))
H1_T$Treatment <- rep("AE", length(H1_T$Treatment))
H1_B$Treatment <- rep("AA", length(H1_B$Treatment))
H2_T$Treatment <- rep("EA", length(H2_T$Treatment))
H2_B$Treatment <- rep("AA", length(H2_B$Treatment))
H3_T$Treatment <- rep("AE", length(H3_T$Treatment))
H3_B$Treatment <- rep("EE", length(H3_B$Treatment))

H0_T$Tank <- rep("H0_T", length(H0_T$Treatment))
H0_B$Tank <- rep("H0_B", length(H0_B$Treatment))
H1_T$Tank <- rep("H1_T", length(H1_T$Treatment))
H1_B$Tank <- rep("H1_B", length(H1_B$Treatment))
H2_T$Tank <- rep("H2_T", length(H2_T$Treatment))
H2_B$Tank <- rep("H2_B", length(H2_B$Treatment))
H3_T$Tank <- rep("H3_T", length(H3_T$Treatment))
H3_B$Tank <- rep("H3_B", length(H3_B$Treatment))

H0_T <- H0_T[c(2:3,7:8)]
H0_B <- H0_B[c(2:3,7:8)]
H1_T <- H1_T[c(2:3,7:8)]
H1_B <- H1_B[c(2:3,7:8)]
H2_T <- H2_T[c(2:3,7:8)]
H2_B <- H2_B[c(2:3,7:8)]
H3_T <- H3_T[c(2:3,7:8)]
H3_B <- H3_B[c(2:3,7:8)]

bind <- rbind(H0_B,H0_T,H1_B,H1_T,H2_B,H2_T,H3_B,H3_T)

#library(tidyr)
#test <- spread(bind, "Tank", "Treatment")

bind$Tank <- as.factor(bind$Tank)
bind$Treatment <- as.factor(bind$Treatment)

library(ggplot2)
library(ggbeeswarm)

attach(bind)
jpeg(filename = "Documents/robertslab/work/FFAR-geoduck/heath-stack-juveniles/length-treatment.jpeg", width = 1000, height = 1000)
ggplot(bind, aes(Treatment,Length, col = Treatment)) + geom_beeswarm(size = 3, cex = 3) + geom_boxplot(aes(Treatment, Length), width = 0.5, alpha = 0.1)
dev.off()

jpeg(filename = "Documents/robertslab/work/FFAR-geoduck/heath-stack-juveniles/length-tank.jpeg", width = 1000, height = 1000)
ggplot(bind, aes(Tank,Length, col = Tank)) + geom_beeswarm(size = 3, cex = 3) + geom_boxplot(aes(Tank, Length), width = 0.5, alpha = 0.1)
dev.off()

#Two-sided test
#t.test(x= Treatment, y= Length, mu = 0, alternative="two.sided", conf.level=0.95, var.equal =F, paired=F)

blockanova <- aov(Length~Treatment+Tank, data=bind)
summary(blockanova)

# Df Sum Sq Mean Sq F value Pr(>F)  
#Treatment    3   44.0  14.667   3.170  0.028 *
#  Tank         4   25.6   6.410   1.385  0.245  
#Residuals   91  421.0   4.626  

#means are not equal so do a post-hoc test

tukey <- TukeyHSD(blockanova, conf.level=0.95)

attach(bind)
results1 <- lm(Length~Treatment+Tank)
anova(results1)

#Evaluate residuals
jpeg(filename = "Documents/robertslab/work/FFAR-geoduck/heath-stack-juveniles/homogenization.jpeg", width = 1000, height = 1000)
par(mfrow=c(1,2))
plot(results1, which = 1)
plot(results1, which = 2)
dev.off()

data.frame(tukey$Treatment[,4])
data.frame(tukey$Tank[,4])


