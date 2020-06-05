library(multcomp)
library(multcompView)
library(Hmisc)
library(gridExtra)
library("viridis")
library(stringi)
library(stringr)
library(plyr)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(betareg)
library(seqinr)
library(vcfR)
library(seqinr)
library(devtools)
library(parallel)
library(iterators)
library(foreach)
library(doParallel)
library(IRanges)
library(GenomicRanges)
library(plyr)
library(stringi)
library(stringr)
library(lme4)
library(lmerTest)
library(betareg)
library(seqinr)
library(vcfR)
library(seqinr)
library(devtools)
library(pROC)
library(Epi)
library(VennDiagram)
library(gridExtra)
library(eulerr)
library(stringdist)
library(reshape2)
library(gtools)
library(stringdist)
library(RecordLinkage)
library(radmixture)
library(FactoMineR)
library(scatterplot3d)
library(RColorBrewer)
library(scales)
library(circlize)
library(DescTools)
library(ggplot2)
library(dplyr)
library(stringr)

options(stringsAsFactors = F)

ld <- read.table("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy_ld/ld_sv347/347.ld",header=T)
str(ld)
head(ld)
dim(ld)
str(ld)

lld <- ld[,c("SNP_A", "CHR_A", "BP_A", "CHR_B", "BP_B", "SNP_B", "R2")]
dim(lld)
head(lld)

lld$dist <- (lld$BP_B - lld$BP_A)
sum(lld$dist == 0)
summary(lld$dist)
range(lld$dist, na.rm = T)

#split chr
lld$CHR_A <- gsub("\\_.*", "", lld$SNP_A)
lld$CHR_B <- gsub("\\_.*", "", lld$SNP_B)

lld <- lld[(lld$CHR_A == lld$CHR_B),]

#find the lowest/highest spot
smooth_vals = predict(loess(R2~dist,lld), lld$dist, method = "gam")
plot(x = lld$dist, y = smooth_vals)
determine <- cbind.data.frame(smooth_vals, lld$dist)
detdiff <- determine#determine[determine$`lld$dist` < 100000,]
lowest <- detdiff[detdiff$smooth_vals == min(detdiff$smooth_vals),]
highest <- detdiff[detdiff$smooth_vals == max(detdiff$smooth_vals),]
halfpoint <- (highest[1,]$smooth_vals+lowest$smooth_vals)/2

halfline <- detdiff[round(detdiff$smooth_vals, 4) == round(halfpoint, 4),]

#plot
LDplot <- ggplot(lld, aes(dist, R2), size=1) + 
  geom_point() + 
  geom_smooth(aes(y = smooth_vals, x = dist), se=F) +
  labs(title = paste("LD of ", "8,170 SVs (3,417 pairs)", sep = "")) + 
  xlab("Distance between SVs (bp)") + ylab("Correlation") + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold")) 
LDplot

pp <- LDplot + geom_vline(xintercept=mean(halfline$`lld$dist`), color = "green", size=0.5)+
  theme_bw() + 
  scale_x_continuous(breaks = seq(100000, 600000, by = 100000)) +
  theme(axis.text.x = element_text(angle=0, size = 9))+
  geom_text(aes(x=mean(halfline$`lld$dist`), label="19kb\n", y=0.3), colour="green", angle=90,text=element_text(size=8))

pp

tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/lumpy_ld.tiff", width=8, height=5, units="in", res=300)
pp
dev.off()


