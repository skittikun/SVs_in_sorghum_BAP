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

right = function (string, char){
  substr(string,nchar(string)-(char-1),nchar(string))
}


left = function (string,char){
  substr(string,1,char)
}


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



# import the data
LDfile <- list.files(path="/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy_ld/ld_sv347", pattern="*.ld", full.names=T, recursive=FALSE)
LDcase <- gsub("*.ld","", gsub("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy_ld/ld_sv347/", "",LDfile))
LDcase <- LDcase[-6]
LDcase
halfset <- c()
lld_all <- read.table(LDfile[1],header=T)
smoothlist <- list()
for(i in 1:length(LDcase)){
#i = 2
      LDcase[i]
      ld <- read.table(LDfile[i],header=T)
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
      #lowest <- detdiff[detdiff$smooth_vals == min(detdiff$smooth_vals),]
      if(LDcase[i] != "oncds"){
        lowest <- detdiff[detdiff$smooth_vals == min(detdiff$smooth_vals),]
      }else{
        lowest <- detdiff[detdiff$smooth_vals == min(detdiff[detdiff$`lld$dist`<=250000,]$smooth_vals),]
      }
      highest <- detdiff[detdiff$smooth_vals == max(detdiff[detdiff$`lld$dist`<=600000,]$smooth_vals),]
      halfpoint <- (highest[1,]$smooth_vals+lowest$smooth_vals)/2
      
      halfline <- detdiff[detdiff$`lld$dist`<=600000,] %>%
        arrange(abs(smooth_vals -  halfpoint)) %>%
        dplyr::slice(1)
      halfset <- rbind.data.frame(halfset, cbind.data.frame(LDcase[i],halfline[1,]))
      
      dfforlist <- cbind.data.frame(lld[,c("SNP_A", "SNP_B", "R2", "dist")], smooth_vals) %>% setNames(c("SNP_A", "SNP_B", paste(c("R2", "dist","smooth"), LDcase[i], sep = "_")))
      
      #filter over shoot over 60k
      tocut <- dfforlist[dfforlist[,4] >= 600000 & dfforlist[,3] > 0.25,]
      dfforlist <- dfforlist[dfforlist[,4] <= 600000,] 
      dfforlist[(dfforlist$SNP_A %in% tocut$SNP_A & dfforlist$SNP_B %in% tocut$SNP_B),ncol(dfforlist)] <- NA
      
      
      listin <- list(dfforlist)
      names(listin) <-  LDcase[i]
      smoothlist[i] <- listin
      names(smoothlist)[i] <-  LDcase[i]
}
str(smoothlist)
halfset
      

#make smooth lines
lldm <- smoothlist[[1]]#lldmm <- join_all(smoothlist, by= c("SNP_A", "SNP_B"), type='left', match = "first")
for(i in 2:length(smoothlist)){
  lldm <- merge(lldm, smoothlist[[i]], by.x = c("SNP_A", "SNP_B"), by.y = c("SNP_A", "SNP_B"), all.x = T, all.y = T)
}

dim(lldm)
head(lldm)

lldm <- lldm %>% mutate(dist = coalesce(dist_347, dist_oncds,dist_ongene,dist_sigcds, dist_siggene),
                        R2 = coalesce(R2_347, R2_oncds,R2_ongene,R2_sigcds, R2_siggene), 
                        smooth = coalesce(smooth_347, smooth_oncds,smooth_ongene,smooth_sigcds, smooth_siggene)) 
lldm$cat <- 1
lldm[!is.na(lldm$smooth_ongene),]$cat <- 2
lldm[!is.na(lldm$smooth_oncds),]$cat <- 3

lldm$cat <- as.factor(lldm$cat)
dim(lldm)
table(lldm$cat)

colorset <- c("black", "blue", "red", "gold", "green")
      LDplot <- ggplot(lldm, aes(x = dist, y = R2, color=cat), size=1) + 
        geom_point() + 
        scale_color_manual(values=c("grey", "blue", "red"), labels = c("all SVs", "SVs on genes", "SVs on CDS"), name = "") + 
        geom_smooth(aes(y = smooth_347, x = lldm$`dist_347`), se=F, color = colorset[1]) +
        geom_smooth(aes(y = smooth_ongene, x = lldm$`dist_ongene`), se=F,colour=colorset[2]) +
        geom_smooth(aes(y = smooth_oncds, x = lldm$`dist_oncds`), se=F, color = colorset[3]) +
        #geom_smooth(aes(y = smooth_sigcds, x = lldm$`dist`), se=F,colour=colorset[4]) +
        #geom_smooth(aes(y = smooth_siggene, x = lldm$`dist`), se=F,colour=colorset[5]) +
        labs(title = paste("LD of ", "8,170 SVs (", format(nrow(lldm),big.mark=",",scientific=FALSE), " pairs)", sep = "")) + 
        xlab("Distance between SVs (bp)") + ylab("Correlation") + 
        theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold")) +
        theme_bw() + theme(legend.position = c(0.8, 0.8),legend.title = element_text(colour="black", size=14, 
                                                                                    face="bold"))
      LDplot

      halfset   
      
#plot with verticle lines
pp <- LDplot +
  geom_vline(xintercept = halfset$`lld$dist`[1], color = colorset[1], size = 0.5, linetype="dashed") +
  geom_text(aes(x=halfset$`lld$dist`[1]), label=paste(left(halfset$`lld$dist`[1], 2), "kb\n", sep = ""), y=0.3, colour=colorset[1], angle=90) +
  
  geom_vline(xintercept = halfset$`lld$dist`[3], color = colorset[2], size = 0.5, linetype="dashed") +
  geom_text(aes(x=halfset$`lld$dist`[3]), label=paste("\n", left(halfset$`lld$dist`[3], 2), "kb", sep = ""), y=0.3, colour=colorset[2], angle=90) +
  
  geom_vline(xintercept = halfset$`lld$dist`[2], color = colorset[3], size = 0.5, linetype="dashed") +
  geom_text(aes(x=halfset$`lld$dist`[2]), label=paste("\n", left(halfset$`lld$dist`[2], 2), "kb", sep = ""), y=0.3, colour=colorset[3], angle=90) +
  
  scale_x_continuous(breaks = seq(100000, 600000, by = 100000)) +
  theme(axis.text.x = element_text(angle=0, size = 10))
pp


tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/lumpy_ld_with_subcat.tiff", width=8, height=5, units="in", res=300)
pp
dev.off()
