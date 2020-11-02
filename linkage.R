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


# import the data
LDfile <- list.files(path="/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy_ld/ld_sv347", pattern="*.ld", full.names=T, recursive=FALSE)
LDcase <- gsub("*.ld","", gsub("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy_ld/ld_sv347/", "",LDfile))
LDcase <- LDcase[1:3]
LDcase
halfset <- c()
lld_all <- read.table(LDfile[1],header=T)
smoothlist <- list()
windowlist <- list()
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
      
      #sliding window 50 bp
      windowsize <- 100
      
      windowdt <- c()
      for(j in seq(min(lld$dist), max(lld$dist)+windowsize, by = windowsize)){
        lld_w <- lld[lld$dist < j & lld$dist > (j - windowsize), ]
        R2_w <- mean(lld_w$R2, na.rm = T)
        dist_w <- sum(j, (j-windowsize))/2
        windowdt <- rbind.data.frame(windowdt,cbind.data.frame(R2_w, dist_w))
        print(j)
      }
      windowdt <- windowdt[is.finite(windowdt$R2_w),]
      #windowdt <- windowdt[windowdt$dist_w <= 600000,]
      
      #find the lowest/highest spot
      smooth_vals = predict(loess(R2_w~dist_w,windowdt), windowdt$dist_w, method = "gam")
      plot(x = windowdt$dist_w, y = smooth_vals)
      determine <- cbind.data.frame(smooth_vals, windowdt$dist_w)
      detdiff <- determine#determine[determine$`lld$dist` < 100000,]
      #lowest <- detdiff[detdiff$smooth_vals == min(detdiff$smooth_vals),]
      if(!(LDcase[i] %in% c("oncds", "347"))){
        lowest <- detdiff[detdiff$smooth_vals == min(detdiff$smooth_vals),]
      }else{
        lowest <- detdiff[detdiff$smooth_vals == min(detdiff[detdiff$`windowdt$dist_w`<=250000,]$smooth_vals),]
      }
      highest <- detdiff[detdiff$smooth_vals == max(detdiff[detdiff$`windowdt$dist_w`<=600000,]$smooth_vals),]
      halfpoint <- (highest[1,]$smooth_vals+lowest$smooth_vals)/2
      
      halfline <- detdiff[detdiff$`windowdt$dist_w`<=600000,] %>%
        arrange(abs(smooth_vals -  halfpoint)) %>%
        dplyr::slice(1)
      halfset <- rbind.data.frame(halfset, cbind.data.frame(LDcase[i],halfline[1,]))

      dfforlist <- cbind.data.frame(lld[,c("SNP_A", "SNP_B", "R2", "dist")]) %>% setNames(c("SNP_A", "SNP_B", paste(c("R2", "dist"), LDcase[i], sep = "_")))
      
      #filter over shoot over 60k
      tocut <- dfforlist[dfforlist[,4] >= 600000 & dfforlist[,3] > 0.25,]
      dfforlist <- dfforlist[dfforlist[,4] <= 600000,] 
      dfforlist[(dfforlist$SNP_A %in% tocut$SNP_A & dfforlist$SNP_B %in% tocut$SNP_B),ncol(dfforlist)] <- NA
      
      
      listin <- list(dfforlist)
      names(listin) <-  LDcase[i]
      smoothlist[i] <- listin
      names(smoothlist)[i] <-  LDcase[i]
      
      windowin <- list(windowdt)
      names(windowin) <- LDcase[i]
      windowlist[i] <- windowin
      names(windowlist)[i]  <-  LDcase[i]
}
str(smoothlist)
halfset
      

colorset <- c("black", "blue", "red", "gold", "green")
LDplot <- ggplot() + 
  geom_point(data = smoothlist$`347`, aes(x = dist_347, y = R2_347), size=1, color="grey") +
  geom_point(data = smoothlist$ongene, aes(x = dist_ongene, y = R2_ongene), size=1, color="blue") +
  geom_point(data = smoothlist$oncds, aes(x = dist_oncds, y = R2_oncds), size=1, color="red") +
  geom_smooth(data = windowlist$`347`, aes(y = R2_w, x = dist_w), se=F, color = colorset[1]) +
  geom_smooth(data = windowlist$ongene, aes(y = R2_w, x = dist_w), se=F, color = colorset[2]) +
  geom_smooth(data = windowlist$oncds, aes(y = R2_w, x = dist_w), se=F, color = colorset[3]) +
  xlab("Distance between SVs (bp)") + ylab("Correlation") + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold")) +
  xlim(0, 600000) +
  theme_bw() + theme(legend.position = c(0.8, 0.8),legend.title = element_text(colour="black", size=14, 
                                                                               face="bold"))
halfset   
      
#plot with verticle lines
pp <- LDplot +
  geom_vline(xintercept = halfset$`windowdt$dist_w`[1], color = colorset[1], size = 0.5, linetype="dashed") +
  geom_text(aes(x=halfset$`windowdt$dist_w`[1]), label=paste(left(halfset$`windowdt$dist_w`[1], 2), "kb\n", sep = ""), y=0.3, colour=colorset[1], angle=90) +
  
  geom_vline(xintercept = halfset$`windowdt$dist_w`[3], color = colorset[2], size = 0.5, linetype="dashed") +
  geom_text(aes(x=halfset$`windowdt$dist_w`[3]), label=paste("\n", left(halfset$`windowdt$dist_w`[3], 2), "kb", sep = ""), y=0.3, colour=colorset[2], angle=90) +
  
  geom_vline(xintercept = halfset$`windowdt$dist_w`[2], color = colorset[3], size = 0.5, linetype="dashed") +
  geom_text(aes(x=halfset$`windowdt$dist_w`[2]), label=paste("\n", left(halfset$`windowdt$dist_w`[2], 2), "kb", sep = ""), y=0.3, colour=colorset[3], angle=90) +
  
  scale_x_continuous(breaks = seq(100000, 600000, by = 100000)) +
  theme(axis.text.x = element_text(angle=0, size = 10)) +
  xlim(0, 600000)
pp


tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/lumpy_ld_with_subcat_windows.tiff", width=8, height=5, units="in", res=300)
pp
dev.off()
