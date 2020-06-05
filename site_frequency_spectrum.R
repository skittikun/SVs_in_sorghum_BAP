library(ggplot2)
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
library(Rfast)
library(lmomco)

options(stringsAsFactors = F)

round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}

gff <- read.table("/users/ksongsom/Ref_Genomes/BTx/gene_anno.txt", quote = "", fill = T, header = T, sep = "\t",na.strings=c("","NA") ,stringsAsFactors=FALSE,check.names=FALSE)
altref <- cbind.data.frame(c("Sobic.004G301500", "Sobic.004G301600", "Sobic.004G301650"), c(64031808, 64036221, 64038396), c(64032785, 64037462, 64039002)) %>% setNames(c("gene", "pos", "end"))
sugargene <- read.csv("/Users/ksongsom/OneDrive/postdoc/VIT/sugargenes_clean.csv", header = T, stringsAsFactors = F, na.strings=c("","NA"), sep = "\t")
background <- read.csv("/Users/ksongsom/OneDrive/postdoc/NAM/Panel/BAP/PanelAccessions-BAP.csv", header = T, stringsAsFactors = F)

#lumpy
#dt <- readRDS("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/zeroone_352_fit_precise.RDS")
#dim(dt)
#dt <- "/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/zeroone_352_fit.txt"
#lumpy <- read.table(dt, quote = "", fill = T, header = T, sep = "\t", na.strings=c("","NA"), stringsAsFactors=FALSE, comment.char = "", check.names = F)
#dim(lumpy)
lumpy <- readRDS("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/zeroone_352_fit.RDS")
dim(lumpy)
head(lumpy)
lumpy <- lumpy[grepl("Chr", lumpy$chr),]
sum(lumpy$maf_col >= 0.05, na.rm = T)
sum(lumpy$missing_col <= 0.2, na.rm = T)
sum(lumpy$maf_col >= 0.05 & lumpy$missing_col <= 0.2, na.rm = T)
lumpy$end <- ifelse(lumpy$svtype == "BND", lumpy$pos, lumpy$end)
lumpy$svlen <- lumpy$end - lumpy$pos

mafcut <- 0.05
missingcut <- 0.2
sizecut <- 100000
sizecut_l <- 50
fillp <- lumpy[lumpy$maf_col >= mafcut & lumpy$missing_col <= missingcut & lumpy$svlen <= sizecut & lumpy$svlen >= sizecut_l, ]
fillp <- fillp[is.na(fillp$imprecise),]
fillp$svlen <- fillp$end -fillp$pos
fillp <- fillp[abs(fillp$svlen) <= 100000,]
dim(fillp)
head(fillp)

##table1 for paper
dim(lumpy)
table(lumpy$svtype)
lumpy_f1 <- lumpy[is.na(lumpy$imprecise),]
dim(lumpy_f1)
table(lumpy_f1$svtype)

lumpy_f2 <- lumpy_f1[lumpy_f1$svtype != "BND",]

lumpy_f3 <- lumpy_f2[lumpy_f2$svlen <= sizecut & lumpy_f2$svlen >= sizecut_l,]
dim(lumpy_f3)
table(lumpy_f3$svtype)

lumpy_f4 <- lumpy_f3[lumpy_f3$missing_col <= missingcut,]
dim(lumpy_f4)
table(lumpy_f4$svtype)

lumpy_f5 <- lumpy_f4[lumpy_f4$maf_col >= mafcut,]
dim(lumpy_f5)
table(lumpy_f5$svtype)
summary(lumpy_f5$svlen)

#Updated table
steps <- cbind(
  t(cbind(NA, t(table(lumpy_f1$svtype)))),
  round(t(cbind(NA, t(table(lumpy_f1$svtype))))*100/nrow(lumpy_f1), 3),
  t(cbind(NA, t(table(lumpy_f2$svtype)), NA)),
  round(t(cbind(NA, t(table(lumpy_f2$svtype)), NA))*100/nrow(lumpy_f2), 3),
  t(cbind(NA, t(table(lumpy_f3$svtype)), NA)),
  round(t(cbind(NA, t(table(lumpy_f3$svtype)), NA))*100/nrow(lumpy_f3), 3),
  t(cbind(NA, t(table(lumpy_f4$svtype)), NA)),
  round(t(cbind(NA, t(table(lumpy_f4$svtype)), NA))*100/nrow(lumpy_f4), 3),
  t(cbind(NA, t(table(lumpy_f5$svtype)), NA)),
  round(t(cbind(NA, t(table(lumpy_f5$svtype)), NA))*100/nrow(lumpy_f5), 3)
)
steps[1,] <- cbind(nrow(lumpy_f1), NA,nrow(lumpy_f2),NA,nrow(lumpy_f3), NA,nrow(lumpy_f4), NA,nrow(lumpy_f5), NA)
steps

#table 2 by chromosome
table(lumpy_f5$svtype)
bychr <- lumpy_f5 %>% group_by(chr, svtype) %>% summarise(n = length(format))
bychr

summary(lumpy_f5$svlen)
lumpy_f5 %>% group_by(svtype) %>% summarise(min = min(svlen), max = max(svlen), median = median(svlen))

#gatk

MAFcal <- function(rr){
  #rr <- rr[41:length(rr)]
  predat <- table(t(rr[rr != -1]))
  minor <- as.numeric(names(which(predat == min(predat))[1]))
  maf <- sum(rr == minor, na.rm = T)/sum(rr != -1, na.rm = T)
  return(maf)
}

NAcal <- function(rr){
  #rr <- rr[41:length(rr)]
  nasum <- sum(rr == -1)/length(rr)
  return(nasum)
}

mafrow <- c()
narow <- c()
for(i in 2:ncol(nofilsv)){
  mafrow <- c(mafrow, MAFcal(nofilsv[,i]))
  narow <- c(narow, NAcal(nofilsv[,i]))
  print(paste(i, i/ncol(nofilsv), sep = "_"))
}

#plot lumpy
          
cluster <- read.csv("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/cluster_list.csv", row.names = 1)
dim(lumpy)
colnames(lumpy)
nrow(cluster)

table(lumpy$svtype)
head(lumpy[lumpy$svtype == "DEL",])
dim(fillp)
str(fillp)
colnames(fillp)
dim(fillp[is.na(fillp$imprecise),])
table(fillp$svtype)
summary(lumpy[lumpy$svtype == "DUP",]$svlen)

dellumpy <- lumpy[lumpy$svtype == "DEL",]
head(dellumpy)
sum(dellumpy$pos <= dellumpy$end)
sum(dellumpy$pos > dellumpy$end)
sum(lumpy$imprecise, na.rm = T)
dim(lumpy)
table(lumpy$svtype)
table(fillp$svtype)
table(fillp2$svtype)

dellumpy[dellumpy$pos > dellumpy$end,]

fillp2 <- lumpy[lumpy$missing_col <= missingcut & lumpy$svlen <= sizecut & lumpy$maf_col != 0, ]
dim(fillp2)
colnames(fillp2)
table(fillp2$svtype)
head(fillp2)

lumpyindividual <- colnames(fillp2)[(which(colnames(fillp2) == "format")+1):which(colnames(fillp2)== "rio")]
#fillp2$min_n <- round2(fillp2$maf_col*length(lumpyindividual), 0)
MINcal <- function(rr){
  #rr <- rr[41:length(rr)]
  predat <- table(t(rr[rr != -1]))
  minor <- as.numeric(names(which(predat == min(predat))[1]))
  maf <- sum(rr == minor, na.rm = T)
  return(maf)
}
min_n_count <- c()
for(i in 1:nrow(fillp2)){
  min_n_count <- c(min_n_count, MINcal(fillp2[i,(which(colnames(fillp2)=="144134"):which(colnames(fillp2)=="rio"))]))
  print(paste(i, i/nrow(fillp2), sep = "_"))
}
fillp2$min_n <- min_n_count

head(fillp2)

#saveRDS(fillp2, "/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/fillp2.RDS")
fillp2 <-fillp4
fillp2_bin <- fillp4_bin
#test and plot lumpy
indilength <- length(lumpyindividual)
w_lp <- (nrow(fillp2)/harsumf((2*indilength-1)))/(genomelength)
w_g <- (nrow(filgatk)/harsumf((2*indilength-1)))/genomelength

size_l <- (seq(1, indilength/2, 1))
size_l <- c(size_l, last(size_l)+1) #unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpy <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpy <- rbind(expect_lumpy, cbind.data.frame(count, (w_lp*genomelength)*(1/count)))  
}
colnames(expect_lumpy) <- c("Var1", "Freq")
expect_lumpy
sum(expect_lumpy$Freq)
nrow(fillp2)
observe_lumpy <- as.data.frame(table(fillp2$min_n_count))
dim(observe_lumpy)
lumpy_test <- merge(observe_lumpy, expect_lumpy, by = "Var1", all.x = T, all.y = T, suffixes = c("_observe", "_expect"))
head(lumpy_test)
lumpy_test <- na.omit(lumpy_test)
chisq.test(lumpy_test[,"Freq_observe"], lumpy_test[,"Freq_expect"])

#test by type del
lpdel <- fillp2_bin[fillp2_bin$svtype == "DEL",]
w_lpdel <- (sum(fillp2_bin[fillp2_bin$svtype == "DEL",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpydel <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpydel <- rbind(expect_lumpydel, cbind.data.frame(count, (w_lpdel*genomelength)*(1/count)))  
}
colnames(expect_lumpydel) <- c("Var1", "Freq")
expect_lumpydel
observe_lumpydel <- fillp2_bin[fillp2_bin$svtype == "DEL",c("min_n_count", "n" )] %>% setNames(c("Var1", "Freq"))
dim(observe_lumpydel)
lumpy_testdel <- merge(observe_lumpydel, expect_lumpydel, by = "Var1", all.x = T, all.y = T, suffixes = c("_observe", "_expect"))
head(lumpy_testdel)
lumpy_testdel <- na.omit(lumpy_testdel)
chisq.test(lumpy_testdel[,"Freq_observe"], lumpy_testdel[,"Freq_expect"])

##bnd
lpbnd <- fillp2_bin[fillp2_bin$svtype == "BND",]
w_lpbnd <- (sum(fillp2_bin[fillp2_bin$svtype == "BND",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpybnd <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpybnd <- rbind(expect_lumpybnd, cbind.data.frame(count, (w_lpbnd*genomelength)*(1/count)))  
}
colnames(expect_lumpybnd) <- c("Var1", "Freq")
expect_lumpybnd
observe_lumpybnd <- fillp2_bin[fillp2_bin$svtype == "BND",c("min_n_count", "n" )] %>% setNames(c("Var1", "Freq"))
dim(observe_lumpybnd)
lumpy_testbnd <- merge(observe_lumpybnd, expect_lumpybnd, by = "Var1", all.x = T, all.y = T, suffixes = c("_observe", "_expect"))
head(lumpy_testbnd)
lumpy_testbnd <- na.omit(lumpy_testbnd)
chisq.test(lumpy_testbnd[,"Freq_observe"], lumpy_testbnd[,"Freq_expect"])

##dup
lpdup <- fillp2_bin[fillp2_bin$svtype == "DUP",]
w_lpdup <- (sum(fillp2_bin[fillp2_bin$svtype == "DUP",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpydup <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpydup <- rbind(expect_lumpydup, cbind.data.frame(count, (w_lpdup*genomelength)*(1/count)))  
}
colnames(expect_lumpydup) <- c("Var1", "Freq")
expect_lumpydup
observe_lumpydup <- fillp2_bin[fillp2_bin$svtype == "DUP",c("min_n", "n" )] %>% setNames(c("Var1", "Freq"))
dim(observe_lumpydup)
lumpy_testdup <- merge(observe_lumpydup, expect_lumpydup, by = "Var1", all.x = T, all.y = T, suffixes = c("_observe", "_expect"))
head(lumpy_testdup)
lumpy_testdup <- na.omit(lumpy_testdup)
chisq.test(lumpy_testdup[,"Freq_observe"], lumpy_testdup[,"Freq_expect"])

##inv
lpinv <- fillp2_bin[fillp2_bin$svtype == "INV",]
w_lpinv <- (sum(fillp2_bin[fillp2_bin$svtype == "INV",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpyinv <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpyinv <- rbind(expect_lumpyinv, cbind.data.frame(count, (w_lpinv*genomelength)*(1/count)))  
}
colnames(expect_lumpyinv) <- c("Var1", "Freq")
expect_lumpyinv
observe_lumpyinv <- fillp2_bin[fillp2_bin$svtype == "INV",c("min_n_count", "n" )] %>% setNames(c("Var1", "Freq"))
dim(observe_lumpyinv)
lumpy_testinv <- merge(observe_lumpyinv, expect_lumpyinv, by = "Var1", all.x = T, all.y = T, suffixes = c("_observe", "_expect"))
head(lumpy_testinv)
lumpy_testinv <- na.omit(lumpy_testinv)
chisq.test(lumpy_testinv[,"Freq_observe"], lumpy_testinv[,"Freq_expect"])

#keepfillp2_bin <- fillp2_bin  
#fillp2_bin <- keepfillp2_bin 
head(fillp2_bin)
head(lumpy_test)
test <- fillp2_bin %>% group_by(min_n_count) %>% summarise(sum(pct_n))

#calculate expected percentage from each type of sv
lumpy_testinv$Freq_expect_percent <- lumpy_testinv$Freq_expect*100/sum(lumpy_testinv$Freq_expect, na.rm = T)
lumpy_testdel$Freq_expect_percent <- lumpy_testdel$Freq_expect*100/sum(lumpy_testdel$Freq_expect, na.rm = T)
lumpy_testdup$Freq_expect_percent <- lumpy_testdup$Freq_expect*100/sum(lumpy_testdup$Freq_expect, na.rm = T)
lumpy_testbnd$Freq_expect_percent <- lumpy_testbnd$Freq_expect*100/sum(lumpy_testbnd$Freq_expect, na.rm = T)

lumpyexpect <- Reduce(function(x, y) merge(x, y, all=TRUE, by = "Var1"), list(lumpy_testinv[,c(1,4)], lumpy_testdel[,c(1,4)], lumpy_testdup[,c(1,4)], lumpy_testbnd[,c(1,4)]))
head(lumpyexpect)
lumpyexpect$sumexpect <- rowSums(lumpyexpect[,-1],na.rm = T)
expect_lumpy$Freq/sum(expect_lumpy$Freq)
#expectation calculation
#Watterson’s estimator (\Theta_w) were
#del
w_lpdel <- (sum(fillp2_bin[fillp2_bin$svtype == "DEL",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpydel <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpydel <- rbind(expect_lumpydel, cbind.data.frame(count, (w_lpdel*genomelength)*(1/count)))  
}
colnames(expect_lumpydel) <- c("Var1", "Freq")

#bnd
lpbnd <- fillp2_bin[fillp2_bin$svtype == "BND",]
w_lpbnd <- (sum(fillp2_bin[fillp2_bin$svtype == "BND",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpybnd <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpybnd <- rbind(expect_lumpybnd, cbind.data.frame(count, (w_lpbnd*genomelength)*(1/count)))  
}
colnames(expect_lumpybnd) <- c("Var1", "Freq")

##dup
lpdup <- fillp2_bin[fillp2_bin$svtype == "DUP",]
w_lpdup <- (sum(fillp2_bin[fillp2_bin$svtype == "DUP",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpydup <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpydup <- rbind(expect_lumpydup, cbind.data.frame(count, (w_lpdup*genomelength)*(1/count)))  
}
colnames(expect_lumpydup) <- c("Var1", "Freq")

##inv
lpinv <- fillp2_bin[fillp2_bin$svtype == "INV",]
w_lpinv <- (sum(fillp2_bin[fillp2_bin$svtype == "INV",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpyinv <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpyinv <- rbind(expect_lumpyinv, cbind.data.frame(count, (w_lpinv*genomelength)*(1/count)))  
}
colnames(expect_lumpyinv) <- c("Var1", "Freq")

#lumpyexpect <- Reduce(function(x, y) merge(x, y, all=TRUE, by = "Var1"), list(expect_lumpyinv[,c(1,3)], expect_lumpydel[,c(1,3)], expect_lumpydup[,c(1,3)], expect_lumpytbnd[,c(1,3)]))
#head(lumpyexpect)
#lumpyexpect$sumexpect <- rowSums(lumpyexpect[,-1],na.rm = T)
colnames(lumpyexpect) <- c("Var1", "Freq.inv", "Freq.del", "Freq.dup", "Freq.bnd", "sumexpect")

expect_lumpydel$per <- expect_lumpydel$Freq*100/sum(expect_lumpydel$Freq, na.rm = T)
expect_lumpydup$per <- expect_lumpydup$Freq*100/sum(expect_lumpydup$Freq, na.rm = T)
expect_lumpyinv$per <- expect_lumpyinv$Freq*100/sum(expect_lumpyinv$Freq, na.rm = T)
lumpyexpect_each <- Reduce(function(x, y) merge(x, y, all=TRUE, by = "Var1"), list(expect_lumpydel, expect_lumpydup, expect_lumpyinv))
colnames(lumpyexpect_each) <- c("Var1", "Freq.del", "per.del", "Freq.dup", "per.dup", "Freq.inv", "per.inv")


    plumpy <- ggplot(fillp2, aes(maf_col, fill = svtype)) + geom_histogram(binwidth = 0.005) + 
      labs(y = "SV counts") +
      scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                               "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
      theme_bw() + 
      ggtitle("Lumpy SV site frequency") +
      scale_x_continuous(breaks = pretty(fillp2$maf_col, n = 10))
    plumpy
    tiff("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy.tiff", width=5, height=5, units="in", res=300)
    plumpy
    dev.off()

      plumpy <- ggplot(fillp2, aes(min_n, fill = svtype)) + geom_histogram(binwidth = 1) + 
        labs(y = "SV counts") +
        scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                                 "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
        theme_bw() + 
        ggtitle("Lumpy SV site frequency") +
        scale_x_continuous(breaks = pretty(fillp2$min_n, n = 10))
      plumpy
      tiff("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy_n.tiff", width=5, height=5, units="in", res=300)
      plumpy
      dev.off()
      
      plumpy_small <- ggplot(fillp2, aes(maf_col, fill = svtype)) + geom_histogram(binwidth = 0.005) + 
        labs(y = "SV counts") +
        scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                                 "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
        theme_bw() + 
        ggtitle("Lumpy SV site frequency") +
        scale_x_continuous(breaks = pretty(fillp$maf_col, n = 10), limits = c(0.02, 0.5)) +
        ylim(0, 4000)
      plumpy_small
      tiff("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy_small.tiff", width=5, height=5, units="in", res=300)
      plumpy_small
      dev.off()

#lumpy percentage 
fillp2$bin <- cut_interval(abs(fillp2$maf_col), length=0.005)
fillp2_bin <- as.data.frame(fillp2 %>%
                              group_by(min_n, svtype) %>%
                              tally  %>%
                              group_by(svtype) %>%
                              mutate(pct_n=round((100*n)/sum(n), 2)))
 #with BND
        plumpy_bin <- ggplot(fillp2_bin, aes(x = min_n, y = pct_n, fill = svtype, label = pct_n)) +
          geom_bar(stat = "identity") + 
          labs(y = "SV abundance percentage", x = "n") +
          scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                                   "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
          theme_bw() + scale_x_continuous(breaks = pretty(fillp$min_n, n = 100)) + theme(axis.text.x = element_text(angle = 90)) +
          ggtitle("Lumpy SV site frequency percentage")
        plumpy_bin
        
        tiff("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy_percent_n.tiff", width=5, height=5, units="in", res=300)
        plumpy_bin
        dev.off()


        plumpy_bin_dodge <- ggplot(fillp2_bin, aes(x = min_n, y = pct_n, fill = svtype, label = pct_n)) +
          geom_bar(stat = "identity", position='dodge') + 
          labs(y = "SV abundance percentage", x = "n") +
          scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                                   "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
          theme_bw() + scale_x_continuous(breaks = pretty(fillp$min_n, n = 100)) + theme(axis.text.x = element_text(angle = 90)) +
          ggtitle("Lumpy SV site frequency percentage")
        plumpy_bin_dodge
        tiff("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy_percent_n_dodge.tiff", width=5, height=5, units="in", res=300)
        plumpy_bin_dodge
        dev.off()

        plumpy_bin_dodge_small <- ggplot(fillp2_bin, aes(x = min_n, y = pct_n, fill = svtype, label = pct_n)) +
          geom_bar(stat = "identity", position='dodge') + 
          labs(y = "SV abundance percentage", x = "n") +
          scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                                   "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
          theme_bw() + scale_x_continuous(breaks = pretty(fillp$min_n, n = 100)) + theme(axis.text.x = element_text(angle = 90)) +
          ggtitle("Lumpy SV site frequency percentage") +
          ylim(0, 50)
        plumpy_bin_dodge_small
        tiff("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy_percent_n_dodge_small.tiff", width=5, height=5, units="in", res=300)
        plumpy_bin_dodge_small
        dev.off()
        


       
#no BND
        fillp4 <- lumpy_f4#[lumpy_f4$maf_col != 0,]#lumpy[lumpy$missing_col <= missingcut & lumpy$svlen <= sizecut & lumpy$maf_col != 0 & lumpy$svtype != "BND", ]
        dim(fillp4)
        dim(lumpy_f5)
        #new min_n_count
        
        MINcal <- function(rr){
          #rr <- rr[41:length(rr)]
          #rr <- fillp4[1,(which(colnames(fillp4)=="144134"):which(colnames(fillp4)=="rio"))]
          predat <- table(t(rr[rr != -1]))
          minor <- as.numeric(names(which(predat == min(predat))[1]))
          if(length(predat) != 1){
          
          maf <- sum(rr == minor, na.rm = T)}else{
            if(length(predat) == 1 & sum(rr == minor, na.rm = T) > sum(predat)/2){
              maf = sum(predat)-sum(rr == minor, na.rm = T)
            }else{if(sum(rr == minor, na.rm = T) == sum(predat)/2){
              maf = sum(predat)/2
            }else{maf = NA}}
            }
          return(maf)
        }
        min_n_count <- c()
        for(i in 1:nrow(fillp4)){
          min_n_count <- c(min_n_count, MINcal(fillp4[i,(which(colnames(fillp4)=="144134"):which(colnames(fillp4)=="rio"))]))
          print(paste(i, i/nrow(fillp4), sep = "_"))
        }
        
        fillp4$min_n_count <- min_n_count
        fillp4$min_n <- round2(fillp4$maf_col*length(lumpyindividual), 0)
        
        #consider cluster
        head(cluster)
        
        min_n_count_k3_1 <- c()
        for(i in 1:nrow(fillp4)){
          min_n_count_k3_1 <- c(min_n_count_k3_1, MINcal(fillp4[i,colnames(fillp4) %in% cluster[cluster$k3_cluster_all == 1,]$ID]))
          print(paste(i, i/nrow(fillp4), sep = "_"))
        }
        length(cluster[cluster$k3_cluster_all == 1,]$ID)
        min_n_count_k3_1
        
        fillp4$min_n_count_k3_1 <- min_n_count_k3_1

        min_n_count_k3_2 <- c()
        for(i in 1:nrow(fillp4)){
          min_n_count_k3_2 <- c(min_n_count_k3_2, MINcal(fillp4[i,colnames(fillp4) %in% cluster[cluster$k3_cluster_all == 2,]$ID]))
          print(paste(i, i/nrow(fillp4), sep = "_"))
        }
        
        fillp4$min_n_count_k3_2 <- min_n_count_k3_2
        
        min_n_count_k3_3 <- c()
        for(i in 1:nrow(fillp4)){
          min_n_count_k3_3 <- c(min_n_count_k3_3, MINcal(fillp4[i,colnames(fillp4) %in% cluster[cluster$k3_cluster_all == 3,]$ID]))
          print(paste(i, i/nrow(fillp4), sep = "_"))
        }
        
        fillp4$min_n_count_k3_3 <- min_n_count_k3_3
        
        #lumpy percentage 
        #fillp4$bin <- cut_interval(abs(fillp4$maf_col), length=0.005)
        fillp4$bin <- cut_interval(abs(fillp4$min_n_count), length=0.005)
        fillp4_bin <- as.data.frame(fillp4 %>%
                                      group_by(min_n_count, svtype) %>%
                                      tally  %>%
                                      group_by(svtype) %>%
                                      mutate(pct_n=round((100*n)/sum(n), 2)))
        fillp2_bin <- merge(fillp2_bin, lumpyexpect, by.x = "min_n_count", by.y = "Var1", all.x = T, all.y = F)
        plumpy_bin <- ggplot(fillp2_bin, aes(x = min_n, y = pct_n, fill = svtype, label = pct_n)) +
          geom_bar(stat = "identity") + 
          labs(y = "SV abundance percentage", x = "n") +
          scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                                   "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
          theme_bw() + scale_x_continuous(breaks = pretty(fillp2_bin$min_n, n = 100)) + theme(axis.text.x = element_text(angle = 90)) +
          ggtitle("Lumpy SV site frequency percentage") +
          geom_line(aes(x = min_n, y = sumexpect), linetype = "dashed") + theme(legend.position="top") 
        plumpy_bin
        lumpyexpect_each$sumexpect <- lumpyexpect_each$per.del*3
        fillp42_bin <- merge(fillp4_bin, lumpyexpect_each, by.x = "min_n_count", by.y = "Var1", all.x = T, all.y = F)
        fillp42_bin$sumexpect_nobnd <- fillp42_bin$sumexpect#-fillp42_bin$Freq.bnd
        fillp42_bin <- fillp42_bin[fillp42_bin$min_n_count != 0,]
        plumpy_bin <- ggplot(fillp42_bin, aes(x = min_n_count, y = pct_n, fill = svtype, label = pct_n)) +
          geom_bar(stat = "identity") + 
          labs(y = "SV abundance percentage", x = "Numbers of SV") +
          scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                                   "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
          theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
          ggtitle("Lumpy SV site frequency percentage (all genotypes)") +
          geom_line(aes(x = min_n_count, y = sumexpect_nobnd), linetype = "dotted", col = "black") +
          geom_text(aes(label = round(n, 1)), color = "black", size = 1, position = position_stack(vjust = 0.5)) +
          scale_x_continuous(breaks = pretty(fillp42_bin$min_n, n = 50)+1) 
        plumpy_bin
        
        tiff("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy_percent_n_nobnd_with_num_allgenotype_fixed.tiff", width=10, height=5, units="in", res=300)
        plumpy_bin
        dev.off()
        
            #group 1
                    fillp4$bin_k3_1 <- cut_interval(abs(fillp4$min_n_count_k3_1), length=0.005)
                    fillp4_bin <- as.data.frame(fillp4 %>%
                                                  group_by(min_n_count_k3_1, svtype) %>%
                                                  tally  %>%
                                                  group_by(svtype) %>%
                                                  mutate(pct_n=round((100*n)/sum(n), 2)))
                    fillp42_bin <- merge(fillp4_bin, lumpyexpect, by.x = "min_n_count_k3_1", by.y = "Var1", all.x = T, all.y = F)
                    fillp42_bin$sumexpect_nobnd <- (fillp42_bin$sumexpect-fillp42_bin$Freq.bnd)#*(length(cluster[cluster$k3_cluster_all == 1,]$ID)/indilength)
                    plumpy_bin <- ggplot(fillp42_bin, aes(x = min_n_count_k3_1, y = pct_n, fill = svtype, label = pct_n)) +
                      geom_bar(stat = "identity") + 
                      labs(y = "SV abundance percentage", x = "Numbers of SV") +
                      scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                                               "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
                      theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
                      ggtitle("Lumpy SV site frequency percentage group k = 1") +
                      geom_line(aes(x = min_n_count_k3_1, y = sumexpect_nobnd), linetype = "dotted", col = "black") +
                      geom_text(aes(label = round(n, 1)), color = "black", size = 1, position = position_stack(vjust = 0.5)) +
                      scale_x_continuous(breaks = pretty(fillp42_bin$min_n, n = 100)) 
                    plumpy_bin
                    
                    tiff("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy_percent_n_nobnd_with_num_k3_1.tiff", width=10, height=5, units="in", res=300)
                    plumpy_bin
                    dev.off()
        
                #group 2
                    fillp4$bin_k3_2 <- cut_interval(abs(fillp4$min_n_count_k3_2), length=0.005)
                    fillp4_bin <- as.data.frame(fillp4 %>%
                                                  group_by(min_n_count_k3_2, svtype) %>%
                                                  tally  %>%
                                                  group_by(svtype) %>%
                                                  mutate(pct_n=round((100*n)/sum(n), 2)))
                    fillp42_bin <- merge(fillp4_bin, lumpyexpect, by.x = "min_n_count_k3_2", by.y = "Var1", all.x = T, all.y = F)
                    fillp42_bin$sumexpect_nobnd <- fillp42_bin$sumexpect-fillp42_bin$Freq.bnd
                    plumpy_bin <- ggplot(fillp42_bin, aes(x = min_n_count_k3_2, y = pct_n, fill = svtype, label = pct_n)) +
                      geom_bar(stat = "identity") + 
                      labs(y = "SV abundance percentage", x = "Numbers of SV") +
                      scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                                               "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
                      theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
                      ggtitle("Lumpy SV site frequency percentage group k = 2") +
                      geom_line(aes(x = min_n_count_k3_2, y = sumexpect_nobnd), linetype = "dotted", col = "black") +
                      geom_text(aes(label = round(n, 1)), color = "black", size = 1, position = position_stack(vjust = 0.5)) +
                      scale_x_continuous(breaks = pretty(fillp42_bin$min_n, n = 100)) 
                    plumpy_bin
                    
                    tiff("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy_percent_n_nobnd_with_num_k3_2.tiff", width=10, height=5, units="in", res=300)
                    plumpy_bin
                    dev.off()
                    
                #group 3
                    fillp4$bin_k3_3 <- cut_interval(abs(fillp4$min_n_count_k3_3), length=0.005)
                    fillp4_bin <- as.data.frame(fillp4 %>%
                                                  group_by(min_n_count_k3_3, svtype) %>%
                                                  tally  %>%
                                                  group_by(svtype) %>%
                                                  mutate(pct_n=round((100*n)/sum(n), 2)))
                    fillp42_bin <- merge(fillp4_bin, lumpyexpect, by.x = "min_n_count_k3_3", by.y = "Var1", all.x = T, all.y = F)
                    fillp42_bin$sumexpect_nobnd <- fillp42_bin$sumexpect-fillp42_bin$Freq.bnd
                    plumpy_bin <- ggplot(fillp42_bin, aes(x = min_n_count_k3_3, y = pct_n, fill = svtype, label = pct_n)) +
                      geom_bar(stat = "identity") + 
                      labs(y = "SV abundance percentage", x = "Numbers of SV") +
                      scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                                               "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
                      theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
                      ggtitle("Lumpy SV site frequency percentage group k = 3") +
                      geom_line(aes(x = min_n_count_k3_3, y = sumexpect_nobnd), linetype = "dotted", col = "black") +
                      geom_text(aes(label = round(n, 1)), color = "black", size = 1, position = position_stack(vjust = 0.5)) +
                      scale_x_continuous(breaks = pretty(fillp42_bin$min_n, n = 100)) 
                    plumpy_bin
                    
                    tiff("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy_percent_n_nobnd_with_num_k3_3.tiff", width=10, height=5, units="in", res=300)
                    plumpy_bin
                    dev.off()
                    
        #side by side
        head(expect_lumpybnd)
        
        expect_lumpydel$per <- expect_lumpydel$Freq*100/sum(expect_lumpydel$Freq, na.rm = T)
        expect_lumpydup$per <- expect_lumpydup$Freq*100/sum(expect_lumpydup$Freq, na.rm = T)
        expect_lumpyinv$per <- expect_lumpyinv$Freq*100/sum(expect_lumpyinv$Freq, na.rm = T)
        lumpyexpect_each <- Reduce(function(x, y) merge(x, y, all=TRUE, by = "Var1"), list(expect_lumpydel, expect_lumpydup, expect_lumpyinv))
        colnames(lumpyexpect_each) <- c("Var1", "Freq.del", "per.del", "Freq.dup", "per.dup", "Freq.inv", "per.inv")
        
        fillp5_bin <- merge(fillp4_bin, lumpyexpect_each, by.x = "min_n_count", by.y = "Var1", all.x = T, all.y = F)
        plumpy_bin_dodge <- ggplot(fillp5_bin, aes(x = min_n_count, y = pct_n, fill = svtype, label = pct_n)) +
          geom_bar(stat = "identity", position="dodge", width=0.9) + 
          labs(y = "SV abundance percentage", x = "numbers of SV") +
          scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                                   "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
          theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 90)) +
          ggtitle("Lumpy SV site frequency percentage") +
          geom_line(aes(x = min_n_count, y = per.dup), linetype = "dotted", col = "black") +
          scale_x_continuous(breaks = pretty(fillp42_bin$min_n, n = 50)+1) 
        plumpy_bin_dodge
        tiff("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy_percent_n_dodge_nobnd_fixed.tiff", width=10, height=5, units="in", res=300)
        plumpy_bin_dodge
        dev.off()
        
            #test chi-sq of each type
        fillp5_bin_del <- fillp5_bin[fillp5_bin$svtype == "DEL", ]
        chisq.test(fillp5_bin_del$pct_n, fillp5_bin_del$per.del)
        
        fillp5_bin_dup <- fillp5_bin[fillp5_bin$svtype == "DUP", ]
        chisq.test(fillp5_bin_dup$pct_n, fillp5_bin_dup$per.dup)
        
        fillp5_bin_inv <- fillp5_bin[fillp5_bin$svtype == "INV", ]
        chisq.test(fillp5_bin_inv$pct_n, fillp5_bin_inv$per.inv)
        

                

#lumpy insertion
dt <- "/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/fil_352_pruned_pe3_dp10_ins.txt"
lumpy_ins <- read.table(dt, quote = "", fill = T, header = T, sep = "\t", na.strings=c("","NA"), stringsAsFactors=FALSE, comment.char = "", check.names = F)
sum(lumpy_ins$maf_col >= 0.05, na.rm = T)
sum(lumpy_ins$missing_col <= 0.2, na.rm = T)
sum(lumpy_ins$maf_col >= 0.05 & lumpy_ins$missing_col <= 0.2, na.rm = T)
lumpy_ins[is.na(lumpy_ins$end),]$end <- lumpy_ins[is.na(lumpy_ins$end),]$pos + 100000
lumpy_ins[is.na(lumpy_ins$pos),]$pos <- lumpy_ins[is.na(lumpy_ins$pos),]$end + 100000
lumpy_ins$svlen <- lumpy_ins$end - lumpy_ins$pos


lumpy_ins$missing_col <- apply(lumpy_ins[, (which(colnames(lumpy_ins)=="format")+1):which(colnames(lumpy_ins)=="rio")], 1, MAFcal)
lumpy_ins$maf_col <- apply(lumpy_ins[, (which(colnames(lumpy_ins)=="format")+1):which(colnames(lumpy_ins)=="rio")], 1, NAcal)
lumpy_ins$min_n <- apply(lumpy_ins[, (which(colnames(lumpy_ins)=="format")+1):which(colnames(lumpy_ins)=="rio")], 1, function(x) sum(x, na.rm = T))
hist(lumpy_ins$min_n)

lumpy_ins$min_n <- ifelse(lumpy_ins$min_n >= length(lumpyindividual)/2, length(lumpyindividual)-lumpy_ins$min_n, lumpy_ins$min_n)


fillp3 <- lumpy_ins[lumpy_ins$missing_col <= missingcut & lumpy_ins$svlen <= sizecut, ]
fillp3$bin <- cut_interval(abs(fillp3$maf_col), length=0.005)
fillp3_bin <- as.data.frame(fillp3 %>%
                              group_by(bin, svtype) %>%
                              tally  %>%
                              group_by(svtype) %>%
                              mutate(pct=round((100*n)/sum(n), 2)))
fillp3_bin$xname <- gsub("]", "",gsub(".*,", "", as.character(fillp3_bin$bin)))

plumpy_bin_small_ins <- ggplot(fillp3_bin, aes(x = bin, y = pct, fill = svtype, label = pct)) +
  geom_bar(stat = "identity") + 
  labs(y = "SV abundance percentage", x = "allele frequency") +
  scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                           "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
  theme_bw() + scale_x_discrete(labels=unique(fillp3_bin$xname)) + theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Lumpy SV site frequency percentage") +
  ylim(0, 50)
plumpy_bin_small_ins

tiff("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy_small_percent_insertion.tiff", width=5, height=5, units="in", res=300)
plumpy_bin_small_ins
dev.off()
#lumpy percentage
fillp3$bin <- cut_interval(abs(fillp3$maf_col), length=0.005)
fillp3_bin <- as.data.frame(fillp3 %>%
                              group_by(bin, svtype) %>%
                              tally  %>%
                              group_by(svtype) %>%
                              mutate(pct=round((100*n)/sum(n), 2)))
fillp3_bin$xname <- gsub("]", "",gsub(".*,", "", as.character(fillp3_bin$bin)))

plumpy_bin <- ggplot(fillp3_bin, aes(x = bin, y = pct, fill = svtype, label = pct)) +
  geom_bar(stat = "identity") + 
  labs(y = "SV abundance percentage", x = "allele frequency") +
  scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                           "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
  theme_bw() + scale_x_discrete(labels=unique(fillp3_bin$xname)) + theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Lumpy SV site frequency percentage")
plumpy_bin

tiff("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy_percent.tiff", width=5, height=5, units="in", res=300)
plumpy_bin
dev.off()

plumpy_bin_small <- ggplot(fillp3_bin, aes(x = bin, y = pct, fill = svtype, label = pct)) +
  geom_bar(stat = "identity") + 
  labs(y = "SV abundance percentage", x = "allele frequency") +
  scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                           "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
  theme_bw() + scale_x_discrete(labels=unique(fillp3_bin$xname)) + theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Lumpy SV site frequency percentage") +
  ylim(0, 50)
plumpy_bin_small

tiff("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy_small_percent.tiff", width=5, height=5, units="in", res=300)
plumpy_bin_small
dev.off()

#compare with theory
#simulate theta
#compare with chi-square


sorghumgenome <- as.data.frame(cbind(c(paste("Chr0", 1:9, sep = ""), "Chr10"), rep(1, 10),c(80884392, 77742459, 74386277, 68658214, 71854669, 61277060, 65505356,62686529,59416394,61233695)))
colnames(sorghumgenome) <- c("Chromosome", "start", "end")
sorghumgenome$Chromosome <- as.factor(sorghumgenome$Chromosome)
sorghumgenome$start <- as.numeric(as.character(sorghumgenome$start))
sorghumgenome$end <- as.numeric(as.character(sorghumgenome$end))
genomelength <- sum(sorghumgenome$end)


harsumf <- function(input){
  harsum <- 0
  for(i in 1:input){
    harsum <- sum(harsum, 1/i)
  }
  return(harsum)
}

indilength <- length(lumpyindividual)
w_lp <- (nrow(fillp2)/harsumf((2*indilength-1)))/(genomelength)
w_g <- (nrow(filgatk)/harsumf((2*indilength-1)))/genomelength

#gatk
head(filgatk_bin)
head(fillp2_bin)
size_g <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(filgatk_bin$xname))
expect_gatk <- c()
for(i in 1:length(size_g)){
  count <- size_g[i]
  expect_gatk <- rbind(expect_gatk, cbind.data.frame(count, (w_g*genomelength)*(1/count)))  
}
colnames(expect_gatk) <- c("Var1", "Freq")
expect_gatk
observe_gatk <- as.data.frame(table(filgatk$min_n))
dim(observe_gatk)
gatk_test <- merge(observe_gatk, expect_gatk, by = "Var1", all.x = T, all.y = T, suffixes = c("_observe", "_expect"))
head(gatk_test)
chisq.test(gatk_test[,"Freq_observe"], gatk_test[,"Freq_expect"])

#lumpy
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpy <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpy <- rbind(expect_lumpy, cbind.data.frame(count, (w_lp*genomelength)*(1/count)))  
}
colnames(expect_lumpy) <- c("Var1", "Freq")
expect_lumpy
sum(expect_lumpy$Freq)
nrow(fillp2)

observe_lumpy <- as.data.frame(table(fillp2$min_n))
dim(observe_lumpy)
lumpy_test <- merge(observe_lumpy, expect_lumpy, by = "Var1", all.x = T, all.y = T, suffixes = c("_observe", "_expect"))
head(lumpy_test)
colSums(lumpy_test[,2:3], na.rm = T)
lumpy_test <- na.omit(lumpy_test)
chisq.test(lumpy_test[,"Freq_observe"], lumpy_test[,"Freq_expect"])

#observe_gatk <- as.data.frame(filgatk %>% group_by(min_n) %>% summarise(count = sum(n, na.rm= T)))

chisq.test(expect_gatk$`(pi_g * genomelength) * (1/count)`, observe_gatk$count)

nrow(filgatk) == sum(expect_gatk$`(pi * genomelength) * (1/count)`)


#test and plot gatk
head(filgatk_bin)
head(fillp2_bin)
size_g <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(filgatk_bin$xname))
expect_gatk <- c()
for(i in 1:length(size_g)){
  count <- size_g[i]
  expect_gatk <- rbind(expect_gatk, cbind.data.frame(count, (w_g*genomelength)*(1/count)))  
}
colnames(expect_gatk) <- c("Var1", "Freq")
expect_gatk
observe_gatk <- as.data.frame(table(filgatk$min_n))
dim(observe_gatk)
gatk_test <- merge(observe_gatk, expect_gatk, by = "Var1", all.x = T, all.y = T, suffixes = c("_observe", "_expect"))
head(gatk_test)
chisq.test(gatk_test[,"Freq_observe"], gatk_test[,"Freq_expect"])

ppgatk_bin <- ggplot(filgatk_bin, aes(x = min_n, y = pct_n, fill = svtype, label = min_n)) +
  geom_bar(data = filgatk_bin, stat = "identity", width = 0.5) + 
  labs(y = "SV abundance percentage", x = "allele frequency") +
  scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                           "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
  theme_bw() + scale_x_discrete(labels=unique(filgatk_bin$xname)) + theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Gatk SV site frequency percentage") +
  scale_x_continuous(breaks = pretty(filgatk$min_n, n = 100))
geom_line(aes(x = c(as.numeric(gatk_test$Var1),as.numeric(gatk_test$Var1)), y = c(gatk_test$Freq_expect, gatk_test$Freq_expect))) 
ppgatk_bin

#tiff("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/gatk_percent_n.tiff", width=5, height=5, units="in", res=300)
ppgatk_bin
dev.off()

ppgatk_bin <- ggplot(filgatk_bin, aes(x = min_n, y = pct_n, fill = svtype, label = pct_n)) +
  geom_bar(stat = "identity", position='dodge', width = 0.5) + 
  labs(y = "SV abundance percentage", x = "allele frequency") +
  scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                           "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
  theme_bw() + scale_x_discrete(labels=unique(filgatk_bin$xname)) + theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Gatk SV site frequency percentage") +
  scale_x_continuous(breaks = pretty(filgatk$min_n, n = 100))
ppgatk_bin

ppgatk_bin_small <- ggplot(filgatk_bin, aes(x = min_n, y = pct_n, fill = svtype, label = pct_n))  +
  geom_bar(stat = "identity", position='dodge', width = 0.5) + 
  labs(y = "SV abundance percentage", x = "allele frequency") +
  scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                           "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
  theme_bw() + scale_x_discrete(labels=unique(filgatk_bin$xname), limits = c(5, 20)) + 
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Gatk SV site frequency percentage") +
  scale_x_continuous(breaks = pretty(filgatk$min_n, n = 100))
ppgatk_bin_small
#tiff("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/gatk_percent.tiff", width=5, height=5, units="in", res=300)
ppgatk_bin_small
dev.off()



      lumpy_test$Freq_expect_percent <- lumpy_test$Freq_expect/sum(lumpy_test$Freq_expect, na.rm = T)
      fillp2_bin <- merge(fillp2_bin, lumpyexpect, by.x = "min_n", by.y = "Var1", all.x = T, all.y = F)
      plumpy_bin <- ggplot(fillp2_bin, aes(x = min_n, y = pct_n, fill = svtype, label = pct_n)) +
        geom_bar(stat = "identity") + 
        labs(y = "SV abundance percentage", x = "numbers of SV") +
        scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                                 "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
        theme_bw(base_size = 12) + scale_x_continuous(breaks = pretty(fillp2_bin$min_n, n = 100)) + theme(axis.text.x = element_text(angle = 90)) +
        ggtitle("Lumpy SV site frequency percentage") +
        geom_text(aes(label = round(n, 1)), color = "black", size = 1, position = position_stack(vjust = 0.5), srt = 90) +
        geom_line(aes(x = min_n, y = sumexpect), linetype = "dotted") 
      
      plumpy_bin
      
      tiff("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy_percent_n_with_BND_withnum.tiff", width=10, height=5, units="in", res=300)
      plumpy_bin
      dev.off()
      
      
      plumpy_bin_dodge <- ggplot(fillp2_bin, aes(x = min_n, y = pct_n, fill = svtype, label = pct_n)) +
        geom_bar(stat = "identity", position="dodge", width=0.9) + 
        labs(y = "SV abundance percentage", x = "n") +
        scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                                 "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
        theme_bw(base_size = 12) + scale_x_continuous(breaks = pretty(fillp2_bin$min_n, n = 100)) + theme(axis.text.x = element_text(angle = 90)) +
        ggtitle("Lumpy SV site frequency percentage")+
        geom_line(aes(x = min_n, y = Freq_expect_percent.x), linetype = "dotted") 
      
      plumpy_bin_dodge
      tiff("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy_percent_n_dodge_withBND.tiff", width=10, height=5, units="in", res=300)
      plumpy_bin_dodge
      dev.off()
      
      
#check with Zach
dim(fillp2)    
fillp2[,1:10]      
table(fillp2$svtype)

sobic.004G301500
Chr04:64031808..64032785
sobic.004G301600
Chr04:64036221..64037462
sobic.004G301650
Chr04:64038396..64039002


