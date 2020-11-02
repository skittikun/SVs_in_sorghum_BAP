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
library(stringi)
library(stringr)
library(lme4)
library(lmerTest)
library(seqinr)
library(vcfR)
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


lumpy <- readRDS("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/zeroone_352_fit.RDS")

#change format and add columns
lumpy <- lumpy[grepl("Chr", lumpy$chr),]
lumpy$end <- ifelse(lumpy$svtype == "BND", lumpy$pos, lumpy$end)
lumpy$svlen <- lumpy$end - lumpy$pos

#remove individual with all NA
lumpy <- lumpy[!(colnames(lumpy) %in% c("455307", "562985"))]

#filter criteria
mafcut <- 0.05
missingcut <- 0.2
sizecut <- 100000
sizecut_l <- 50
fillp <- lumpy[lumpy$maf_col >= mafcut & lumpy$missing_col <= missingcut & lumpy$svlen <= sizecut & lumpy$svlen >= sizecut_l, ]
fillp <- fillp[is.na(fillp$imprecise),]
fillp$svlen <- fillp$end -fillp$pos
fillp <- fillp[abs(fillp$svlen) <= sizecut,]
fillp <- fillp[abs(fillp$svlen) >= sizecut_l,]

##table1 for paper

#remove imprecise
lumpy_f1 <- lumpy[is.na(lumpy$imprecise),]
dim(lumpy_f1)
table(lumpy_f1$svtype)

#remove BND
lumpy_f2 <- lumpy_f1[lumpy_f1$svtype != "BND",]

#cut of size between 50 to 100k bp 
lumpy_f3 <- lumpy_f2[lumpy_f2$svlen <= sizecut & lumpy_f2$svlen >= sizecut_l,]
dim(lumpy_f3)
table(lumpy_f3$svtype)

#remove rows with missing value higher than .20
lumpy_f4 <- lumpy_f3[lumpy_f3$missing_col <= missingcut,]
dim(lumpy_f4)
table(lumpy_f4$svtype)

#remove rows with maf less than .05
lumpy_f5 <- lumpy_f4[lumpy_f4$maf_col >= mafcut,]
dim(lumpy_f5)
table(lumpy_f5$svtype)
summary(lumpy_f5$svlen)

#Updated table 1
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

#plot lumpy
lumpyindividual <- colnames(lumpy_f4)[(which(colnames(lumpy_f4)=="144134"):which(colnames(lumpy_f4)=="rio"))]

#harmonic sum for Watterson
harsumf <- function(input){
  harsum <- 0
  for(i in 1:input){
    harsum <- sum(harsum, 1/i)
  }
  return(harsum)
}
#sorghum chromosome size
sorghumgenome <- as.data.frame(cbind(c(paste("Chr0", 1:9, sep = ""), "Chr10"), rep(1, 10),c(80884392, 77742459, 74386277, 68658214, 71854669, 61277060, 65505356,62686529,59416394,61233695)))
colnames(sorghumgenome) <- c("Chromosome", "start", "end")
sorghumgenome$Chromosome <- as.factor(sorghumgenome$Chromosome)
sorghumgenome$start <- as.numeric(as.character(sorghumgenome$start))
sorghumgenome$end <- as.numeric(as.character(sorghumgenome$end))
genomelength <- sum(sorghumgenome$end)




#test and plot lumpy
fillp4 <- lumpy_f4
dim(fillp4)

#remove sv with whole row 
onlygeno <- fillp4[,(which(colnames(fillp4)=="144134"):which(colnames(fillp4)=="rio"))]
head(onlygeno)
noalternative <- c()
for(i in 1:nrow(onlygeno)){
  rrr <- onlygeno[i, ]
  noalternative[i] <- ifelse(length(apply(rrr, 1, table))==1, 1, NA)
  print(paste(i, i/nrow(onlygeno), sep= "_"))
}

fillp4 <- fillp4[is.na(noalternative),]
dim(fillp4) #after removing no variation 18711 SV

indilength <- length(lumpyindividual)



#calculate minor allele amount in each site  
#min_n_count

MINcal <- function(rr){
  #rr <- fillp4[901,(which(colnames(fillp4)=="144134"):which(colnames(fillp4)=="rio"))] for testing
  #make missing value -9
  #rr[is.na(rr)] <- -9
  predat <- table(t(rr[rr != -1]))
  minor <- as.numeric(names(which(predat == min(predat))[1]))
  halfway <- round(indilength/2, 0)
  if(length(predat) != 1){
    
    #maf <- sum(rr == minor, na.rm = T)}else{
    #  if(length(predat) == 1 & sum(rr == minor, na.rm = T) > halfway){
    #    maf = sum(rr == minor, na.rm = T)-(sum(predat)/2)
    #      #indilength-sum(rr == minor, na.rm = T)
    #      #halfway-(sum(rr == minor, na.rm = T)-halfway)
    #  }else{if(sum(rr == minor, na.rm = T) == halfway){
    #    maf = (sum(predat)/2)#halfway
    #  }else{maf = NA}}
    #}
  
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
table(min_n_count)

fillp4$min_n_count <- min_n_count

observe_lumpy <- as.data.frame(table(fillp4$min_n_count))
observe_lumpy$Var1 <- as.numeric(observe_lumpy$Var1)
dim(observe_lumpy)

#Watterson’s estimator of Theta
w_lp <- (nrow(fillp4)/harsumf((2*indilength-1)))/(genomelength)

size_l <- (seq(1, indilength/2, 1))
size_l <- c(size_l, last(size_l)+1)
size_l <- as.numeric(observe_lumpy$Var1)
expect_lumpy <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpy <- rbind(expect_lumpy, cbind.data.frame(count, (w_lp*genomelength)*(1/count)))  
}
colnames(expect_lumpy) <- c("Var1", "Freq")
expect_lumpy
sum(expect_lumpy$Freq)
nrow(fillp4)


lumpy_test <- merge(y = observe_lumpy, x = expect_lumpy, by = "Var1", all.x = T, all.y = T, suffixes = c("_expect","_observe"))
head(lumpy_test)
tail(lumpy_test)
lumpy_test[is.na(lumpy_test$Freq_expect),]
lumpy_test[is.na(lumpy_test)] <- 0
#lumpy_test <- na.omit(lumpy_test)
chisq.test(lumpy_test[,"Freq_observe"]/10, p = lumpy_test[,"Freq_expect"]/sum(lumpy_test[,"Freq_expect"]))

ks.test(x= lumpy_test[,"Freq_observe"], y= lumpy_test[,"Freq_expect"])$p.val
ks.test(x= lumpy_test[,"Freq_observe"], y= lumpy_test[,"Freq_expect"]/sum(lumpy_test[,"Freq_expect"]))

observeinput <- lumpy_test[,"Freq_observe"]
expectinput <- lumpy_test[,"Freq_expect"]
chi2 <- sum(((observeinput-expectinput)^2)/expectinput)
degreef <- length(observeinput)-1
qchisq(.95, df=7)
chi2
qchisq(.95, df=degreef)

#lumpy percentage 
#fillp4$bin <- cut_interval(abs(fillp4$maf_col), length=0.005)
fillp4_bin <- as.data.frame(fillp4 %>%
                              group_by(min_n_count, svtype) %>%
                              tally  %>%
                              group_by(svtype) %>%
                              mutate(pct_n=round((100*n)/sum(n), 2)))
dim(fillp4)
dim(fillp4_bin)

#test by type del
#separate each type
#calculate expectation by each type
lpdel <- fillp4_bin[fillp4_bin$svtype == "DEL",]
w_lpdel <- (sum(fillp4_bin[fillp4_bin$svtype == "DEL",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpydel <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpydel <- rbind(expect_lumpydel, cbind.data.frame(count, (w_lpdel*genomelength)*(1/count)))  
}
colnames(expect_lumpydel) <- c("Var1", "Freq")
expect_lumpydel
observe_lumpydel <- fillp4_bin[fillp4_bin$svtype == "DEL",c("min_n_count", "n" )] %>% setNames(c("Var1", "Freq"))
dim(observe_lumpydel)
lumpy_testdel <- merge(observe_lumpydel, expect_lumpydel, by = "Var1", all.x = T, all.y = T, suffixes = c("_observe", "_expect"))
head(lumpy_testdel)
lumpy_testdel <- na.omit(lumpy_testdel)
chisq.test(lumpy_testdel[,"Freq_observe"], p = lumpy_testdel[,"Freq_expect"]/sum(lumpy_testdel[,"Freq_expect"]))

##bnd
lpbnd <- fillp4_bin[fillp4_bin$svtype == "BND",]
w_lpbnd <- (sum(fillp4_bin[fillp4_bin$svtype == "BND",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpybnd <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpybnd <- rbind(expect_lumpybnd, cbind.data.frame(count, (w_lpbnd*genomelength)*(1/count)))  
}
colnames(expect_lumpybnd) <- c("Var1", "Freq")
expect_lumpybnd
observe_lumpybnd <- fillp4_bin[fillp4_bin$svtype == "BND",c("min_n_count", "n" )] %>% setNames(c("Var1", "Freq"))
dim(observe_lumpybnd)
lumpy_testbnd <- merge(observe_lumpybnd, expect_lumpybnd, by = "Var1", all.x = T, all.y = T, suffixes = c("_observe", "_expect"))
head(lumpy_testbnd)
lumpy_testbnd <- na.omit(lumpy_testbnd)
chisq.test(lumpy_testbnd[,"Freq_observe"], p = lumpy_testbnd[,"Freq_expect"]/sum(lumpy_testbnd[,"Freq_expect"]))

##dup
lpdup <- fillp4_bin[fillp4_bin$svtype == "DUP",]
w_lpdup <- (sum(fillp4_bin[fillp4_bin$svtype == "DUP",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpydup <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpydup <- rbind(expect_lumpydup, cbind.data.frame(count, (w_lpdup*genomelength)*(1/count)))  
}
colnames(expect_lumpydup) <- c("Var1", "Freq")
expect_lumpydup
observe_lumpydup <- fillp4_bin[fillp4_bin$svtype == "DUP",c("min_n_count", "n" )] %>% setNames(c("Var1", "Freq"))
dim(observe_lumpydup)
lumpy_testdup <- merge(observe_lumpydup, expect_lumpydup, by = "Var1", all.x = T, all.y = T, suffixes = c("_observe", "_expect"))
head(lumpy_testdup)
lumpy_testdup <- na.omit(lumpy_testdup)
chisq.test(lumpy_testdup[,"Freq_observe"], p = lumpy_testdup[,"Freq_expect"]/sum(lumpy_testdup[,"Freq_expect"]))

##inv
lpinv <- fillp4_bin[fillp4_bin$svtype == "INV",]
w_lpinv <- (sum(fillp4_bin[fillp4_bin$svtype == "INV",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpyinv <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpyinv <- rbind(expect_lumpyinv, cbind.data.frame(count, (w_lpinv*genomelength)*(1/count)))  
}
colnames(expect_lumpyinv) <- c("Var1", "Freq")
expect_lumpyinv
observe_lumpyinv <- fillp4_bin[fillp4_bin$svtype == "INV",c("min_n_count", "n" )] %>% setNames(c("Var1", "Freq"))
dim(observe_lumpyinv)
lumpy_testinv <- merge(observe_lumpyinv, expect_lumpyinv, by = "Var1", all.x = T, all.y = T, suffixes = c("_observe", "_expect"))
head(lumpy_testinv)
lumpy_testinv <- na.omit(lumpy_testinv)
chisq.test(lumpy_testinv[,"Freq_observe"], p = lumpy_testinv[,"Freq_expect"]/sum(lumpy_testinv[,"Freq_expect"]))

#calculate expected percentage from each type of sv
lumpy_testinv$Freq_expect_percent <- lumpy_testinv$Freq_expect*100/sum(lumpy_testinv$Freq_expect, na.rm = T)
lumpy_testdel$Freq_expect_percent <- lumpy_testdel$Freq_expect*100/sum(lumpy_testdel$Freq_expect, na.rm = T)
lumpy_testdup$Freq_expect_percent <- lumpy_testdup$Freq_expect*100/sum(lumpy_testdup$Freq_expect, na.rm = T)

#test on each type
chisq.test(lumpy_testdel$Freq_observe, p = lumpy_testdel$Freq_expect_percent/100)$p.value #sig
chisq.test(lumpy_testdup$Freq_observe, p = lumpy_testdup$Freq_expect_percent/100)$p.value #sig
chisq.test(lumpy_testinv$Freq_observe, p = lumpy_testinv$Freq_expect_percent/100)$p.value #not sig

#lumpy_testbnd$Freq_expect_percent <- lumpy_testbnd$Freq_expect*100/sum(lumpy_testbnd$Freq_expect, na.rm = T)

lumpyexpect <- Reduce(function(x, y) merge(x, y, all=TRUE, by = "Var1"), list(lumpy_testinv[,c(1,4)], lumpy_testdel[,c(1,4)], lumpy_testdup[,c(1,4)]))
head(lumpyexpect)
lumpyexpect$sumexpect <- rowSums(lumpyexpect[,-1],na.rm = T)
colnames(lumpyexpect) <- c("Var1", "Freq.inv", "Freq.del", "Freq.dup", "sumexpect")

#expectation calculation
#Watterson’s estimator (\Theta_w) 
#del
w_lpdel <- (sum(fillp4_bin[fillp4_bin$svtype == "DEL",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpydel <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpydel <- rbind(expect_lumpydel, cbind.data.frame(count, (w_lpdel*genomelength)*(1/count)))  
}
colnames(expect_lumpydel) <- c("Var1", "Freq")

#bnd
lpbnd <- fillp4_bin[fillp4_bin$svtype == "BND",]
w_lpbnd <- (sum(fillp4_bin[fillp4_bin$svtype == "BND",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpybnd <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpybnd <- rbind(expect_lumpybnd, cbind.data.frame(count, (w_lpbnd*genomelength)*(1/count)))  
}
colnames(expect_lumpybnd) <- c("Var1", "Freq")

##dup
lpdup <- fillp4_bin[fillp4_bin$svtype == "DUP",]
w_lpdup <- (sum(fillp4_bin[fillp4_bin$svtype == "DUP",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpydup <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpydup <- rbind(expect_lumpydup, cbind.data.frame(count, (w_lpdup*genomelength)*(1/count)))  
}
colnames(expect_lumpydup) <- c("Var1", "Freq")

##inv
lpinv <- fillp4_bin[fillp4_bin$svtype == "INV",]
w_lpinv <- (sum(fillp4_bin[fillp4_bin$svtype == "INV",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpyinv <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpyinv <- rbind(expect_lumpyinv, cbind.data.frame(count, (w_lpinv*genomelength)*(1/count)))  
}
colnames(expect_lumpyinv) <- c("Var1", "Freq")

expect_lumpydel$per <- expect_lumpydel$Freq*100/sum(expect_lumpydel$Freq, na.rm = T)
expect_lumpydup$per <- expect_lumpydup$Freq*100/sum(expect_lumpydup$Freq, na.rm = T)
expect_lumpyinv$per <- expect_lumpyinv$Freq*100/sum(expect_lumpyinv$Freq, na.rm = T)
lumpyexpect_each <- Reduce(function(x, y) merge(x, y, all=TRUE, by = "Var1"), list(expect_lumpydel, expect_lumpydup, expect_lumpyinv))
colnames(lumpyexpect_each) <- c("Var1", "Freq.del", "per.del", "Freq.dup", "per.dup", "Freq.inv", "per.inv")

       
#no BND
       
        #lumpy percentage 
        #fillp4$bin <- cut_interval(abs(fillp4$maf_col), length=0.005)
        fillp4$bin <- cut_interval(abs(fillp4$min_n_count), length=0.005)
        fillp4_bin <- as.data.frame(fillp4 %>%
                                      group_by(min_n_count, svtype) %>%
                                      tally  %>%
                                      group_by(svtype) %>%
                                      mutate(pct_n=round((100*n)/sum(n), 2)))
        lumpyexpect_each$sumexpect <- lumpyexpect_each$per.del*3
        fillp42_bin <- merge(fillp4_bin, lumpyexpect_each, by.x = "min_n_count", by.y = "Var1", all.x = T, all.y = T)
        fillp42_bin$sumexpect_nobnd <- fillp42_bin$sumexpect
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
          scale_x_continuous(breaks = pretty(fillp42_bin$min_n_count, n = 50)+1) 
        plumpy_bin
        
        tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/draft2/Figures/Fig3_lumpy_percent.tiff", width=10, height=5, units="in", res=300)
        plumpy_bin
        dev.off()
        
        #test grand 
        lumpy_test <- na.omit(lumpy_test)
        chisq.test(lumpy_test[,"Freq_observe"], p = lumpy_test[,"Freq_expect"]/sum(lumpy_test[,"Freq_expect"]))$p.value #sig
        chisq.test(lumpy_test[,"Freq_observe"]/10, p = lumpy_test[,"Freq_expect"]/sum(lumpy_test[,"Freq_expect"]))$p.value
        #test on each type
        chisq.test(lumpy_testdel$Freq_observe/10, p = lumpy_testdel$Freq_expect_percent/100)$p.value #sig
        chisq.test(lumpy_testdup$Freq_observe/10, p = lumpy_testdup$Freq_expect_percent/100)$p.value #sig
        chisq.test(lumpy_testinv$Freq_observe/10, p = lumpy_testinv$Freq_expect_percent/100)$p.value #not sig
        
        chisq.test(fillp42_bin[fillp42_bin$svtype== "DEL",]$n, 
                   p = fillp42_bin[fillp42_bin$svtype== "DEL",]$per.del)
        table(fillp4$svtype)

        ks.test(x = lumpy_test[,"Freq_observe"], y = lumpy_test[,"Freq_expect"])$p.value
        ks.test(x = lumpy_testdel$Freq_observe, y = lumpy_testdel$Freq_expect)$p.value
        ks.test(x = lumpy_testdup$Freq_observe, y = lumpy_testdup$Freq_expect)$p.value
        ks.test(x = lumpy_testinv$Freq_observe, y = lumpy_testinv$Freq_expect)$p.value
        
        #by type
        plumpy_bin_del <- ggplot(fillp42_bin[fillp42_bin$svtype == "DEL",], aes(x = min_n_count, y = pct_n, fill = svtype, label = pct_n)) +
          geom_bar(stat = "identity") + 
          labs(y = "SV abundance percentage", x = "Numbers of SV") +
          scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                                   "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
          theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
          ggtitle("Lumpy SV site frequency percentage (all genotypes)") +
          geom_line(aes(x = min_n_count, y = sumexpect_nobnd/3), linetype = "dotted", col = "black") +
          geom_text(aes(label = round(n, 1)), color = "black", size = 1, position = position_stack(vjust = 0.5)) +
          scale_x_continuous(breaks = pretty(fillp42_bin$min_n_count, n = 50)+1) 
        plumpy_bin_del
        
        tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/draft2/Figures/Supplement1_del_lumpy_percent.tiff", width=10, height=5, units="in", res=300)
        plumpy_bin_del
        dev.off()
        
        
        plumpy_bin_dup <- ggplot(fillp42_bin[fillp42_bin$svtype == "DUP",], aes(x = min_n_count, y = pct_n, fill = svtype, label = pct_n)) +
          geom_bar(stat = "identity") + 
          labs(y = "SV abundance percentage", x = "Numbers of SV") +
          scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                                   "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
          theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
          ggtitle("Lumpy SV site frequency percentage (all genotypes)") +
          geom_line(aes(x = min_n_count, y = sumexpect_nobnd/3), linetype = "dotted", col = "black") +
          geom_text(aes(label = round(n, 1)), color = "black", size = 1, position = position_stack(vjust = 0.5)) +
          scale_x_continuous(breaks = pretty(fillp42_bin$min_n_count, n = 50)+1) 
        plumpy_bin_dup
        
        tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/draft2/Figures/Supplement1_dup_lumpy_percent.tiff", width=10, height=5, units="in", res=300)
        plumpy_bin_dup
        dev.off()

        plumpy_bin_inv <- ggplot(fillp42_bin[fillp42_bin$svtype == "INV",], aes(x = min_n_count, y = pct_n, fill = svtype, label = pct_n)) +
          geom_bar(stat = "identity") + 
          labs(y = "SV abundance percentage", x = "Numbers of SV") +
          scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                                   "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
          theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
          ggtitle("Lumpy SV site frequency percentage (all genotypes)") +
          geom_line(aes(x = min_n_count, y = sumexpect_nobnd/3), linetype = "dotted", col = "black") +
          geom_text(aes(label = round(n, 1)), color = "black", size = 1, position = position_stack(vjust = 0.5)) +
          scale_x_continuous(breaks = pretty(fillp42_bin$min_n_count, n = 50)+1) 
        plumpy_bin_inv
        
        tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/draft2/Figures/Supplement1_inv_lumpy_percent.tiff", width=10, height=5, units="in", res=300)
        plumpy_bin_inv
        dev.off()        
        
        #need to run after "site_frequency_spectrum"

library(GenomicRanges)
library(IRanges)
#frequency spectrum gene
#load gff
gff <- readRDS("/Users/ksongsom/OneDrive/postdoc/SV/illumina/gene_anno.rds")
gene <- gff[gff$feature == "gene", ]
cds <- gff[gff$feature == "CDS", ]

genegrange <- GRanges(gene$Chr, IRanges(start=as.numeric(gene$start),end=as.numeric(gene$end)))
cdsgrange <- GRanges(cds$Chr, IRanges(start=as.numeric(cds$start),end=as.numeric(cds$end)))
                     
lumpy_f4$svid <- paste(lumpy_f4$chr, lumpy_f4$pos, lumpy_f4$end, lumpy_f4$svtype, sep = "_")

ovlp <- lumpy_f4[,c("svid", "chr", "pos", "end")] %>% purrr::set_names(c("svid", "Chr", "start", "end"))
gr <- makeGRangesFromDataFrame(ovlp, keep.extra.columns = TRUE)

geneoverlapnames <- subsetByOverlaps(gr, genegrange)$svid
cdsoverlapnames <- subsetByOverlaps(gr, cdsgrange)$svid

genetab <- lumpy_f4[lumpy_f4$svid %in% geneoverlapnames,]
cdstab <- lumpy_f4[lumpy_f4$svid %in% cdsoverlapnames,]
dim(genetab)
dim(cdstab)


#remove sv with whole row 
onlygeno <- genetab[,(which(colnames(genetab)=="144134"):which(colnames(genetab)=="rio"))]
head(onlygeno)
noalternative <- c()
for(i in 1:nrow(onlygeno)){
  rrr <- onlygeno[i, ]
  noalternative[i] <- ifelse(length(apply(rrr, 1, table))==1, 1, NA)
  print(paste(i, i/nrow(onlygeno), sep= "_"))
}

genetab <- genetab[is.na(noalternative),]
#remove sv with whole row 
onlygeno <- cdstab[,(which(colnames(cdstab)=="144134"):which(colnames(cdstab)=="rio"))]
head(onlygeno)
noalternative <- c()
for(i in 1:nrow(onlygeno)){
  rrr <- onlygeno[i, ]
  noalternative[i] <- ifelse(length(apply(rrr, 1, table))==1, 1, NA)
  print(paste(i, i/nrow(onlygeno), sep= "_"))
}

cdstab <- cdstab[is.na(noalternative),]

#test and plot lumpy GENES
fillp4 <- genetab
dim(fillp4)

indilength <- length(lumpyindividual)

#Watterson’s estimator of Theta
w_lp <- (nrow(fillp4)/harsumf((2*indilength-1)))/(genomelength)

size_l <- (seq(1, indilength/2, 1))
size_l <- c(size_l, last(size_l)+1)
expect_lumpy <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpy <- rbind(expect_lumpy, cbind.data.frame(count, (w_lp*genomelength)*(1/count)))  
}
colnames(expect_lumpy) <- c("Var1", "Freq")
expect_lumpy
sum(expect_lumpy$Freq)
nrow(fillp4)

#calculate minor allele amount in each site  
#min_n_count

MINcal <- function(rr){
  #rr <- fillp4[901,(which(colnames(fillp4)=="144134"):which(colnames(fillp4)=="rio"))] for testing
  #make missing value -9
  #rr[is.na(rr)] <- -9
  predat <- table(t(rr[rr != -1]))
  minor <- as.numeric(names(which(predat == min(predat))[1]))
  halfway <- round(indilength/2, 0)
  if(length(predat) != 1){
    
    #maf <- sum(rr == minor, na.rm = T)}else{
    #  if(length(predat) == 1 & sum(rr == minor, na.rm = T) > halfway){
    #    maf = sum(rr == minor, na.rm = T)-(sum(predat)/2)
    #      #indilength-sum(rr == minor, na.rm = T)
    #      #halfway-(sum(rr == minor, na.rm = T)-halfway)
    #  }else{if(sum(rr == minor, na.rm = T) == halfway){
    #    maf = (sum(predat)/2)#halfway
    #  }else{maf = NA}}
    #}
    
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
table(min_n_count)

fillp4$min_n_count <- min_n_count
#fillp4$min_n <- round2(fillp4$maf_col*length(lumpyindividual), 0)


observe_lumpy <- as.data.frame(table(fillp4$min_n_count))
dim(observe_lumpy)
lumpy_test <- merge(y = observe_lumpy, x = expect_lumpy, by = "Var1", all.x = T, all.y = F, suffixes = c("_expect","_observe"))
head(lumpy_test)
tail(lumpy_test)
lumpy_test[is.na(lumpy_test$Freq_expect),]

lumpy_test <- na.omit(lumpy_test)
chisq.test(lumpy_test[,"Freq_observe"], p = lumpy_test[,"Freq_expect"]/sum(lumpy_test[,"Freq_expect"]))


#lumpy percentage 
#fillp4$bin <- cut_interval(abs(fillp4$maf_col), length=0.005)
fillp4_bin <- as.data.frame(fillp4 %>%
                              group_by(min_n_count, svtype) %>%
                              tally  %>%
                              group_by(svtype) %>%
                              mutate(pct_n=round((100*n)/sum(n), 2)))
dim(fillp4)
dim(fillp4_bin)

#test by type del
#separate each type
#calculate expectation by each type
lpdel <- fillp4_bin[fillp4_bin$svtype == "DEL",]
w_lpdel <- (sum(fillp4_bin[fillp4_bin$svtype == "DEL",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpydel <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpydel <- rbind(expect_lumpydel, cbind.data.frame(count, (w_lpdel*genomelength)*(1/count)))  
}
colnames(expect_lumpydel) <- c("Var1", "Freq")
expect_lumpydel
observe_lumpydel <- fillp4_bin[fillp4_bin$svtype == "DEL",c("min_n_count", "n" )] %>% setNames(c("Var1", "Freq"))
dim(observe_lumpydel)
lumpy_testdel <- merge(observe_lumpydel, expect_lumpydel, by = "Var1", all.x = T, all.y = T, suffixes = c("_observe", "_expect"))
head(lumpy_testdel)
lumpy_testdel <- na.omit(lumpy_testdel)
chisq.test(lumpy_testdel[,"Freq_observe"], p = lumpy_testdel[,"Freq_expect"]/sum(lumpy_testdel[,"Freq_expect"]))

##bnd
lpbnd <- fillp4_bin[fillp4_bin$svtype == "BND",]
w_lpbnd <- (sum(fillp4_bin[fillp4_bin$svtype == "BND",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpybnd <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpybnd <- rbind(expect_lumpybnd, cbind.data.frame(count, (w_lpbnd*genomelength)*(1/count)))  
}
colnames(expect_lumpybnd) <- c("Var1", "Freq")
expect_lumpybnd
observe_lumpybnd <- fillp4_bin[fillp4_bin$svtype == "BND",c("min_n_count", "n" )] %>% setNames(c("Var1", "Freq"))
dim(observe_lumpybnd)
lumpy_testbnd <- merge(observe_lumpybnd, expect_lumpybnd, by = "Var1", all.x = T, all.y = T, suffixes = c("_observe", "_expect"))
head(lumpy_testbnd)
lumpy_testbnd <- na.omit(lumpy_testbnd)
chisq.test(lumpy_testbnd[,"Freq_observe"], p = lumpy_testbnd[,"Freq_expect"]/sum(lumpy_testbnd[,"Freq_expect"]))

##dup
lpdup <- fillp4_bin[fillp4_bin$svtype == "DUP",]
w_lpdup <- (sum(fillp4_bin[fillp4_bin$svtype == "DUP",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpydup <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpydup <- rbind(expect_lumpydup, cbind.data.frame(count, (w_lpdup*genomelength)*(1/count)))  
}
colnames(expect_lumpydup) <- c("Var1", "Freq")
expect_lumpydup
observe_lumpydup <- fillp4_bin[fillp4_bin$svtype == "DUP",c("min_n_count", "n" )] %>% setNames(c("Var1", "Freq"))
dim(observe_lumpydup)
lumpy_testdup <- merge(observe_lumpydup, expect_lumpydup, by = "Var1", all.x = T, all.y = T, suffixes = c("_observe", "_expect"))
head(lumpy_testdup)
lumpy_testdup <- na.omit(lumpy_testdup)
chisq.test(lumpy_testdup[,"Freq_observe"], p = lumpy_testdup[,"Freq_expect"]/sum(lumpy_testdup[,"Freq_expect"]))

##inv
lpinv <- fillp4_bin[fillp4_bin$svtype == "INV",]
w_lpinv <- (sum(fillp4_bin[fillp4_bin$svtype == "INV",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpyinv <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpyinv <- rbind(expect_lumpyinv, cbind.data.frame(count, (w_lpinv*genomelength)*(1/count)))  
}
colnames(expect_lumpyinv) <- c("Var1", "Freq")
expect_lumpyinv
observe_lumpyinv <- fillp4_bin[fillp4_bin$svtype == "INV",c("min_n_count", "n" )] %>% setNames(c("Var1", "Freq"))
dim(observe_lumpyinv)
lumpy_testinv <- merge(observe_lumpyinv, expect_lumpyinv, by = "Var1", all.x = T, all.y = T, suffixes = c("_observe", "_expect"))
head(lumpy_testinv)
lumpy_testinv <- na.omit(lumpy_testinv)
chisq.test(lumpy_testinv[,"Freq_observe"], p = lumpy_testinv[,"Freq_expect"]/sum(lumpy_testinv[,"Freq_expect"]))

#calculate expected percentage from each type of sv
lumpy_testinv$Freq_expect_percent <- lumpy_testinv$Freq_expect*100/sum(lumpy_testinv$Freq_expect, na.rm = T)
lumpy_testdel$Freq_expect_percent <- lumpy_testdel$Freq_expect*100/sum(lumpy_testdel$Freq_expect, na.rm = T)
lumpy_testdup$Freq_expect_percent <- lumpy_testdup$Freq_expect*100/sum(lumpy_testdup$Freq_expect, na.rm = T)

#test on each type
chisq.test(lumpy_testdel$Freq_observe, p = lumpy_testdel$Freq_expect_percent/100)$p.value #sig
chisq.test(lumpy_testdup$Freq_observe, p = lumpy_testdup$Freq_expect_percent/100)$p.value #not sig
chisq.test(lumpy_testinv$Freq_observe, p = lumpy_testinv$Freq_expect_percent/100)$p.value #not sig

#lumpy_testbnd$Freq_expect_percent <- lumpy_testbnd$Freq_expect*100/sum(lumpy_testbnd$Freq_expect, na.rm = T)

lumpyexpect <- Reduce(function(x, y) merge(x, y, all=TRUE, by = "Var1"), list(lumpy_testinv[,c(1,4)], lumpy_testdel[,c(1,4)], lumpy_testdup[,c(1,4)]))
head(lumpyexpect)
lumpyexpect$sumexpect <- rowSums(lumpyexpect[,-1],na.rm = T)
colnames(lumpyexpect) <- c("Var1", "Freq.inv", "Freq.del", "Freq.dup", "sumexpect")

#expectation calculation
#Watterson’s estimator (\Theta_w) 
#del
w_lpdel <- (sum(fillp4_bin[fillp4_bin$svtype == "DEL",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpydel <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpydel <- rbind(expect_lumpydel, cbind.data.frame(count, (w_lpdel*genomelength)*(1/count)))  
}
colnames(expect_lumpydel) <- c("Var1", "Freq")

#bnd
lpbnd <- fillp4_bin[fillp4_bin$svtype == "BND",]
w_lpbnd <- (sum(fillp4_bin[fillp4_bin$svtype == "BND",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpybnd <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpybnd <- rbind(expect_lumpybnd, cbind.data.frame(count, (w_lpbnd*genomelength)*(1/count)))  
}
colnames(expect_lumpybnd) <- c("Var1", "Freq")

##dup
lpdup <- fillp4_bin[fillp4_bin$svtype == "DUP",]
w_lpdup <- (sum(fillp4_bin[fillp4_bin$svtype == "DUP",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpydup <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpydup <- rbind(expect_lumpydup, cbind.data.frame(count, (w_lpdup*genomelength)*(1/count)))  
}
colnames(expect_lumpydup) <- c("Var1", "Freq")

##inv
lpinv <- fillp4_bin[fillp4_bin$svtype == "INV",]
w_lpinv <- (sum(fillp4_bin[fillp4_bin$svtype == "INV",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpyinv <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpyinv <- rbind(expect_lumpyinv, cbind.data.frame(count, (w_lpinv*genomelength)*(1/count)))  
}
colnames(expect_lumpyinv) <- c("Var1", "Freq")

expect_lumpydel$per <- expect_lumpydel$Freq*100/sum(expect_lumpydel$Freq, na.rm = T)
expect_lumpydup$per <- expect_lumpydup$Freq*100/sum(expect_lumpydup$Freq, na.rm = T)
expect_lumpyinv$per <- expect_lumpyinv$Freq*100/sum(expect_lumpyinv$Freq, na.rm = T)
lumpyexpect_each <- Reduce(function(x, y) merge(x, y, all=TRUE, by = "Var1"), list(expect_lumpydel, expect_lumpydup, expect_lumpyinv))
colnames(lumpyexpect_each) <- c("Var1", "Freq.del", "per.del", "Freq.dup", "per.dup", "Freq.inv", "per.inv")


#no BND

#lumpy percentage 
#fillp4$bin <- cut_interval(abs(fillp4$maf_col), length=0.005)
fillp4$bin <- cut_interval(abs(fillp4$min_n_count), length=0.005)
fillp4_bin <- as.data.frame(fillp4 %>%
                              group_by(min_n_count, svtype) %>%
                              tally  %>%
                              group_by(svtype) %>%
                              mutate(pct_n=round((100*n)/sum(n), 2)))
lumpyexpect_each$sumexpect <- lumpyexpect_each$per.del*3
fillp42_bin <- merge(fillp4_bin, lumpyexpect_each, by.x = "min_n_count", by.y = "Var1", all.x = T, all.y = T)
fillp42_bin$sumexpect_nobnd <- fillp42_bin$sumexpect
fillp42_bin <- fillp42_bin[fillp42_bin$min_n_count != 0,]
plumpy_bin_gene <- ggplot(fillp42_bin, aes(x = min_n_count, y = pct_n, fill = svtype, label = pct_n)) +
  geom_bar(stat = "identity") + 
  labs(y = "SV abundance percentage", x = "Numbers of SV") +
  scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                           "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Lumpy SV (genes) site frequency percentage (all genotypes)") +
  geom_line(aes(x = fillp42_bin$min_n_count, y = fillp42_bin$sumexpect_nobnd), linetype = "dotted", col = "black") +
  geom_text(aes(label = round(n, 1)), color = "black", size = 1, position = position_stack(vjust = 0.5)) +
  scale_x_continuous(breaks = pretty(fillp42_bin$min_n_count, n = 50)+1) 
plumpy_bin_gene

tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/draft2/Figures/Supplement2_lumpy_percent_gene.tiff", width=10, height=5, units="in", res=300)
plumpy_bin_gene
dev.off()

#test grand 
chisq.test(lumpy_test[,"Freq_observe"]/5, p = lumpy_test[,"Freq_expect"]/sum(lumpy_test[,"Freq_expect"]))$p.value #sig

#test on each type
chisq.test(lumpy_testdel$Freq_observe/5, p = lumpy_testdel$Freq_expect_percent/100)$p.value #sig
chisq.test(lumpy_testdup$Freq_observe/5, p = lumpy_testdup$Freq_expect_percent/100)$p.value #sig
chisq.test(lumpy_testinv$Freq_observe/5, p = lumpy_testinv$Freq_expect_percent/100)$p.value #not sig

chisq.test(fillp42_bin[fillp42_bin$svtype== "DEL",]$n, 
           p = fillp42_bin[fillp42_bin$svtype== "DEL",]$per.del)
table(fillp4$svtype)


ks.test(x = lumpy_test[,"Freq_observe"], y = lumpy_test[,"Freq_expect"])$p.value
ks.test(x = lumpy_testdel$Freq_observe, y = lumpy_testdel$Freq_expect)$p.value
ks.test(x = lumpy_testdup$Freq_observe, y = lumpy_testdup$Freq_expect)$p.value
ks.test(x = lumpy_testinv$Freq_observe, y = lumpy_testinv$Freq_expect)$p.value

















































































#test and plot lumpy CDS
fillp4 <- cdstab
dim(fillp4)

indilength <- length(lumpyindividual)

#Watterson’s estimator of Theta
w_lp <- (nrow(fillp4)/harsumf((2*indilength-1)))/(genomelength)

size_l <- (seq(1, indilength/2, 1))
size_l <- c(size_l, last(size_l)+1)
expect_lumpy <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpy <- rbind(expect_lumpy, cbind.data.frame(count, (w_lp*genomelength)*(1/count)))  
}
colnames(expect_lumpy) <- c("Var1", "Freq")
expect_lumpy
sum(expect_lumpy$Freq)
nrow(fillp4)

#calculate minor allele amount in each site  
#min_n_count
min_n_count <- c()
for(i in 1:nrow(fillp4)){
  min_n_count <- c(min_n_count, MINcal(fillp4[i,(which(colnames(fillp4)=="144134"):which(colnames(fillp4)=="rio"))]))
  print(paste(i, i/nrow(fillp4), sep = "_"))
}
table(min_n_count)

fillp4$min_n_count <- min_n_count
#fillp4$min_n <- round2(fillp4$maf_col*length(lumpyindividual), 0)


observe_lumpy <- as.data.frame(table(fillp4$min_n_count))
dim(observe_lumpy)
lumpy_test <- merge(y = observe_lumpy, x = expect_lumpy, by = "Var1", all.x = T, all.y = F, suffixes = c("_expect","_observe"))
head(lumpy_test)
tail(lumpy_test)
lumpy_test[is.na(lumpy_test$Freq_expect),]

lumpy_test <- na.omit(lumpy_test)
chisq.test(lumpy_test[,"Freq_observe"], p = lumpy_test[,"Freq_expect"]/sum(lumpy_test[,"Freq_expect"]))


#lumpy percentage 
#fillp4$bin <- cut_interval(abs(fillp4$maf_col), length=0.005)
fillp4_bin <- as.data.frame(fillp4 %>%
                              group_by(min_n_count, svtype) %>%
                              tally  %>%
                              group_by(svtype) %>%
                              mutate(pct_n=round((100*n)/sum(n), 2)))
dim(fillp4)
dim(fillp4_bin)

#test by type del
#separate each type
#calculate expectation by each type
lpdel <- fillp4_bin[fillp4_bin$svtype == "DEL",]
w_lpdel <- (sum(fillp4_bin[fillp4_bin$svtype == "DEL",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpydel <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpydel <- rbind(expect_lumpydel, cbind.data.frame(count, (w_lpdel*genomelength)*(1/count)))  
}
colnames(expect_lumpydel) <- c("Var1", "Freq")
expect_lumpydel
observe_lumpydel <- fillp4_bin[fillp4_bin$svtype == "DEL",c("min_n_count", "n" )] %>% setNames(c("Var1", "Freq"))
dim(observe_lumpydel)
lumpy_testdel <- merge(observe_lumpydel, expect_lumpydel, by = "Var1", all.x = T, all.y = T, suffixes = c("_observe", "_expect"))
head(lumpy_testdel)
lumpy_testdel <- na.omit(lumpy_testdel)
chisq.test(lumpy_testdel[,"Freq_observe"], p = lumpy_testdel[,"Freq_expect"]/sum(lumpy_testdel[,"Freq_expect"]))

##bnd
lpbnd <- fillp4_bin[fillp4_bin$svtype == "BND",]
w_lpbnd <- (sum(fillp4_bin[fillp4_bin$svtype == "BND",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpybnd <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpybnd <- rbind(expect_lumpybnd, cbind.data.frame(count, (w_lpbnd*genomelength)*(1/count)))  
}
colnames(expect_lumpybnd) <- c("Var1", "Freq")
expect_lumpybnd
observe_lumpybnd <- fillp4_bin[fillp4_bin$svtype == "BND",c("min_n_count", "n" )] %>% setNames(c("Var1", "Freq"))
dim(observe_lumpybnd)
lumpy_testbnd <- merge(observe_lumpybnd, expect_lumpybnd, by = "Var1", all.x = T, all.y = T, suffixes = c("_observe", "_expect"))
head(lumpy_testbnd)
lumpy_testbnd <- na.omit(lumpy_testbnd)
chisq.test(lumpy_testbnd[,"Freq_observe"], p = lumpy_testbnd[,"Freq_expect"]/sum(lumpy_testbnd[,"Freq_expect"]))

##dup
lpdup <- fillp4_bin[fillp4_bin$svtype == "DUP",]
w_lpdup <- (sum(fillp4_bin[fillp4_bin$svtype == "DUP",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpydup <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpydup <- rbind(expect_lumpydup, cbind.data.frame(count, (w_lpdup*genomelength)*(1/count)))  
}
colnames(expect_lumpydup) <- c("Var1", "Freq")
expect_lumpydup
observe_lumpydup <- fillp4_bin[fillp4_bin$svtype == "DUP",c("min_n_count", "n" )] %>% setNames(c("Var1", "Freq"))
dim(observe_lumpydup)
lumpy_testdup <- merge(observe_lumpydup, expect_lumpydup, by = "Var1", all.x = T, all.y = T, suffixes = c("_observe", "_expect"))
head(lumpy_testdup)
lumpy_testdup <- na.omit(lumpy_testdup)
chisq.test(lumpy_testdup[,"Freq_observe"], p = lumpy_testdup[,"Freq_expect"]/sum(lumpy_testdup[,"Freq_expect"]))

##inv
lpinv <- fillp4_bin[fillp4_bin$svtype == "INV",]
w_lpinv <- (sum(fillp4_bin[fillp4_bin$svtype == "INV",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpyinv <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpyinv <- rbind(expect_lumpyinv, cbind.data.frame(count, (w_lpinv*genomelength)*(1/count)))  
}
colnames(expect_lumpyinv) <- c("Var1", "Freq")
expect_lumpyinv
observe_lumpyinv <- fillp4_bin[fillp4_bin$svtype == "INV",c("min_n_count", "n" )] %>% setNames(c("Var1", "Freq"))
dim(observe_lumpyinv)
lumpy_testinv <- merge(observe_lumpyinv, expect_lumpyinv, by = "Var1", all.x = T, all.y = T, suffixes = c("_observe", "_expect"))
head(lumpy_testinv)
lumpy_testinv <- na.omit(lumpy_testinv)
chisq.test(lumpy_testinv[,"Freq_observe"], p = lumpy_testinv[,"Freq_expect"]/sum(lumpy_testinv[,"Freq_expect"]))

#calculate expected percentage from each type of sv
lumpy_testinv$Freq_expect_percent <- lumpy_testinv$Freq_expect*100/sum(lumpy_testinv$Freq_expect, na.rm = T)
lumpy_testdel$Freq_expect_percent <- lumpy_testdel$Freq_expect*100/sum(lumpy_testdel$Freq_expect, na.rm = T)
lumpy_testdup$Freq_expect_percent <- lumpy_testdup$Freq_expect*100/sum(lumpy_testdup$Freq_expect, na.rm = T)

#test on each type
chisq.test(lumpy_testdel$Freq_observe, p = lumpy_testdel$Freq_expect_percent/100)$p.value #sig
chisq.test(lumpy_testdup$Freq_observe, p = lumpy_testdup$Freq_expect_percent/100)$p.value #not sig
chisq.test(lumpy_testinv$Freq_observe, p = lumpy_testinv$Freq_expect_percent/100)$p.value #not sig

#lumpy_testbnd$Freq_expect_percent <- lumpy_testbnd$Freq_expect*100/sum(lumpy_testbnd$Freq_expect, na.rm = T)

lumpyexpect <- Reduce(function(x, y) merge(x, y, all=TRUE, by = "Var1"), list(lumpy_testinv[,c(1,4)], lumpy_testdel[,c(1,4)], lumpy_testdup[,c(1,4)]))
head(lumpyexpect)
lumpyexpect$sumexpect <- rowSums(lumpyexpect[,-1],na.rm = T)
colnames(lumpyexpect) <- c("Var1", "Freq.inv", "Freq.del", "Freq.dup", "sumexpect")

#expectation calculation
#Watterson’s estimator (\Theta_w) 
#del
w_lpdel <- (sum(fillp4_bin[fillp4_bin$svtype == "DEL",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpydel <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpydel <- rbind(expect_lumpydel, cbind.data.frame(count, (w_lpdel*genomelength)*(1/count)))  
}
colnames(expect_lumpydel) <- c("Var1", "Freq")

#bnd
lpbnd <- fillp4_bin[fillp4_bin$svtype == "BND",]
w_lpbnd <- (sum(fillp4_bin[fillp4_bin$svtype == "BND",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpybnd <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpybnd <- rbind(expect_lumpybnd, cbind.data.frame(count, (w_lpbnd*genomelength)*(1/count)))  
}
colnames(expect_lumpybnd) <- c("Var1", "Freq")

##dup
lpdup <- fillp4_bin[fillp4_bin$svtype == "DUP",]
w_lpdup <- (sum(fillp4_bin[fillp4_bin$svtype == "DUP",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpydup <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpydup <- rbind(expect_lumpydup, cbind.data.frame(count, (w_lpdup*genomelength)*(1/count)))  
}
colnames(expect_lumpydup) <- c("Var1", "Freq")

##inv
lpinv <- fillp4_bin[fillp4_bin$svtype == "INV",]
w_lpinv <- (sum(fillp4_bin[fillp4_bin$svtype == "INV",]$n, na.rm = T)/harsumf((2*indilength-1)))/(genomelength)
size_l <- (seq(1, indilength/2, 1))#unique(indilength*as.numeric(fillumpy_bin$xname))
expect_lumpyinv <- c()
for(i in 1:length(size_l)){
  count <- size_l[i]
  expect_lumpyinv <- rbind(expect_lumpyinv, cbind.data.frame(count, (w_lpinv*genomelength)*(1/count)))  
}
colnames(expect_lumpyinv) <- c("Var1", "Freq")

expect_lumpydel$per <- expect_lumpydel$Freq*100/sum(expect_lumpydel$Freq, na.rm = T)
expect_lumpydup$per <- expect_lumpydup$Freq*100/sum(expect_lumpydup$Freq, na.rm = T)
expect_lumpyinv$per <- expect_lumpyinv$Freq*100/sum(expect_lumpyinv$Freq, na.rm = T)
lumpyexpect_each <- Reduce(function(x, y) merge(x, y, all=TRUE, by = "Var1"), list(expect_lumpydel, expect_lumpydup, expect_lumpyinv))
colnames(lumpyexpect_each) <- c("Var1", "Freq.del", "per.del", "Freq.dup", "per.dup", "Freq.inv", "per.inv")


#no BND

#lumpy percentage 
#fillp4$bin <- cut_interval(abs(fillp4$maf_col), length=0.005)
fillp4$bin <- cut_interval(abs(fillp4$min_n_count), length=0.005)
fillp4_bin <- as.data.frame(fillp4 %>%
                              group_by(min_n_count, svtype) %>%
                              tally  %>%
                              group_by(svtype) %>%
                              mutate(pct_n=round((100*n)/sum(n), 2)))
lumpyexpect_each$sumexpect <- lumpyexpect_each$per.del*3
fillp42_bin <- merge(fillp4_bin, lumpyexpect_each, by.x = "min_n_count", by.y = "Var1", all.x = T, all.y = T)
fillp42_bin$sumexpect_nobnd <- fillp42_bin$sumexpect
fillp42_bin <- fillp42_bin[fillp42_bin$min_n_count != 0,]
plumpy_bin_cds <- ggplot(fillp42_bin, aes(x = min_n_count, y = pct_n, fill = svtype, label = pct_n)) +
  geom_bar(stat = "identity") + 
  labs(y = "SV abundance percentage", x = "Numbers of SV") +
  scale_fill_manual("SV types", values = c("DEL" = "red", "INS" = "blue",  "BND" = "purple",
                                           "DUP" = "green", "RPL" = "grey", "INV" = "gold")) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Lumpy SV (CDS) site frequency percentage (all genotypes)") +
  geom_line(aes(x = fillp42_bin$min_n_count, y = fillp42_bin$sumexpect_nobnd), linetype = "dotted", col = "black") +
  geom_text(aes(label = round(n, 1)), color = "black", size = 1, position = position_stack(vjust = 0.5)) +
  scale_x_continuous(breaks = pretty(fillp42_bin$min_n_count, n = 50)+1) 
plumpy_bin_cds

tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/draft2/Figures/Supplement2_lumpy_percent_cds.tiff", width=10, height=5, units="in", res=300)
plumpy_bin_cds
dev.off()

#test grand 
chisq.test(lumpy_test[,"Freq_observe"]/5, p = lumpy_test[,"Freq_expect"]/sum(lumpy_test[,"Freq_expect"]))$p.value #sig

#test on each type
chisq.test(lumpy_testdel$Freq_observe/5, p = lumpy_testdel$Freq_expect_percent/100)$p.value #sig
chisq.test(lumpy_testdup$Freq_observe/5, p = lumpy_testdup$Freq_expect_percent/100)$p.value #sig
chisq.test(lumpy_testinv$Freq_observe/5, p = lumpy_testinv$Freq_expect_percent/100)$p.value #not sig

chisq.test(fillp42_bin[fillp42_bin$svtype== "DEL",]$n, 
           p = fillp42_bin[fillp42_bin$svtype== "DEL",]$per.del)
table(fillp4$svtype)

ks.test(x = lumpy_test[,"Freq_observe"], y = "dexp")$p.value
ks.test(x = lumpy_test[,"Freq_observe"], y = lumpy_test[,"Freq_expect"])$p.value
ks.test(x = lumpy_testdel$Freq_observe, y = lumpy_testdel$Freq_expect)$p.value
ks.test(x = lumpy_testdel$Freq_observe, y = "dexp")$p.value
ks.test(x = lumpy_testdup$Freq_observe, y = lumpy_testdup$Freq_expect)$p.value
ks.test(x = lumpy_testinv$Freq_observe, y = lumpy_testinv$Freq_expect)$p.value
