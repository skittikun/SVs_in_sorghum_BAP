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

options(stringsAsFactors = F)
split = detectCores()
cl = makeCluster(split)
registerDoParallel(cl)

startline <- function(filepath){
  con = file(filepath, "r")
  readstart <- 1
  while ( T ) {
    line = readLines(con, n = 1)
    if(grepl("#CHROM", line) ) {
      
      break
    }
    readstart <- readstart+1
  }
  close(con)
  return(readstart)
}

dt <- "/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/zeroone_352.txt"
dat <- read.table(dt, quote = "", fill = T, header = T, sep = "\t", na.strings=c("","NA"), stringsAsFactors=FALSE, comment.char = "", check.names = F)#, skip = (startline(dt)-1))
dat$INFO <- gsub(";IMPRECISE", ";IMPRECISE=1", gsub(";SECONDARY", ";SECONDARY=1", dat$INFO))

          #checkdat
          dat$prpos <- NULL
          dat$prend <- NULL
          dat[1:10,]
          sum(dat$pe >= 3)
          hist(dat$pe)
          dupdat <- dat[dat$svtype == "DUP",]
          dupdat[1:10,]
          mafcut <- 0.05
          missingcut <- 0.2 #0.05
          sizecut <- 100000
          pecut <- 3
          
          dat$end <- ifelse(dat$svtype == "BND", dat$pos, dat$end)
          dat$svlen <- dat$end -dat$pos
    
          fildt <- dat[dat$maf_col >= mafcut & dat$missing_col <= missingcut & dat$svlen <= sizecut & dat$pe >= pecut, ]
          dim(fildt)
          fildupdat <- fildt[fildt$svtype == "DUP",]
          fildupdat[1:10,]
          hist(fildupdat$svlen)
          #

#split backbone from d1

saminfo <- gregexpr(";", sample(dat$INFO, 10))
lsam <- c()
for(j in 1:length(saminfo)){
  lsam <- c(lsam, length(saminfo[[j]]))
}
max(lsam)

strinput <- dat$INFO[2]

forsep_df <- c()
forname_df <- c()
for(i in 1:length(gregexpr(";", strinput)[[1]])){
  allpos <- gregexpr(";", strinput)[[1]]
  allend <- (gregexpr("=", strinput)[[1]][-1])
  allpos > allend
  strpos <- allpos[i]
  strend <- allend[i]
  forsep <- substr(strinput, strpos, strend)
  forname <- tolower(substr(strinput, strpos+1, strend-1))
  forsep_df <- c(forsep_df, forsep)
  forname_df <- c(forname_df, forname)
}
forsep_df
forname_df

sdat <- dat
for(i in length(forsep_df):1){
  sdat <- sdat %>% separate(INFO, c("INFO", forname_df[i]), sep = forsep_df[i])
  print(paste(i, i/length(forsep_df), sep = "_"))
  }

sdat <- sdat %>% separate(INFO, c("INFO", "end"), sep = ";END=") %>% separate(INFO, c("svtype", "svlen"), sep = ";SVLEN=") 
sdat$svtype <- gsub("SVTYPE=", "", sdat$svtype)
colnames(sdat)[c(1:7, 27)] <- c("chr", "pos", "ID", "ref", "alt", "qual", "filter", "format")
sdat[1:20, 1:29]

#k = (which(colnames(sdat) == "format")+1)
#sdat[,k]

write.table(sdat, "/projects/cooper_research1/ksongsom/speedseq_result/TERRA/vcf/merge/terra_sort_merge_CNV_352_add.vcf", row.names = F, sep = "\t", quote = F)

#convert to 0/1
depthcutoff <- 10
for(k in (which(colnames(sdat) == "format")+1):ncol(sdat)){
  homalttf <- grepl("1/1", sdat[,colnames(sdat)[k]])
  tmptab <- (as.data.frame(do.call(rbind, str_split(sdat[,colnames(sdat)[k]], ":" )),stringsAsFactors=FALSE))
  deptheach <- tmptab[,5]
  deptheachtf <- deptheach >= depthcutoff
  sdat[,colnames(sdat)[k]] <- ifelse(homalttf & deptheachtf, 1, ifelse(!homalttf & deptheachtf, 0, NA))
  print(paste(k, k/ncol(sdat), sep = "_"))
}

head(sdat)
sum(sdat[,(which(colnames(sdat) == "format")+1):ncol(sdat)] == 1, na.rm = T)

MAF_function <- function(input){
  aa <- as.numeric(input)
  tab_aa <- as.data.frame(table(aa))
  maf <- round(min(tab_aa$Freq)/sum(tab_aa$Freq), digits = 2)
  maf <- ifelse(maf == Inf, NA, ifelse(maf == 1, 0, maf))
  return(maf)
}

maf_col <- apply(sdat[,(which(colnames(sdat) == "format")+1):ncol(sdat)], 1, MAF_function)
sum(maf_col >= 0.05, na.rm = T)

missing_col <- round(rowSums(is.na(sdat[,(which(colnames(sdat) == "format")+1):ncol(sdat)]))/length((which(colnames(sdat) == "format")+1):ncol(sdat)), digits = 2)
sum(missing_col <= 0.20, na.rm = T)
sum(maf_col >= 0.05 & missing_col <= 0.20, na.rm = T)

alldat <- cbind.data.frame(sdat, maf_col, missing_col)
dim(alldat)
fil_sdat <- sdat[maf_col >= 0.05 & missing_col <= 0.20,]
dim(fil_sdat)

table(fil_sdat$svtype)
#write.table(fil_sdat, "/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/fil_352.txt", row.names = F, sep = "\t", quote = F)
write.table(fil_sdat, "/projects/cooper_research1/ksongsom/speedseq_result/TERRA/vcf/merge/fil_352.txt", row.names = F, sep = "\t", quote = F)
write.table(alldat, "/projects/cooper_research1/ksongsom/speedseq_result/TERRA/vcf/merge/zeroone_352.txt", row.names = F, sep = "\t", quote = F)


dt <- "/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/zeroone_352.txt"
dat <- read.table(dt, quote = "", fill = T, header = T, sep = "\t", na.strings=c("","NA"), stringsAsFactors=FALSE, comment.char = "", check.names = F)
sum(dat$maf_col >= 0.05, na.rm = T)
sum(dat$missing_col <= 0.2, na.rm = T)
sum(dat$maf_col >= 0.05 & dat$missing_col <= 0.2, na.rm = T)
dat[is.na(dat$end),]$end <- dat[is.na(dat$end),]$pos + 100000
dat[is.na(dat$pos),]$pos <- dat[is.na(dat$pos),]$end + 100000
dat$svlen <- dat$end - dat$pos

mafcut <- 0.05
missingcut <- 0.2 #0.05
sizecut <- 100000
fildt <- dat[dat$maf_col >= mafcut & dat$missing_col <= missingcut, ]
fildt$end <- ifelse(fildt$svtype == "BND", fildt$pos, fildt$end)
fildt$svlen <- fildt$end -fildt$pos
fildt <- fildt[is.na(fildt$imprecise),]
fildt <- fildt[fildt$pe >= 3,]
fildt <- fildt[fildt$svlen <= sizecut, ]
dim(fildt)
head(fildt)
sum(is.na(fildt$end))
sum(is.na(fildt$pos))

head(dat)
altref <- cbind.data.frame(c("Chr04", "Chr04", "Chr04"), c("Sobic.004G301500", "Sobic.004G301600", "Sobic.004G301650"), c(64031808, 64036221, 64038396), c(64032785, 64037462, 64039002)) %>% setNames(c("chr", "gene", "pos", "end"))

#editrange <- dat[dat$end >= dat$pos, ]
#rangedat <- GRanges(seqnames = editrange$chr, ranges = IRanges(editrange$pos, editrange$end))
#rangefildat <- GRanges(seqnames = fildt$chr, ranges = IRanges(fildt$pos, fildt$end))
#rangealtref <- GRanges(seqnames = altref$chr, ranges = IRanges(altref$pos, altref$end))

#altoverdat <- findOverlapPairs(rangedat, rangealtref)
#altoverdatfil <- findOverlapPairs(rangefildat, rangealtref)

vitoverfildt <- c()
for(r in 1:nrow(altref)){
  #vitoverfildt <- rbind.data.frame(vitoverfildt, fildt[fildt$chr == altref$chr[r] & fildt$pos <= altref$end[r] & fildt$end >= altref$pos[r], ])
  vitoverfildt <- rbind.data.frame(vitoverfildt, dat[dat$chr == altref$chr[r] & dat$pos <= altref$end[r] & dat$end >= altref$pos[r], ])
  vitoverfildt <- vitoverfildt[!duplicated(vitoverfildt),]
}

vitfil <- (vitoverfildt[vitoverfildt$maf_col >= mafcut & abs(vitoverfildt$svlen) <= sizecut & vitoverfildt$missing_col <= missingcut,])
dim(vitfil)
table(vitfil$svtype)
sumvit <- colSums(vitfil[, (which(colnames(vitfil) == "format")+1):(which(colnames(vitfil) == "rio"))], na.rm = T)
sum(sumvit == 0)
sum(sumvit != 0)

vitfil[vitfil$svtyp != "BND",]
sumvit2 <- colSums(vitfil[vitfil$svtyp != "BND",(which(colnames(vitfil) == "format")+1):(which(colnames(vitfil) == "rio"))], na.rm = T)
sum(sumvit2 == 0)
sum(sumvit2 != 0)

summary(vitoverfildt$maf_col)
summary(vitoverfildt$missing_col)

sumvitdt <- cbind.data.frame(gsub("", "", names(sumvit)), as.data.frame(sumvit)) %>% setNames(c("PI", "sumvit"))

fil_sdat <- read.table("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/fil_352.txt", header = T, sep = "\t", na.strings=c("","NA"), stringsAsFactors=FALSE, comment.char = "", check.names = F)
dim(fil_sdat)
str(fil_sdat)
table(fil_sdat$svtype)



#add function
#gff <- read.table("/Users/ksongsom/OneDrive/postdoc/SV/illumina/gene_anno.txt", quote = "", fill = T, header = T, sep = "\t",na.strings=c("","NA") ,stringsAsFactors=FALSE,check.names=FALSE)
gff <- read.table("/users/ksongsom/Ref_Genomes/BTx/gene_anno.txt", quote = "", fill = T, header = T, sep = "\t",na.strings=c("","NA") ,stringsAsFactors=FALSE,check.names=FALSE)
pn <- fil_sdat

#with BND has no end add the length for 100kb
pn[is.na(pn$end),]$end <- pn[is.na(pn$end),]$pos + 100000
pn[is.na(pn$pos),]$pos <- pn[is.na(pn$pos),]$end + 100000

crops <-c()

mrna <- gff[gff$feature == "mRNA",]
cds <- gff[gff$feature == "CDS",]
five <- gff[gff$feature =="five-prime_UTR",]
gene <- gff[gff$feature =="gene",]
three <- gff[gff$feature =="three_prime_UTR",]

mRNA <- c()
for(i in 1:dim(pn)[1]){
  mRNA_name <- paste(mrna[mrna$Chr == pn[i,]$chr & (mrna$end >= pn[i,]$pos & mrna$start <= pn[i,]$end),]$Name, collapse= "& ")
  
  panther <- paste(unique(na.omit(mrna[mrna$Chr == pn[i,]$chr & (mrna$end >= pn[i,]$pos & mrna$start <= pn[i,]$end),]$Panther)), collapse= "& ")
  kog <- paste(unique(na.omit(mrna[mrna$Chr == pn[i,]$chr & (mrna$end >= pn[i,]$pos & mrna$start <= pn[i,]$end),]$KOG)), collapse= "& ")
  ec <- paste(unique(na.omit(mrna[mrna$Chr == pn[i,]$chr & (mrna$end >= pn[i,]$pos & mrna$start <= pn[i,]$end),]$ec)), collapse= "& ")
  go <- paste(unique(na.omit(mrna[mrna$Chr == pn[i,]$chr & (mrna$end >= pn[i,]$pos & mrna$start <= pn[i,]$end),]$GO)), collapse= "& ")
  
  arabi_name <- paste(unique(na.omit(mrna[mrna$Chr == pn[i,]$chr & (mrna$end >= pn[i,]$pos & mrna$start <= pn[i,]$end),]$Best.hit.arabi.name)), collapse= "& ")
  arabi_defline <- paste(unique(na.omit(mrna[mrna$Chr == pn[i,]$chr & (mrna$end >= pn[i,]$pos & mrna$start <= pn[i,]$end),]$arabi.defline)), collapse= "& ")
  rice_name <- paste(unique(na.omit(mrna[mrna$Chr == pn[i,]$chr & (mrna$end >= pn[i,]$pos & mrna$start <= pn[i,]$end),]$Best.hit.rice.name)), collapse= "& ")
  rice_defline <- paste(unique(na.omit(mrna[mrna$Chr == pn[i,]$chr & (mrna$end >= pn[i,]$pos & mrna$start <= pn[i,]$end),]$rice.defline)), collapse= "& ")
  mRNA <- rbind(mRNA, cbind(mRNA_name, panther, kog, ec, go, arabi_name, arabi_defline, rice_name, rice_defline))
  print(paste(crops, "mrna", i/dim(pn)[1], sep = "_"))    
}


CDS <- c()
for(i in 1:dim(pn)[1]){
  meff <- paste(cds[cds$Chr == pn[i,]$chr & (cds$end >= pn[i,]$pos & cds$start <= pn[i,]$end),]$ID, collapse= "& ")
  CDS <- append(CDS, meff)
  print(paste(crops, "CDS", i/dim(pn)[1], sep = "_"))  
}

fiveUTR <- c()
for(i in 1:dim(pn)[1]){
  meff <- paste(five[five$Chr == pn[i,]$chr & (five$end >= pn[i,]$pos & five$start <= pn[i,]$end),]$ID, collapse= "& ")
  fiveUTR <- append(fiveUTR, meff)
  print(paste(crops, "fiveUTR", i/dim(pn)[1], sep = "_"))  
}

genes <- c()
for(i in 1:dim(pn)[1]){
  meff <- paste(gene[gene$Chr == pn[i,]$chr & (gene$end >= pn[i,]$pos & gene$start <= pn[i,]$end),]$Name, collapse= "& ")
  genes <- append(genes, meff)
  print(paste(crops, "genes", i/dim(pn)[1], sep = "_"))  
}

threeUTR <- c()
for(i in 1:dim(pn)[1]){
  meff <- paste(three[three$Chr == pn[i,]$chr & (three$end >= pn[i,]$pos & three$start <= pn[i,]$end),]$ID, collapse= "& ")
  threeUTR <- append(threeUTR, meff)
  print(paste(crops, "three", i/dim(pn)[1], sep = "_"))  
}

#cover 80%

coverage <- 0.8

mRNA80 <- c()
for(i in 1:dim(pn)[1]){
  mRNA_name80 <- paste(mrna[mrna$Chr == pn[i,]$chr & (mrna$end >= pn[i,]$pos & mrna$start <= pn[i,]$end),]$Name, collapse= "& ")
  
  panther80 <- paste(unique(na.omit(mrna[mrna$Chr == pn[i,]$chr & (mrna$end >= pn[i,]$pos & mrna$start <= pn[i,]$end) & (pn[i,]$end-mrna$start >= coverage*mrna$len | mrna$end - pn[i,]$pos >= coverage*mrna$len),]$Panther)), collapse= "& ")
  kog80 <- paste(unique(na.omit(mrna[mrna$Chr == pn[i,]$chr & (mrna$end >= pn[i,]$pos & mrna$start <= pn[i,]$end) & (pn[i,]$end-mrna$start >= coverage*mrna$len | mrna$end - pn[i,]$pos >= coverage*mrna$len),]$KOG)), collapse= "& ")
  ec80 <- paste(unique(na.omit(mrna[mrna$Chr == pn[i,]$chr & (mrna$end >= pn[i,]$pos & mrna$start <= pn[i,]$end) & (pn[i,]$end-mrna$start >= coverage*mrna$len | mrna$end - pn[i,]$pos >= coverage*mrna$len),]$ec)), collapse= "& ")
  go80 <- paste(unique(na.omit(mrna[mrna$Chr == pn[i,]$chr & (mrna$end >= pn[i,]$pos & mrna$start <= pn[i,]$end) & (pn[i,]$end-mrna$start >= coverage*mrna$len | mrna$end - pn[i,]$pos >= coverage*mrna$len),]$GO)), collapse= "& ")
  
  arabi_name80 <- paste(unique(na.omit(mrna[mrna$Chr == pn[i,]$chr & (mrna$end >= pn[i,]$pos & mrna$start <= pn[i,]$end) & (pn[i,]$end-mrna$start >= coverage*mrna$len | mrna$end - pn[i,]$pos >= coverage*mrna$len),]$Best.hit.arabi.name)), collapse= "& ")
  arabi_defline80 <- paste(unique(na.omit(mrna[mrna$Chr == pn[i,]$chr & (mrna$end >= pn[i,]$pos & mrna$start <= pn[i,]$end & (pn[i,]$end-mrna$start >= coverage*mrna$len | mrna$end - pn[i,]$pos >= coverage*mrna$len)),]$arabi.defline)), collapse= "& ")
  rice_name80 <- paste(unique(na.omit(mrna[mrna$Chr == pn[i,]$chr & (mrna$end >= pn[i,]$pos & mrna$start <= pn[i,]$end) & (pn[i,]$end-mrna$start >= coverage*mrna$len | mrna$end - pn[i,]$pos >= coverage*mrna$len),]$Best.hit.rice.name)), collapse= "& ")
  rice_defline80 <- paste(unique(na.omit(mrna[mrna$Chr == pn[i,]$chr & (mrna$end >= pn[i,]$pos & mrna$start <= pn[i,]$end & (pn[i,]$end-mrna$start >= coverage*mrna$len | mrna$end - pn[i,]$pos >= coverage*mrna$len)),]$rice.defline)), collapse= "& ")
  mRNA80 <- rbind(mRNA80, cbind(mRNA_name80, panther80, kog80, ec80, go80, arabi_name80, arabi_defline80, rice_name80, rice_defline80))
  print(paste(crops, "mrna80", i/dim(pn)[1], sep = "_"))  
}

CDS80 <- c()
for(i in 1:dim(pn)[1]){
  meff <- paste(cds[cds$Chr == pn[i,]$chr & (cds$end >= pn[i,]$pos & cds$start <= pn[i,]$end) & (pn[i,]$end-cds$start >= coverage*cds$len | cds$end - pn[i,]$pos >= coverage*cds$len),]$ID, collapse= "& ")
  CDS80 <- append(CDS80, meff)
  print(paste(crops, "CDS80", i/dim(pn)[1], sep = "_"))  
}

fiveUTR80 <- c()
for(i in 1:dim(pn)[1]){
  meff <- paste(five[five$Chr == pn[i,]$chr & (five$end >= pn[i,]$pos & five$start <= pn[i,]$end) & (pn[i,]$end-five$start >= coverage*five$len | five$end - pn[i,]$pos >= coverage*five$len),]$ID, collapse= "& ")
  fiveUTR80 <- append(fiveUTR80, meff)
  print(paste(crops, "fiveUTR80", i/dim(pn)[1], sep = "_"))  
}

genes80 <- c()
for(i in 1:dim(pn)[1]){
  meff <- paste(gene[gene$Chr == pn[i,]$chr & (gene$end >= pn[i,]$pos & gene$start <= pn[i,]$end) & (pn[i,]$end-gene$start >= coverage*gene$len | gene$end - pn[i,]$pos >= coverage*gene$len),]$Name, collapse= "& ")
  genes80 <- append(genes80, meff)
  print(paste(crops, "genes80", i/dim(pn)[1], sep = "_"))  
}

threeUTR80 <- c()
for(i in 1:dim(pn)[1]){
  meff <- paste(three[three$Chr == pn[i,]$chr & (three$end >= pn[i,]$pos & three$start <= pn[i,]$end) & (pn[i,]$end-three$start >= coverage*three$len | three$end - pn[i,]$pos >= coverage*three$len),]$ID, collapse= "& ")
  threeUTR80 <- append(threeUTR80, meff)
  print(paste(crops, "threeUTR80", i/dim(pn)[1], sep = "_"))  
}

tab180 <- cbind(pn, genes, CDS, fiveUTR, threeUTR, mRNA, genes80, CDS80, fiveUTR80, threeUTR80, mRNA80)

head(tab180)
dim(tab180)
tab180$chr
#tab180 <- tab180[!duplicated(tab180),]
#save overlap
#write.table(tab180,paste("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/352_go180.txt", sep = ""), row.names = F, sep = "\t", quote = F)
write.table(tab180,paste("/projects/cooper_research1/ksongsom/speedseq_result/TERRA/vcf/merge/352_go180.txt", sep = ""), row.names = F, sep = "\t", quote = F)

print("DONE")

tab180 <- read.csv("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/352_go180.txt", header = T, stringsAsFactors = F, na.strings=c("","NA"), sep = "\t", check.names = F)
tab180f <- tab180[abs(tab180$svlen) <= 100000,]
#check sugar
#sugargene <- read.csv("/Users/ksongsom/OneDrive/postdoc/VIT/sugargenes_clean.csv", header = T, stringsAsFactors = F, na.strings=c("","NA"), sep = "\t")
sugargene <- read.csv("/projects/cooper_research1/ksongsom/speedseq_result/TERRA/vcf/merge/sugargenes_clean.csv", header = T, stringsAsFactors = F, na.strings=c("","NA"), sep = "\t")
sugargene$genes
tab180$genes
sum(grepl(paste(sugargene$genes, collapse = "|"), tab180$genes))

tab180sugar <- tab180f[grepl(paste(sugargene$genes, collapse = "|"), tab180f$genes),]
dim(tab180sugar)
table(tab180sugar$svtype)
hist(tab180sugar$svlen)
