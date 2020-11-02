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

dt <- "/projects/cooper_research1/ksongsom/speedseq_result/TERRA/vcf/merge/352nocheck.vcf"

#dt <- "/projects/cooper_research1/ksongsom/speedseq_result/TERRA/vcf/merge/352precise.vcf"
dat <- read.table(dt, quote = "", fill = T, header = T, sep = "\t", na.strings=c("","NA"), stringsAsFactors=FALSE, comment.char = "", check.names = F, skip = (startline(dt)-1))
dat$INFO <- gsub(";IMPRECISE", ";IMPRECISE=1", gsub(";SECONDARY", ";SECONDARY=1", dat$INFO))

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
}

sdat <- sdat %>% separate(INFO, c("INFO", "end"), sep = ";END=") %>% separate(INFO, c("svtype", "svlen"), sep = ";SVLEN=") 
sdat$svtype <- gsub("SVTYPE=", "", sdat$svtype)
colnames(sdat)[c(1:7)] <- c("chr", "pos", "ID", "ref", "alt", "qual", "filter")
colnames(sdat) <- tolower(colnames(sdat))

#pe cutoff
pecutoff <- 3
sdat <- sdat[sdat$pe >= pecutoff,]

k = (which(colnames(sdat) == "format")+1)
k

sdat <- sdat[,!(colnames(sdat) %in% c("prpos", "prend"))]

#convert to 0/1 and -1
depthcutoff <- 10
for(k in (which(colnames(sdat) == "format")+1):ncol(sdat)){
oneone <- grepl("1/1", sdat[,colnames(sdat)[k]])
zeroone <- grepl("0/1", sdat[,colnames(sdat)[k]])
zerozero <- grepl("0/0", sdat[,colnames(sdat)[k]])
tmptab <- (as.data.frame(do.call(rbind, str_split(sdat[,colnames(sdat)[k]], ":" )),stringsAsFactors=FALSE))
deptheach <- tmptab[,5]
deptheachtf <- deptheach >= depthcutoff
sdat[,colnames(sdat)[k]] <- ifelse(oneone & deptheachtf, 1, ifelse(zerozero & deptheachtf, -1, ifelse(zeroone & deptheachtf, 0, NA)))
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

#add end for BND
sdat$end <- ifelse(sdat$svtype == "BND", sdat$pos, sdat$end)

alldat <- cbind.data.frame(sdat, maf_col, missing_col)
write.table(alldat, "/projects/cooper_research1/ksongsom/speedseq_result/TERRA/vcf/merge/zeroone_352_pruned_pe3_dp10.txt", row.names = F, sep = "\t", quote = F)

fil_sdat <- sdat[maf_col >= 0.05 & missing_col <= 0.20,]
#write.table(fil_sdat, "/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/fil_352.txt", row.names = F, sep = "\t", quote = F)
write.table(fil_sdat, "/projects/cooper_research1/ksongsom/speedseq_result/TERRA/vcf/merge/fil_352_nocheck.txt", row.names = F, sep = "\t", quote = F)



#add function
#gff <- read.table("/Users/ksongsom/OneDrive/postdoc/SV/illumina/gene_anno.txt", quote = "", fill = T, header = T, sep = "\t",na.strings=c("","NA") ,stringsAsFactors=FALSE,check.names=FALSE)
gff <- read.table("/users/ksongsom/Ref_Genomes/BTx/gene_anno.txt", quote = "", fill = T, header = T, sep = "\t",na.strings=c("","NA") ,stringsAsFactors=FALSE,check.names=FALSE)
pn <- fil_sdat

#with BND has no end add the length for 100kb
#pn[is.na(pn$end),]$end <- pn[is.na(pn$end),]$pos + 100000
#pn[is.na(pn$pos),]$pos <- pn[is.na(pn$pos),]$end + 100000
pn$end <- ifelse(pn$svtype == "BND", pn$pos, pn$end)

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
write.table(tab180,paste("/projects/cooper_research1/ksongsom/speedseq_result/TERRA/vcf/merge/352_go180_nocheck.txt", sep = ""), row.names = F, sep = "\t", quote = F)

print("DONE")


#check sugar
#sugargene <- read.csv("/Users/ksongsom/OneDrive/postdoc/VIT/sugargenes_clean.csv", header = T, stringsAsFactors = F, na.strings=c("","NA"), sep = "\t")
sugargene <- read.csv("/projects/cooper_research1/ksongsom/speedseq_result/TERRA/vcf/merge/sugargenes_clean.csv", header = T, stringsAsFactors = F, na.strings=c("","NA"), sep = "\t")
sugargene$genes
tab180$genes
sum(grepl(paste(sugargene$genes, collapse = "|"), tab180$genes))

stopCluster(cl)

