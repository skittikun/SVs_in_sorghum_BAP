library(dplyr)
library(tidyr)
library(readxl)
library(rlang)
options(stringsAsFactors = F)



#scp ksongsom@hpc.uncc.edu://projects/cooper_research1/ksongsom/speedseq_result/TERRA/vcf/merge/admixture/* /Users/ksongsom/OneDrive/postdoc/publication/terra/result/admixture_del/

background <- read.csv("/Users/ksongsom/OneDrive/postdoc/NAM/Panel/BAP/PanelAccessions-BAP.csv", header = T, stringsAsFactors = F)
background$PI <- gsub("LEOTI", "leoti", gsub("RIO","rio", gsub("-dark","", gsub("WRAY", "653616", gsub("ATLAS", "586537", background$Taxa)))))
background$PI <- gsub("PI_", "", background$PI)
altref <- cbind.data.frame(c("Sobic.004G301500", "Sobic.004G301600", "Sobic.004G301650"), c(64031808, 64036221, 64038396), c(64032785, 64037462, 64039002)) %>% setNames(c("gene", "pos", "end"))
vitoverfildt <- c()
for(r in 1:nrow(altref)){
  #vitoverfildt <- rbind.data.frame(vitoverfildt, fildt[fildt$chr == altref$chr[r] & fildt$pos <= altref$end[r] & fildt$end >= altref$pos[r], ])
  vitoverfildt <- rbind.data.frame(vitoverfildt, dat[dat$chr == altref$chr[r] & dat$pos <= altref$end[r] & dat$end >= altref$pos[r], ])
  vitoverfildt <- vitoverfildt[!duplicated(vitoverfildt),]
}


#cross validation error for determining k
  eachk <- list.files(path="/Users/ksongsom/OneDrive/postdoc/publication/terra/result/admixture_del", pattern=".out", full.names=T, recursive=FALSE)
  eachk
  eachk <- eachk[c(1,12, 14:20,2:11, 13)]
  #eachk <- eachk[-length(eachk)]
  kname <- gsub("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/admixture_del/log", "", gsub(".out", "", eachk))
  kname
#plot curve k
  k <- sort(as.numeric(kname))
  k
  cbtab <- c()
  for(i in k){
  text <- readLines(eachk[i])
  cv <- text[grep("CV", text)]
  cverror <- as.numeric(gsub(".*: ", "", cv))
  cbtab <- append(cbtab, cverror)
  print(i)
  }
  
  cbtab <- cbind.data.frame(k, cbtab)
  cbtab <- cbtab[order(cbtab$k),]
  plot(x = cbtab$k, y = cbtab$cbtab, xlab = "k", ylab = "cross-validation error", xaxt='n')
  axis(side = 1, at = cbtab$k,labels = T)
  lines(x = cbtab$k, y = cbtab$cbtab)
  text(x = cbtab$k, y = cbtab$cbtab+(0.01*max(cbtab$cbtab)), labels = round(cbtab$cbtab, 2))

#admixture plot

tbl <- read.table("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/admixture_del/onlydel.8.Q")
head(tbl)
barplot(t(as.matrix(tbl)), col=rainbow(3),
        xlab="Individual #", ylab="Ancestry", border=NA)


# install dependencies and devtools
#install.packages(c("Cairo","ggplot2","gridExtra","gtable","tidyr","devtools"),dependencies=T)

# install pophelper package from GitHub
#devtools::install_github('royfrancis/pophelper')

# load library for use
library(pophelper)
sfiles <- list.files(path="/Users/ksongsom/OneDrive/postdoc/publication/terra/result/admixture_del", pattern=".Q", full.names=T, recursive=FALSE)
slist <- readQ(sfiles)
names(slist) <- gsub(".Q", "", gsub("347.", "", names(slist)))
slist <- slist[str_sort(names(slist), numeric = T)]
str(slist)

setwd("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/admixture_del/")
plotQ(slist,imgoutput="join", sortind = "all", sharedindlab = F,
      outputfilename="347_del",exportpath=getwd())

fam <- read.table("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/admixture_del/onlydel.fam", header = F, sep = " ", col.names = c("fam", "ID", NA, NA, NA, NA))
fam
str(fam)

#load map
map <- (read_excel("/Users/ksongsom/OneDrive/postdoc/NAM/Panel/BAP/update_sorghum_location.xlsx", sheet = 1)) #from google sheet
#update elevation
map$ID <- gsub("PI ", "", map$Genotype)
head(map)
map <- map %>% separate(Coordinate, c("lat", "lon"), sep = ", ")
map <- as.data.frame(map)
#load cluster
#update cluster from SNP
updated1 <- read.xlsx("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/updated_sup_tab1.xlsx")
updated1$ID <- gsub("PI ", "", updated1$Genotype)
all(updated1$ID %in% colnames(pnin))
cluster <- updated1
colnames(cluster)[ncol(cluster)] <- "k8_cluster_all"
cluster$k8_cluster_all <- as.character(cluster$Cluster)
head(cluster)
cluster$ID <- gsub("PI ", "", cluster$Genotype)
fam
cluster <- cluster[match(fam$ID, cluster$ID),]

#try k = 8
setwd("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/admixture_del/")
k8 <- slist$onlydel.8 
str(k8)
attr(x = k8, which = "ind")
attr(x = k8, which = "k")
#colnames(k8) <- paste("Cluster", c(8, 3, 2, 5, 6, 1, 4, 7), sep = "")

k8 <- k8[,order(paste("Cluster", c(8, 3, 2, 5, 6, 1, 4, 7), sep = ""))]
ks8 <- list(k8)
names(ks8) <- "k8"
dim(ks8$k8)
str(ks8)
cluster$k8_cluster_all
dim(cluster)
cluster$Genotype

cluster$ID <- gsub("PI ", "", cluster$Genotype)
sum(cluster$ID %in% fam$fam)
reorder_set <- match(fam$fam, cluster$ID)

cluster$k8_cluster_all
cluster$Cluster

#this is the publication plot
#ks8$k8 <- ks8$k8[reorder_set,]
plotQ(ks8,  grplab=cluster[, "k8_cluster_all", drop = F], ordergrp=T, showindlab=F, sortind = "all",
      outputfilename="347_k8_group", grplabsize = 1, showlegend=T,
      #clustercol=c("red", "deepskyblue1", "dark green", "green", "purple", "orange", "gold", "pink"), 
      clustercol=c(
        "orange",
        "deepskyblue1", 
        "dark green", 
        "gold",
        "green",
        "purple",
        "pink",
        "red" ), 
      exportpath=getwd())

for(CLU in 2:length(slist)){
  k8 <- slist[[CLU]]
  str(k8)
  ks8 <- list(k8)
  names(ks8) <- "k8"
  
  #cluster$ID <- gsub("PI ", "", cluster$Genotype)
  #sum(cluster$ID %in% fam$fam)
  #reorder_set <- match(fam$fam, cluster$ID)
  
  ks8$k8 <- ks8$k8[reorder_set,]
  plotQ(ks8,  grplab=cluster[, "k8_cluster_all", drop = F], ordergrp=T, showindlab=F, sortind = "all",
        outputfilename=paste("347_group", CLU, sep = "_"), grplabsize = 1, showlegend=T,
        clustercol=c("red", "deepskyblue1", "dark green", "green", "purple", "orange", "gold", "pink"), 
        exportpath=getwd())
  print(CLU)
}

#red = group 1 
#deepskyblue1 = group 2 
#dark green = group 3 
#green = group 4 
#purple = group 5
#orange = group 6
#gold = group 7
#pink = group 8


