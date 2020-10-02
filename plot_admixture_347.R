library(dplyr)
library(tidyr)
library(readxl)
library(rlang)
options(stringsAsFactors = F)



#scp ksongsom@hpc.uncc.edu://projects/cooper_research1/ksongsom/speedseq_result/TERRA/vcf/merge/admixture/* /Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy347/

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

vitfil <- (vitoverfildt[vitoverfildt$maf_col >= mafcut & abs(vitoverfildt$svlen) <= sizecut & vitoverfildt$missing_col <= missingcut,])
dim(vitfil)
table(vitfil$svtype)
sumvit <- colSums(vitfil[, (which(colnames(vitfil) == "format")+1):(which(colnames(vitfil) == "rio"))], na.rm = T)
background <-  merge(background, sumvitdt, by.y = "PI", by.x = "PI", all.x = T, all.y = F)
unique(background$Race)
head(background)
background$vit01 <- background$sumvit
background$vit01 <- ifelse(background$vit01 != 0, "alt", "ref")
background$vit01[is.na(background$vit01)] <- "Unknown"
background$vit01 
background$sumvit[is.na(background$vit01)] <- "Unknown"

tbl <- read.table("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy347/347.3.Q")
head(tbl)
barplot(t(as.matrix(tbl)), col=rainbow(3),
        xlab="Individual #", ylab="Ancestry", border=NA)

#cross validation error for determining k
  eachk <- list.files(path="/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy347", pattern=".out", full.names=T, recursive=FALSE)
  #eachk <- eachk[-length(eachk)]
  kname <- gsub("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy347/log", "", gsub(".out", "", eachk))
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

tbl <- read.table("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy347/347.3.Q")
head(tbl)
barplot(t(as.matrix(tbl)), col=rainbow(3),
        xlab="Individual #", ylab="Ancestry", border=NA)


# install dependencies and devtools
#install.packages(c("Cairo","ggplot2","gridExtra","gtable","tidyr","devtools"),dependencies=T)

# install pophelper package from GitHub
#devtools::install_github('royfrancis/pophelper')

# load library for use
library(pophelper)
sfiles <- list.files(path="/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy347", pattern=".Q", full.names=T, recursive=FALSE)
slist <- readQ(sfiles)
names(slist) <- gsub(".Q", "", gsub("347.", "", names(slist)))
slist <- slist[str_sort(names(slist), numeric = T)]
str(slist)

setwd("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy347/")
plotQ(slist,imgoutput="join", sortind = "all", sharedindlab = F,
      outputfilename="347_untouch")

fam <- read.table("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy347/347.fam", header = F, sep = " ", col.names = c("fam", "ID", NA, NA, NA, NA))
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
cluster1 <- read.csv("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/cluster_list.csv")
cluster2 <- read.csv("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/cluster_list_update.csv") %>% set_names(c("PI", "k5_cluster_all"))
cluster <- cbind.data.frame(cluster1, cluster2)
colnames(cluster)[ncol(cluster)] <- "k5_cluster_all"
cluster$k5_cluster_all <- as.character(cluster$k5_cluster_all)
head(cluster)

#try k = 5
setwd("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/lumpy347/")
k5 <- slist[[5]]
str(k5)
ks5 <- list(k5)
names(ks5) <- "k5"
plotQ(ks5,  grplab=cluster[, "k5_cluster_all", drop = F], ordergrp=T, showindlab=F, sortind = "all",
      outputfilename="347_K5_group", grplabsize = 1, showlegend=T,
      clustercol=c("green", "blue", "purple", "gold", "red"))

map[map$Race == "NA",]$Race <- "unknown"
map[grepl("-",map$Race),]$Race <- "Mixed"
table(map$Race)
nmap <- map[map$ID %in% cluster$PI,]
dim(nmap)

plotQ(ks5,  grplab=nmap[, "Race", drop = F], ordergrp=T, showindlab=F, sortind = "all",
      outputfilename="347_K5_group_race", grplabsize = 1, showlegend=T)

#try k = 4
k4 <- slist$347.4.Q
str(k4)
ks4 <- list(k4)
names(ks4) <- "k4"
plotQ(ks4,  grplab=background[, "race5", drop = F], ordergrp=T, showindlab=F, sortind = "all",
      outputfilename="K4_group", grplabsize = 0.8, showlegend=T)
plotQ(ks4, ordergrp=T, showindlab=T, sortind = "all", indlabsize = 0.6,
      outputfilename="K4", grplabsize = 0.8, dpi=1000, imgtype="tiff", indlabangle = 90)

#remove ungrooup
rownames(k5) <- fam[,2]
back <- background[background$race5 != "ungroup",]
newk5 <- k5[rownames(k5) %in% back$PI,]
head(newk5)
dim(newk5)
head(back)
newks5 <- list(newk5)
names(newks5) <- "newk5"
plotQ(newks5,  grplab=back[, "race5", drop = F], ordergrp=T, showindlab=F, sortind = "all",
      outputfilename="K5_group_sub", grplabsize = 0.8, showlegend=T)

plotQ(newks5, ordergrp=F, showindlab=F, sortind = "all",
      outputfilename="K5_group_no_group", grplabsize = 0.8, showlegend=T)

tr1 <- tabulateQ(qlist=slist)

#all
plotQ(slist[3:5],imgoutput="join", sortind = "all", sharedindlab = F,
      outputfilename="347_name_all_group",grplab=background[, "race5", drop = F], ordergrp=T, grplabangle = 30, grplabsize = 0.5, dpi = 3000)
plotQ(slist[c(5,17)],imgoutput="join", sortind = "all", sharedindlab = F,
      outputfilename="347_name_all_race",grplab=background[, "Race", drop = F], ordergrp=T, grplabangle = 30, grplabsize = 0.5, dpi = 3000)
plotQ(slist[4:5],imgoutput="join", sortind = "all", sharedindlab = F,
      outputfilename="347_name_all_type",grplab=background[, "Type", drop = F], ordergrp=T, grplabangle = 30, grplabsize = 0.5, dpi = 3000)
plotQ(slist[1:2],imgoutput="join", sortind = "all", sharedindlab = F,
      outputfilename="347_name_all_vit01",grplab=background[, "vit01", drop = F], ordergrp=T, grplabangle = 30, grplabsize = 1.5, dpi = 3000)
plotQ(slist[1:8],imgoutput="join", sortind = "all", sharedindlab = F,
      outputfilename="347_name_all_sumvit",grplab=background[, "sumvit", drop = F], ordergrp=T, grplabangle = 30, grplabsize = 1, dpi = 3000)

head(background)
head(k5)
head(fam)

#k = 3 by origin coordinate

coor <- read.csv("/Users/ksongsom/OneDrive/postdoc/programs/country-capitals.csv", header = T)
head(coor)
sum(background$Origin %in% coor$CountryName)

#seperate in 4 panel
cutlat <- coor[coor$CountryName == "Central African Republic",]$CapitalLatitude
cutlong <- coor[coor$CountryName == "Central African Republic",]$CapitalLongitude

coor$news <- NA
head(coor)

coor$news <- ifelse(coor$CapitalLatitude >= cutlat & coor$CapitalLongitude >= cutlong, "NE",
                    ifelse(coor$CapitalLatitude >= cutlat & coor$CapitalLongitude <= cutlong, "NW",
                           ifelse(coor$CapitalLatitude <= cutlat & coor$CapitalLongitude >= cutlong, "SE", "SW")))

bc <- merge(x = background, y = coor, by.x = "Origin", by.y = "CountryName", all.x = F, all.y = F)
k3 <- slist$347.3.Q
str(k3)
rownames(k3) <- fam[,2]
dim(k3)
k3 <- k3[rownames(k3) %in% bc$PI,]
head(k3)
head(bc)
bc3 <- bc[match(rownames(k3), bc$PI),]
head(bc3)
ks3 <- list(k3)
names(ks3) <- "k3"
plotQ(ks3,  grplab=bc[, "news", drop = F], ordergrp=T, showindlab=F, sortind = "all",
      outputfilename="k3_group_news", grplabsize = 0.8, showlegend=T)
plotQ(ks3, ordergrp=T, showindlab=T, sortind = "all", indlabsize = 0.6,
      outputfilename="k3_news", grplabsize = 0.8, dpi=1000, imgtype="tiff", indlabangle = 90)
