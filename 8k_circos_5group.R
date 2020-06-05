library(circlize)
library(RColorBrewer)
library(RCircos)
library(plyr)
library(stringi)
library(stringr)
library(dplyr)
library(tidyr)
detach(package:plyr)

########################################circos
#######################################
sorghumgenome <- as.data.frame(cbind(c(paste("Chr0", 1:9, sep = ""), "Chr10"), rep(1, 10),c(80884392, 77742459, 74386277, 68658214, 71854669, 61277060, 65505356,62686529,59416394,61233695)))
colnames(sorghumgenome) <- c("Chromosome", "start", "end")
sorghumgenome$Chromosome <- as.factor(sorghumgenome$Chromosome)
sorghumgenome$start <- as.numeric(as.character(sorghumgenome$start))
sorghumgenome$end <- as.numeric(as.character(sorghumgenome$end))

tab <- readRDS("/Users/ksongsom/OneDrive/postdoc/publication/terra/code/files/347_go.rds")
dim(tab)
head(tab)

table(tab$svtype)
table(tab$svtype)*100/nrow(tab)
tab %>% group_by(chr) %>% summarise(min=min(svlen, na.rm = T), max = max(svlen, na.rm = T), count = n(), percent = n()*100/nrow(tab))

tab2 <- tab[!is.na(tab$genes),]
tab2 %>% group_by(chr) %>% summarise(min=min(svlen, na.rm = T), max = max(svlen, na.rm = T), count = n(), percent = n()*100/nrow(tab2))

head(tab)

table(tab$genes)
sum(is.na(tab$genes))
sum(!is.na(tab$genes))
sum(is.na(tab$genes))*100/nrow(tab)
sum(!is.na(tab$genes))*100/nrow(tab)

sum(is.na(tab$CDS))
sum(!is.na(tab$CDS))
#


tab$ID <- paste(tab$chr, tab$pos, tab$svtype, tab$end, sep = "_")

pnin <- tab[order(tab$pos),]
pnin <- pnin[grepl("Chr", pnin$chr),] 
table(pnin$svtype)
head(pnin)

#background
background <- read.csv("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/dat_back.cvs")
head(background)
cluster1 <- read.csv("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/cluster_list.csv")
cluster2 <- read.csv("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/cluster_list_update.csv") %>% set_names(c("PI", "k5_cluster_all"))
cluster <- cbind.data.frame(cluster1, cluster2)
colnames(cluster)[ncol(cluster)] <- "k5_cluster_all"
head(cluster)
all(cluster$ID == cluster[,6])
table(cluster$k5_cluster_all)

sum(cluster$ID %in% background$Taxa)

backclust <- merge(cluster, background, by.x = "ID", by.y = "Taxa", all.x = T, all.y = F)
head(backclust)
bc_sub <- backclust[,c("PI", "k3_cluster_all", "k7_cluster_all", "k5_cluster_all", "Origin", "Race", "PHOTOPERIOD", "Type")]
bc_sub$PHOTOPERIOD <- tolower(gsub("PHOTOPERIOD_", "", bc_sub$PHOTOPERIOD))
head(bc_sub)
bc_sub$PI[(bc_sub$PI %in% colnames(pnin))]
sum(colnames(pnin) %in% cluster$PI)






#tiff("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/8k_pval.tiff", width=20, height=15, units="in", res=300)
tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/8k_pval_5groups.tiff", width=20, height=15, units="in", res=300)
tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/8k_pval_5groups_sugar.tiff", width=20, height=15, units="in", res=300)
#setup the template
gapdegree = 1
circos.par("track.height" = 0.2, "clock.wise" = TRUE, start.degree = 90, gap.degree=gapdegree)
circos.initialize(factors = sorghumgenome$Chromosome, xlim=as.matrix(sorghumgenome[,2:3]))
circos.trackPlotRegion(ylim=c(0,1), bg.col=rep("white", 10), bg.border=rep("white",10))
str(sorghumgenome)

### Label the sorghum chromosomes
chr.center = as.integer(sorghumgenome$end/2)
y = rep(1, 10) + uy(5, "mm")
for (i in 1:10) {
  circos.text(x=chr.center[i], y=y[i], paste("Chr", i, sep=" "), sector.index=sorghumgenome$Chromosome[i])
}

### Add rectangles for the alignment blocks for all sv (track 1)
my.colors=brewer.pal(9, "Set1")
my.colors = c(my.colors, "hotpink4")

coords <- cbind(sorghumgenome, my.colors)
for (i in 1:nrow(coords)) {
  circos.rect(xleft=coords$start[i], ybottom=0,
              xright=coords$end[i], ytop=1,
              sector.index=coords$Chromosome[i], col="white", border="black",track.index=1)
}

coords <- pnin

ybottom <- (0:dim(coords)[1]%%10)/10
ytop <- ybottom+0.07

colmat <- ifelse(coords$svtype == "DEL", "red", 
                 ifelse(coords$svtype == "INS", "blue", 
                        ifelse(coords$svtype == "INV", "yellow", 
                               ifelse(coords$svtype == "RPL", "grey", 
                                      ifelse(coords$svtype == "DUP", "green", "purple")))))



for(i in 1:nrow(coords)){
  circos.rect(xleft=coords$pos[i], ybottom=ybottom[i], 
              xright=coords$end[i], ytop=ytop[i],
              sector.index=coords$chr[i], track.index=1,
              border=colmat[i], col = colmat[i])
  print(paste(i/nrow(coords)))
}

#add track 2 for number individual all type (track 2)
head(pnin)
stackpn <- data.frame(pnin[,1:(which(colnames(pnin)=="format"))], stack(pnin[,(which(colnames(pnin)=="format")+1):(which(colnames(pnin)=="Rio"))]))
head(stackpn)
length(unique(stackpn$ind))
#stackpn <- data.frame(pn[,1:(which(colnames(pn)=="format"))], stack(pn[,(which(colnames(pn)=="format")+1):ncol(pn)]))

totalindi <- length((which(colnames(pnin)=="format")+1):which(colnames(pnin)=="Rio"))

head(stackpn)
dim(stackpn)

sumstackpn <- stackpn %>% group_by(chr, pos, end, svtype) %>% summarise(numindi = sum(values, na.rm = T))
head(sumstackpn)
table(sumstackpn$numindi)
sumstackpn$percent <- sumstackpn$numindi/totalindi

#add background
head(bc_sub)
#stackpn_back <- merge(stackpn, bc_sub, by.y = "PI", by.x = "ind", all.x = F, all.y = F)
stackpn_back <- merge(stackpn, cluster, by.y = "PI", by.x = "ind", all.x = F, all.y = F)
head(stackpn_back)
dim(stackpn)
dim(stackpn_back)

length(unique(stackpn_back$ind))
length(unique(stackpn_back$PI))

circos.trackPlotRegion(ylim=c(0,1), bg.col=rep("white", 10), bg.border=rep("white",10))
coords <- cbind(sorghumgenome, my.colors)
for (i in 1:nrow(coords)) {
  circos.rect(xleft=coords$start[i], ybottom=0,
              xright=coords$end[i], ytop=1,
              #sector.index=coords$Chromosome[i], col=my.colors[coords$my.colors[i]], border=NA,track.index=2)
              sector.index=coords$Chromosome[i], col="white", border="white",track.index=2)
}


coords <- as.data.frame(sumstackpn)
intv <- 500000
allpeaktab <- c()
for(i in 1:nrow(sorghumgenome)){
  peaktab <- c()
  for(j in seq(sorghumgenome$start[i], sorghumgenome$end[i], by = intv)){
    tab <- coords[grepl(sorghumgenome$Chromosome[i], coords$chr) & coords$pos >= j & coords$pos < (j+intv),]
    if(nrow(tab) == 0){tab[1,] <- 0}
    #peaktab <- rbind.data.frame(peaktab, cbind.data.frame(as.character(sorghumgenome$Chromosome[i]), (j+(j+intv))/2, mean(tab$percent, na.rm = T)))
    peaktab <- rbind.data.frame(peaktab, cbind.data.frame(as.character(sorghumgenome$Chromosome[i]), (j+(j+intv))/2, mean(tab$percent, na.rm = T)))
    print(paste(i, j, sep = "_"))
  }
  allpeaktab <- rbind.data.frame(allpeaktab, peaktab)
}
colnames(allpeaktab) <- c("Chromosome", "mid_x", "peak_y")
str(allpeaktab)
head(allpeaktab)
allpeaktab$peak_y <- allpeaktab$peak_y/max(allpeaktab$peak_y, na.rm = T)
#allpeaktab$peak_y <- replace_na(allpeaktab$peak_y, 0)
coords <- allpeaktab
for(i in 1:nrow(sorghumgenome)){
  circos.lines(x=coords[coords$Chromosome==sorghumgenome$Chromosome[i],]$mid_x, 
               y=coords[coords$Chromosome==sorghumgenome$Chromosome[i],]$peak_y, 
               sector.index=sorghumgenome$Chromosome[i], 
               track.index = 1,
               type = "l", col = "dark blue")
  print(paste(i/nrow(sorghumgenome)))
}

#table of 5 groups
table(bc_sub$k5_cluster_all)

#peak frequency by 5 groups
sumstackpn1 <- stackpn_back[stackpn_back$k5_cluster_all == 1,] %>% group_by(chr, pos, end, svtype) %>% summarise(numindi = sum(values, na.rm = T))
head(sumstackpn1)
head(sumstackpn)
all(sumstackpn1$numindi <= table(cluster$k5_cluster_all)[1])
sumstackpn1$percent <- sumstackpn1$numindi/table(cluster$k5_cluster_all)[1]

sumstackpn2 <- stackpn_back[stackpn_back$k5_cluster_all == 2,] %>% group_by(chr, pos, end, svtype) %>% summarise(numindi = sum(values, na.rm = T))
head(sumstackpn2)
head(sumstackpn)
all(sumstackpn2$numindi <= table(cluster$k5_cluster_all)[2])
sumstackpn2$percent <- sumstackpn2$numindi/table(cluster$k5_cluster_all)[2]

sumstackpn3 <- stackpn_back[stackpn_back$k5_cluster_all == 3,] %>% group_by(chr, pos, end, svtype) %>% summarise(numindi = sum(values, na.rm = T))
head(sumstackpn3)
head(sumstackpn)
all(sumstackpn3$numindi <= table(cluster$k5_cluster_all)[3])
sumstackpn3$percent <- sumstackpn3$numindi/table(cluster$k5_cluster_all)[3]

sumstackpn4 <- stackpn_back[stackpn_back$k5_cluster_all == 4,] %>% group_by(chr, pos, end, svtype) %>% summarise(numindi = sum(values, na.rm = T))
head(sumstackpn4)
head(sumstackpn)
all(sumstackpn4$numindi <= table(cluster$k5_cluster_all)[4])
sumstackpn4$percent <- sumstackpn4$numindi/table(cluster$k5_cluster_all)[4]

sumstackpn5 <- stackpn_back[stackpn_back$k5_cluster_all == 5,] %>% group_by(chr, pos, end, svtype) %>% summarise(numindi = sum(values, na.rm = T))
head(sumstackpn5)
head(sumstackpn)
all(sumstackpn5$numindi <= table(cluster$k5_cluster_all)[5])
sumstackpn5$percent <- sumstackpn5$numindi/table(cluster$k5_cluster_all)[5]

head(sumstackpn1)
head(sumstackpn2)
head(sumstackpn3)
head(sumstackpn4)
head(sumstackpn5)

sumstackpn1$svid <- paste(sumstackpn1$chr, sumstackpn1$pos, sumstackpn1$end, sumstackpn1$svtype, sep = "_")
sumstackpn2$svid <- paste(sumstackpn2$chr, sumstackpn2$pos, sumstackpn2$end, sumstackpn2$svtype, sep = "_")
sumstackpn3$svid <- paste(sumstackpn3$chr, sumstackpn3$pos, sumstackpn3$end, sumstackpn3$svtype, sep = "_")
sumstackpn4$svid <- paste(sumstackpn4$chr, sumstackpn4$pos, sumstackpn4$end, sumstackpn4$svtype, sep = "_")
sumstackpn5$svid <- paste(sumstackpn5$chr, sumstackpn3$pos, sumstackpn5$end, sumstackpn5$svtype, sep = "_")

sumstackpn1$group <- 1
sumstackpn2$group <- 2
sumstackpn3$group <- 3
sumstackpn4$group <- 4
sumstackpn5$group <- 5

dtsum <- cbind.data.frame(sumstackpn1, sumstackpn2, sumstackpn3, sumstackpn4, sumstackpn5)#rbind.data.frame(sumstackpn1, sumstackpn2, sumstackpn3)#
head(dtsum)
tail(dtsum)


#pvalue column
pvalcol <- c()
for(i in 1:nrow(dtsum)){
eachrow <- dtsum[i, grepl("percent", colnames(dtsum))]*100#nrow(cluster)
che <- chisq.test(eachrow, p = rep(1/5,5), correct = T)
pval <- che[[3]]
pvalcol <- c(pvalcol, pval)
print(paste(i, i/nrow(dtsum)))
#need to add which group has highest
}
dtsum <- cbind.data.frame(dtsum, pvalcol)
dtsum$pvalsig <- ifelse(pvalcol <= 0.01/totalindi, "**", ifelse(pvalcol <= 0.05/totalindi, "*", NA))
table(dtsum$pvalsig)
dim(dtsum)
head(dtsum)

#cluster1

coords <- as.data.frame(sumstackpn1)

allpeaktab <- c()
for(i in 1:nrow(sorghumgenome)){
  peaktab <- c()
  for(j in seq(sorghumgenome$start[i], sorghumgenome$end[i], by = intv)){
    tab <- coords[grepl(sorghumgenome$Chromosome[i], coords$chr) & coords$pos >= j & coords$pos < (j+intv),]
    if(nrow(tab) == 0){tab[1,] <- 0}
    #peaktab <- rbind.data.frame(peaktab, cbind.data.frame(as.character(sorghumgenome$Chromosome[i]), (j+(j+intv))/2, mean(tab$percent, na.rm = T)))
    peaktab <- rbind.data.frame(peaktab, cbind.data.frame(as.character(sorghumgenome$Chromosome[i]), (j+(j+intv))/2, mean(tab$percent, na.rm = T)))
    print(paste(i, j, sep = "_"))
  }
  allpeaktab <- rbind.data.frame(allpeaktab, peaktab)
}
colnames(allpeaktab) <- c("Chromosome", "mid_x", "peak_y")
str(allpeaktab)
head(allpeaktab)
allpeaktab$peak_y <- allpeaktab$peak_y/max(allpeaktab$peak_y, na.rm = T)
#allpeaktab$peak_y <- replace_na(allpeaktab$peak_y, 0)
coords <- allpeaktab
for(i in 1:nrow(sorghumgenome)){
  circos.lines(x=coords[coords$Chromosome==sorghumgenome$Chromosome[i],]$mid_x, 
               y=coords[coords$Chromosome==sorghumgenome$Chromosome[i],]$peak_y, 
               sector.index=sorghumgenome$Chromosome[i], 
               track.index = 2,
               type = "l", col = "green")
  print(paste(i/nrow(sorghumgenome)))
}

#cluster2
coords <- as.data.frame(sumstackpn2)

allpeaktab <- c()
for(i in 1:nrow(sorghumgenome)){
  peaktab <- c()
  for(j in seq(sorghumgenome$start[i], sorghumgenome$end[i], by = intv)){
    tab <- coords[grepl(sorghumgenome$Chromosome[i], coords$chr) & coords$pos >= j & coords$pos < (j+intv),]
    if(nrow(tab) == 0){tab[1,] <- 0}
    #peaktab <- rbind.data.frame(peaktab, cbind.data.frame(as.character(sorghumgenome$Chromosome[i]), (j+(j+intv))/2, mean(tab$percent, na.rm = T)))
    peaktab <- rbind.data.frame(peaktab, cbind.data.frame(as.character(sorghumgenome$Chromosome[i]), (j+(j+intv))/2, mean(tab$percent, na.rm = T)))
    print(paste(i, j, sep = "_"))
  }
  allpeaktab <- rbind.data.frame(allpeaktab, peaktab)
}
colnames(allpeaktab) <- c("Chromosome", "mid_x", "peak_y")
str(allpeaktab)
head(allpeaktab)
allpeaktab$peak_y <- allpeaktab$peak_y/max(allpeaktab$peak_y, na.rm = T)
#allpeaktab$peak_y <- replace_na(allpeaktab$peak_y, 0)
coords <- allpeaktab
for(i in 1:nrow(sorghumgenome)){
  circos.lines(x=coords[coords$Chromosome==sorghumgenome$Chromosome[i],]$mid_x, 
               y=coords[coords$Chromosome==sorghumgenome$Chromosome[i],]$peak_y, 
               sector.index=sorghumgenome$Chromosome[i], 
               track.index = 2,
               type = "l", col = "gold")
  print(paste(i/nrow(sorghumgenome)))
}

#cluster3
coords <- as.data.frame(sumstackpn3)

allpeaktab <- c()
for(i in 1:nrow(sorghumgenome)){
  peaktab <- c()
  for(j in seq(sorghumgenome$start[i], sorghumgenome$end[i], by = intv)){
    tab <- coords[grepl(sorghumgenome$Chromosome[i], coords$chr) & coords$pos >= j & coords$pos < (j+intv),]
    if(nrow(tab) == 0){tab[1,] <- 0}
    #peaktab <- rbind.data.frame(peaktab, cbind.data.frame(as.character(sorghumgenome$Chromosome[i]), (j+(j+intv))/2, mean(tab$percent, na.rm = T)))
    peaktab <- rbind.data.frame(peaktab, cbind.data.frame(as.character(sorghumgenome$Chromosome[i]), (j+(j+intv))/2, mean(tab$percent, na.rm = T)))
    print(paste(i, j, sep = "_"))
  }
  allpeaktab <- rbind.data.frame(allpeaktab, peaktab)
}
colnames(allpeaktab) <- c("Chromosome", "mid_x", "peak_y")
str(allpeaktab)
head(allpeaktab)
allpeaktab$peak_y <- allpeaktab$peak_y/max(allpeaktab$peak_y, na.rm = T)
#allpeaktab$peak_y <- replace_na(allpeaktab$peak_y, 0)
coords <- allpeaktab
for(i in 1:nrow(sorghumgenome)){
  circos.lines(x=coords[coords$Chromosome==sorghumgenome$Chromosome[i],]$mid_x, 
               y=coords[coords$Chromosome==sorghumgenome$Chromosome[i],]$peak_y, 
               sector.index=sorghumgenome$Chromosome[i], 
               track.index = 2,
               type = "l", col = "red")
  print(paste(i/nrow(sorghumgenome)))
}

#cluster4
coords <- as.data.frame(sumstackpn4)

allpeaktab <- c()
for(i in 1:nrow(sorghumgenome)){
  peaktab <- c()
  for(j in seq(sorghumgenome$start[i], sorghumgenome$end[i], by = intv)){
    tab <- coords[grepl(sorghumgenome$Chromosome[i], coords$chr) & coords$pos >= j & coords$pos < (j+intv),]
    if(nrow(tab) == 0){tab[1,] <- 0}
    #peaktab <- rbind.data.frame(peaktab, cbind.data.frame(as.character(sorghumgenome$Chromosome[i]), (j+(j+intv))/2, mean(tab$percent, na.rm = T)))
    peaktab <- rbind.data.frame(peaktab, cbind.data.frame(as.character(sorghumgenome$Chromosome[i]), (j+(j+intv))/2, mean(tab$percent, na.rm = T)))
    print(paste(i, j, sep = "_"))
  }
  allpeaktab <- rbind.data.frame(allpeaktab, peaktab)
}
colnames(allpeaktab) <- c("Chromosome", "mid_x", "peak_y")
str(allpeaktab)
head(allpeaktab)
allpeaktab$peak_y <- allpeaktab$peak_y/max(allpeaktab$peak_y, na.rm = T)
#allpeaktab$peak_y <- replace_na(allpeaktab$peak_y, 0)
coords <- allpeaktab
for(i in 1:nrow(sorghumgenome)){
  circos.lines(x=coords[coords$Chromosome==sorghumgenome$Chromosome[i],]$mid_x, 
               y=coords[coords$Chromosome==sorghumgenome$Chromosome[i],]$peak_y, 
               sector.index=sorghumgenome$Chromosome[i], 
               track.index = 2,
               type = "l", col = "blue")
  print(paste(i/nrow(sorghumgenome)))
}

#cluster5
coords <- as.data.frame(sumstackpn5)

allpeaktab <- c()
for(i in 1:nrow(sorghumgenome)){
  peaktab <- c()
  for(j in seq(sorghumgenome$start[i], sorghumgenome$end[i], by = intv)){
    tab <- coords[grepl(sorghumgenome$Chromosome[i], coords$chr) & coords$pos >= j & coords$pos < (j+intv),]
    if(nrow(tab) == 0){tab[1,] <- 0}
    #peaktab <- rbind.data.frame(peaktab, cbind.data.frame(as.character(sorghumgenome$Chromosome[i]), (j+(j+intv))/2, mean(tab$percent, na.rm = T)))
    peaktab <- rbind.data.frame(peaktab, cbind.data.frame(as.character(sorghumgenome$Chromosome[i]), (j+(j+intv))/2, mean(tab$percent, na.rm = T)))
    print(paste(i, j, sep = "_"))
  }
  allpeaktab <- rbind.data.frame(allpeaktab, peaktab)
}
colnames(allpeaktab) <- c("Chromosome", "mid_x", "peak_y")
str(allpeaktab)
head(allpeaktab)
allpeaktab$peak_y <- allpeaktab$peak_y/max(allpeaktab$peak_y, na.rm = T)
#allpeaktab$peak_y <- replace_na(allpeaktab$peak_y, 0)
coords <- allpeaktab
for(i in 1:nrow(sorghumgenome)){
  circos.lines(x=coords[coords$Chromosome==sorghumgenome$Chromosome[i],]$mid_x, 
               y=coords[coords$Chromosome==sorghumgenome$Chromosome[i],]$peak_y, 
               sector.index=sorghumgenome$Chromosome[i], 
               track.index = 2,
               type = "l", col = "purple")
  print(paste(i/nrow(sorghumgenome)))
}

#present second ring in combination
head(sumstackpn2)
head(sumstackpn3)
head(sumstackpn)
dim(sumstackpn1)
dim(sumstackpn2)
dim(sumstackpn3)
dim(sumstackpn4)
dim(sumstackpn5)

coords <- rbind.fill(sumstackpn1, sumstackpn2, sumstackpn3, sumstackpn4, sumstackpn5)
table(coords$group)
dim(coords)

allrect <- c()
for(i in 1:nrow(sorghumgenome)){
  recttab <- c()
  for(j in 1:length(seq(sorghumgenome$start[i], sorghumgenome$end[i], by = intv))){
    xleftj <- seq(sorghumgenome$start[i], sorghumgenome$end[i], by = intv)
    xrightj <- c(seq(sorghumgenome$start[i], sorghumgenome$end[i], by = intv)[-1]-1,sorghumgenome$end[i])
    tab <- coords[grepl(sorghumgenome$Chromosome[i], coords$chr) & coords$pos >= xleftj[j] & coords$pos < xrightj[j],]
    if(nrow(tab) == 0){tab[1,] <- 0}
    rset <- length(unique(tab$group))
    namegroup <- unique(tab$group)
    sdf <- cbind.data.frame(as.character(sorghumgenome$Chromosome[i]), xleftj[j], xrightj[j])
    headdf <- cbind.data.frame(rbind.data.frame(sdf, sdf[rep(1, rset-1), ]), namegroup)
    #ytopj <- as.data.frame(tab %>% group_by(group) %>% summarise(sumeach = (sum(numindi, na.rm = T)/sum(tab$numindi, na.rm = T))) %>% mutate(cumsum = cumsum(sumeach)))
    ytopj <- as.data.frame(tab %>% group_by(group) %>% summarise(sumeach = (sum(percent, na.rm = T)/sum(tab$percent, na.rm = T))) %>% mutate(cumsum = cumsum(sumeach)))
    ybottomj <- c(0, (ytopj$cumsum[-nrow(ytopj)]-0.001))
    #recttab <- rbind.data.frame(recttab, cbind.data.frame(as.character(sorghumgenome$Chromosome[i]), (j+(j+intv))/2, mean(tab$percent, na.rm = T)))
    recttab <- rbind.data.frame(recttab, cbind.data.frame(headdf, ybottomj, ytopj$cumsum))
    print(paste(i, j, sep = "_"))
  }
  allrect <- rbind.data.frame(allrect, recttab)
}
allrect <- allrect %>% set_names(c("chr", "pos", "end", "namegroup", "ybottom", "ytop"))
coords <- allrect
colmat <- ifelse(coords$namegroup == 1, "yellow", 
                 ifelse(coords$namegroup == 2, "purple", 
                        ifelse(coords$namegroup == 3, "green", 
                               ifelse(coords$namegroup == 4, "red", 
                                      ifelse(coords$namegroup == 5, "blue", NA)))))
table(colmat)

for(i in 1:nrow(coords)){
  circos.rect(xleft=coords$pos[i], ybottom=coords$ybottom[i], 
              xright=coords$end[i], ytop=coords$ytop[i],
              sector.index=coords$chr[i], track.index=2,
              border=colmat[i], col = colmat[i])
  print(paste(i/nrow(coords)))
}
#blue = all
#yellow = group 1 yellow
#green = group 2 purple
#red = group 3 green
#blue = group 4 red
#purple = gropu 5 blue


###add pval stars track 3
head(dtsum)
dtsum$logp <- -log(dtsum$pvalcol)
dim(dtsum)


circos.trackPlotRegion(ylim=c(0,1), bg.col=rep("white", 10), bg.border=rep("white",10))
coords <- cbind(sorghumgenome, my.colors)
for (i in 1:nrow(coords)) {
  circos.rect(xleft=coords$start[i], ybottom=0,
              xright=coords$end[i], ytop=1,
              #sector.index=coords$Chromosome[i], col=my.colors[coords$my.colors[i]], border=NA,track.index=2)
              sector.index=coords$Chromosome[i], col="white", border="white",track.index=3)
}

coords <- as.data.frame(dtsum)
allpeaktab <- c()
for(i in 1:nrow(sorghumgenome)){
  peaktab <- c()
  for(j in seq(sorghumgenome$start[i], sorghumgenome$end[i], by = intv)){
    tab <- coords[grepl(sorghumgenome$Chromosome[i], coords$chr) & coords$pos >= j & coords$pos < (j+intv),]
    if(nrow(tab) == 0){tab[1,] <- 0}
    #peaktab <- rbind.data.frame(peaktab, cbind.data.frame(as.character(sorghumgenome$Chromosome[i]), (j+(j+intv))/2, mean(tab$percent, na.rm = T)))
    peaktab <- rbind.data.frame(peaktab, cbind.data.frame(as.character(sorghumgenome$Chromosome[i]), (j+(j+intv))/2, mean(tab$logp, na.rm = T)))
    print(paste(i, j, sep = "_"))
  }
  allpeaktab <- rbind.data.frame(allpeaktab, peaktab)
}
colnames(allpeaktab) <- c("Chromosome", "mid_x", "peak_y")
str(allpeaktab)
head(allpeaktab)
maxy = max(allpeaktab$peak_y, na.rm = T)
allpeaktab$peak_y <- allpeaktab$peak_y/maxy
#allpeaktab$peak_y <- replace_na(allpeaktab$peak_y, 0)
pval_cutoff <- (-log(0.05/(nrow(allpeaktab)*intv)))/maxy

coords <- allpeaktab
for(i in 1:nrow(sorghumgenome)){
  circos.lines(x=coords[coords$Chromosome==sorghumgenome$Chromosome[i],]$mid_x, 
               y=coords[coords$Chromosome==sorghumgenome$Chromosome[i],]$peak_y, 
               sector.index=sorghumgenome$Chromosome[i], 
               track.index = 3,
               type = "l", col = "blue")
  print(paste(i/nrow(sorghumgenome)))
}

for(i in 1:nrow(sorghumgenome)){
  circos.lines(y = rep(pval_cutoff, 2), x = as.numeric(sorghumgenome[i,2:3]), track.index = 3, 
               type = "s", sector.index = sorghumgenome$Chromosome[i], col = "red")
  print(i)
  }


text(0, 0, "Structure variations \n among \n 347 sorghums", cex = 2)
text(0, 0.90, "Overall SV", cex = 1, srt = 90)
text(0, 0.67, "Abundance", cex = 1, srt = 90)
text(0, 0.40, "Chi-square", cex = 1, srt = 90)

dev.off()































#add track X for number individual sweet 
dim(background)
head(stackpn)
length(unique(stackpn$ind))
stackpn_back <- merge(stackpn, background, by.y = "PI", by.x = "ind", all.x = T, all.y = F)
sumstackpn_back_sweet <- stackpn_back[stackpn_back$Type == "Sweet", ] %>% group_by(chr, pos, end, svtype) %>% summarise(numindi = sum(values, na.rm = T))
head(sumstackpn_back_sweet)
table(stackpn_back$Type)

totalindi_sweet <- length(na.omit(unique(stackpn_back[stackpn_back$Type == "Sweet", ]$ind)))

sumstackpn_back_sweet$percent <- sumstackpn_back_sweet$numindi/totalindi_sweet

circos.trackPlotRegion(ylim=c(0,1), bg.col=rep("white", 10), bg.border=rep("white",10))
coords <- cbind(sorghumgenome, my.colors)
for (i in 1:nrow(coords)) {
  circos.rect(xleft=coords$start[i], ybottom=0,
              xright=coords$end[i], ytop=1,
              #sector.index=coords$Chromosome[i], col=my.colors[coords$my.colors[i]], border=NA,track.index=2)
              sector.index=coords$Chromosome[i], col="white", border="white",track.index=3)
}

coords <- sumstackpn_back_sweet
#coords <- dat5in[missing <= 0.2,]
intv <- 100000
allpeaktab <- c()
for(i in 1:nrow(sorghumgenome)){
  peaktab <- c()
  for(j in seq(sorghumgenome$start[i], sorghumgenome$end[i], by = intv)){
    tab <- coords[coords$chr == sorghumgenome$Chromosome[i] & coords$pos >= j & coords$pos < (j+intv),]
    peaktab <- rbind.data.frame(peaktab, cbind.data.frame(as.character(sorghumgenome$Chromosome[i]), (j+(j+intv))/2, mean(tab$percent, na.rm = T)))
    #peaktab <- rbind.data.frame(peaktab, cbind.data.frame(as.character(sorghumgenome$Chromosome[i]), (j+(j+intv))/2, mean(tab$sumpresent_sweet, na.rm = T)))
    print(paste(i, j, sep = "_"))
  }
  allpeaktab <- rbind.data.frame(allpeaktab, peaktab)
}
colnames(allpeaktab) <- c("Chromosome", "mid_x", "peak_y")
str(allpeaktab)
#allpeaktab$peak_y <- allpeaktab$peak_y/totalindi
allpeaktab$peak_y <- allpeaktab$peak_y/max(allpeaktab$peak_y, na.rm= T)

coords <- allpeaktab
for(i in 1:nrow(sorghumgenome)){
  circos.lines(x=coords[coords$Chromosome==sorghumgenome$Chromosome[i],]$mid_x, 
               y=coords[coords$Chromosome==sorghumgenome$Chromosome[i],]$peak_y, 
               sector.index=sorghumgenome$Chromosome[i], 
               track.index = 3,
               type = "l", col = "gold")
  print(paste(i/nrow(sorghumgenome)))
}

#add track 4 for number individual sweet and non-sweet
#non-sweet
circos.trackPlotRegion(ylim=c(0,1), bg.col=rep("white", 10), bg.border=rep("white",10))
coords <- cbind(sorghumgenome, my.colors)
for (i in 1:nrow(coords)) {
  circos.rect(xleft=coords$start[i], ybottom=0,
              xright=coords$end[i], ytop=1,
              #sector.index=coords$Chromosome[i], col=my.colors[coords$my.colors[i]], border=NA,track.index=2)
              sector.index=coords$Chromosome[i], col="white", border="white",track.index=4)
}

#sumstackpn_back_nonsweet <- stackpn_back[stackpn_back$Type != "Sweet", ] %>% group_by(chr, pos, end, svtype) %>% summarise(numindi = sum(values, na.rm = T))
head(sumstackpn_back_nonsweet)
table(stackpn_back$Type)

totalindi_nonsweet <- totalindi- totalindi_sweet

sumstackpn_back_nonsweet$percent <- sumstackpn_back_nonsweet$numindi/totalindi_nonsweet

coords <- sumstackpn_back_nonsweet
#coords <- dat5in[missing <= 0.2,]
intv <- 100000
allpeaktab <- c()
for(i in 1:nrow(sorghumgenome)){
  peaktab <- c()
  for(j in seq(sorghumgenome$start[i], sorghumgenome$end[i], by = intv)){
    tab <- coords[coords$chr == sorghumgenome$Chromosome[i] & coords$pos >= j & coords$pos < (j+intv),]
    peaktab <- rbind.data.frame(peaktab, cbind.data.frame(as.character(sorghumgenome$Chromosome[i]), (j+(j+intv))/2, mean(tab$percent, na.rm = T)))
    #peaktab <- rbind.data.frame(peaktab, cbind.data.frame(as.character(sorghumgenome$Chromosome[i]), (j+(j+intv))/2, mean(tab$sumpresent_nonsweet, na.rm = T)))
    print(paste(i, j, sep = "_"))
  }
  allpeaktab <- rbind.data.frame(allpeaktab, peaktab)
}
colnames(allpeaktab) <- c("Chromosome", "mid_x", "peak_y")
str(allpeaktab)
allpeaktab$peak_y <- allpeaktab$peak_y/max(allpeaktab$peak_y, na.rm= T)
#allpeaktab$peak_y <- allpeaktab$peak_y/totalindi
coords <- allpeaktab
for(i in 1:nrow(sorghumgenome)){
  circos.lines(x=coords[coords$Chromosome==sorghumgenome$Chromosome[i],]$mid_x, 
               y=coords[coords$Chromosome==sorghumgenome$Chromosome[i],]$peak_y, 
               sector.index=sorghumgenome$Chromosome[i], 
               track.index = 4,
               type = "l", col = "green4")
  print(paste(i/nrow(sorghumgenome)))
}



