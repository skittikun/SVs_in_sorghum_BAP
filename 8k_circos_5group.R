library("scales")
library(openxlsx)
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
tab <- tab[tab$svlen>= 50,]
dim(tab)
head(tab)
tab[grepl("GO:0051015",tab$go),]

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

tab$ID <- paste(tab$chr, tab$pos, tab$svtype, tab$end, sep = "_")

pnin <- tab[order(tab$pos),]
pnin <- pnin[grepl("Chr", pnin$chr),] 
table(pnin$svtype)
head(pnin)

#background
background <- read.csv("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/dat_back.cvs")
head(background)
cluster1 <- read.csv("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/cluster_list.csv")
cluster2 <- read.csv("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/cluster_list_update.csv") 
colnames(cluster2) <- (c("PI", "k5_cluster_all"))
cluster <- cbind.data.frame(cluster1, cluster2)
colnames(cluster)[ncol(cluster)] <- "k5_cluster_all"
head(cluster)
head(updated1$ID)
all(cluster$ID == cluster[,6])
table(cluster$Cluster)


#update cluster from SNP
updated1 <- read.xlsx("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/updated_sup_tab1.xlsx")
updated1$ID <- gsub("PI ", "", updated1$Genotype)
sum(updated1$ID %in% background$Taxa)
all(updated1$ID %in% colnames(pnin))
cluster <- updated1

#load snp and sv diversity
windowdat <- readRDS("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/pi_window_snp_820.rds")
windowdat_sv <- readRDS("/Users/ksongsom/OneDrive/postdoc/publication/terra/result//pi_window_sv_deletions.rds")

#remove row with now snp
windowsize <- 500000
windowdat[windowdat$size == 1,]$sum_pi <- NA
windowdat[windowdat$size == 1,]$size <- NA

windowdat_sv[windowdat_sv$size == 1,]$sum_pi <- NA
windowdat_sv[windowdat_sv$size == 1,]$size <- NA

windowdat_sv$pi_per_window <- (windowdat_sv$sum_pi/windowsize)
windowdat$pi_per_window <- (windowdat$sum_pi/windowsize)

windowdat_sv$pi_per_window[is.na(windowdat_sv$pi_per_window)] <- 0
windowdat$pi_per_window[is.na(windowdat$pi_per_window)] <- 0

#plot circos
tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/8k_pval_8groups_fixedpval_sv_noname_with_pi_big_green.tiff", width=20, height=15, units="in", res=300)

#setup the template
gapdegree = 1
circos.par("track.height" = 0.22, "clock.wise" = TRUE, start.degree = 90, gap.degree=gapdegree)
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
                        ifelse(coords$svtype == "INV", "gold", 
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
#only deletions
pnin <- pnin[grepl("DEL", pnin$svid),]
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
stackpn_back <- merge(stackpn, cluster, by.y = "ID", by.x = "ind", all.x = F, all.y = F)
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

#SV abundance
#coords <- allpeaktab
#for(i in 1:nrow(sorghumgenome)){
#  circos.lines(x=coords[coords$Chromosome==sorghumgenome$Chromosome[i],]$mid_x, 
#               y=coords[coords$Chromosome==sorghumgenome$Chromosome[i],]$peak_y, 
#               sector.index=sorghumgenome$Chromosome[i], 
#               track.index = 1,
#               type = "l", col = "dark blue", cex = 2)
#  print(paste(i/nrow(sorghumgenome)))
#}

#add pi of both snp and sv
windowdat$mid_x <- (windowdat$upper + windowdat$lower)/2
windowdat$peak_y <- windowdat$pi_per_window/max(windowdat$pi_per_window, na.rm = T)
coords <- windowdat
for(i in 1:nrow(sorghumgenome)){
  circos.lines(x=coords[coords$windowchr==sorghumgenome$Chromosome[i],]$mid_x, 
               y=coords[coords$windowchr==sorghumgenome$Chromosome[i],]$peak_y, 
               sector.index=sorghumgenome$Chromosome[i], 
               track.index = 1,
               type = "l", col = "green", cex = 2, lwd = 2)
  print(paste(i/nrow(sorghumgenome)))
}

windowdat_sv$mid_x <- (windowdat_sv$upper + windowdat_sv$lower)/2
windowdat_sv$peak_y <- windowdat_sv$pi_per_window/max(windowdat_sv$pi_per_window, na.rm = T)
coords <- windowdat_sv
for(i in 1:nrow(sorghumgenome)){
  circos.lines(x=coords[coords$windowchr==sorghumgenome$Chromosome[i],]$mid_x, 
               y=coords[coords$windowchr==sorghumgenome$Chromosome[i],]$peak_y, 
               sector.index=sorghumgenome$Chromosome[i], 
               track.index = 1,
               type = "l", col = "dark blue", cex = 2, lwd = 2)
  print(paste(i/nrow(sorghumgenome)))
}


#table of 8 groups
table(bc_sub  )

#peak frequency by 8 groups
sumstackpn1 <- stackpn_back[stackpn_back$Cluster == 1,] %>% group_by(chr, pos, end, svtype) %>% summarise(numindi = sum(values, na.rm = T))
head(sumstackpn1)
head(sumstackpn)
all(sumstackpn1$numindi <= table(cluster$Cluster)[1])
sumstackpn1$percent <- sumstackpn1$numindi/table(cluster$Cluster)[1]

sumstackpn2 <- stackpn_back[stackpn_back$Cluster == 2,] %>% group_by(chr, pos, end, svtype) %>% summarise(numindi = sum(values, na.rm = T))
head(sumstackpn2)
head(sumstackpn)
all(sumstackpn2$numindi <= table(cluster$Cluster)[2])
sumstackpn2$percent <- sumstackpn2$numindi/table(cluster$Cluster)[2]

sumstackpn3 <- stackpn_back[stackpn_back$Cluster == 3,] %>% group_by(chr, pos, end, svtype) %>% summarise(numindi = sum(values, na.rm = T))
head(sumstackpn3)
head(sumstackpn)
all(sumstackpn3$numindi <= table(cluster$Cluster)[3])
sumstackpn3$percent <- sumstackpn3$numindi/table(cluster$Cluster)[3]

sumstackpn4 <- stackpn_back[stackpn_back$Cluster == 4,] %>% group_by(chr, pos, end, svtype) %>% summarise(numindi = sum(values, na.rm = T))
head(sumstackpn4)
head(sumstackpn)
all(sumstackpn4$numindi <= table(cluster$Cluster)[4])
sumstackpn4$percent <- sumstackpn4$numindi/table(cluster$Cluster)[4]

sumstackpn5 <- stackpn_back[stackpn_back$Cluster == 5,] %>% group_by(chr, pos, end, svtype) %>% summarise(numindi = sum(values, na.rm = T))
head(sumstackpn5)
head(sumstackpn)
all(sumstackpn5$numindi <= table(cluster$Cluster)[5])
sumstackpn5$percent <- sumstackpn5$numindi/table(cluster$Cluster)[5] 

sumstackpn6 <- stackpn_back[stackpn_back$Cluster == 6,] %>% group_by(chr, pos, end, svtype) %>% summarise(numindi = sum(values, na.rm = T))
head(sumstackpn6)
head(sumstackpn)
all(sumstackpn6$numindi <= table(cluster$Cluster)[6])
sumstackpn6$percent <- sumstackpn6$numindi/table(cluster$Cluster)[6] 

sumstackpn7 <- stackpn_back[stackpn_back$Cluster == 7,] %>% group_by(chr, pos, end, svtype) %>% summarise(numindi = sum(values, na.rm = T))
head(sumstackpn7)
head(sumstackpn)
all(sumstackpn7$numindi <= table(cluster$Cluster)[7])
sumstackpn7$percent <- sumstackpn7$numindi/table(cluster$Cluster)[7] 

sumstackpn8 <- stackpn_back[stackpn_back$Cluster == 8,] %>% group_by(chr, pos, end, svtype) %>% summarise(numindi = sum(values, na.rm = T))
head(sumstackpn8)
head(sumstackpn)
all(sumstackpn8$numindi <= table(cluster$Cluster)[8])
sumstackpn8$percent <- sumstackpn8$numindi/table(cluster$Cluster)[8] 

head(sumstackpn1)
head(sumstackpn2)
head(sumstackpn3)
head(sumstackpn4)
head(sumstackpn5)
head(sumstackpn6)
head(sumstackpn7)
head(sumstackpn8)




sumstackpn1$svid <- paste(sumstackpn1$chr, sumstackpn1$pos, sumstackpn1$end, sumstackpn1$svtype, sep = "_")
sumstackpn2$svid <- paste(sumstackpn2$chr, sumstackpn2$pos, sumstackpn2$end, sumstackpn2$svtype, sep = "_")
sumstackpn3$svid <- paste(sumstackpn3$chr, sumstackpn3$pos, sumstackpn3$end, sumstackpn3$svtype, sep = "_")
sumstackpn4$svid <- paste(sumstackpn4$chr, sumstackpn4$pos, sumstackpn4$end, sumstackpn4$svtype, sep = "_")
sumstackpn5$svid <- paste(sumstackpn5$chr, sumstackpn3$pos, sumstackpn5$end, sumstackpn5$svtype, sep = "_")
sumstackpn6$svid <- paste(sumstackpn6$chr, sumstackpn3$pos, sumstackpn6$end, sumstackpn6$svtype, sep = "_")
sumstackpn7$svid <- paste(sumstackpn7$chr, sumstackpn3$pos, sumstackpn7$end, sumstackpn7$svtype, sep = "_")
sumstackpn8$svid <- paste(sumstackpn8$chr, sumstackpn3$pos, sumstackpn8$end, sumstackpn8$svtype, sep = "_")



sumstackpn1$group <- 1
sumstackpn2$group <- 2
sumstackpn3$group <- 3
sumstackpn4$group <- 4
sumstackpn5$group <- 5
sumstackpn6$group <- 6
sumstackpn7$group <- 7
sumstackpn8$group <- 8

dtsum <- cbind.data.frame(sumstackpn1, sumstackpn2, sumstackpn3, sumstackpn4, sumstackpn5,
                          sumstackpn6, sumstackpn7, sumstackpn8)#rbind.data.frame(sumstackpn1, sumstackpn2, sumstackpn3)#
head(dtsum)
tail(dtsum)
dim(dtsum)
table(dtsum$group)

#pvalue column
#load observed sum
observetable <- read.csv("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/expected_percentage_withoutNA.csv")
observetable <- observetable[grepl("DEL", observetable$eachsv),]
dim(observetable[grepl("DEL", observetable$eachsv),])
head(observetable[,c(1,7)])
table(updated1$Cluster)
head(pnin)
dim(observetable)
head(updated1)
all(updated1$ID %in% colnames(pnin))
#calculate expected percentage regarding missing value

#expectedtable <- c()
#for(i in 1:nrow(pnin)){
  #only name without na
#  stillin <- row.names(na.omit(t(pnin[i,colnames(pnin) %in% updated1$ID])))
  #count each cluster
  
#  eachexpected <- (c(round(table(updated1[updated1$ID %in% stillin, c("ID", "Cluster")]$Cluster)/length(stillin), 2), length(stillin)))
#  expectedtable <- rbind.data.frame(expectedtable, eachexpected)
#  print(paste(i/nrow(pnin)))
#}

#expectedtable <- cbind.data.frame(pnin$svid, expectedtable)
#colnames(expectedtable) <- c("eachsv", paste("Cluster_", 1:8, sep = ""), "obssum")

#write.xlsx(expectedtable, "/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/expected_percentage_from_SNP.xlsx")
expectedtable <- read.xlsx("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/expected_percentage_from_SNP.xlsx")
expectedtable <- expectedtable[grepl("DEL", expectedtable$eachsv),]

pvalcol <- c()
for(i in 1:nrow(dtsum)){
  SVID <- dtsum$svid[i]
  expectsv <- expectedtable[grepl(SVID, expectedtable$eachsv),2:9]
  expectsv <- as.data.frame(t(expectsv))[,1]
  eachrow <- dtsum[i, grepl("percent", colnames(dtsum))]*100#nrow(cluster)
  expectsv <- expectsv/sum(expectsv)
  che <- chisq.test(eachrow, p = (expectsv), correct = T)
  pval <- che$p.value
  pvalcol <- c(pvalcol, pval)
  print(paste(i, i/nrow(dtsum)))
#need to add which group has highest
}

dtsum <- cbind.data.frame(dtsum, pvalcol)
dtsum$pvalsig <- ifelse(pvalcol <= 0.01/totalindi, "**", ifelse(pvalcol <= 0.05/totalindi, "*", NA))
table(dtsum$pvalsig)
dim(dtsum)
head(dtsum)
dtsum[grep("Chr01_1011344", dtsum$svid),]


#present second ring in combination
head(sumstackpn2)
head(sumstackpn3)
head(sumstackpn)
dim(sumstackpn1)
dim(sumstackpn2)
dim(sumstackpn3)
dim(sumstackpn4)
dim(sumstackpn5)
dim(sumstackpn6)
dim(sumstackpn7)
dim(sumstackpn8)


coords <- plyr::rbind.fill(sumstackpn1, sumstackpn2, sumstackpn3, sumstackpn4, sumstackpn5,
                           sumstackpn6, sumstackpn7, sumstackpn8)
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
allrect <- allrect %>% purrr::set_names(c("chr", "pos", "end", "namegroup", "ybottom", "ytop"))
coords <- allrect
colmat <- ifelse(coords$namegroup == 1, "red", 
                 ifelse(coords$namegroup == 2, "deepskyblue1", 
                        ifelse(coords$namegroup == 3, "dark green", 
                               ifelse(coords$namegroup == 4, "green", 
                                      ifelse(coords$namegroup == 5, "purple", 
                                             ifelse(coords$namegroup == 6, "orange", 
                                                    ifelse(coords$namegroup == 7, "gold", 
                                                           ifelse(coords$namegroup == 8, "pink", 
                                             NA))))))))
                  
table(colmat)

for(i in 1:nrow(coords)){
  circos.rect(xleft=coords$pos[i], ybottom=coords$ybottom[i], 
              xright=coords$end[i], ytop=coords$ytop[i],
              sector.index=coords$chr[i], track.index=2,
              border=colmat[i], col = colmat[i])
  print(paste(i/nrow(coords)))
}
#blue = all
#red = group 1 
#deepskyblue1 = group 2 
#dark green = group 3 
#green = group 4 
#purple = group 5
#orange = group 6
#gold = group 7
#pink = group 8


###add pval track 3
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

#Bonferoni's cutoff
pval_cutoff <- (-log(0.05/(nrow(allpeaktab)*intv)))/maxy

coords <- allpeaktab
for(i in 1:nrow(sorghumgenome)){
  circos.lines(x=coords[coords$Chromosome==sorghumgenome$Chromosome[i],]$mid_x, 
               y=coords[coords$Chromosome==sorghumgenome$Chromosome[i],]$peak_y, 
               sector.index=sorghumgenome$Chromosome[i], 
               track.index = 3,
               type = "l", col = "black")
  print(paste(i/nrow(sorghumgenome)))
}

for(i in 1:nrow(sorghumgenome)){
  circos.lines(y = rep(pval_cutoff, 2), x = as.numeric(sorghumgenome[i,2:3]), track.index = 3, 
               type = "s", sector.index = sorghumgenome$Chromosome[i], col = "red")
  print(i)
  }




#add cluster specific
#add track 4
circos.trackPlotRegion(ylim=c(0,1), bg.col=rep("white", 10), bg.border=rep("white",10))
coords <- cbind(sorghumgenome, my.colors)
for (i in 1:nrow(coords)) {
  circos.rect(xleft=coords$start[i], ybottom=0,
              xright=coords$end[i], ytop=1,
              #sector.index=coords$Chromosome[i], col=my.colors[coords$my.colors[i]], border=NA,track.index=2)
              sector.index=coords$Chromosome[i], col="white", border="white",track.index=4)
}

functab <- read.xlsx("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/svlist_sigGO_1570_del.xlsx")
dim(functab)
for(i in 1:nrow(functab)){
  functab$Chr[i] <- as.character(str_split(functab$eachsv[i], "_")[[1]][1])
  functab$start[i] <- as.numeric(str_split(functab$eachsv[i], "_")[[1]][2])
  functab$end[i] <-   as.numeric(str_split(functab$eachsv[i], "_")[[1]][3])
  functab$svtype[i] <- as.character(str_split(functab$eachsv[i], "_")[[1]][4])
  functab$mid[i] <- (functab$start[i] + functab$end[i])/2
  functab$funcdiff[i] <- abs(functab$mid[i+1]-functab$mid[i])
}
functab$funcdiff[nrow(functab)] <- functab$start[nrow(functab)]
head(functab)
sum(functab$overGO != "")
sum(functab$eachpvalwithoutNA < 0.05 & functab$overthreshold_list %in% c(1:8) & !is.na(functab$genes.x))

#add color

functab$functabcol <- ifelse(functab$overthreshold_list == 1, "red", 
                 ifelse(functab$overthreshold_list == 2, "deepskyblue", 
                        ifelse(functab$overthreshold_list == 3, "dark green", 
                               ifelse(functab$overthreshold_list == 4, "green", 
                                      ifelse(functab$overthreshold_list == 5, "purple", 
                                             ifelse(functab$overthreshold_list == 6, "orange", 
                                             ifelse(functab$overthreshold_list == 7, "gold", 
                                             ifelse(functab$overthreshold_list == 8, "pink", NA))))))))

#37 unique effect function
#functab <- functab[functab$eachpvalwithoutNA < 0.05 & functab$overthreshold_list %in% c(1:5) & !is.na(functab$genes.x),]
dim(functab)
table(functab$Chr)
table(functab$group)
(table(functab$overthreshold_list))

functab$cleanfunction <- (gsub(", unclassified", "", gsub(", putative", "", gsub(", expressed", "", functab$rice_defline))))
functab$cleanfunction <- ifelse(is.na(functab$cleanfunction) | grepl("expressed protein", functab$cleanfunction), "", functab$cleanfunction)
functab$cleanfunction <- functab$genes

functab$eachsvshort <- paste("C. ", functab$overthreshold_list, " ", #"\n",
                             functab$cleanfunction, sep = "")
functab$eachsvshort <- gsub("& Sobic.", "\n", functab$eachsvshort)
functab <- functab[order(functab$overthreshold_list, decreasing = F),]
#label affected and significant genes
for(i in 1:nrow(functab)){
  circos.text(#x=functab$start[i], 
    x=functab$mid[i], 
    y=0.5, labels=functab$eachsvshort[i], 
    sector.index=functab$Chr[i], track.index= 4, facing="clockwise", 
    adj=c(0,0.), cex = 0.5, niceFacing = T, col = functab$functabcol[i])
}




labsorghum <- sorghumgenome
labsorghum$cumend <- cumsum(labsorghum$end)
labsorghum$cumstart <- labsorghum$cumend-labsorghum$end+1
labsorghum
gapdegree = 1

degreetab <- c()
for(i in 1:nrow(sorghumgenome)){
  cellstart <- get.cell.meta.data("cell.start.degree",sector.index = sorghumgenome$Chromosome[i])-(gapdegree)
  cellend <- get.cell.meta.data("cell.end.degree",sector.index = sorghumgenome$Chromosome[i])+(gapdegree)
  degreetab <- rbind.data.frame(degreetab, 
                                cbind.data.frame(sorghumgenome$Chromosome[i], cellstart, cellend, labsorghum$end[i]/(cellstart-cellend)))
}
degreetab
colnames(degreetab) <- c("Chromosome", "cellstart", "cellend", "onedegree")
degreetab[degreetab$Chromosome=="Chr03",]$onedegree <- labsorghum$end[which(labsorghum$Chromosome == "Chr03")]/((360-degreetab[degreetab$Chromosome=="Chr03",]$cellend)+degreetab[degreetab$Chromosome=="Chr03",]$cellstart)

newfunctab <- merge(x= functab, y = labsorghum, by.x = "Chr", by.y = "Chromosome", all.x = T, suffixes = c("", "_chrom")) %>% mutate(newstart = start+cumstart, newend = end+cumstart)
newfunctab <- merge(x= newfunctab, y= degreetab, by.x = "Chr", by.y = "Chromosome", all.x = T, suffixes = c("", "_degree"))
newfunctab  <- newfunctab %>% mutate(startdegree = (cellstart-(start/(onedegree)))) %>% 
  mutate(enddegree = (cellstart-(end/(onedegree))))


for(i in 1:nrow(newfunctab)){
  #for(i in 1:8){
  input <- newfunctab$startdegree[i]
  inputend <- newfunctab$enddegree[i]
  scalediff <- rescale(abs(abs(newfunctab$startdegree) -  abs(newfunctab$enddegree)), to=c(0,0.3))[i]
  draw.sector(start.degree=input, 
              end.degree= (input+0.3),#input+scalediff, 
              #end.degree = inputend,
              col =  newfunctab$functabcol[i], border = NA, clock.wise = F,
              rou1 = 0.55, rou2 = 0.30)
}

#text(0, 0, "Genomic \n structural variations \n among \n 347 sorghums", cex = 2)
#text(0, 0.90, "Overall SV", cex = 1, srt = 90)
#text(0, 0.67, "Abundance", cex = 1, srt = 90)
#text(0, 0.40, "Chi-square", cex = 1, srt = 90)

dev.off()
