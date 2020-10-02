library(openxlsx)
library(diveRsity)
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
library(Rgraphviz)

setwd("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/go/cluster")

right = function (string, char){
  substr(string,nchar(string)-(char-1),nchar(string))
}


left = function (string,char){
  substr(string,1,char)
}

#load gene
tab19 <- readRDS("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/347_go19_filtered.rds")
tab19 <- tab19[tab19$svlen >= 50,]
dim(tab19)
dim(tab19[!is.na(tab19$genes),])
dim(tab19[!is.na(tab19$CDS),])
tab19[!is.na(tab19$CDS),]%>% dplyr::group_by(chr) %>% summarise(n = n())
tab19[tab19$svtype == "DEL",] %>% dplyr::group_by(chr) %>% summarise(n = n())
tab19$svid <- paste(tab19$chr, tab19$pos, tab19$end, tab19$svtype, sep = "_")

#load gff
gff <- read.table("/Users/ksongsom/OneDrive/postdoc/SV/illumina/gene_anno.txt", quote = "", fill = T, header = T, sep = "\t",na.strings=c("","NA") ,stringsAsFactors=FALSE,check.names=FALSE)

#1 filter only in gene
effgene <- tab19[!is.na(tab19$genes),]
head(effgene)
dim(effgene)

effgene_sv <- effgene$svid
effgene_go <- effgene$go
unlist(strsplit(effgene_go, "& "))
effgene_gene <- effgene$genes
effgene_gene_interest <- unique(unlist(strsplit(effgene_gene, "& ")))
str(effgene_gene_interest)
length(effgene_gene_interest)

#2 filter only protein coding
#2.1 filter tab to only cluster-specific SVs
data_k5_cluster_all <- read.xlsx(paste("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/data_k5_cluster_allnoNA_weight_withtotal.xlsx", sep = ""))
head(data_k5_cluster_all)
clustersv <- (data_k5_cluster_all[!is.na(data_k5_cluster_all$overthreshold_list) & data_k5_cluster_all$eachpvalwithoutNA < 0.05,]$eachsv)

effprot <- tab19[!is.na(tab19$CDS) & tab19$svid %in% clustersv,]
dim(effprot)
effprot$CDS
effprot_prot <- effprot$CDS
effprot_prot_interest <- unique(unlist(strsplit(effprot_prot, "& ")))

#3 GO enrichmen analysis
#https://bioconductor.org/packages/release/bioc/html/topGO.html
#https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf

library(topGO)
library(ALL)
library(genefilter)
library(ViSEAGO)

#my case
gogff <- gff[!is.na(gff$GO),c("Name", "GO")]
head(gogff)
mygeneID2GO <- strsplit(gogff$GO, ",")

names(mygeneID2GO) <- gsub("Sobic_", "Sobic.", gsub("\\..*", "", gsub("Sobic.", "Sobic_", gogff$Name)))
mygeneID2GO <- mygeneID2GO[unique(names(mygeneID2GO))]
myGO2geneID <- inverseList(mygeneID2GO)

effprot_prot_interest <- unique(gsub("Sobic_", "Sobic.", gsub("\\..*", "", gsub("Sobic.", "Sobic_", effprot_prot_interest))))
length(effprot_prot_interest)

#gene names
mygeneList_effgene_gene <- factor(as.integer(names(mygeneID2GO) %in% effgene_gene_interest))
names(mygeneList_effgene_gene) <- names(mygeneID2GO)

mygeneList_effprot_prot <- factor(as.integer(names(mygeneID2GO) %in% effprot_prot_interest))
names(mygeneList_effprot_prot) <- names(mygeneID2GO)


myGOdata_gene <- new("topGOdata", ontology = "MF", allGenes = mygeneList_effgene_gene,
                     annot = annFUN.gene2GO, gene2GO = mygeneID2GO)

myGOdata_prot <- new("topGOdata", ontology = "MF", allGenes = mygeneList_effprot_prot,
                     annot = annFUN.gene2GO, gene2GO = mygeneID2GO)

topDiffGenes <- function(allScore){
  return(allScore < 0.01)
}

#stat summary
#CDS
numGenes(myGOdata_prot)
gs_prot <- geneScore(myGOdata_prot, whichGenes = effprot_prot_interest)

sg_prot <- sigGenes(myGOdata_prot)

length(effprot_prot_interest)
length(sg_prot)
numSigGenes(myGOdata_prot)

usedGO(myGOdata_prot)
termStat(myGOdata_prot, usedGO(myGOdata_prot))

#gene
numGenes(myGOdata_gene)
gs_gene <- geneScore(myGOdata_gene, whichGenes = effgene_gene_interest)

sg_gene <- sigGenes(myGOdata_gene)

length(effgene_gene_interest)
length(sg_gene)
numSigGenes(myGOdata_gene)

usedGO(myGOdata_gene)
termStat(myGOdata_gene, usedGO(myGOdata_gene))


#classic count test CDS
gene.universe_prot <- genes(myGOdata_prot)
go.genes_prot <- genesInTerm(myGOdata_prot, usedGO(myGOdata_prot))[[1]]
sig.genes_prot <- sigGenes(myGOdata_prot)

my.group_prot <- new("classicCount", testStatistic = GOFisherTest, 
                     name = "fisher", allMembers = gene.universe_prot, 
                     groupMembers = go.genes_prot,
                     sigMembers = sig.genes_prot)
contTable(my.group_prot)
runTest(my.group_prot)

test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher_prot <- getSigGroups(myGOdata_prot, test.stat)

test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
resultKS_prot <- getSigGroups(myGOdata_prot, test.stat)

test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight_prot <- getSigGroups(myGOdata_prot, test.stat)

resultFisher_prot
resultKS_prot
resultWeight_prot
score(resultWeight_prot)
hist(score(resultWeight_prot), 50, xlab = "p-values")
geneData(resultWeight_prot)

top10_prot <- as.integer(geneData(resultWeight_prot)[2]*0.1)
allRes_prot <- GenTable(myGOdata_prot, classic = resultFisher_prot, 
                        KS = resultKS_prot, weight = resultWeight_prot,
                        orderBy = "weight", ranksOf = "classic", topNodes = top10_prot)
goID_prot <- allRes_prot[2, "GO.ID"]

print(showGroupDensity(myGOdata_prot, goID_prot, ranks = TRUE))

showSigOfNodes(myGOdata_prot, score(resultWeight_prot), firstSigNodes = 5, useInfo = 'all')
printGraph(myGOdata_prot, resultFisher_prot, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = F)
printGraph(myGOdata_prot, resultFisher_prot, firstSigNodes = 1, fn.prefix = "tGO", useInfo = "all", pdfSW = F)
printGraph(myGOdata_prot, resultWeight_prot, firstSigNodes = 10, resultFisher_prot, fn.prefix = "tGO", useInfo = "def")
printGraph(myGOdata_prot, resultWeight_prot, firstSigNodes = 20, resultFisher_prot, fn.prefix = "tGO", useInfo = "def")

#write.csv(allRes_prot, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/go_CDS_cluster.csv")







#classic count test gene
gene.universe_gene <- genes(myGOdata_gene)
go.genes_gene <- genesInTerm(myGOdata_gene, usedGO(myGOdata_gene))[[1]]
sig.genes_gene <- sigGenes(myGOdata_gene)

my.group_gene <- new("classicCount", testStatistic = GOFisherTest, 
                     name = "fisher", allMembers = gene.universe_gene, 
                     groupMembers = go.genes_gene,
                     sigMembers = sig.genes_gene)
contTable(my.group_gene)
runTest(my.group_gene)

test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher_gene <- getSigGroups(myGOdata_gene, test.stat)

test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
resultKS_gene <- getSigGroups(myGOdata_gene, test.stat)

resultFisher_gene
resultKS_gene

test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight_gene <- getSigGroups(myGOdata_gene, test.stat)

resultWeight_gene
score(resultWeight_gene)
hist(score(resultWeight_gene), 50, xlab = "p-values")
geneData(resultWeight_gene)

top10_gene <- as.integer(geneData(resultWeight_gene)[2]*0.1)
allRes_gene <- GenTable(myGOdata_gene, classic = resultFisher_gene, 
                        KS = resultKS_gene, weight = resultWeight_gene,
                        orderBy = "weight", ranksOf = "classic", topNodes = top10_gene)
goID_gene <- allRes_gene[2, "GO.ID"]

print(showGroupDensity(myGOdata_gene, goID_gene, ranks = TRUE))

showSigOfNodes(myGOdata_gene, score(resultWeight_gene), firstSigNodes = 5, useInfo = 'all')
printGraph(myGOdata_gene, resultFisher_gene, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

printGraph(myGOdata_gene, resultWeight_gene, firstSigNodes = 10, resultFisher_gene, fn.prefix = "tGO", useInfo = "def")


#write.csv(allRes_gene, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/go_gene_cluster.csv")


#test fisher between stress related gene vs 
#extract cds stress related go
myGOdata_prot
gene.universe_prot <- genes(myGOdata_prot)
go.genes_prot <- genesInTerm(myGOdata_prot, usedGO(myGOdata_prot))[[1]]
test <- GenTable(myGOdata_prot, usedGO(myGOdata_prot))
test
all_count_prot <- as.integer(geneData(resultWeight_prot)[2])
allRes_count_prot <- GenTable(myGOdata_prot, classic = resultFisher_prot, 
                              KS = resultKS_prot, weight = resultWeight_prot,
                              orderBy = "weight", ranksOf = "classic", topNodes = all_count_prot)

saveRDS(allRes_count_prot, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/go_216_cds_cluster.rds")

#compare between all and cluster specific
allRes_count_prot_all <- readRDS("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/go_216_cds.rds")
allRes_count_prot_cluster <- readRDS("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/go_216_cds_cluster.rds")

write.xlsx(allRes_count_prot_all, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/go_all.xlsx")
write.xlsx(allRes_count_prot_cluster, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/go_cluster.xlsx")


mprot <- merge(x = allRes_count_prot_all, y = allRes_count_prot_cluster[,!(colnames(allRes_count_prot_cluster) %in% c("Term", "Annotated"))], 
               by = "GO.ID", all.x = T, all.y = T, suffixes = c(" all", " cluster specific"))
head(mprot)
mprot$`weight all` <- as.numeric(mprot$`weight all`)
mprot$`weight cluster specific` <- as.numeric(mprot$`weight cluster specific`)
mprot <- mprot[order(mprot$`weight all`, decreasing = F),]

write.xlsx(mprot, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/go_216_cluster.xlsx")


mprot2 <- merge(x = allRes_count_prot_cluster, y = allRes_count_prot_all[,!(colnames(allRes_count_prot_all) %in% c("Term", "Annotated"))], 
               by = "GO.ID", all.x = T, all.y = T, suffixes = c(" cluster specific", " all"))
head(mprot2)
mprot2$`weight all` <- as.numeric(mprot2$`weight all`)
mprot2$`weight cluster specific` <- as.numeric(mprot2$`weight cluster specific`)
mprot2 <- mprot2[order(mprot2$`weight cluster`, decreasing = F),]


write.xlsx(mprot2, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/go_216_cluster_sv.xlsx")