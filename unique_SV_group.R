library(tidyr)
library(plyr)
library(stringi)
library(stringr)
library(dplyr)
library(lme4)
library(lmerTest)
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
library(readxl)
library(rlang)

#load map
map <- (read_excel("/Users/ksongsom/OneDrive/postdoc/NAM/Panel/BAP/update_sorghum_location.xlsx", sheet = 1)) #from google sheet
#update elevation
map$ID <- gsub("PI ", "", map$Genotype)
head(map)
map <- map %>% separate(Coordinate, c("lat", "lon"), sep = ", ")

#load cluster
cluster <- read.csv("/Users/ksongsom/OneDrive/postdoc/publication/terra/code/files/five_cluster.csv")
head(cluster)
dim(cluster)

#background
background <- read.csv("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/dat_back.cvs")
head(background)

#load geno
tab <- readRDS("/Users/ksongsom/OneDrive/postdoc/publication/terra/code/files/347_go.rds")
dim(tab)
head(tab)

#merge map cluster gene
mapclus <- merge(map, cluster, by.x = "ID", by.y = "PI")
dim(map)
dim(cluster)
dim(mapclus)
str(mapclus)
mapclus$Race <- ifelse(grepl("-", mapclus$Race), "Mixed", mapclus$Race)
mapclus2 <- mapclus %>% group_by(k5_cluster_all, Race) %>% summarise(number_genotypes = n())
uncutoff <- 0.7

criteria <- c(quo(k3_cluster_all), quo(k5_cluster_all), quo(k7_cluster_all), quo(k3_cluster_Ethiopia), quo(Race), quo(Type))
datlist <- list()
for(j in seq_along(criteria)){
  j = 2
  critdat <- c()
      for(i in 1:nrow(tab)){
        eachrow <- tab[i,]
        eachsv <- eachrow$svid
        eachtab <- as.data.frame(t(eachrow[,(which(colnames(eachrow)=="format")+1):(which(colnames(eachrow)=="Rio"))])) %>% rlang::set_names("individual")
        eachtab$ID <- rownames(eachtab)
        mtab <- merge(eachtab, mapclus, by.x = "ID", by.y = "ID")
      
        eachmapclus <- mapclus %>% group_by(!!criteria[[j]]) %>% summarise(len = n())
        eachmtab <- mtab %>% group_by(!!criteria[[j]]) %>% summarise(sumindi = sum(individual, na.rm = T)) %>% 
          mutate(percent = sumindi/eachmapclus$len)
        expectsum <- sum(eachmtab$sumindi)
        expectport <- expectsum/sum(eachmapclus$len)
        eachmtab$expectpercent <- eachmapclus$len/sum(eachmapclus$len)
        eachmtab$expectindi <- eachmapclus$len*expectport
        eachmtab$percentweight <- eachmtab$sumindi/expectsum#*sum(eachmapclus$len)
        if(length(unique(eachmtab$sumindi)) == 1){
          eachmtab$sumindi[length(eachmtab$sumindi)] <- eachmtab$sumindi[length(eachmtab$sumindi)]+1
        }
        eachmtab <- cbind.data.frame(eachmtab, eachmapclus$len)
        eachpval <- chisq.test(eachmtab$sumindi, eachmtab$expectindi)$p.value
        #eachpval <- chisq.test(eachmtab$percentweight, eachmtab$expectpercent)$p.value
        eachmtab <- eachmtab %>% mutate(number = paste(round(percentweight, 2), " (", sumindi, ")", sep = ""))#, "/", eachmapclus$len
        overthreshold <- ifelse(eachmtab$percentweight>=uncutoff, 1, NA)
        overthreshold_list <- ifelse(all(is.na(overthreshold)), NA, paste(eachmtab[,1][!is.na(overthreshold)], sep = ", "))

        eachresult <- t(eachmtab[,"number"])
        colnames(eachresult) <- paste(as.character(criteria[[j]][2]), eachmtab[,1], sep = "_")
        eachresult <- cbind.data.frame(eachsv, eachresult, eachpval, overthreshold_list)
        critdat <- rbind.data.frame(critdat, eachresult)
        print(paste(j, i, j/length(criteria), i/nrow(tab), sep = "_"))
      }
  datlist <- list(datlist,assign(paste("data", as.character(criteria[[j]][2]), sep = "_"), critdat ))
  write.csv(critdat, paste("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/","data_", as.character(criteria[[j]][2], ".csv"), sep = ""), row.names = F)
}

#uniqueness among 5 k-groups
data_k5_cluster_all <- read.csv("/Users/ksongsom/OneDrive/postdoc/publication/terra/code/files/data_k5_cluster_all.csv")
head(data_k5_cluster_all)
sum(data_k5_cluster_all$eachpval < 0.05) #non of them pass the Chi-square test
dim(data_k5_cluster_all)
sum(!is.na(data_k5_cluster_all$overthreshold_list)) #1548 pass threshold of .7
table(data_k5_cluster_all[(!is.na(data_k5_cluster_all$overthreshold_list)),]$overthreshold_list)


#uniqueness in function from 5 k-groups
tab200 <- readRDS("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/347_go200_filtered.rds")
head(tab200)
tab200$svid
summary(tab200$maf_col)

k = 2   
dtinput <- data_k5_cluster_all #or datlist[[k]]
as.character(criteria[[k]][2])
sum(!is.na(dtinput$overthreshold_list))
sigdt <- dtinput[(!is.na(dtinput$overthreshold_list)),]
table(dtinput[(!is.na(dtinput$overthreshold_list)),]$overthreshold_list)
table(gsub("\\_.*", "", sigdt$eachsv))

sum(sigdt$eachsv %in% tab$svid)
functab <- merge(sigdt, tab200[, c("svid", "genes", "rice_defline")], by.x = "eachsv", by.y = "svid", all.x = T, all.y = F) #, "genes200", "rice_defline200" 
dim(functab)
head(functab)
functab <- functab[!is.na(functab$genes),]
table(functab$rice_defline)
table(functab$rice_defline)[order(table(functab$rice_defline), decreasing = F)]

#five group

k = 2   
dtinput <- datlist[[k]]
as.character(criteria[[k]][2])
sum(!is.na(dtinput$overthreshold_list))
sigdt <- dtinput[(!is.na(dtinput$overthreshold_list)),]
head(sigdt)
table(dtinput[(!is.na(dtinput$overthreshold_list)),]$overthreshold_list)
table(gsub("\\_.*", "", sigdt$eachsv))

sum(sigdt$eachsv %in% tab$svid)
functab <- merge(sigdt, tab200[, c("svid", "genes", "rice_defline")], by.x = "eachsv", by.y = "svid", all.x = T, all.y = F) #, "genes200", "rice_defline200" 
dim(functab)
head(functab)
functab <- functab[!is.na(functab$genes),]
table(functab$overthreshold_list)
table(gsub("\\_.*", "", functab$eachsv))
write.csv(functab, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/function_5group.csv", row.names <- F)
table(functab$rice_defline)
table(functab$rice_defline)[order(table(functab$rice_defline), decreasing = F)]

f1 <- table(functab[functab$overthreshold_list == 1,]$rice_defline)
f2 <- table(functab[functab$overthreshold_list == 2,]$rice_defline)
f3 <- table(functab[functab$overthreshold_list == 3,]$rice_defline)
f4 <- table(functab[functab$overthreshold_list == 4,]$rice_defline)
f5 <- table(functab[functab$overthreshold_list == 5,]$rice_defline)

uf1 <- setdiff(names(f1), c(names(f2), names(f3), names(f4), names(f5)))
uf2 <- setdiff(names(f2), c(names(f1), names(f3), names(f4), names(f5)))
uf3 <- setdiff(names(f3), c(names(f2), names(f1), names(f4), names(f5)))
uf4 <- setdiff(names(f4), c(names(f2), names(f3), names(f1), names(f5)))
uf5 <- setdiff(names(f5), c(names(f2), names(f3), names(f4), names(f1)))

write.csv(uf1, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/function_1.csv", row.names <- F)
write.csv(uf2, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/function_2.csv", row.names <- F)
write.csv(uf3, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/function_3.csv", row.names <- F)
write.csv(uf4, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/function_4.csv", row.names <- F)
write.csv(uf5, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/function_5.csv", row.names <- F)

#three group Ethiopia

k = 4   
dtinput <- datlist[[k]]
as.character(criteria[[k]][2])
sum(!is.na(dtinput$overthreshold_list))
sigdt <- dtinput[(!is.na(dtinput$overthreshold_list)),]
table(dtinput[(!is.na(dtinput$overthreshold_list)),]$overthreshold_list)
table(gsub("\\_.*", "", sigdt$eachsv))

sum(sigdt$eachsv %in% tab$svid)
functab <- merge(sigdt, tab200[, c("svid", "genes", "rice_defline")], by.x = "eachsv", by.y = "svid", all.x = T, all.y = F) #, "genes200", "rice_defline200" 
dim(functab)
head(functab)
functab <- functab[!is.na(functab$genes),]
table(functab$rice_defline)
table(functab$rice_defline)[order(table(functab$rice_defline), decreasing = F)]

f1 <- table(functab[functab$overthreshold_list == 1,]$rice_defline)
f2 <- table(functab[functab$overthreshold_list == 2,]$rice_defline)
f3 <- table(functab[functab$overthreshold_list == 3,]$rice_defline)
f4 <- table(functab[functab$overthreshold_list == 4,]$rice_defline)
f5 <- table(functab[functab$overthreshold_list == 5,]$rice_defline)

uf1 <- setdiff(names(f1), c(names(f2), names(f3), names(f4), names(f5)))
uf2 <- setdiff(names(f2), c(names(f1), names(f3), names(f4), names(f5)))
uf3 <- setdiff(names(f3), c(names(f2), names(f1), names(f4), names(f5)))
uf4 <- setdiff(names(f4), c(names(f2), names(f3), names(f1), names(f5)))
uf5 <- setdiff(names(f5), c(names(f2), names(f3), names(f4), names(f1)))

#Race

k = 5   
dtinput <- datlist[[k]]
as.character(criteria[[k]][2])
sum(!is.na(dtinput$overthreshold_list))

sigdt <- dtinput[(!is.na(dtinput$overthreshold_list)),]
table(dtinput[(!is.na(dtinput$overthreshold_list)),]$overthreshold_list)
table(gsub("\\_.*", "", sigdt$eachsv))

sum(sigdt$eachsv %in% tab$svid)
functab <- merge(sigdt, tab200[, c("svid", "genes", "rice_defline")], by.x = "eachsv", by.y = "svid", all.x = T, all.y = F) #, "genes200", "rice_defline200" 
dim(functab)
head(functab)
functab <- functab[!is.na(functab$genes),]
table(functab$overthreshold_list)
table(functab$rice_defline)
table(functab$rice_defline)[order(table(functab$rice_defline), decreasing = F)]

f1 <- table(functab[functab$overthreshold_list == "Guinea",]$rice_defline)
names(f1)
f2 <- table(functab[functab$overthreshold_list == "NA",]$rice_defline)
names(f2)
write.csv(names(f1), "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/function_race_g.csv", row.names <- F)


#Type

k = 6
dtinput <- datlist[[k]]
dim(dtinput)
as.character(criteria[[k]][2])
sum(!is.na(dtinput$overthreshold_list))
dtinput[dtinput$eachpval < 0.05,]
sigdt <- dtinput[(!is.na(dtinput$overthreshold_list)),]
table(dtinput$overthreshold_list)
table(dtinput[(!is.na(dtinput$overthreshold_list)),]$overthreshold_list)
table(gsub("\\_.*", "", sigdt$eachsv))


sum(sigdt$eachsv %in% tab$svid)
functab <- merge(sigdt, tab200[, c("svid", "genes", "rice_defline")], by.x = "eachsv", by.y = "svid", all.x = T, all.y = F) #, "genes200", "rice_defline200" 
dim(functab)
head(functab)
functab <- functab[!is.na(functab$genes),]
table(functab$overthreshold_list)
table(functab$rice_defline)
table(functab$rice_defline)[order(table(functab$rice_defline), decreasing = F)]

f1 <- table(functab[functab$overthreshold_list == "Cellulosic",]$rice_defline)
f2 <- table(functab[functab$overthreshold_list == "Sweet",]$rice_defline)

uf1 <- setdiff(names(f1), c(names(f2)))
uf2 <- setdiff(names(f2), c(names(f1)))

write.csv(uf1, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/function_type_cel.csv", row.names <- F)
write.csv(uf2, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/function_type_sweet.csv", row.names <- F)


#race VS cluster
head(background)
head(cluster)
sum(cluster$ID %in% background$Taxa)
cluster_back <- merge(cluster, background, by.x = "ID", by.y = "Taxa", all.x = T, all.y = F)
head(cluster_back)
cluster_back[,c("k5_cluster_all", "Race")]
reshape(cluster_back[,c("k5_cluster_all", "Race")], direction = "wide", idvar = "Race", timevar = "k5_cluster_all")
cluster_race_tab <-as.data.frame(cluster_back[,c("k5_cluster_all", "Race")] %>% group_by(k5_cluster_all, Race) %>% tally())
write.csv(cluster_race_tab, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/cluster_race_tab.csv", row.names <- F)
