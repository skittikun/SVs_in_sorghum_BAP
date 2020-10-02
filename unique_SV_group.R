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
library(tidyverse)
library(readxl)
library(openxlsx)

options(warn=-1)
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
tab <- tab[tab$svlen >= 50,]
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
locutoff <- 0.3

#only 5 cluster seperation
criteria <- c(quo(k3_cluster_all), quo(k5_cluster_all), quo(k7_cluster_all), quo(k3_cluster_Ethiopia), quo(Race), quo(Type))
datlist <- list()
expectedlist <- list()
bootlist <- list()
svlist <- list()
svpval <- c()
#for(j in seq_along(criteria)){
  j = 2
  critdat <- c()
  expectedtab <- c()
  
      for(i in 1:nrow(tab)){
      
        eachrow <- tab[i,]
        
        #eachrow <- tab[tab$chr == "Chr01" & tab$pos == 10091090 & tab$end == 10092908,] for testing Chr01_10091090_10092908_DEL
        eachsv <- eachrow$svid
        eachtab <- as.data.frame(t(eachrow[,(which(colnames(eachrow)=="format")+1):(which(colnames(eachrow)=="Rio"))])) %>% rlang::set_names("individual")
        eachtab$ID <- rownames(eachtab)
        mtab <- merge(eachtab, mapclus, by.x = "ID", by.y = "ID")
      
        #calculate expected and percentage of the sv
        eachmapclus <- mapclus %>% group_by(!!criteria[[j]]) %>% summarise(len = n())
        eachmtab <- mtab %>% group_by(!!criteria[[j]]) %>% summarise(obsindi = sum(individual, na.rm = T), 
                                                                     totalwithoutNA = sum(!is.na(individual))) %>% 
          mutate(percentfrom347 = obsindi/eachmapclus$len) %>% 
          mutate(percentwithoutNA = obsindi/totalwithoutNA) 
        obssum <- sum(eachmtab$obsindi)
       
        eachmtab$expectpercent347 <- eachmapclus$len/sum(eachmapclus$len)
        eachmtab$expectpercentwithoutNA <- eachmtab$totalwithoutNA/sum(eachmtab$totalwithoutNA)
        
        eachmtab$percentabundance <- eachmtab$obsindi/obssum
        
        #calculate weight observation to compare
        (prod(eachmtab$totalwithoutNA)/eachmtab$totalwithoutNA)
        eachmtab$weightobsindi <- eachmtab$obsindi*prod(eachmtab$totalwithoutNA)/eachmtab$totalwithoutNA
        eachmtab$weightobsindi/prod(eachmtab$totalwithoutNA)
        eachmtab$perweightobsindi <- eachmtab$weightobsindi/sum(eachmtab$weightobsindi)
        
        if(length(unique(eachmtab$obsindi)) == 1){
          eachmtab$obsindi[length(eachmtab$obsindi)] <- eachmtab$obsindi[length(eachmtab$obsindi)]+1
        }
        eachmtab <- cbind.data.frame(eachmtab, eachmapclus$len)
        eachpval347 <- chisq.test(eachmtab$obsindi, p =  eachmtab$expectpercent347)$p.value
        eachpvalwithoutNA <- chisq.test(eachmtab$obsindi, p= eachmtab$expectpercentwithoutNA)$p.value
        
        
        #make table by showing the percentwwithoutNA and in parenthesis with actual counts
        eachmtab <- eachmtab %>% mutate(number = paste(round(percentabundance, 2), " (", obsindi, "/", totalwithoutNA, ")", sep = ""))#, "/", eachmapclus$len
        
        #criteria of unique in clusters have more than 70% within group and other group less than 30%
        overthreshold <- ifelse(round(eachmtab$percentabundance, digits = 2)>=uncutoff, 1, 2)
        overthreshold_list <- ifelse(all(overthreshold == 2) | sum(overthreshold==1)>1, NA, 
                                     paste(eachmtab[,1][overthreshold==1], collapse = ", "))
        
        eachresult <- t(eachmtab[,"number"])
        colnames(eachresult) <- paste(as.character(criteria[[j]][2]), eachmtab[,1], sep = "_")
        eachresult <- cbind.data.frame(eachsv, eachresult, eachpvalwithoutNA, overthreshold_list)
        critdat <- rbind.data.frame(critdat, eachresult)
        print(paste(j, i, j/length(criteria), i/nrow(tab), "over", sum(critdat$eachpvalwithoutNA < 0.05 & !is.na(critdat$overthreshold_list)), sep = "_"))
        expectedtab <- rbind.data.frame(expectedtab, cbind.data.frame(eachsv,as.data.frame(t(eachmtab[,c("expectpercentwithoutNA")])),obssum))
      
        
        #if chisquare pval < 0.05 and percent abundance is higher than .7 need to determine if it is by chance via permutation
        if(eachpvalwithoutNA < 0.05 & !(all(is.na(overthreshold)) | sum(overthreshold == 2, na.rm = T) == 5 | sum(overthreshold==1, na.rm = T)>1)){
          overgroup <- which(eachmtab$percentabundance >= uncutoff)

          #permutation
          set.seed(40)
          bootn = 1000
          bootlist_k <- list()
          pval_k <-c()
          for(k in 1:(bootn)){
            pereachtab <- cbind.data.frame(eachtab$ID, sample(eachtab$individual))
            colnames(pereachtab) <- c("ID", "individual")
            permtab <- merge(pereachtab, mapclus, by.x = "ID", by.y = "ID")
            pereachmtab <- permtab %>% group_by(!!criteria[[j]]) %>% summarise(perobsindi = sum(individual, na.rm = T), 
                                                                               pertotalwithoutNA = sum(!is.na(individual))) %>% 
              mutate(perpercentfrom347 = perobsindi/eachmapclus$len) 
            pereachmtab$perpercentweightwithouNA <- pereachmtab$perobsindi/pereachmtab$pertotalwithoutNA
            pereachmtab$perexpectperwithoutNA <- pereachmtab$pertotalwithoutNA/sum(pereachmtab$pertotalwithoutNA)
            
            pereachmtab$percentabundance <- pereachmtab$perobsindi/sum(pereachmtab$perobsindi)
            
            bootlist_k[[k]] <-  cbind.data.frame(eachmtab, pereachmtab[,-1])
           
            #if any of percent abundance over cutoff (.7), other less than .3 and that .7 happens in the same cluster count the indicent as one
            pval_k[k] <- ifelse(sum(pereachmtab$percentabundance >= uncutoff)==0, 0,
                                ifelse(sum(pereachmtab$percentabundance >= uncutoff)==1 & 
                                         sum(pereachmtab$percentabundance < locutoff) == 4 &
                                         which(pereachmtab$percentabundance >= uncutoff)==overgroup, 1, 0))
            print(paste("i=", i, "k=", k, "pvalper=", sum(svpval$`sum(pval_k == 1)`==1), sep = "_")) 
          }
          
          bootlist[[i]] <- bootlist_k
          names(bootlist)[i] <- eachsv
          
          #
          unk <- (unlist(bootlist_k))
          mean_boot <- c(mean(as.numeric(unk[grepl("percentabundance1", names(unk))])),
                         mean(as.numeric(unk[grepl("percentabundance2", names(unk))])),
                         mean(as.numeric(unk[grepl("percentabundance3", names(unk))])),
                         mean(as.numeric(unk[grepl("percentabundance4", names(unk))])),
                         mean(as.numeric(unk[grepl("percentabundance5", names(unk))])))
          
          sd_boot <- c(sd(as.numeric(unk[grepl("percentabundance1", names(unk))])),
                       sd(as.numeric(unk[grepl("percentabundance2", names(unk))])),
                       sd(as.numeric(unk[grepl("percentabundance3", names(unk))])),
                       sd(as.numeric(unk[grepl("percentabundance4", names(unk))])),
                       sd(as.numeric(unk[grepl("percentabundance5", names(unk))])))
          
          eachmtab$mean_boot <- mean_boot
          eachmtab$sd_boot <- sd_boot
          svpval <- rbind.data.frame(svpval, cbind.data.frame(eachsv, chisq.test(eachmtab$obsindi, p= eachmtab$expectpercentwithoutNA)$p.value, sum(pval_k==1)))
          
          svlist[[i]] <- eachmtab
          names(svlist)[i] <- eachsv
          
          
        }
        }
  expectedlist[[j]] <- expectedtab 
  names(expectedlist[[j]]) <- paste("data", as.character(criteria[[j]][2]), sep = "_")
  datlist[[j]] <- critdat
  names(datlist)[j] <- paste("data", as.character(criteria[[j]][2]), sep = "_")
  
  colnames(svpval) <- c("eachsv", "pval_chisquare", "pval_permutation")
  #colnames(svpval) <- colnames(cbind.data.frame(eachsv, chisq.test(eachmtab$obsindi, p= eachmtab$expectpercentwithoutNA)$p.value, sum(pval_k==1)))
  
  #write.csv(critdat, paste("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/","data_", as.character(criteria[[j]][2]),"noNA_weight_withtotal", ".csv", sep = ""), row.names = F)
  #write.xlsx(critdat, paste("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/","data_", as.character(criteria[[j]][2]),"noNA_weight_withtotal", ".xlsx", sep = ""))
  #write.csv(expectedtab, "/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/expected_percentage_withoutNA.csv", row.names = F, quote = F)
  #saveRDS(expectedtab,"/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/expected_percentage_withoutNA.rds")
  #write.xlsx(expectedtab, "/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/expected_percentage_withoutNA.xlsx")

# 
  #saveRDS(datlist, "/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/datlist.rds")
  saveRDS(expectedlist, "/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/expectedlist.rds")
  
saveRDS(svlist, "/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/bootstrap_list_sumeachsv_permutation.rds")
saveRDS(bootlist, "/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/bootstrap_list_permutation.rds")
write.xlsx(svpval, "/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/pval_chisquare_bootstrap_permutation.xlsx")

#check permutation
svpval <- read.xlsx("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/pval_chisquare_bootstrap_permutation.xlsx")
svlist <- readRDS("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/bootstrap_list_sumeachsv_permutation.rds")
bootlist <- readRDS("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/bootstrap_list_permutation.rds")
critdat <- read.xlsx(paste("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/","data_", as.character(criteria[[j]][2]),"noNA_weight_withtotal", ".xlsx", sep = ""))

str(critdat)
dim(critdat)
dim(critdat[critdat$eachpvalwithoutNA < 0.05,])
dim(critdat[critdat$eachpvalwithoutNA < 0.05 & !is.na(critdat$overthreshold_list),])

cluster1 <- as.numeric(gsub("\\(.*", "", critdat$k5_cluster_all_1))
cluster2 <- as.numeric(gsub("\\(.*", "", critdat$k5_cluster_all_2))
cluster3 <- as.numeric(gsub("\\(.*", "", critdat$k5_cluster_all_3))
cluster4 <- as.numeric(gsub("\\(.*", "", critdat$k5_cluster_all_4))
cluster5 <- as.numeric(gsub("\\(.*", "", critdat$k5_cluster_all_5))

critdat <- cbind.data.frame(critdat, cluster1, cluster2, cluster3, cluster4, cluster5)
meetthreshold <- ifelse((cluster1>=uncutoff | cluster2>=uncutoff | cluster3>=uncutoff |
                          cluster4>=uncutoff | cluster5>=uncutoff) &
                          (cluster1>=locutoff | cluster2>=locutoff | cluster3>=locutoff |
                             cluster4>=locutoff | cluster5>=locutoff), 1, NA)
sum(cluster1>=uncutoff)
sum(cluster2>=uncutoff)
sum(cluster3>=uncutoff)
sum(cluster4>=uncutoff)
sum(cluster5>=uncutoff)



#find unique genes among grouping method
datlist <- readRDS("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/datlist.rds")
#uniqueness among 5 k-groups
data_k5_cluster_all <- read.xlsx(paste("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/","data_", as.character(criteria[[j]][2]),"noNA_weight_withtotal", ".xlsx", sep = ""))
head(data_k5_cluster_all)
sum(data_k5_cluster_all$eachpvalwithoutNA < 0.05) #7208 significantly different
length(svlist)
dim(data_k5_cluster_all)
#by clusters
dim(svpval)
table(svpval$pval_permutation)
data_k5_cluster_all[!is.na(data_k5_cluster_all$overthreshold_list),]
dim(data_k5_cluster_all[data_k5_cluster_all$eachpvalwithoutNA < 0.05 & !is.na(data_k5_cluster_all$overthreshold_list),])  #1570 pass threshold of .7 and on other more than .3
table(critdat[critdat$eachpvalwithoutNA < 0.05 & !is.na(critdat$overthreshold_list),]$overthreshold_list) #number of SV unique to each group


#uniqueness in functions from 5 k-groups
tab200 <- readRDS("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/347_go200_filtered.rds")
tab200 <- tab200[tab200$svlen >= 50, ]

#only CDS affected
tab200 <- tab200[!is.na(tab200$CDS),]
colnames(tab200)
tab200$svid
summary(tab200$maf_col)

#five groups
k = 2   
dtinput <- critdat
as.character(criteria[[k]][2])
sum(!is.na(dtinput$overthreshold_list))
sigdt <- dtinput[(!is.na(dtinput$overthreshold_list)) & dtinput$eachpvalwithoutNA < 0.05,]
head(sigdt)
table(sigdt$overthreshold_list)
table(gsub("\\_.*", "", sigdt$eachsv))

sum(sigdt$eachsv %in% tab$svid)
functab <- merge(sigdt, tab200[, c("svid", "genes", "go", "rice_defline")], by.x = "eachsv", by.y = "svid", all.x = T, all.y = F) #, "genes200", "rice_defline200" 
dim(tab200)
dim(functab)
head(functab)
functab <- functab[!is.na(functab$genes),]
table(functab$overthreshold_list)
table(gsub("\\_.*", "", functab$eachsv))
#write.csv(functab, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/function_5group_noNA_weight.csv", row.names <- F)
table(functab$rice_defline)
table(functab$rice_defline)[order(table(functab$rice_defline), decreasing = F)]

#unique functions
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

length(f1)
length(f2)
length(f3)
length(f4)
length(f5)
length(uf1)
length(uf2)
length(uf3)
length(uf4)
length(uf5)


write.csv(uf1, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/function_1_noNA_weight.csv", row.names <- F)
write.csv(uf2, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/function_2_noNA_weight.csv", row.names <- F)
write.csv(uf3, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/function_3_noNA_weight.csv", row.names <- F)
write.csv(uf4, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/function_4_noNA_weight.csv", row.names <- F)
write.csv(uf5, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/function_5_noNA_weight.csv", row.names <- F)

#unique gene affected by CDS SV
f1 <- table(unlist(str_split(functab[functab$overthreshold_list == 1,]$genes, "& ")))
f2 <- table(unlist(str_split(functab[functab$overthreshold_list == 2,]$genes, "& ")))
f3 <- table(unlist(str_split(functab[functab$overthreshold_list == 3,]$genes, "& ")))
f4 <- table(unlist(str_split(functab[functab$overthreshold_list == 4,]$genes, "& ")))
f5 <- table(unlist(str_split(functab[functab$overthreshold_list == 5,]$genes, "& ")))

uf1 <- setdiff(names(f1), c(names(f2), names(f3), names(f4), names(f5)))
uf2 <- setdiff(names(f2), c(names(f1), names(f3), names(f4), names(f5)))
uf3 <- setdiff(names(f3), c(names(f2), names(f1), names(f4), names(f5)))
uf4 <- setdiff(names(f4), c(names(f2), names(f3), names(f1), names(f5)))
uf5 <- setdiff(names(f5), c(names(f2), names(f3), names(f4), names(f1)))

length(f1)
length(f2)
length(f3)
length(f4)
length(f5)
length(uf1)
length(uf2)
length(uf3)
length(uf4)
length(uf5)
write.csv(uf1, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/gene_1_noNA_weight.csv", row.names <- F)
write.csv(uf2, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/gene_2_noNA_weight.csv", row.names <- F)
write.csv(uf3, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/gene_3_noNA_weight.csv", row.names <- F)
write.csv(uf4, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/gene_4_noNA_weight.csv", row.names <- F)
write.csv(uf5, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/gene_5_noNA_weight.csv", row.names <- F)

suptab2 <- rbind.data.frame(cbind.data.frame(as.character(criteria[[k]][2]), "Cluster 1", uf1) %>% purrr::set_names("Grouping method", "Groups", "Genes"),
                 cbind.data.frame(as.character(criteria[[k]][2]), "Cluster 2", uf2) %>% purrr::set_names("Grouping method", "Groups", "Genes"),
                 cbind.data.frame(as.character(criteria[[k]][2]), "Cluster 3", uf3) %>% purrr::set_names("Grouping method", "Groups", "Genes"),
                 cbind.data.frame(as.character(criteria[[k]][2]), "Cluster 4", uf4) %>% purrr::set_names("Grouping method", "Groups", "Genes"),
                 cbind.data.frame(as.character(criteria[[k]][2]), "Cluster 5", uf5) %>% purrr::set_names("Grouping method", "Groups", "Genes"))
saveRDS(functab, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/sigfive_noNA_weight.rds")

table(suptab2$Groups)

#add function to unique genes
#load gff
gff <- readRDS("/Users/ksongsom/OneDrive/postdoc/SV/illumina/gene_anno.rds")
head(gff)

gene_function <- c()
for(i in 1:nrow(suptab2)){
  gene_function[i] <- paste(na.omit(gff[grepl(suptab2$Genes[i], gff$Name),]$rice.defline), collapse = "& ")  
  print(paste(i, i/nrow(suptab2)))
 }

suptab2$`Gene annotationes` <- gene_function

#write.xlsx(suptab2, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/supplement2.xlsx")


#race VS cluster table
head(background)
head(cluster)
sum(cluster$ID %in% background$Taxa)
cluster_back <- merge(cluster, background, by.x = "ID", by.y = "Taxa", all.x = T, all.y = F)
head(cluster_back)
cluster_back[,c("k5_cluster_all", "Race")]
reshape(cluster_back[,c("k5_cluster_all", "Race")], direction = "wide", idvar = "Race", timevar = "k5_cluster_all")
cluster_race_tab <-as.data.frame(cluster_back[,c("k5_cluster_all", "Race")] %>% group_by(k5_cluster_all, Race) %>% tally())
write.csv(cluster_race_tab, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/cluster_race_tab_noNA.csv", row.names <- F)


#make supplement percentage group
sup2_2 <- data_k5_cluster_all[data_k5_cluster_all$eachpvalwithoutNA < 0.05,]
sup2_2 <- sup2_2[!is.na(sup2_2$overthreshold_list),]
dim(sup2_2)
head(sup2_2)
dim(tab200)
colnames(tab200)
dim(svpval)
sup2_2 <- merge(sup2_2, svpval[,-2], by.x = "eachsv", by.y = "eachsv", all.x = T, all.y = F)
sup2_2m <- merge(sup2_2, tab200[, c("svid", "genes", "rice_defline", "go")], by.x = "eachsv", by.y = "svid", all.x = T, all.y = F )
head(sup2_2m)
dim(sup2_2m)

write.xlsx(sup2_2m, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/sup2.2.xlsx")
