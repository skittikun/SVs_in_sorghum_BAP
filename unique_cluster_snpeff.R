library(tidyverse)
library(openxlsx)
library(tidyr)
library(dplyr)

tab19 <- readRDS("/Users/ksongsom/OneDrive/postdoc/SV/illumina/lumpy/347_go19_filtered.rds")
tab19 <- tab19[tab19$svlen >= 50 & tab19$svtype == "DEL",]
dim(tab19)

updatedt <- tab19

untab <- read.xlsx("/Users/ksongsom/OneDrive/postdoc/publication/terra/to_submit_revised/Songsomboon_table_revised1.xlsx", sheet = "Supplement table 3", startRow = 2)


#update new name of deletions for easy tracking
updatedt$newID <- paste(updatedt$chr, updatedt$pos, updatedt$end, updatedt$svtype, sep = "_")

uniqueID <- c()
for(i in 1:nrow(updatedt)){
  eachrow <- updatedt[i,]
  uniqueID <- c(uniqueID, ifelse(eachrow$newID %in% untab$SV.ID, 
                                 paste(eachrow$newID, "unique", 
                                       untab[untab$SV.ID %in% eachrow$newID,]$Unique.to.cluster, sep = "_"), 
                                 paste(eachrow$newID, "general", sep = "_"))
  )
  print(i)
}

updatedt$uniqueID <- uniqueID


updatedt[1,1:10]

#make geno tab
gnt <- updatedt[,(which(colnames(updatedt) == "144134")):(which(colnames(updatedt) == "Rio"))]

gnt <- t(gnt)

library(dplyr)

#col.names = c("chr", "pos", "ID", "ref", "alt", "qual", "filter", "info", "format", baseped[,1]))
snpheader <- read.csv("/Users/ksongsom/OneDrive/postdoc/SV/illumina/TERRA/header_fil.txt", header = F)

library(vcfR)
#vcf <- read.vcfR("/Users/ksongsom/OneDrive/postdoc/SV/illumina/TERRA/fil_del_eff_v3.vcf")
vcf <- readRDS("/Users/ksongsom/OneDrive/postdoc/SV/illumina/TERRA/fil_del_eff_readvcf.vcf.rds")
#vcf_dt <- cbind(as.data.frame(getFIX(vcf)), INFO2df(vcf))
#saveRDS(vcf, "/Users/ksongsom/OneDrive/postdoc/SV/illumina/TERRA/fil_del_eff_readvcf.vcf.rds")
vcf_gt <- as.data.frame(extract.gt(vcf, element = "GT"))
vcf_depth <- as.data.frame(extract.gt(vcf, element = "DP"))
depth_th <- 10
depth_tf <- vcf_depth >= depth_th

test1 <- vcf_gt
vcf_01 <- ifelse(test1 == "0/0" & !is.na(test1), 0, 
                 ifelse((test1 == "1/1" |
                           test1 == "0/1" |
                           test1 == "1/0")  & !is.na(test1), 1, NA))

vcf_01_dp <- ifelse(depth_tf, vcf_01, NA)

headcol <- colnames(vcf_01_dp)
templatecol <- colnames(tab19)[which(colnames(tab19) == ("144134")):which(colnames(tab19) == ("Rio"))]
headcol[!(headcol  %in% templatecol)]

#edit head column
headcol <- gsub("620157","534157",headcol)
headcol <- gsub("655996","655999",headcol)
headcol <- gsub("656020","656021",headcol)
headcol <- gsub("656026","656025",headcol)
headcol <- gsub("656051","656050",headcol)

headcol <- gsub("BTx623genome","BTx623",headcol)
headcol <- gsub("rio","Rio",headcol)

#update colnames
colnames(vcf_01_dp) <- headcol

#remove
#"455307"       "562985"
vcf_01_dp <- vcf_01_dp[,!(colnames(vcf_01_dp) %in% c("455307","562985"))]

#saveRDS(vcf_01_dp, "/Users/ksongsom/OneDrive/postdoc/SV/illumina/TERRA/fil_del_eff_readvcf_01.vcf.rds")

vcf_dt <- vcf_dt[,!(colnames(vcf_dt) %in% c("PRPOS", "PREND", "SNAME"))]
str(vcf_dt )
vcf_dt_full <- cbind.data.frame(vcf_dt, vcf_01_dp)
vcf_dt_full <- vcf_dt_full[vcf_dt_full$SVTYPE == "DEL",]
#saveRDS(vcf_dt_full, "/Users/ksongsom/OneDrive/postdoc/SV/illumina/TERRA/fil_del_eff_readvcf_01_full.vcf.rds")


#only 7766
vcf_dt_full$svid <- paste(vcf_dt_full$CHROM, vcf_dt_full$POS, vcf_dt_full$END, sep = "_")
tab19$svid <- paste(tab19$chr, tab19$pos, tab19$end, sep = "_")
sum(vcf_dt_full$svid %in% tab19$svid) #from all 7766 deletions have 7591 showed in snpeff
dim(vcf_dt_full)

vcf_dt_fill <- vcf_dt_full[vcf_dt_full$svid %in% tab19$svid,]
#saveRDS(vcf_dt_fill, "/Users/ksongsom/OneDrive/postdoc/SV/illumina/TERRA/fil_del_eff_readvcf_01_fill.vcf.rds")

vcf_dt_fill <- readRDS("/Users/ksongsom/OneDrive/postdoc/SV/illumina/TERRA/fil_del_eff_readvcf_01_fill.vcf.rds")

tosplit <- vcf_dt_fill[,c("ANN")]
#tosplit <- gsub("\\|\\|\\|\\|\\|\\|", "\\|", tosplit)
#tosplit <- gsub("\\|\\|\\|\\|\\|", "\\|", tosplit)
#tosplit <- gsub("\\|\\|", "\\|", tosplit)

vcf_dt_fill$ANN <- tosplit

vcf_dt_fill <- vcf_dt_fill %>% separate(ANN, sep = "\\|", 
                                        c("info", "variant_type", "mod1", "gene1", "gene2", "transcript", "mrna", "prot", "detail_sv", 
                                          "eff_num", "eff1",
                                          #"eff_num", "seq", 
                                          "intergenic_region1", "mod2", "gene3", "gene4", "intergenic_region2", "mrna2", "eff2")) %>%
  separate(eff1, sep = "\\,", c("intergenic_region1", "svtype2")) 
head(vcf_dt_fill)

vcf_dt_full[, which(colnames(vcf_dt_full) == "144134"):which(colnames(vcf_dt_full) == "Rio")]
tab19[, which(colnames(tab19) == "144134"):which(colnames(tab19) == "Rio")]

##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' 

vcf_dt_fill <- vcf_dt_fill %>% separate(ANN, sep = "\\|", c("Type", "Functional annotations", "Allele",
                                                            "Annotation","Annotation_Impact","Gene_Name", "Gene_ID", "Feature_Type", "Feature_ID",
                                                            "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos", "cDNA.length",
                                                            "CDS.pos", "CDS.length", "AA.pos", "AA.length", "Distance", "ERRORS", "WARNINGS", "INFO")) 
table(vcf_dt_fill$cDNA.pos)
table(vcf_dt_fill$`Functional annotations`)
head(vcf_dt_fill)

#check known sv Chr09_ 59145680_ 59153091_DEL
dim(vcf_dt_fill)
vcf_dt_fill[vcf_dt_fill$CHROM == "Chr09" & vcf_dt_fill$POS == 59145680,]


#update gene names
#vcf_dt_fill
#genenew <- read.table("/Users/ksongsom/OneDrive/postdoc/reference_Sbicolor_JGI/Sbicolor/Sbicolor/v3.1.1/annotation/Sbicolor_454_v3.1.1.synonym.txt", sep = "\t", header = F)
#colnames(genenew) <- c("Gene_names", "oldgene")
#genenew$Gene_names <- gsub(paste(c("\\.1", "\\.2", "\\.3", "\\.4", "\\.5", "\\.6", "\\.7", "\\.8", "\\.9"), collapse = "|"), "", genenew$Gene_names)
#all(vcf_dt_fill$Annotation %in% genenew$oldgene)

#gsub(paste(genenew$oldgene, collapse = "|"), genenew$Gene_names, vcf_dt_fill$Annotation)

#alleachrow <- c()
#for(i in 1:nrow(vcf_dt_fill)){
#  eachrow <- vcf_dt_fill[i,]$Annotation
#  for(j in 1:nrow(genenew)){
#    eachgene <- genenew[j,]
#    eachrow <- ifelse(grepl(eachgene$oldgene, eachrow), 
#                      gsub(eachgene$oldgene, eachgene$Gene_names, eachrow), eachrow)
#  }
#  alleachrow <- c(alleachrow, eachrow)
#  print(paste(i, i/nrow(vcf_dt_fill)))
#}


#saveRDS(vcf_dt_fill, "/Users/ksongsom/OneDrive/postdoc/SV/illumina/TERRA/fil_del_eff_readvcf_01_fill_split.vcf.rds")

#making crimaric

#which(colnames(tab19) == "144134"):which(colnames(tab19) == "Rio")
filtab19 <- tab19[tab19$svid %in% vcf_dt_fill$svid, ]
all(filtab19$svid %in% vcf_dt_fill$svid)

reorder_dt_fill <- vcf_dt_fill[match(vcf_dt_fill$svid, filtab19$svid),]

left <- reorder_dt_fill[,c(1:which(colnames(reorder_dt_fill) == "NMD"), ncol(reorder_dt_fill))]
right <- filtab19[which(colnames(filtab19) == "144134"):ncol(filtab19)]


update19 <- merge(left, right, by = "svid", all.x = T, all.y = T)
dim(update19)
head(update19)

str(tab19)
head(vcf_dt_fill)
vcf_dt_fill_mini <- vcf_dt_fill[,c("svid", "variant_type", "gene1", "mod1")]
tab19_withsnpeff <- merge(tab19, vcf_dt_fill_mini, by = "svid", all.x = T, all.y = T)
dim(tab19_withsnpeff)

#saveRDS(update19, "/Users/ksongsom/OneDrive/postdoc/SV/illumina/TERRA/fil_del_eff_readvcf_01_fill_split_sveff.vcf.rds")
#saveRDS(tab19_withsnpeff, "/Users/ksongsom/OneDrive/postdoc/SV/illumina/TERRA/fil_del_eff_readvcf_01_fill_split_sveff_all.vcf.rds")

tab19_withsnpeff <- readRDS("/Users/ksongsom/OneDrive/postdoc/SV/illumina/TERRA/fil_del_eff_readvcf_01_fill_split_sveff_all.vcf.rds")
head(tab19_withsnpeff)

#unique to cluster
clusteru <- read.xlsx("/Users/ksongsom/OneDrive/postdoc/publication/terra/to_submit_revised/Songsomboon_table_revised1.xlsx", sheet = "Supplement table 3", startRow = 3)
dim(clusteru )
clusteru$svid <- gsub("\\_DEL", "", clusteru$SV.ID)
all(update19$svid %in% clusteru$svid)
sum(clusteru$svid %in% update19$svid)

sup3 <- merge(clusteru, left[,c("svid", "variant_type", "gene1", "mod1")], by = "svid", all.x = T, all.y = F)

sup3[sup3$svid == "Chr01_46010099_46010163",]

#each type
table(sup3$variant_type)
strsplit(sup3$variant_type, split = "&")
table(unlist(strsplit(sup3$variant_type, split = "&")))
typeofeff <- unique(unlist(strsplit(sup3$variant_type, split = "&")))
typeofeff1 <- typeofeff
typeofeff1[is.na(typeofeff1)] <- "undetermined"
sup3 %>% group_by(Unique.to.cluster, variant_type) %>% summarise(NNN = n())

summarytab <- c()
summarytab_uns <- typeofeff1 
ncluster <- 1:8
for(k in 1:length(ncluster)){
  eachcluster <- sup3[sup3$Unique.to.cluster == k,]
  sumtype <- c()
  for(l in 1:length(typeofeff)){
    sumtype <- c(sumtype, 
                 ifelse(!is.na(typeofeff[l]),
                        sum(grepl(typeofeff[l], eachcluster$variant_type), na.rm = T),
                        sum(is.na(eachcluster$variant_type)))
    )
  }
  
  summarytab <- rbind.data.frame(summarytab, cbind.data.frame(ncluster[k], typeofeff1, sumtype))
  summarytab_uns <- cbind.data.frame(summarytab_uns, sumtype)
}


summarytab_uns <- cbind.data.frame(summarytab_uns, rowSums(summarytab_uns[,-1]))

colnames(summarytab_uns) <- c("Type of mutation", paste("Cluster", 1:8, sep = ""), "Total")

summarytab_uns <- summarytab_uns[order(summarytab_uns$Total, decreasing = T),]

summarytab_uns <- rbind.data.frame(summarytab_uns, c("total", colSums(summarytab_uns[,-1])))

#####include impact
summarytab_imp <- c()
summarytab_uns_imp <- typeofeff1 
eachtype_cluster_all <- c()
ncluster <- 1:8
for(k in 1:length(ncluster)){
  eachcluster <- sup3[sup3$Unique.to.cluster == k,]
  sumtype <- c()
  eachtype_cluster_each <- c()
  for(l in 1:length(typeofeff)){
    eachsum <- ifelse(!is.na(typeofeff[l]),
                      sum(grepl(typeofeff[l], eachcluster$variant_type), na.rm = T),
                      sum(is.na(eachcluster$variant_type)))
    sumtype <- c(sumtype, eachsum)
                 
    if(!is.na(typeofeff[l])){
      eachtype <- eachcluster[grepl(typeofeff[l], eachcluster$variant_type),]  
    }else{
      eachtype <- eachcluster[is.na(eachcluster$variant_type),]
    }
    eachtype_1 <- as.data.frame(eachtype %>% group_by(mod1) %>% summarise(impactongene = n()))
    if(eachsum != 0){
      eachtype_cluster_each <- (rbind.data.frame(eachtype_cluster_each, cbind.data.frame(typeofeff[l], eachtype_1, eachsum )))  
    }
    
  }
  if(is.null(dim(eachtype_cluster_all))){
    eachtype_cluster_all <- eachtype_cluster_each
  }else{
    eachtype_cluster_all <- merge(eachtype_cluster_all, 
                                  eachtype_cluster_each, by = c("typeofeff[l]", "mod1"), 
                                  all.x = T, all.y = T, suffixes = c(k-1, k))
    
  }
  #summarytab_imp <- rbind.data.frame(summarytab_imp, cbind.data.frame(ncluster[k], typeofeff1, sumtype))
  
}

eachtype_cluster_all$`typeofeff[l]`[is.na(eachtype_cluster_all$`typeofeff[l]`)] <- "undetermined"
eachtype_cluster_all$mod1[is.na(eachtype_cluster_all$mod1)] <- "NO IMPACT"
eachtype_cluster_all[is.na(eachtype_cluster_all)] <- 0

eachtype_cluster_all$Total <- rowSums(eachtype_cluster_all[,grepl("impactongene", colnames(eachtype_cluster_all))])

eachtype_cluster_all$mod1 <- factor(eachtype_cluster_all$mod1, levels = c("NO IMPACT", "MODIFIER", "LOW", 
                                                 "MODERATE","HIGH"))

eachtype_cluster_all <- eachtype_cluster_all[order(eachtype_cluster_all$Total, eachtype_cluster_all$mod1, decreasing = T),]

eachtype_cluster_all <- eachtype_cluster_all[,!grepl("eachsum", colnames(eachtype_cluster_all))]

#write.xlsx(eachtype_cluster_all, "/Users/ksongsom/OneDrive/postdoc/publication/terra/snpeff_unique_cluster_sup5.xlsx")
#write.xlsx(summarytab, "/Users/ksongsom/OneDrive/postdoc/publication/terra/snpeff_unique_cluster.xlsx")
#write.xlsx(summarytab_uns, "/Users/ksongsom/OneDrive/postdoc/publication/terra/snpeff_unique_cluster_uns.xlsx")

summary(lm(sumtype ~ `ncluster[k]`*typeofeff, data = summarytab)) #intergenic_region & upstream_gene_variant 

sup3



#write.xlsx(sup3, "/Users/ksongsom/OneDrive/postdoc/publication/terra/sup3.xlsx")
#synonymous and non-synonymous https://pcingola.github.io/SnpEff/se_inputoutput/


nonsyn <- vcf_dt_fill[grepl(paste(c("missense_variant", "initiator_codon_variant", "stop_retained_variant"), collapse = "|"), vcf_dt_fill$variant_type),]
synony <- vcf_dt_fill[grepl(paste(c("synonymous_variant", "start_retained", "stop_retained_variant"), collapse = "|"), vcf_dt_fill$variant_type),]
dim(nonsyn)
dim(synony)

sum(grepl("syno", vcf_dt_fill$variant_type))



#compare between cluster-specific and non-specific

tab19_withsnpeff <- readRDS("/Users/ksongsom/OneDrive/postdoc/SV/illumina/TERRA/fil_del_eff_readvcf_01_fill_split_sveff_all.vcf.rds")

cluspec <- tab19_withsnpeff[grepl(paste(clusteru$svid, collapse = "|"), tab19_withsnpeff$svid),]
nonspec <- tab19_withsnpeff[!grepl(paste(clusteru$svid, collapse = "|"), tab19_withsnpeff$svid),]

effimpact <- unique(nonspec$mod1)

summarytab_non <- c()
eachtype_non_all <- c()
sumtype <- c()
for(l in 1:length(typeofeff)){
  eachsumtype <- ifelse(!is.na(typeofeff[l]),
                        sum(grepl(typeofeff[l], nonspec$variant_type), na.rm = T),
                        sum(is.na(nonspec$variant_type)))
  sumtype <- c(sumtype, eachsumtype)
  if(!is.na(typeofeff[l])){
    eachtype <- nonspec[grepl(typeofeff[l], nonspec$variant_type),]  
  }else{
    eachtype <- nonspec[is.na(nonspec$variant_type),]  
  }
  eachtype_1 <- as.data.frame(eachtype %>% group_by(mod1) %>% summarise(impactongene = n()))
  eachtype_non_all <- rbind.data.frame(eachtype_non_all, cbind.data.frame("nonster", typeofeff[l], eachtype_1, eachsumtype))
}
summarytab_non <- cbind.data.frame("nonster", typeofeff1, sumtype)


eachtype_non_all$mod1 <- factor(eachtype_non_all$mod1, levels = c("MODIFIER", "LOW", 
                                                                  "MODERATE","HIGH"))

eachtype_non_all <- eachtype_non_all[order(eachtype_non_all$eachsumtype,
                                           
                                           eachtype_non_all$mod1, 
                                           
                                           decreasing = T),]

summarytab_clu <- c()
eachtype_clu_all <- c()
sumtype <- c()
for(l in 1:length(typeofeff)){
  eachsumtype <- ifelse(!is.na(typeofeff[l]),
                        sum(grepl(typeofeff[l], cluspec$variant_type), na.rm = T),
                        sum(is.na(cluspec$variant_type)))
  sumtype <- c(sumtype, eachsumtype)
  if(!is.na(typeofeff[l])){
    eachtype <- cluspec[grepl(typeofeff[l], cluspec$variant_type),]  
  }else{
    eachtype <- cluspec[is.na(cluspec$variant_type),]  
  }
  eachtype_1 <- as.data.frame(eachtype %>% group_by(mod1) %>% summarise(impactongene = n()))
  eachtype_clu_all <- rbind.data.frame(eachtype_clu_all, cbind.data.frame("cluster", typeofeff[l], eachtype_1, eachsumtype))
}
summarytab_clu <- cbind.data.frame("cluster", typeofeff1, sumtype)


eachtype_clu_all$mod1 <- factor(eachtype_clu_all$mod1, levels = c("MODIFIER", "LOW", 
                                                                  "MODERATE","HIGH"))

eachtype_clu_all <- eachtype_clu_all[order(eachtype_clu_all$eachsumtype, 
                       eachtype_clu_all$mod1, 
                       decreasing = T),]
summarytab_clu
summarytab_non

eachtype_non_all
eachtype_clu_all
colnames(eachtype_non_all) <- colnames(eachtype_clu_all)
test3 <- rbind.data.frame(eachtype_non_all, eachtype_clu_all)

sup4 <- merge(eachtype_non_all[,-1], eachtype_clu_all[,-1], by = c("typeofeff[l]", "mod1"), suffixes = c("_none", "_cluster"), all.x = T, all.y = T)
sup4 <- sup4[order(sup4$eachsumtype_none, sup4$mod1, decreasing = T),]

#write.xlsx(sup4, "/Users/ksongsom/OneDrive/postdoc/publication/terra/sup4.xlsx")

summary(aov(lm(impactongene ~ `"cluster"`*`typeofeff[l]`+mod1, data = test3)))

#test Chisquare

expect_p <- (summarytab_clu$sumtype + summarytab_non$sumtype)/sum(summarytab_clu$sumtype, summarytab_non$sumtype)
chisq.test(x = summarytab_clu$sumtype, p = expect_p)
chisq.test(x = summarytab_non$sumtype, p = expect_p)

chisq.test(x = summarytab_clu$sumtype, p = rep(nrow(summarytab_clu)))
chisq.test(x = summarytab_non$sumtype, p = expect_p)


colnames(summarytab_non) <- colnames(summarytab_clu)
test3 <- rbind.data.frame(summarytab_clu, summarytab_non)

summary((lm(sumtype ~ `"cluster"`+ typeofeff1 + (`"cluster"`* typeofeff1), data = test3)))
summary((lm(sumtype ~  (`"cluster"`+ typeofeff1), data = test3)))


#known SV
tab19_withsnpeff[tab19_withsnpeff$chr == "Chr09" & tab19_withsnpeff$pos == 59145680,]

##merge table
sup4[is.na(sup4$`typeofeff[l]`),]$`typeofeff[l]` <- "undetermined" 
sup4[is.na(sup4$mod1),]$mod1 <- as.character(sup4[is.na(sup4$mod1),]$mod1)
sup4$mod1 <- as.character(sup4$mod1)
sup4[is.na(sup4$mod1),]$mod1 <- "NO IMPACT"

finaltab <- merge(sup4, eachtype_cluster_all, by = c("typeofeff[l]", "mod1"), all.x = T, all.y = T)
finaltab$mod1 <- factor(finaltab$mod1, levels = c("NO IMPACT", "MODIFIER", "LOW", 
                                                              "MODERATE","HIGH"))

finaltab <- finaltab[order(finaltab$eachsumtype_none,
                       finaltab$mod1, 
                       decreasing = T),]
finaltab
finaltab[is.na(finaltab)] <- 0

eachtype_cluster_all
reshape(eachtype_cluster_all, idvar=`typeofeff[l]`, varying="impactongene1", sep="_", 
        direction = "long", times = "impactongene2")

eachtype_cluster_all_stack <- reshape(eachtype_cluster_all, varying=colnames(eachtype_cluster_all)[grepl("impactongene", colnames(eachtype_cluster_all))],
        direction="long", idvar=c("typeofeff[l]","mod1"),
        v.names="number", timevar="cluster")

summary(aov(lm(number ~ `typeofeff[l]`*mod1*cluster, data = eachtype_cluster_all_stack)))

#write.xlsx(finaltab, "/Users/ksongsom/OneDrive/postdoc/publication/terra/sup_combine.xlsx")


#add gene functions
suptab <- read.xlsx("/Users/ksongsom/OneDrive/postdoc/publication/terra/to_submit_revised/Songsomboon_table_revised1.xlsx", sheet = "Supplement table 3", startRow = 2)
head(suptab)
dim(suptab)
gff <- readRDS("/Users/ksongsom/OneDrive/postdoc/SV/illumina/gene_anno.rds")
head(gff)

mrna <- gff[gff$feature == "mRNA", ]
shortname <- gsub(paste(paste("\\.", 1:30, sep = ""), collapse = "|"), "", mrna$Name)
mrna <- mrna[!duplicated(shortname),]

gene_function_all <- c()
for(m in 1:nrow(suptab)){
#m = 55
affected <- suptab[m,]$Affected.genes
if(!grepl(paste(c("-", "&"), collapse = "|"), affected)){
  
  gene_function <- mrna[grepl(affected, mrna$Name),]$rice.defline
  gene_function <- gsub(", putative", "", gsub(", expressed", "", gene_function))
  #gene_function_all <- c(gene_function_all, gene_function)
}else{
  if(grepl(paste(c("-"), collapse = "|"), affected)){
    
    aff1 <- strsplit(affected, "-")[[1]][1]
    aff2 <- strsplit(affected, "-")[[1]][2]
    allaffected <- c()
    startrow <- which(grepl(aff1, gff$Name))
    endrow <- which(grepl(aff2, mrna$Name))
    gene_function <- paste(na.omit(mrna[startrow:endrow,]$rice.defline), collapse = "& ")
    gene_function <- gsub(", putative", "", gsub(", expressed", "", gene_function))
    #gene_function_all <- c(gene_function_all, gene_function)
  }else{
    if(grepl(paste(c("&"), collapse = "|"), affected)){
      affected
      allaffected <- unlist(strsplit(affected, "&"))
      gene_function <- paste(na.omit(mrna[grepl(paste(allaffected, collapse = "|"), mrna$Name),]$rice.defline), collapse = "& ")
      gene_function <- gsub(", putative", "", gsub(", expressed", "", gene_function))
      #gene_function_all <- c(gene_function_all, gene_function)
    }
  }
}
gene_function_all[m] <- gene_function
print(paste(m, m/nrow(suptab), sep = "_"))
}
head(gene_function_all)
length(gene_function_all)
suptab2 <- cbind.data.frame(suptab, gene_function_all)
head(suptab2)

#write.xlsx(suptab2, "/Users/ksongsom/OneDrive/postdoc/publication/terra/sup3_update.xlsx")

#gene_function_all <- gene_function_all[!is.na(gene_function_all)]
