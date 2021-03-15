library(stats)
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
library(ggbiplot)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra)
library(readxl)
library(openxlsx)

options(stringsAsFactors = F)

fildt <- readRDS("/Users/ksongsom/OneDrive/postdoc/publication/terra/code/files/347_go.rds")
fildt <- fildt[fildt$svlen>= 50,]

fildt <- readRDS("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/AllChr_snp012_3.vcf.012_imputed.rds")
fildt <- readRDS("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/AllChr_snp012_2.vcf.012_imputed.rds")
dim(fildt)

#load map
map <- (read_excel("/Users/ksongsom/OneDrive/postdoc/NAM/Panel/BAP/update_sorghum_location_subregions.xlsx", sheet = 1)) #from google sheet
#update elevation
map$ID <- gsub("PI ", "", map$Genotype)
head(map)
#map <- map %>% separate(Coordinate, c("lat", "lon"), sep = ", ")
map$ID <- gsub("Rio", "rio", map$ID)


#head(fildt)
#dim(fildt)
mafcut <- 0.05
missingcut <- 0.2 #0.05
sizecut <- 100000
#table(fildt$svtype)
#fildt <- fildt[fildt$maf_col >= mafcut & fildt$missing_col <= missingcut, ]
#fildt$end <- ifelse(fildt$svtype == "BND", fildt$pos, fildt$end)
#fildt$svlen <- fildt$end -fildt$pos
#colnames(fildt) <- gsub("X", "", colnames(fildt))
#fildtinfo <- fildt[, (which(colnames(fildt) == "format")+1):which(colnames(fildt) == "Rio")]
#tfildtinfo <- t(fildtinfo)
#rownames(tfildtinfo) <- gsub("X", "", rownames(tfildtinfo))
#colnames(tfildtinfo) <- fildt$id
df <- fildt
#df <- df[!(rownames(df) %in% names(which(rowSums(df, na.rm = T) == 0))),]

      df <- na.omit(df)
      #df <- scale(df)
 
      #define number of k
      set.seed(123)
      
      # function to compute total within-cluster sum of square 
      wss <- function(k) {
        kmeans(df, k, nstart = 10)$tot.withinss
      }
      # Compute and plot wss for k = 1 to k = 15
      k.values <- 1:15
      # extract wss for 2-15 clusters
      wss_values <- map_dbl(k.values, wss)
      
      plot(k.values, wss_values,
           type="b", pch = 19, frame = FALSE, 
           xlab="Number of clusters K",
           ylab="Total within-clusters sum of squares")
      
      #factoextra::fviz_nbclust(df, kmeans, method = "wss")
      factoextra::fviz_nbclust(df, kmeans, method = "silhouette") #the optimal number of clusters is 5
      
      
      k5 <- stats::kmeans(df, centers = 5, nstart = 25)
      k8 <- stats::kmeans(df, centers = 8, nstart = 25)
      #save clusters
      #write.csv(k5, "/Users/ksongsom/OneDrive/postdoc/publication/terra/code/files/five_cluster.csv")
           
      # plots to compare
      p4 <- factoextra::fviz_cluster(k5, geom = c("point", "text"),  data = df, repel = TRUE, labelsize = 8) + ggtitle("k = 5")
      p4
      head(p4)
      head(p4$data)
      
      p8 <- factoextra::fviz_cluster(k9, geom = c("point", "text"),  data = df, repel = TRUE, labelsize = 8) + ggtitle("k = 9")
      p8
      
      #different PC
      #load cluster
      precluster <- read.csv("/Users/ksongsom/OneDrive/postdoc/publication/terra/code/files/five_cluster.csv")
      k5cluster <- cbind.data.frame(names(k5$cluster), k5$cluster)
      colnames(k5cluster) <- c("PI", "k5cluster")
      k8cluster <- cbind.data.frame(names(k8$cluster), k8$cluster)
      colnames(k8cluster) <- c("PI", "k8cluster")
      
      dim(precluster)
      dim(k5cluster)
      dim(k8cluster)
      
      clusterdat <- merge(precluster[,c("PI", "k5_cluster_all")], y = k5cluster, by = "PI")
      clusterdat <- merge(clusterdat, y = k8cluster, by = "PI")
      head(clusterdat)
      
      clusterdat[clusterdat$k5_cluster_all != clusterdat$k5cluster,]
      dim(clusterdat[clusterdat$k5_cluster_all != clusterdat$k5cluster,])
      table(clusterdat$k8cluster)
      table(clusterdat$k5cluster)
      #head(cluster)
      #dim(cluster)
      cluster <- k5$cluster
      cluster <- cbind.data.frame(names(cluster), cluster)
      colnames(cluster) <- c("PI", "cluster")
      #run PCA  5 clusters
            pca_res <- PCA(df, graph = F)
            PCdat <- as.data.frame(pca_res$ind$coord[,1:5])
            
            pca_res <- PCA(df, graph = F, ncp = 8)
            PCdat <- as.data.frame(pca_res$ind$coord[,1:8])
            
            PCdat$id <- rownames(PCdat)
            PCdat <- merge(x = PCdat, y = cluster[,c("PI", "cluster")], by.x = "id", by.y = "PI", all.x = T, all.y = F)
            PCdat <- merge(x = PCdat, y = map[,c("ID", "subregions", "Country")], by.x = "id", by.y = "ID", all.x = T, all.y = F)
            nlevels(as.factor(PCdat$Country))
            nlevels(as.factor(PCdat$subregions))
            PCdat$clusters <- as.character(PCdat$cluster)
            
            PCdat$subregions <- factor(PCdat$subregions, levels = 
                                         c("Asia", "Central Africa",
                                           "Eastern Africa", "Northern Africa",
                                           "Southern Africa", "Western Africa",
                                           "Australia", "Europe", "North America"
                                         ))
            
            #making the text for %variance for axis labels
            eigen <- as.data.frame(pca_res$eig)
            
            pc1 <- round(eigen[1,]$`percentage of variance`, digits = 1)
            pc2 <- round(eigen[2,]$`percentage of variance`, digits = 1)
            pc3 <- round(eigen[3,]$`percentage of variance`, digits = 1)
            pc4 <- round(eigen[4,]$`percentage of variance`, digits = 1)
            pc5 <- round(eigen[5,]$`percentage of variance`, digits = 1)
            pc6 <- round(eigen[6,]$`percentage of variance`, digits = 1)
            pc7 <- round(eigen[7,]$`percentage of variance`, digits = 1)
            pc8 <- round(eigen[8,]$`percentage of variance`, digits = 1)
            
            pc1lab = paste("PC1 (", pc1,"%)", sep = "" )
            pc1lab
            pc2lab = paste("PC2 (", pc2,"%)", sep = "")
            pc2lab
            pc3lab = paste("PC3 (", pc3,"%)", sep = "")
            pc3lab
            pc4lab = paste("PC4 (", pc4,"%)", sep = "")
            pc4lab
            pc5lab = paste("PC5 (", pc5,"%)", sep = "")
            pc5lab
            
            pc8lab = paste("PC8 (", pc8,"%)", sep = "")
            pc8lab
            
            
            #define the circle line
            find_hull12 <- function(PCdat) PCdat[chull(PCdat[,2], PCdat[,3]), ]
            hulls12 <- ddply(PCdat, "clusters", find_hull12)
            
            
            find_hull15 <- function(PCdat) PCdat[chull(PCdat[,2], PCdat[,6]), ]
            hulls15 <- ddply(PCdat, "clusters", find_hull15)
            
            find_hull18 <- function(PCdat) PCdat[chull(PCdat[,2], PCdat[,9]), ]
            hulls18 <- ddply(PCdat, "clusters", find_hull18)
            
            
            find_hull35 <- function(PCdat) PCdat[chull(PCdat[,4], PCdat[,6]), ]
            hulls35 <- ddply(PCdat, "clusters", find_hull35)
            
            plot12 <- ggplot(data = PCdat, aes(x = Dim.1, y = -1*Dim.2, col=clusters, 
                                               fill=clusters, group=clusters)) +
              geom_point(aes(shape=subregions), size = 3) + labs(x = pc1lab, y = pc2lab) +
              scale_color_manual(breaks = c("1", "2", "3", "4", "5"),
                                values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              scale_fill_manual(breaks = c("1", "2", "3", "4", "5"),
                                 values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              geom_polygon(data = hulls12, alpha = 0.3) + theme_bw() +
              theme(text = element_text(size=20)) +
              scale_shape_manual(values=1:nlevels(as.factor(PCdat$subregions)))
            
            #c("gold", "purple", "green", "red", "blue")
            #c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink")
            plot12
            
            plot15 <- ggplot(data = PCdat, aes(x = Dim.1, y = -1*Dim.5, col=clusters, 
                                               fill=clusters, group=clusters)) +
              geom_point(aes(shape=subregions), size = 3) + labs(x = pc1lab, y = pc5lab) +
              scale_color_manual(breaks = c("1", "2", "3", "4", "5"),
                                 values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink")
                                )  + 
              scale_fill_manual(breaks = c("1", "2", "3", "4", "5"),
                                values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink")
                                )  + 
              geom_polygon(data = hulls15, alpha = 0.3) + theme_bw() +
              theme(text = element_text(size=20)) +
              scale_shape_manual(values=1:nlevels(as.factor(PCdat$subregions)))
            
            plot15
            
            plot18 <- ggplot(data = PCdat, aes(x = Dim.1, y = -1*Dim.8, col=clusters, 
                                               fill=clusters, group=clusters)) +
              geom_point(aes(shape=subregions), size = 3) + labs(x = pc1lab, y = pc8lab) +
              scale_color_manual(breaks = c("1", "2", "3", "4", "5"),
                                 values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink")
              )  + 
              scale_fill_manual(breaks = c("1", "2", "3", "4", "5"),
                                values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink")
              )  + 
              geom_polygon(data = hulls18, alpha = 0.3) + theme_bw() +
              theme(text = element_text(size=20)) +
              scale_shape_manual(values=1:nlevels(as.factor(PCdat$subregions)))
            
            plot18
            
            
            plot35 <- ggplot(data = PCdat, aes(x = Dim.3, y = -1*Dim.5, col=clusters, 
                                               fill=clusters, group=clusters)) +
              geom_point(aes(shape=subregions), size = 3) + labs(x = pc3lab, y = pc5lab) +
              scale_color_manual(breaks = c("1", "2", "3", "4", "5"),
                                 values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink")
              )  + 
              scale_fill_manual(breaks = c("1", "2", "3", "4", "5"),
                                values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink")
              )  + 
              geom_polygon(data = hulls35, alpha = 0.3) + theme_bw() +
              theme(text = element_text(size=20)) +
              scale_shape_manual(values=1:nlevels(as.factor(PCdat$subregions)))
            
            plot35
            
            #test combination 5 clusters
            
            find_hull28 <- function(PCdat)PCdat[chull(PCdat[,3], PCdat[,9]), ]
            hulls28 <- ddply(PCdat, "clusters", find_hull28)
            plot28 <- ggplot(data = PCdat, aes(x = Dim.2, y = -1*Dim.8, col=clusters, 
                                               fill=clusters, group=clusters)) +
              geom_point(aes(shape=subregions), size = 3) + labs(x = pc2lab, y = pc8lab) +
              scale_color_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                 values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              scale_fill_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              geom_polygon(data = hulls28, alpha = 0.3) + theme_bw() +
              theme(text = element_text(size=20)) +
              scale_shape_manual(values=1:nlevels(as.factor(PCdat$subregions)))
            
            plot28
            
            
            find_hull38 <- function(PCdat)PCdat[chull(PCdat[,4], PCdat[,9]), ]
            hulls38 <- ddply(PCdat, "clusters", find_hull38)
            plot38 <- ggplot(data = PCdat, aes(x = Dim.3, y = -1*Dim.8, col=clusters, 
                                               fill=clusters, group=clusters)) +
              geom_point(aes(shape=subregions), size = 3) + labs(x = pc3lab, y = pc8lab) +
              scale_color_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                 values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              scale_fill_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              geom_polygon(data = hulls38, alpha = 0.3) + theme_bw() +
              theme(text = element_text(size=20)) +
              scale_shape_manual(values=1:nlevels(as.factor(PCdat$subregions)))
            
            plot38
            
            #save tiff
            tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/pc12.tiff", width=10, height=10, units="in", res=300)
            plot12
            dev.off()
            
            tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/pc15.tiff", width=10, height=10, units="in", res=300)
            plot15
            dev.off()
            
            tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/silhouette.tiff", width=10, height=10, units="in", res=300)
            factoextra::fviz_nbclust(df, kmeans, method = "silhouette") +
              theme(text = element_text(size=20))
            dev.off()
            
            find_hull13 <- function(PCdat) PCdat[chull(PCdat[,2], PCdat[,4]), ]
            hulls13 <- ddply(PCdat, "clusters", find_hull13)
            plot13 <- ggplot(data = PCdat, aes(x = Dim.1, y = -1*Dim.3, col=clusters, 
                                               fill=clusters, group=clusters)) +
              geom_point(aes(shape=clusters), size = 3) + labs(x = pc1lab, y = pc3lab) +
              scale_color_manual(breaks = c("1", "2", "3", "4", "5"),
                                 values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink")
              )  + 
              scale_fill_manual(breaks = c("1", "2", "3", "4", "5"),
                                values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink")
              )  + 
              geom_polygon(data = hulls13, alpha = 0.3) + theme_bw() +
              theme(text = element_text(size=20))
            
            plot13
            
            find_hull14 <- function(PCdat) PCdat[chull(PCdat[,2], PCdat[,5]), ]
            hulls14 <- ddply(PCdat, "clusters", find_hull14)
            plot14 <- ggplot(data = PCdat, aes(x = Dim.1, y = -1*Dim.4, col=clusters, 
                                               fill=clusters, group=clusters)) +
              geom_point(aes(shape=clusters), size = 3) + labs(x = pc1lab, y = pc4lab) +
              scale_color_manual(breaks = c("1", "2", "3", "4", "5"),
                                 values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink")
              )  + 
              scale_fill_manual(breaks = c("1", "2", "3", "4", "5"),
                                values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink")
              )  + 
              geom_polygon(data = hulls14, alpha = 0.3) + theme_bw() +
              theme(text = element_text(size=20))
            
            plot14
            
            tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/pc13.tiff", width=10, height=10, units="in", res=300)
            plot13
            dev.off()
            tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/pc14.tiff", width=10, height=10, units="in", res=300)
            plot14
            dev.off()
            
            
            
      #run PCA  8 clusters
            
            cluster <- k8$cluster
            cluster <- cbind.data.frame(names(cluster), cluster)
            colnames(cluster) <- c("PI", "cluster")
            
            pca_res <- PCA(df, graph = F, ncp = 8)
            PCdat <- as.data.frame(pca_res$ind$coord[,1:8])
            PCdat$id <- rownames(PCdat)
            PCdat <- merge(x = PCdat, y = cluster[,c("PI", "cluster")], by.x = "id", by.y = "PI", all.x = T, all.y = F)
            PCdat <- merge(x = PCdat, y = map[,c("ID", "subregions", "Country")], by.x = "id", by.y = "ID", all.x = T, all.y = F)
             
            nlevels(as.factor(PCdat$Country))
            nlevels(as.factor(PCdat$subregions))
            PCdat$clusters <- as.character(PCdat$cluster)
            
            PCdat$subregions <- factor(PCdat$subregions, levels = 
                                           c("Asia", "Central Africa",
                                             "Eastern Africa", "Northern Africa",
                                             "Southern Africa", "Western Africa",
                                             "Australia", "Europe", "North America"
                                             ))
            
            #making the text for %variance for axis labels
            eigen <- as.data.frame(pca_res$eig)
            
            pc1 <- round(eigen[1,]$`percentage of variance`, digits = 1)
            pc2 <- round(eigen[2,]$`percentage of variance`, digits = 1)
            pc3 <- round(eigen[3,]$`percentage of variance`, digits = 1)
            pc4 <- round(eigen[4,]$`percentage of variance`, digits = 1)
            pc5 <- round(eigen[5,]$`percentage of variance`, digits = 1)
            pc6 <- round(eigen[6,]$`percentage of variance`, digits = 1)
            pc7 <- round(eigen[7,]$`percentage of variance`, digits = 1)
            pc8 <- round(eigen[8,]$`percentage of variance`, digits = 1)
            
            pc1lab = paste("PC1 (", pc1,"%)", sep = "" )
            pc1lab
            pc2lab = paste("PC2 (", pc2,"%)", sep = "")
            pc2lab
            pc3lab = paste("PC3 (", pc3,"%)", sep = "")
            pc3lab
            pc4lab = paste("PC4 (", pc4,"%)", sep = "")
            pc4lab
            pc5lab = paste("PC5 (", pc5,"%)", sep = "")
            pc5lab
            pc6lab = paste("PC6 (", pc6,"%)", sep = "")
            pc6lab
            pc7lab = paste("PC7 (", pc7,"%)", sep = "")
            pc7lab
            pc8lab = paste("PC8 (", pc8,"%)", sep = "")
            pc8lab
            
            
            #define the circle line
            find_hull12 <- function(PCdat) PCdat[chull(PCdat[,2], PCdat[,3]), ]
            hulls12 <- ddply(PCdat, "clusters", find_hull12)
            
            
            find_hull15 <- function(PCdat) PCdat[chull(PCdat[,2], PCdat[,6]), ]
            hulls15 <- ddply(PCdat, "clusters", find_hull15)
            
            find_hull18 <- function(PCdat)PCdat[chull(PCdat[,2], PCdat[,9]), ]
            hulls18 <- ddply(PCdat, "clusters", find_hull18)
            
            plot12 <- ggplot(data = PCdat, aes(x = Dim.1, y = -1*Dim.2, col=clusters, 
                                               fill=clusters, group=clusters)) +
              geom_point(aes(shape=subregions), size = 3) + labs(x = pc1lab, y = pc2lab) +
              scale_color_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                 values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              scale_fill_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              geom_polygon(data = hulls12, alpha = 0.3) + theme_bw() +
              theme(text = element_text(size=20)) +
              scale_shape_manual(values=1:nlevels(as.factor(PCdat$subregions)))
            
            plot12
            
            plot15 <- ggplot(data = PCdat, aes(x = Dim.1, y = -1*Dim.5, col=clusters, 
                                               fill=clusters, group=clusters)) +
              geom_point(aes(shape=subregions), size = 3) + labs(x = pc1lab, y = pc5lab) +
              scale_color_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                 values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              scale_fill_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              geom_polygon(data = hulls15, alpha = 0.3) + theme_bw() +
              theme(text = element_text(size=20)) +
              scale_shape_manual(values=1:nlevels(as.factor(PCdat$subregions)))
            
            plot15
            
            plot18 <- ggplot(data = PCdat, aes(x = Dim.1, y = -1*Dim.8, col=clusters, 
                                               fill=clusters, group=clusters)) +
              geom_point(aes(shape=subregions), size = 3) + labs(x = pc1lab, y = pc8lab) +
              scale_color_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                 values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              scale_fill_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              geom_polygon(data = hulls18, alpha = 0.3) + theme_bw() +
              theme(text = element_text(size=20)) +
              scale_shape_manual(values=1:nlevels(as.factor(PCdat$subregions)))
            
            plot18
            
            
        #test combination 8 cluster
            
            find_hull28 <- function(PCdat)PCdat[chull(PCdat[,3], PCdat[,9]), ]
            hulls28 <- ddply(PCdat, "clusters", find_hull28)
            plot28 <- ggplot(data = PCdat, aes(x = Dim.2, y = -1*Dim.8, col=clusters, 
                                               fill=clusters, group=clusters)) +
              geom_point(aes(shape=subregions), size = 3) + labs(x = pc2lab, y = pc8lab) +
              scale_color_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                 values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              scale_fill_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              geom_polygon(data = hulls28, alpha = 0.3) + theme_bw() +
              theme(text = element_text(size=20)) +
              scale_shape_manual(values=1:nlevels(as.factor(PCdat$subregions)))
            
            plot28
            
            
            find_hull38 <- function(PCdat)PCdat[chull(PCdat[,4], PCdat[,9]), ]
            hulls38 <- ddply(PCdat, "clusters", find_hull38)
            plot38 <- ggplot(data = PCdat, aes(x = Dim.3, y = -1*Dim.8, col=clusters, 
                                               fill=clusters, group=clusters)) +
              geom_point(aes(shape=subregions), size = 3) + labs(x = pc3lab, y = pc8lab) +
              scale_color_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                 values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              scale_fill_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              geom_polygon(data = hulls38, alpha = 0.3) + theme_bw() +
              theme(text = element_text(size=20)) +
              scale_shape_manual(values=1:nlevels(as.factor(PCdat$subregions)))
            
            plot38
            
            
          #race*subregions
            cluster <- k8$cluster
            cluster <- cbind.data.frame(names(cluster), cluster)
            colnames(cluster) <- c("PI", "cluster")
            
            pca_res <- PCA(df, graph = F, ncp = 8)
            PCdat <- as.data.frame(pca_res$ind$coord[,1:8])
            PCdat$id <- rownames(PCdat)
            
            map$Race[(map$Race == "NA")] <- "Unknown"
            map$Race[(map$Race == "Mixed")] <- "Unknown"
            map$subregions_torace <- ifelse(grepl("Africa", map$subregions), map$subregions, "outside Africa")
            map$racereg <- paste(map$Race, map$subregions_torace, sep = "_")
            PCdat <- merge(x = PCdat, y = cluster[,c("PI", "cluster")], by.x = "id", by.y = "PI", all.x = T, all.y = F)
            
            PCdat <- merge(x = PCdat, y = map[,c("ID", "subregions", "Country", "racereg")], by.x = "id", by.y = "ID", all.x = T, all.y = F)
            
            
            nlevels(as.factor(PCdat$Country))
            nlevels(as.factor(PCdat$subregions))
            PCdat$clusters <- as.character(PCdat$cluster)
            
            PCdat$subregions <- factor(PCdat$subregions, levels = 
                                         c("Asia", "Central Africa",
                                           "Eastern Africa", "Northern Africa",
                                           "Southern Africa", "Western Africa",
                                           "Australia", "Europe", "North America"
                                         ))
            
            #making the text for %variance for axis labels
            eigen <- as.data.frame(pca_res$eig)
            
            pc1 <- round(eigen[1,]$`percentage of variance`, digits = 1)
            pc2 <- round(eigen[2,]$`percentage of variance`, digits = 1)
            pc3 <- round(eigen[3,]$`percentage of variance`, digits = 1)
            pc4 <- round(eigen[4,]$`percentage of variance`, digits = 1)
            pc5 <- round(eigen[5,]$`percentage of variance`, digits = 1)
            pc6 <- round(eigen[6,]$`percentage of variance`, digits = 1)
            pc7 <- round(eigen[7,]$`percentage of variance`, digits = 1)
            pc8 <- round(eigen[8,]$`percentage of variance`, digits = 1)
            
            pc1lab = paste("PC1 (", pc1,"%)", sep = "" )
            pc1lab
            pc2lab = paste("PC2 (", pc2,"%)", sep = "")
            pc2lab
            pc3lab = paste("PC3 (", pc3,"%)", sep = "")
            pc3lab
            pc4lab = paste("PC4 (", pc4,"%)", sep = "")
            pc4lab
            pc5lab = paste("PC5 (", pc5,"%)", sep = "")
            pc5lab
            pc6lab = paste("PC6 (", pc6,"%)", sep = "")
            pc6lab
            pc7lab = paste("PC7 (", pc7,"%)", sep = "")
            pc7lab
            pc8lab = paste("PC8 (", pc8,"%)", sep = "")
            pc8lab
            
            
            #define the circle line
            find_hull12 <- function(PCdat) PCdat[chull(PCdat[,2], PCdat[,3]), ]
            hulls12 <- ddply(PCdat, "clusters", find_hull12)
            
            selectedsymbol_n <- c(0:18, 33:38, 42, 43, 47, 60:72, 74:78, 80:90)
            
            table(PCdat$racereg)[order(table(PCdat$racereg), decreasing = T)]
            raceregorder <- names(table(PCdat$racereg)[order(table(PCdat$racereg), decreasing = T)])
            
            selectedsymbol_n1 <- cbind.data.frame(raceregorder, selectedsymbol_n[1:length(raceregorder)])
            
            selectedsymbol <- selectedsymbol_n1$`selectedsymbol_n[1:length(raceregorder)]`[order(selectedsymbol_n1$raceregorder)]
            
            PCdat$racereg <- PCdat$racereg#factor(PCdat$racereg, levels = raceregorder)
            
            plot12 <- ggplot(data = PCdat, aes(x = Dim.1, y = -1*Dim.2, col=clusters, 
                                               fill=clusters, group=clusters)) +
              geom_point(aes(shape=racereg), size = 3) + labs(x = pc1lab, y = pc2lab) +
              scale_color_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                 values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              scale_fill_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              geom_polygon(data = hulls12, alpha = 0.3) + theme_bw() +
              theme(text = element_text(size=10),
                    legend.title = element_text(color = "black", size = 10),
                    legend.text = element_text(color = "black", size = 8)) +
              scale_shape_manual(name = "Races and regions", values=selectedsymbol[1:nlevels(as.factor(PCdat$racereg))])
            
            
            plot12
            
            find_hull15 <- function(PCdat) PCdat[chull(PCdat[,2], PCdat[,6]), ]
            hulls15 <- ddply(PCdat, "clusters", find_hull15)
            
            plot15 <- ggplot(data = PCdat, aes(x = Dim.1, y = -1*Dim.5, col=clusters, 
                                               fill=clusters, group=clusters)) +
              geom_point(aes(shape=racereg), size = 3) + labs(x = pc1lab, y = pc5lab) +
              scale_color_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                 values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              scale_fill_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              geom_polygon(data = hulls15, alpha = 0.3) + theme_bw() +
              theme(text = element_text(size=10),
                    legend.title = element_text(color = "black", size = 10),
                    legend.text = element_text(color = "black", size = 8)) +
              scale_shape_manual(name = "Races and regions", values=selectedsymbol[1:nlevels(as.factor(PCdat$racereg))])
            
            plot15
            
            find_hull18 <- function(PCdat) PCdat[chull(PCdat[,2], PCdat[,9]), ]
            hulls18 <- ddply(PCdat, "clusters", find_hull18)
            
            plot18 <- ggplot(data = PCdat, aes(x = Dim.1, y = -1*Dim.8, col=clusters, 
                                               fill=clusters, group=clusters)) +
              geom_point(aes(shape=racereg), size = 3) + labs(x = pc1lab, y = pc8lab) +
              scale_color_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                 values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              scale_fill_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              geom_polygon(data = hulls18, alpha = 0.3) + theme_bw() +
              theme(text = element_text(size=10),
                    legend.title = element_text(color = "black", size = 10),
                    legend.text = element_text(color = "black", size = 8)) +
              scale_shape_manual(name = "Races and regions", values=selectedsymbol[1:nlevels(as.factor(PCdat$racereg))])
            
            plot18
            
            find_hull28 <- function(PCdat) PCdat[chull(PCdat[,3], PCdat[,9]), ]
            hulls28 <- ddply(PCdat, "clusters", find_hull28)
            
            plot28 <- ggplot(data = PCdat, aes(x = Dim.2, y = -1*Dim.8, col=clusters, 
                                               fill=clusters, group=clusters)) +
              geom_point(aes(shape=racereg), size = 3) + labs(x = pc2lab, y = pc8lab) +
              scale_color_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                 values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              scale_fill_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              geom_polygon(data = hulls28, alpha = 0.3) + theme_bw() +
              theme(text = element_text(size=10),
                    legend.title = element_text(color = "black", size = 10),
                    legend.text = element_text(color = "black", size = 8)) +
              scale_shape_manual(name = "Races and regions", values=selectedsymbol[1:nlevels(as.factor(PCdat$racereg))])
            
            plot28
            
            find_hull38 <- function(PCdat) PCdat[chull(PCdat[,4], PCdat[,9]), ]
            hulls38 <- ddply(PCdat, "clusters", find_hull38)
            
            plot38 <- ggplot(data = PCdat, aes(x = Dim.3, y = -1*Dim.8, col=clusters, 
                                               fill=clusters, group=clusters)) +
              geom_point(aes(shape=racereg), size = 3) + labs(x = pc3lab, y = pc8lab) +
              scale_color_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                 values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              scale_fill_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
                                values=c("red", "deepskyblue1","dark green","green","purple","orange","gold","pink"))  + 
              geom_polygon(data = hulls38, alpha = 0.3) + theme_bw() +
              theme(text = element_text(size=10),
                    legend.title = element_text(color = "black", size = 10),
                    legend.text = element_text(color = "black", size = 8)) +
              scale_shape_manual(name = "Races and regions", values=selectedsymbol[1:nlevels(as.factor(PCdat$racereg))])
            
            plot38
            
            #save tiff
            tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/pc12.tiff", width=10, height=10, units="in", res=300)
            plot12
            dev.off()
            
            tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/pc15.tiff", width=10, height=10, units="in", res=300)
            plot15
            dev.off()
            
            tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/pc18.tiff", width=10, height=10, units="in", res=300)
            plot18
            dev.off()
            
            
      #summary table of race and cluster from 8 clusters
            head(PCdat)
            
            allraces <- unique(gsub("\\_.*", "", PCdat$racereg))
            allraces <- allraces[order(allraces)]

            allraces_dat <- as.data.frame(allraces)
            for(i in 1:length(unique(PCdat$clusters)[order(unique(PCdat$clusters))])){
              eachcluster <-PCdat[PCdat$clusters == unique(PCdat$clusters)[order(unique(PCdat$clusters))][i],]
              sumeachcluster <- eachcluster %>% group_by(racereg) %>% summarise(numindi = n())
              sumeachcluster$race <- gsub("\\_.*", "", sumeachcluster$racereg)
              sumeachcluster_more <- sumeachcluster %>% group_by(race) %>% summarise(sumnumindi = sum(numindi, na.rm = T))
              colnames(sumeachcluster) <- c("racereg", paste("Cluster ", i, sep = ""), "race")
              colnames(sumeachcluster_more) <- c("race", paste("Cluster ", i, sep = ""))
              allraces_dat <- merge(allraces_dat, sumeachcluster_more, by.x = "allraces", by.y = "race", all.x = T, all.y = T)
              }
            
            allraces_dat
            Total <- c(rowSums(allraces_dat[,-1], na.rm = T))#, sum(rowSums(allraces_dat[,-1], na.rm = T)))
            allraces_dat <- cbind.data.frame(allraces_dat, Total)
            allraces_dat_tosave <- rbind.data.frame(allraces_dat, c("total", colSums(allraces_dat[,-1], na.rm = T)))
            
            head(allraces_dat_tosave)
write.xlsx(allraces_dat_tosave, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/cluster_number_race.xlsx")            
            
            #origin country

allbico <- c("156018", "176766", "197542", "213900", "226096", "22913", "521281", "562717", "562730", "563032")
sum(PCdat$id %in% allbico)
sum(allbico %in% map$ID)            
intersect(intersect(PCdat$id, allbico), map$ID)

head(PCdat)


        allcountries <- unique(gsub("\\_.*", "", PCdat$Country))
        allcountries <- allcountries[order(allcountries)]
        
        allcountries_dat <- as.data.frame(allcountries)
        for(i in 1:length(unique(PCdat$clusters)[order(unique(PCdat$clusters))])){
          eachcluster <-PCdat[PCdat$clusters == unique(PCdat$clusters)[order(unique(PCdat$clusters))][i],]
          sumeachcluster <- eachcluster %>% group_by(Country) %>% summarise(numindi = n())
          #sumeachcluster$race <- gsub("\\_.*", "", sumeachcluster$Country)
          sumeachcluster_more <- sumeachcluster #%>% group_by(Country) %>% summarise(sumnumindi = sum(numindi, na.rm = T))
          #colnames(sumeachcluster) <- c("Country", paste("Cluster ", i, sep = ""))
          colnames(sumeachcluster_more) <- c("Country", paste("Cluster ", i, sep = ""))
          allcountries_dat <- merge(allcountries_dat, sumeachcluster_more, by.x = "allcountries", by.y = "Country", all.x = T, all.y = T)
        }
        
        allcountries_dat
        Total <- c(rowSums(allcountries_dat[,-1], na.rm = T))#, sum(rowSums(allcountries_dat[,-1], na.rm = T)))
        allcountries_dat <- cbind.data.frame(allcountries_dat, Total)
        allcountries_dat <- allcountries_dat[order(allcountries_dat$Total, decreasing = T),]
        allcountries_dat_tosave <- rbind.data.frame(allcountries_dat, c("total", colSums(allcountries_dat[,-1], na.rm = T)))
        
        colnames(allcountries_dat_tosave)[1] <- "Countries of origin"
        
        head(allcountries_dat_tosave)
        write.xlsx(allcountries_dat_tosave, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/cluster_number_country.xlsx")            

        #update sup tab 1
        toupdate1 <- read.xlsx("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/to_update_sup_tab1.xlsx")
        
        head(toupdate1)
        
        toupdate1 <- toupdate1[,!(colnames(toupdate1) %in%c("k-mean.clustered.group", "Race"))]
        toupdate1$id <- gsub("Rio", "rio", gsub("PI ", "", toupdate1$Genotype))
        head(PCdat)
        all(toupdate1$id %in% PCdat$id)
        PCdat2 <- PCdat
        PCdat2$race <- gsub("\\_.*", "", PCdat2$racereg)
        updated1 <- (merge(toupdate1, PCdat2[,c("id", "race", "clusters")], by = "id", all.x = T, all.y = T, ))
        colnames(updated1)
        
        updated1 <- updated1[,c(2,3,4,5,6,7,8,10, 11, 9)]
        updated1 <- updated1[order(updated1$Genotype),]

        write.xlsx(updated1, "/Users/ksongsom/OneDrive/postdoc/publication/terra/result/updated_sup_tab1.xlsx")        
        
