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

#load map
map <- (read_excel("/Users/ksongsom/OneDrive/postdoc/NAM/Panel/BAP/update_sorghum_location.xlsx", sheet = 1)) #from google sheet
#update elevation
map$ID <- gsub("PI ", "", map$Genotype)
head(map)
map <- map %>% separate(Coordinate, c("lat", "lon"), sep = ", ")

table(fildt[is.na(fildt$end), ]$svtype)
head(fildt)
dim(fildt)
mafcut <- 0.05
missingcut <- 0.2 #0.05
sizecut <- 100000
table(fildt$svtype)
fildt <- fildt[fildt$maf_col >= mafcut & fildt$missing_col <= missingcut, ]
fildt$end <- ifelse(fildt$svtype == "BND", fildt$pos, fildt$end)
fildt$svlen <- fildt$end -fildt$pos
colnames(fildt) <- gsub("X", "", colnames(fildt))
fildtinfo <- fildt[, (which(colnames(fildt) == "format")+1):which(colnames(fildt) == "Rio")]
tfildtinfo <- t(fildtinfo)
rownames(tfildtinfo) <- gsub("X", "", rownames(tfildtinfo))
colnames(tfildtinfo) <- fildt$id
df <- tfildtinfo
df <- df[!(rownames(df) %in% names(which(rowSums(df, na.rm = T) == 0))),]
dim(df)

      df[is.na(df)] <- 0
      df <- na.omit(df)
      df <- scale(df)
      distance <- factoextra::get_dist(df)
 
      #define number of k
      set.seed(123)
      
      # function to compute total within-cluster sum of square 
      wss <- function(k) {
        kmeans(df, k, nstart = 10 )$tot.withinss
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
      
      
      k5 <- kmeans(df, centers = 5, nstart = 25)

      #save clusters
      #write.csv(k5, "/Users/ksongsom/OneDrive/postdoc/publication/terra/code/files/five_cluster.csv")
           
      # plots to compare
      p4 <- factoextra::fviz_cluster(k5, geom = c("point", "text"),  data = df, repel = TRUE, labelsize = 8) + ggtitle("k = 5")
      p4
      
      #different PC
      #load cluster
      cluster <- read.csv("/Users/ksongsom/OneDrive/postdoc/publication/terra/code/files/five_cluster.csv")
      head(cluster)
      dim(cluster)
            
      #run PCA      
            pca_res <- PCA(df, graph = F)
            PCdat <- as.data.frame(pca_res$ind$coord[,1:5])
            PCdat$id <- rownames(PCdat)
            PCdat <- merge(x = PCdat, y = cluster[,c("PI", "k5_cluster_all")], by.x = "id", by.y = "PI", all.x = T, all.y = F)
            PCdat$clusters <- as.character(PCdat$k5_cluster_all)
            
            #making the text for %variance for axis labels
            eigen <- as.data.frame(pca_res$eig)
            
            pc1 <- round(eigen[1,]$`percentage of variance`, digits = 1)
            pc2 <- round(eigen[2,]$`percentage of variance`, digits = 1)
            pc3 <- round(eigen[3,]$`percentage of variance`, digits = 1)
            pc4 <- round(eigen[4,]$`percentage of variance`, digits = 1)
            pc5 <- round(eigen[5,]$`percentage of variance`, digits = 1)
            
            pc1lab = paste("PC1 (", pc1,"%)", sep = "" )
            pc1lab
            pc2lab = paste("PC2 (", pc2,"%)", sep = "")
            pc2lab
            pc3lab = paste("PC3 (", pc3,"%)", sep = "")
            pc3lab
            pc4lab = paste("PC4 (", pc3,"%)", sep = "")
            pc4lab
            pc5lab = paste("PC5 (", pc5,"%)", sep = "")
            pc5lab
            
            
            #define the circle line
            find_hull12 <- function(PCdat) PCdat[chull(PCdat[,2], PCdat[,3]), ]
            hulls12 <- ddply(PCdat, "clusters", find_hull12)
            
            
            find_hull15 <- function(PCdat) PCdat[chull(PCdat[,2], PCdat[,6]), ]
            hulls15 <- ddply(PCdat, "clusters", find_hull15)
            
            plot12 <- ggplot(data = PCdat, aes(x = Dim.1, y = -1*Dim.2, col=clusters, 
                                               fill=clusters, group=clusters)) +
              geom_point(aes(shape=clusters), size = 3) + labs(x = pc1lab, y = pc2lab) +
              scale_color_manual(breaks = c("1", "2", "3", "4", "5"),
                                values=c("gold", "purple", "green", "red", "blue"))  + 
              scale_fill_manual(breaks = c("1", "2", "3", "4", "5"),
                                 values=c("gold", "purple", "green", "red", "blue"))  + 
              geom_polygon(data = hulls12, alpha = 0.3) + theme_bw() +
              theme(text = element_text(size=20))
            
            plot12
            
            plot15 <- ggplot(data = PCdat, aes(x = Dim.1, y = -1*Dim.5, col=clusters, 
                                               fill=clusters, group=clusters)) +
              geom_point(aes(shape=clusters), size = 3) + labs(x = pc1lab, y = pc5lab) +
              scale_color_manual(breaks = c("1", "2", "3", "4", "5"),
                                 values=c("gold", "purple", "green", "red", "blue")
                                )  + 
              scale_fill_manual(breaks = c("1", "2", "3", "4", "5"),
                                values=c("gold", "purple", "green", "red", "blue")
                                )  + 
              geom_polygon(data = hulls15, alpha = 0.3) + theme_bw() +
              theme(text = element_text(size=20))
            
            plot15
            
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
                                 values=c("gold", "purple", "green", "red", "blue")
              )  + 
              scale_fill_manual(breaks = c("1", "2", "3", "4", "5"),
                                values=c("gold", "purple", "green", "red", "blue")
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
                                 values=c("gold", "purple", "green", "red", "blue")
              )  + 
              scale_fill_manual(breaks = c("1", "2", "3", "4", "5"),
                                values=c("gold", "purple", "green", "red", "blue")
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
            
