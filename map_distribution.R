library(sf)
library(raster)
library(dplyr)
library(tmap)
library(leaflet)
library(cartogram)
library(rworldmap)
library(rworldxtra)
library(plyr)
library(gridExtra)
#library(ggtree)
library(reshape2)
library(scatterpie)
library(colorspace)
library(DescTools)
library(ggplot2)
library(scales)
library(spData)

options(stringsAsFactors = F)
background <- read.csv("/Users/ksongsom/OneDrive/postdoc/NAM/Panel/BAP/PanelAccessions-BAP.csv", header = T, stringsAsFactors = F)
background$PI <- gsub("RIO","rio", gsub(paste(c("PI_", "-dark"), collapse = "|"),"", background$Taxa))
table(background$Origin)
coor <- read.csv("/Users/ksongsom/OneDrive/postdoc/programs/country-capitals.csv", header = T)
head(coor)
sum(background$Origin %in% coor$CountryName)

#update map center of country if province is NA
mapdata <- read.csv("/Users/ksongsom/OneDrive/postdoc/publication/terra/code/files/map_plot347_middle.csv")
mapdata$coordinate <- paste(mapdata$lat, mapdata$lon, sep = ", ")

mapdata$Country[which(!(mapdata$Country %in% coor$CountryName))]
head(mapdata)
dim(mapdata)
head(coor)
dim(mapdata)

#capital
cmapdata <- merge(mapdata, coor, by.x = "Country", by.y = "CountryName", all.x = T, all.y = F)
dim(cmapdata)
cmapdata[1:3,]

#middle of all countries
mid <- read.csv("/Users/ksongsom/OneDrive/postdoc/publication/terra/code/files/middle_database.csv")

world <- map_data('world')
p <- ggplot(world, aes(long, lat)) +
  geom_map(map=world, aes(map_id=region), fill="white", color=NA) +
  coord_quickmap() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))
  
p

mapdata[1:3,]
head(mapdata)

#aggregate clusters with countries
mapdata$cluster <- paste("cluster", mapdata$k5g, sep = "_")
newmap <- dcast(mapdata, Country ~ cluster, value.var = "k5g")

#size of each pie plot
newmap$radius <- rescale(scale(rowSums(newmap[,-1])), to = c(3, 10))
head(newmap)

newmap <- merge(newmap, mid, by.x = "Country", by.y = "name", all.x = T, all.y = F)
tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/world_dis_mid.tiff", width=7, height=5, units="in", res=300)
p + geom_scatterpie(aes(x=longitude, y=latitude, r = radius),
                    data=newmap, alpha=.8, cols = paste("cluster", 1:5, sep = "_"),  color=NA) +
  labs(y = "Latitude", x = "Longitude") +
  scale_fill_manual("5 clusters", breaks = c("cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_5"),
  labels = c("cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_5"),
  values = c("cluster_1" = "gold",
             "cluster_2" = "purple",
             "cluster_3" = "green",
             "cluster_4" = "red",
             "cluster_5" = "blue"))
dev.off()

#zoom Africa
africamap <- cmapdata[cmapdata$ContinentName == "Africa",]
mapdata$coordinate <- paste(mapdata$lat, mapdata$lon, sep = "_")
newmap2 <- dcast(mapdata, coordinate ~ cluster, value.var = "k5g")
head(newmap2)
#newmap$radius <- rowSums(newmap[,-1])*5/(max(rowSums(newmap[,-1])))

newmap2$radius <- rescale(scale(rowSums(newmap2[,-1])), to = c(1,3))#scale(rowSums(newmap[,-1]))

newmap2 <- merge(newmap2, mapdata, by.x = "coordinate", by.y = "coordinate", all.x = T, all.y = F)


mapaf<- get_map(location = c(lon = mean(africamap$CapitalLongitude, na.rm = T), lat = mean(africamap$CapitalLatitude, na.rm = T)), zoom = 3,
                       maptype = "terrain", scale = 4)
mapaf2<- get_map(location = c(lon = mean(africamap$CapitalLongitude, na.rm = T), lat = mean(africamap$CapitalLatitude, na.rm = T)), zoom = 3,
                maptype = "hybrid", scale = 4)
afmap <- ggmap(mapaf) +
  scale_x_continuous(limits = c(-20, 60)) +
  scale_y_continuous(limits = c(-35, 40))


tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/africa_dis_mid.tiff", width=7, height=5, units="in", res=300)
afmap + geom_scatterpie(aes(x=lon, y=lat, r = radius),
                        data=newmap2, alpha=.8, cols = paste("cluster", 1:5, sep = "_"),  color=NA) +
  labs(y = "Latitude", x = "Longitude") +
  scale_fill_manual("5 clusters", breaks = c("cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_5"),
                    labels = c("cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_5"),
                    values = c("cluster_1" = "gold",
                               "cluster_2" = "purple",
                               "cluster_3" = "green",
                               "cluster_4" = "red",
                               "cluster_5" = "blue"))
dev.off()




gb1 <- p + geom_scatterpie(aes(x=longitude, y=latitude, r = radius),
                    data=newmap, alpha=.8, cols = paste("cluster", 1:5, sep = "_"),  color=NA) +
  labs(y = "Latitude", x = "Longitude") +
  scale_fill_manual("5 clusters", breaks = c("cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_5"),
                    labels = c("cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_5"),
                    values = c("cluster_1" = "gold",
                               "cluster_2" = "purple",
                               "cluster_3" = "green",
                               "cluster_4" = "red",
                               "cluster_5" = "blue"), guide = FALSE)
gb1
gb2 <- afmap + geom_scatterpie(aes(x=lon, y=lat, r = radius),
                        data=newmap2, alpha=.8, cols = paste("cluster", 1:5, sep = "_"),  color=NA) +
  labs(y = "Latitude", x = "Longitude") +
  scale_fill_manual("5 clusters", breaks = c("cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_5"),
                    labels = c("cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_5"),
                    values = c("cluster_1" = "gold",
                               "cluster_2" = "purple",
                               "cluster_3" = "green",
                               "cluster_4" = "red",
                               "cluster_5" = "blue")) +
  theme(legend.position = c(0.2,0.25))
gb2

#combine the plots

tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/world_africa_mid.tiff", width=7, height=14, units="in", res=300)
grid.arrange(gb1, gb2, nrow = 2)
dev.off()

#Elevation
mapdata$coordinate <- paste(mapdata$lat, mapdata$lon, sep = "_")
group_name <- c("<623", "623-1200", "1200-1820", "1820-2410", ">2410")
mapdata$elegroup <- cut(mapdata$Elevation, breaks = 5, labels = group_name )

newmap_ele <- dcast(mapdata, coordinate ~ elegroup, value.var = "k5g")
head(newmap_ele)
#newmap$radius <- rowSums(newmap[,-1])*5/(max(rowSums(newmap[,-1])))

newmap_ele$radius <- rescale(scale(rowSums(newmap_ele[,-1])), to = c(1,3))#scale(rowSums(newmap[,-1]))

newmap_ele <- merge(newmap_ele, mapdata, by.x = "coordinate", by.y = "coordinate", all.x = T, all.y = F)

africamap <- newmap[newmap$ContinentName == "Africa",]

mapaf<- get_map(location = c(lon = mean(africamap$CapitalLongitude, na.rm = T), lat = mean(africamap$CapitalLatitude, na.rm = T)), zoom = 3,
                maptype = "terrain", scale = 4)
afmap <- ggmap(mapaf) +
  scale_x_continuous(limits = c(-20, 60)) +
  scale_y_continuous(limits = c(-35, 40))


tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/africa_dis_mid.tiff", width=7, height=5, units="in", res=300)
afmap + geom_scatterpie(aes(x=lon, y=lat, r = radius),
                        data=newmap_ele, alpha=.8, cols =group_name ,  color=NA) +
  labs(y = "Latitude", x = "Longitude") +
  scale_fill_manual("Elevation", breaks = group_name ,
                    labels = group_name,
                    values = c("<623" = "gold",
                               "623-1200" = "purple",
                               "1200-1820" = "green",
                               "1820-2410" = "red",
                               ">2410" = "blue"))
dev.off()


#phenotype correlations
pheno16 <- read.csv("/Users/ksongsom/OneDrive/postdoc/NAM/Panel/BAP/pheno16.csv", header = T, stringsAsFactors = F, sep = ",")
head(pheno16)
pheno17 <- read.csv("/Users/ksongsom/OneDrive/postdoc/NAM/Panel/BAP/pheno17.csv", header = T, stringsAsFactors = F, sep = "\t")
head(pheno17)
background <- read.csv("/Users/ksongsom/OneDrive/postdoc/NAM/Panel/BAP/PanelAccessions-BAP.csv", header = T, stringsAsFactors = F)
dim(pheno16)

background$PI <- gsub("RIO","rio", gsub("-dark","", background$Taxa))
sum(pheno16$PI %in% background$Taxa)
sum(pheno16$PI %in% background$PI)
pheno16$PI[which(!(pheno16$PI %in% background$PI))]

dat_back <- merge(pheno16, background, by.y = "PI", by.x = "PI", all.x = T, all.y = F)
dim(dat_back)

#write.csv(dat_back, "/Users/ksongsom/OneDrive/postdoc/NAM/Panel/BAP/ID_oldmap.csv", row.names = F)
#write.csv(dat_back[!duplicated(dat_back$PI),], "/Users/ksongsom/OneDrive/postdoc/NAM/Panel/BAP/ID_oldmap_unique.csv", row.names = F)


mapdata[1:3,]
fivetest <- c()
for(i in which(colnames(mapdata)=="ADF"):which(colnames(mapdata)=="Hemicellulose")){
inputdt <- na.omit(mapdata[,c(2,i)])
tt <- colnames(mapdata)[i]
pval <- kruskal.test(inputdt[,2] ~ inputdt[,1])$p.value
sig <- ifelse(pval < 0.05/5, "*", NA)
eachtest <- cbind.data.frame(tt, pval, sig)
fivetest <- rbind.data.frame(fivetest, eachtest)
}
fivetest
