library(ggrepel)
library(sf)
library(ggmap)
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
library(openxlsx)
library(readxl)

options(stringsAsFactors = F)
background <- read.csv("/Users/ksongsom/OneDrive/postdoc/NAM/Panel/BAP/PanelAccessions-BAP.csv", header = T, stringsAsFactors = F)
background$PI <- gsub("RIO","rio", gsub(paste(c("PI_", "-dark"), collapse = "|"),"", background$Taxa))
table(background$Origin)
coor <- read.csv("/Users/ksongsom/OneDrive/postdoc/programs/country-capitals.csv", header = T)
head(coor)
sum(background$Origin %in% coor$CountryName)

subregions_dat <- read_excel("/Users/ksongsom/OneDrive/postdoc/NAM/Panel/BAP/subregions.xlsx", sheet = 3)

#update sub table
sub1 <- read_excel("/Users/ksongsom/OneDrive/postdoc/publication/terra/to_submit_revised/Songsomboon_table_revised1.xlsx", sheet = "Supplement table 1", skip = 2)
all(sub1$Country %in% subregions_dat$country)
upsub1 <- merge(sub1, subregions_dat, by.x = "Country", by.y = "country", all.x = T, all.y = F)
head(upsub1)
upsub1 <- upsub1[order(upsub1$Genotype),]
#write.xlsx(upsub1, "/Users/ksongsom/OneDrive/postdoc/publication/terra/to_submit_revised/Songsomboon_table_sup1.xlsx")

#update table 3
tab3 <- read_excel("/Users/ksongsom/OneDrive/postdoc/publication/terra/to_submit_revised/Songsomboon_table_revised1.xlsx", sheet = "Table 3", skip = 2)
tab3$`Countries of origin`[which(!(tab3$`Countries of origin` %in% subregions_dat$country))]
tab3 <- tab3[-nrow(tab3),]
uptab3 <- merge(y = tab3, x = subregions_dat, by.y = "Countries of origin", by.x = "country", all.x = F, all.y = T)
head(uptab3)
uptab3$Total <- as.numeric(uptab3$Total)
uptab3 <- uptab3[order(uptab3$Total, decreasing = T),]
#write.xlsx(uptab3, "/Users/ksongsom/OneDrive/postdoc/publication/terra/to_submit_revised/Songsomboon_table_3.xlsx")


#load map
map <- (read_excel("/Users/ksongsom/OneDrive/postdoc/NAM/Panel/BAP/update_sorghum_location_subregions.xlsx", sheet = 1)) #from google sheet
#update elevation
map$ID <- gsub("PI ", "", map$Genotype)
head(map)
#map <- map %>% separate(Coordinate, c("lat", "lon"), sep = ", ")
map$ID <- gsub("Rio", "rio", map$ID)
map$Race[(map$Race == "NA")] <- "Unknown"
map$Race[(map$Race == "Mixed")] <- "Unknown"
map$subregions_torace <- ifelse(grepl("Africa", map$subregions), map$subregions, "outside Africa")
map$racereg <- paste(map$Race, map$subregions_torace, sep = "_")

all(map$Country %in% subregions$country)

#update map center of country if province is NA
mapdata <- read.csv("/Users/ksongsom/OneDrive/postdoc/publication/terra/code/files/map_plot347_middle.csv")
mapdata$coordinate <- paste(mapdata$lat, mapdata$lon, sep = ", ")

#update cluster from SNP
updated1 <- read.xlsx("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/updated_sup_tab1.xlsx")
updated1$ID <- gsub("PI ", "", updated1$Genotype)
all(updated1$ID %in% mapdata$ID)

mapdata <- (merge(mapdata[, !(colnames(mapdata) %in% c("Race"))], updated1[,c("ID", "Race", "Cluster")], by = "ID"))

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
mapdata$cluster <- paste("cluster", mapdata$Cluster, sep = "_")
newmap <- dcast(mapdata, Country ~ cluster, value.var = "Cluster")

#size of each pie plot
newmap$radius <- rescale(scale(rowSums(newmap[,-1])), to = c(3, 8))
head(newmap)

newmap <- merge(newmap, mid, by.x = "Country", by.y = "name", all.x = T, all.y = F)
#tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/world_dis_mid_8_cluster.#tiff", width=7, height=5, units="in", res=300)
p + geom_scatterpie(aes(x=longitude, y=latitude, r = radius),
                    data=newmap, alpha=.8, cols = paste("cluster", 1:8, sep = "_"),  color=NA) +
  labs(y = "Latitude", x = "Longitude") +
  scale_fill_manual("8 clusters", breaks = c("cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_5", "cluster_6", "cluster_7", "cluster_8"),
  labels = c("cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_5", "cluster_6", "cluster_7", "cluster_8"),
  values = c("cluster_1" = "red",
             "cluster_2" = "deepskyblue1",
             "cluster_3" = "dark green",
             "cluster_4" = "green",
             "cluster_5" = "purple",
             "cluster_6" = "orange", 
             "cluster_7" = "gold",
             "cluster_8" = "pink")) 
#dev.off()

#sub regions
newmap
all(newmap$Country %in% subregions$country)
newmap_sub <- merge(newmap, subregions_dat, by.x = "Country", by.y = "country", all.x = T, all.y = F)
newmap_sub <- newmap_sub %>% group_by(subregions) %>% dplyr::summarise(
  cluster_1 = sum(cluster_1, na.rm = T),
  cluster_2 = sum(cluster_2, na.rm = T),
  cluster_3 = sum(cluster_3, na.rm = T),
  cluster_4 = sum(cluster_4, na.rm = T),
  cluster_5 = sum(cluster_5, na.rm = T),
  cluster_6 = sum(cluster_6, na.rm = T),
  cluster_7 = sum(cluster_7, na.rm = T),
  cluster_8 = sum(cluster_8, na.rm = T),
  latitude = mean(latitude, na.rm = T),
  longitude = mean(longitude, na.rm = T)
)
newmap_sub

newmap_sub$radius <- rescale(scale(rowSums(newmap_sub[,2:9])), to = c(3, 10))

p + geom_scatterpie(aes(x=longitude, y=latitude, r = radius),
                    data=newmap_sub, alpha=.8, cols = paste("cluster", 1:8, sep = "_"),  color=NA) +
  labs(y = "Latitude", x = "Longitude") +
  scale_fill_manual("8 clusters", breaks = c("cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_5", "cluster_6", "cluster_7", "cluster_8"),
                    labels = c("cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_5", "cluster_6", "cluster_7", "cluster_8"),
                    values = c("cluster_1" = "red",
                               "cluster_2" = "deepskyblue1",
                               "cluster_3" = "dark green",
                               "cluster_4" = "green",
                               "cluster_5" = "purple",
                               "cluster_6" = "orange", 
                               "cluster_7" = "gold",
                               "cluster_8" = "pink")) +
  geom_label_repel(data = newmap_sub2, aes(x = lon, y = lat, label = subregions), 
                   fill = "white", box.padding = unit(.4, "lines"),
                   label.padding = unit(.15, "lines"),
                   segment.color = NA, segment.size = 1,
                   direction = "x")

#dev.off()

#zoom Africa
africamap <- cmapdata[cmapdata$ContinentName == "Africa",]
mapdata$coordinate <- paste(mapdata$lat, mapdata$lon, sep = "_")
newmap2 <- dcast(mapdata, coordinate ~ Cluster, value.var = "Cluster")
head(newmap2)
#newmap$radius <- rowSums(newmap[,-1])*5/(max(rowSums(newmap[,-1])))

newmap2$radius <- rescale(scale(rowSums(newmap2[,-1])), to = c(1,2))#scale(rowSums(newmap[,-1]))

newmap2 <- merge(newmap2, mapdata, by.x = "coordinate", by.y = "coordinate", all.x = T, all.y = F)

#API from google

register_google(key = "AIzaSyDVIwuIHmOBGFHH2zGQYrkUXMI--EzZdo0")

mapaf <- ggmap::get_map(location = c(lon = mean(africamap$CapitalLongitude, na.rm = T), 
                              lat = mean(africamap$CapitalLatitude, na.rm = T)), zoom = 3,
                       maptype = "terrain", scale = 4)

mapaf2 <- get_map(location = c(lon = mean(africamap$CapitalLongitude, na.rm = T), lat = mean(africamap$CapitalLatitude, na.rm = T)), zoom = 3,
                maptype = "hybrid", scale = 4)
afmap <- ggmap(mapaf) +
  scale_x_continuous(limits = c(-20, 60)) +
  scale_y_continuous(limits = c(-35, 40))

#newmap2$Cluster <- paste("cluster ", newmap2$Cluster, sep = "")
colnames(newmap2)[2:9] <- paste("cluster_", colnames(newmap2)[2:9], sep = "")
#tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/africa_dis_mid_8_cluster.#tiff", width=7, height=5, units="in", res=300)
afmap + geom_scatterpie(aes(x=lon, y=lat, r = radius),
                        data=newmap2, alpha=.8, cols = paste("cluster_", 1:8, sep = ""),  color=NA) +
  labs(y = "Latitude", x = "Longitude") +
  scale_fill_manual("8 clusters", breaks = c("cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_5", "cluster_6", "cluster_7", "cluster_8"),
                    labels = c("cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_5", "cluster_6", "cluster_7", "cluster_8"),
                    values = c("cluster_1" = "red",
                               "cluster_2" = "deepskyblue1",
                               "cluster_3" = "dark green",
                               "cluster_4" = "green",
                               "cluster_5" = "purple",
                               "cluster_6" = "orange", 
                               "cluster_7" = "gold",
                               "cluster_8" = "pink"))    
#dev.off()


#sub regions
newmap
all(newmap2$Country %in% subregions_dat$country)
newmap_sub2 <- merge(newmap2, subregions_dat, by.x = "Country", by.y = "country", all.x = T, all.y = F)
newmap_sub2 <- newmap_sub2 %>% group_by(subregions) %>% dplyr::summarise(
  cluster_1 = sum(cluster_1, na.rm = T),
  cluster_2 = sum(cluster_2, na.rm = T),
  cluster_3 = sum(cluster_3, na.rm = T),
  cluster_4 = sum(cluster_4, na.rm = T),
  cluster_5 = sum(cluster_5, na.rm = T),
  cluster_6 = sum(cluster_6, na.rm = T),
  cluster_7 = sum(cluster_7, na.rm = T),
  cluster_8 = sum(cluster_8, na.rm = T),
  lat = mean(lat, na.rm = T),
  lon = mean(lon, na.rm = T)
)
newmap_sub2

newmap_sub2$radius <- rescale(scale(rowSums(newmap_sub2[,2:9])), to = c(3, 8))

afmap + geom_scatterpie(aes(x=lon, y=lat, r = radius, label = subregions),
                        data=newmap_sub2, alpha=.8, cols = paste("cluster_", 1:8, sep = ""),  color=NA) +
  labs(y = "Latitude", x = "Longitude") +
  scale_fill_manual("8 clusters", breaks = c("cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_5", "cluster_6", "cluster_7", "cluster_8"),
                    labels = c("cluster_1", "cluster_2", "cluster_3", "cluster_4", "cluster_5", "cluster_6", "cluster_7", "cluster_8"),
                    values = c("cluster_1" = "red",
                               "cluster_2" = "deepskyblue1",
                               "cluster_3" = "dark green",
                               "cluster_4" = "green",
                               "cluster_5" = "purple",
                               "cluster_6" = "orange", 
                               "cluster_7" = "gold",
                               "cluster_8" = "pink")) +
  geom_label_repel(data = newmap_sub2, aes(x = lon, y = lat, label = subregions), 
                   fill = "white", box.padding = unit(.4, "lines"),
                   label.padding = unit(.15, "lines"),
                   segment.color = NA, segment.size = 1,
                   direction = "x")





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

#tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/world_africa_mid.#tiff", width=7, height=14, units="in", res=300)
grid.arrange(gb1, gb2, nrow = 2)
#dev.off()

#Elevation
mapdata$coordinate <- paste(mapdata$lat, mapdata$lon, sep = "_")
group_name <- c("<623", "623-1200", "1200-1820", "1820-2410", ">2410")
mapdata$elegroup <- cut(mapdata$Elevation, breaks = 5, labels = group_name )

newmap_ele <- dcast(mapdata, coordinate ~ elegroup, value.var = "Cluster")
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


#tiff("/Users/ksongsom/OneDrive/postdoc/publication/terra/result/africa_dis_mid.#tiff", width=7, height=5, units="in", res=300)
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
#dev.off()

