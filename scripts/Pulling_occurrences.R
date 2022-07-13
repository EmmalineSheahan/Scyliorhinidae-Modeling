library(dismo)
library(robis)
library(sp)
library(sf)
library(dplyr)
library(rnaturalearth)
library(voluModel)
library(spatialEco)

# pulling occurrence for Poroderma africanum
africanum_gbif_raw <- gbif(genus = 'Poroderma', species = 'africanum', geo = T,
                           removeZeros = T, download = T)

africanum_obis_raw <- occurrence(scientificname = "Poroderma africanum", qcfields = T)
class(africanum_obis_raw)

# cleaning occurrences from each dataset
# gbif
head(africanum_gbif_raw)
africanum_gbif_raw <- africanum_gbif_raw[complete.cases(africanum_gbif_raw$depth),]
africanum_gbif_raw <- africanum_gbif_raw[africanum_gbif_raw$depth > 0.0,]
africanum_gbif_raw <- africanum_gbif_raw[africanum_gbif_raw$depth < 4000.0,]
africanum_gbif_raw <- africanum_gbif_raw %>% select(lat, lon, depth, eventDate)
dim(africanum_gbif_raw)
africanum_gbif <- distinct(africanum_gbif_raw)
dim(africanum_gbif)
unique(africanum_gbif$lat)
africanum_gbif_coords <- africanum_gbif %>% select(lat, lon)
unique(africanum_gbif_coords)

# obis
dim(africanum_obis_raw)
africanum_obis_raw <- africanum_obis_raw[complete.cases(africanum_obis_raw$decimalLatitude),]
africanum_obis_raw <- africanum_obis_raw[complete.cases(africanum_obis_raw$decimalLongitude),]
africanum_obis_raw <- africanum_obis_raw[complete.cases(africanum_obis_raw$depth),]
africanum_obis_raw <- africanum_obis_raw[africanum_obis_raw$depth > 0.0,]
africanum_obis_raw <- africanum_obis_raw[africanum_obis_raw$depth < 4000.0,]
africanum_obis_raw <- africanum_obis_raw %>% select(decimalLatitude, decimalLongitude, depth, eventDate)
africanum_obis <- distinct(africanum_obis_raw)
colnames(africanum_obis) <- colnames(africanum_gbif)

# combining data and performing final cleaning
africanum_raw <- rbind(africanum_gbif, africanum_obis)
dim(africanum_raw)
africanum_raw <- africanum_raw[,1:3]
africanum_clean <- distinct(africanum_raw)
dim(africanum_clean)
unique(africanum_clean)
which(is.na(africanum_clean))
africanum_clean[106,]
africanum_clean <- na.omit(africanum_clean)
colnam <- c("latitude", "longitude", "depth")
colnames(africanum_clean) <- colnam
africanum_clean_sp <- africanum_clean

# create spatial points
coordinates(africanum_clean_sp) <- ~longitude+latitude
proj4string(africanum_clean_sp) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]
plot(land)
pointMap(occs = africanum_clean, ptCol = "orange", landCol = "black",
         spName = "Poroderma africanum", ptSize = 3,
         land = land)

# Removing points which intersect with land
land_sp <- as_Spatial(land)
land_sp <- spTransform(land_sp, CRSobj = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
africanum_clean_sp2 <- erase.point(africanum_clean_sp, land_sp, inside = T)
plot(africanum_clean_sp2, col = "red")
plot(land_sp, add = T)
dim(africanum_clean_sp2@coords)
class(africanum_clean_sp2@coords)
africanum_clean_sp2@data
which(!unique(africanum_clean$latitude))
which(!unique(africanum_clean$longitude))
ref <- data.frame(africanum_clean_sp2@coords)
lonremove <- which(!(africanum_clean$longitude %in% ref$longitude))
test <- africanum_clean[-lonremove,]
latremove <- which(!(test$latitude %in% ref$latitude))
test <- test[-latremove,]
africanum_final <- test
save(africanum_final, file = './data/africanum_final.Rdata')

# final plot
png('./figures/Poroderma_africanum_points.png')
plot(land)
pointMap(occs = africanum_final, ptCol = "orange", landCol = "black",
         spName = "Poroderma africanum", ptSize = 3,
         land = land)
dev.off()

# Pulling occurrences for Poroderma pantherinum
pantherinum_gbif_raw <- gbif(genus = 'Poroderma', species = 'pantherinum', geo = T,
                           removeZeros = T, download = T)

pantherinum_obis_raw <- occurrence(scientificname = "Poroderma pantherinum", qcfields = T)

# cleaning occurrences from each dataset
# gbif
dim(pantherinum_gbif_raw)
pantherinum_gbif_raw <- pantherinum_gbif_raw[complete.cases(pantherinum_gbif_raw$depth),]
pantherinum_gbif_raw <- pantherinum_gbif_raw[pantherinum_gbif_raw$depth > 0.0,]
pantherinum_gbif_raw <- pantherinum_gbif_raw[pantherinum_gbif_raw$depth < 4000.0,]
pantherinum_gbif_raw <- pantherinum_gbif_raw %>% select(lat, lon, depth, eventDate)
dim(pantherinum_gbif_raw)
pantherinum_gbif <- distinct(pantherinum_gbif_raw)
dim(pantherinum_gbif)
pantherinum_gbif_coords <- pantherinum_gbif %>% select(lat, lon)
unique(pantherinum_gbif_coords)

# obis
dim(pantherinum_obis_raw)
pantherinum_obis_raw <- pantherinum_obis_raw[complete.cases(pantherinum_obis_raw$decimalLatitude),]
pantherinum_obis_raw <- pantherinum_obis_raw[complete.cases(pantherinum_obis_raw$decimalLongitude),]
pantherinum_obis_raw <- pantherinum_obis_raw[complete.cases(pantherinum_obis_raw$depth),]
pantherinum_obis_raw <- pantherinum_obis_raw[pantherinum_obis_raw$depth > 0.0,]
pantherinum_obis_raw <- pantherinum_obis_raw[pantherinum_obis_raw$depth < 4000.0,]
pantherinum_obis_raw <- pantherinum_obis_raw %>% select(decimalLatitude, decimalLongitude, depth, eventDate)
pantherinum_obis <- distinct(pantherinum_obis_raw)
dim(pantherinum_obis)
colnames(pantherinum_obis) <- colnames(pantherinum_gbif)

# combining data
pantherinum_raw <- rbind(pantherinum_gbif, pantherinum_obis)
dim(pantherinum_raw)
pantherinum_raw <- pantherinum_raw[,1:3]
pantherinum_clean <- distinct(pantherinum_raw)
dim(pantherinum_clean)
unique(pantherinum_clean)
which(is.na(pantherinum_clean))
pantherinum_clean[31,]
pantherinum_clean <- na.omit(pantherinum_clean)
colnam <- c("latitude", "longitude", "depth")
colnames(pantherinum_clean) <- colnam
pantherinum_clean_sp <- pantherinum_clean

# create spatial points
coordinates(pantherinum_clean_sp) <- ~longitude+latitude
proj4string(pantherinum_clean_sp) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]
plot(land)
pointMap(occs = pantherinum_clean, ptCol = "orange", landCol = "black",
         spName = "Poroderma pantherinum", ptSize = 3,
         land = land)

# Removing points which intersect with land
land_sp <- as_Spatial(land)
land_sp <- spTransform(land_sp, CRSobj = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
pantherinum_clean_sp2 <- erase.point(pantherinum_clean_sp, land_sp, inside = T)
plot(pantherinum_clean_sp2, col = "red")
plot(land_sp, add = T)
which(!unique(pantherinum_clean$latitude))
which(!unique(pantherinum_clean$longitude))
ref <- data.frame(pantherinum_clean_sp2@coords)
lonremove <- which(!(pantherinum_clean$longitude %in% ref$longitude))
test <- pantherinum_clean[-lonremove,]
which(!(test$latitude %in% ref$latitude))
pantherinum_final <- test
save(pantherinum_final, file = './data/pantherinum_final.Rdata')

# final plot
png('./figures/Poroderma_pantherinum_points.png')
plot(land)
pointMap(occs = pantherinum_final, ptCol = "green", landCol = "black",
         spName = "Poroderma pantherinum", ptSize = 3,
         land = land)
dev.off()

# both species plot
png('./figures/both_species_points.png')
pointCompMap(occs1 = africanum_final, occs2 = pantherinum_final, spName = "Poroderma africanum and Poroderma pantherinum",
             land = land, occs1Col = "orange", occs2Col = "green", occs1Name = "Poroderma africanum",
             occs2Name = "Poroderma pantherinum", landCol = "black", ptSize = 3)
dev.off()

# investigate data
min(africanum_final$depth)
max(africanum_final$depth)

min(pantherinum_final$depth)
max(pantherinum_final$depth)

africanum_final[which(africanum_final$depth > 110),]

pantherinum_final[which(pantherinum_final$depth > 256),]
