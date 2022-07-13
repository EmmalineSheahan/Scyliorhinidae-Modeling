library(voluModel)
library(dplyr)
library(sp)
library(raster)
library(rnaturalearth)

load('./data/africanum_final.Rdata')
load('./data/pantherinum_final.Rdata')
load('./data/africanum_accarea.Rdata')
load('./data/pantherinum_accarea.Rdata')
load('./data/rasters/aou_3d.Rdata')
load('./data/rasters/aou_Bottom.Rdata')
load('./data/rasters/aou_Surface.Rdata')
load('./data/rasters/nitrate_3d.Rdata')
load('./data/rasters/nitrate_Bottom.Rdata')
load('./data/rasters/nitrate_Surface.Rdata')
load('./data/rasters/oxygen_3d.Rdata')
load('./data/rasters/oxygen_Bottom.Rdata')
load('./data/rasters/oxygen_Surface.Rdata')
load('./data/rasters/phosphate_3d.Rdata')
load('./data/rasters/phosphate_Bottom.Rdata')
load('./data/rasters/phosphate_Surface.Rdata')
load('./data/rasters/salinity_3d.Rdata')
load('./data/rasters/salinity_Bottom.Rdata')
load('./data/rasters/salinity_Surface.Rdata')
load('./data/rasters/temperature_3d.Rdata')
load('./data/rasters/temperature_Bottom.Rdata')
load('./data/rasters/temperature_Surface.Rdata')

# pleminary glm models for poroderma africanum
# the max depth value for africanum is 424m, investigating max depth in accessible area
bathy <- raster('./data/rasters/ETOPO1_Bed_c_geotiff.tif')
africanum_sp <- africanum_final
coordinates(africanum_sp) <- ~longitude+latitude
proj4string(africanum_sp) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]
bathy1 <- intersect(bathy, africanum_accarea)
plot(bathy1)
plot(africanum_accarea, add = T)
plot(land, col = "black", add = T)
bathy2 <- mask(bathy1, mask = africanum_accarea)
min(bathy2@data@values, na.rm = T)
plot(bathy2)
plot(africanum_sp, col = "red", add = T)
plot(land, col = "black", add = T)

# I'm going to remove all depth layers beyond a 100 m buffer beyond the deepest point at which africanum is found
# in order to avoid sampling at depths which are inaccessible to it

max(africanum_final$depth)
names(aou_3d)
# cutting off depth sampling for africanum at 550 m
aou_3d_africanum <- aou_3d[[1:38]]
nitrate_3d_africanum <- nitrate_3d[[1:38]]
oxygen_3d_africanum <- oxygen_3d[[1:38]]
phosphate_3d_africanum <- phosphate_3d[[1:38]]
salinity_3d_africanum <- salinity_3d[[1:38]]
temperature_3d_africanum <- temperature_3d[[1:38]]

# Gets the layer index for each occurrence by matching to depth
layerNames <- as.numeric(gsub("[X]", "", names(temperature_3d_africanum)))
africanum_final$index <- unlist(lapply(africanum_final$depth, FUN = function(x) which.min(abs(layerNames - x))))
indices <- unique(africanum_final$index)

# downsampling occurrences
africanum_down <- data.frame()
for(i in indices){
  tempPoints <- africanum_final[africanum_final$index==i,]
  tempPoints <- downsample(tempPoints, temperature_3d_africanum[[1]])
  tempPoints$depth <- rep(layerNames[[i]], times = nrow(tempPoints))
  africanum_down <- rbind(africanum_down, tempPoints)
}

save(africanum_down, file = './data/africanum_down.Rdata')

print(paste0("Original number of points: ", nrow(africanum_final), "; number of downsampled occs: ", 
             nrow(africanum_down)))
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]
pointCompMap(occs1 = africanum_final, occs2 = africanum_down, 
             occs1Name = "Original", occs2Name = "Downsampled", 
             spName = "Poroderma africanum", 
             land = land)

aouVals <- xyzSample(occs = africanum_down, aou_3d_africanum)
tempVals <- xyzSample(occs = africanum_down, temperature_3d_africanum)
nitVals <- xyzSample(occs = africanum_down, nitrate_3d_africanum)
phosVals <- xyzSample(occs = africanum_down, phosphate_3d_africanum)
salVals <- xyzSample(occs = africanum_down, salinity_3d_africanum)
oxyVals <- xyzSample(occs = africanum_down, oxygen_3d_africanum)
vals <- cbind(africanum_down, oxyVals, tempVals, aouVals, nitVals, phosVals, salVals)
colnames(vals) <- c("decimalLongitude", "decimalLatitude", "depth", "Oxygen", "Temperature", "AOU", 
                    "Nitrate", "Phosphate", "Salinity")
vals <- vals[complete.cases(vals),]
row.names(vals) <- NULL
occsWdata <- vals

# Add response as a column
occsWdata$response <- rep(1, times = nrow(occsWdata))

# Background
backgroundVals <- mSampling3D(occs = africanum_down, 
                              envBrick = temperature_3d_africanum, 
                              mShp = africanum_accarea, 
                              depthLimit = "all")
aouVals <- xyzSample(occs = backgroundVals, aou_3d_africanum)
tempVals <- xyzSample(occs = backgroundVals, temperature_3d_africanum)
nitVals <- xyzSample(occs = backgroundVals, nitrate_3d_africanum)
phosVals <- xyzSample(occs = backgroundVals, phosphate_3d_africanum)
salVals <- xyzSample(occs = backgroundVals, salinity_3d_africanum)
oxyVals <- xyzSample(occs = backgroundVals, oxygen_3d_africanum)
vals <- cbind(backgroundVals, oxyVals, tempVals, aouVals, nitVals, phosVals, salVals)
colnames(vals) <- c("decimalLongitude", "decimalLatitude", "depth", "Oxygen", "Temperature", "AOU", 
                    "Nitrate", "Phosphate", "Salinity")
vals <- vals[complete.cases(vals),]
row.names(vals) <- NULL
backgroundWdata <- vals
# Add response as a column
backgroundWdata$response <- rep(0, times = nrow(backgroundWdata))

# Sample background points weighted by distance from centroid of occurrence environments
suitableCentroid <- apply(occsWdata[,c("Oxygen", "Temperature", "AOU", "Nitrate", "Phosphate", "Salinity")], 
                          MARGIN = 2, FUN = mean)
backgroundWdata$distance <- apply(backgroundWdata[,c("Oxygen", "Temperature", "AOU", 
                                                     "Nitrate", "Phosphate", "Salinity")], MARGIN = 1, 
                                  FUN = function(x) dist(rbind(suitableCentroid, x)))
backgroundWdata$sampleWeight <- (backgroundWdata$distance - 
                                   min(backgroundWdata$distance))/(max(backgroundWdata$distance)-
                                                                     min(backgroundWdata$distance))
sampleForAbsence <- sample(x = rownames(backgroundWdata), 
                           size = nrow(occsWdata) * 4, 
                           prob = backgroundWdata$sampleWeight)
backgroundWdata <- backgroundWdata[match(sampleForAbsence, 
                                         rownames(backgroundWdata)),]

# Unite datasets
datForMod <- rbind(occsWdata, backgroundWdata[,colnames(occsWdata)])
rm(suitableCentroid, sampleForAbsence, backgroundWdata, occsWdata)

# generate models with all 6 variables
glmModel_africanum <- glm(formula = response ~ Temperature *  AOU * Oxygen * Nitrate * Phosphate * Salinity, 
                family = binomial(link = "logit"),  data = datForMod)
summary(glmModel_africanum)
glmModel_africanum$coefficients
save(glmModel_africanum, file = './results/glmModel_africanum.Rdata')

layerNames <- as.numeric(gsub("[X]", "", names(temperature_3d_africanum)))
index <- seq(from = match(min(datForMod$depth), layerNames), 
             to = match(max(datForMod$depth), layerNames), by = 1)
depthPred <- NULL
for(j in index){
  depthPreds <- stack(temperature_3d_africanum[[j]], aou_3d_africanum[[j]], oxygen_3d_africanum[[j]],
                      nitrate_3d_africanum[[j]], phosphate_3d_africanum[[j]], salinity_3d_africanum[[j]])
  crs(depthPreds) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  names(depthPreds) <- c("Temperature", "AOU", "Oxygen", "Nitrate", "Phosphate", "Salinity")
  depthPred[[j]] <- mask(predict(depthPreds, glmModel_africanum), africanum_accarea)
  depthPred[[j]] <- crop(depthPred[[j]], africanum_accarea)
  names(depthPred[[j]]) <- layerNames[[j]]
}
glmPred <- stack(depthPred[!unlist(lapply(depthPred, FUN = function(x) is.null(x)))])

glmThreshold <- quantile(xyzSample(datForMod[datForMod$response == 1,
                                             c("depth", "decimalLongitude", "decimalLatitude")], 
                                   brick(glmPred)), .1, na.rm = T)[[1]] # MS90
glmThresholded <- glmPred > glmThreshold
glmThresholded <- reclassify(glmThresholded, c(NA, NA, 0), include.lowest = T)

layerNames <- as.numeric(gsub("[X]", "", names(glmThresholded)))
datForMod$index <- unlist(lapply(datForMod$depth, FUN = function(x) which.min(abs(layerNames - x))))
indices <- unique(datForMod$index)
pdf('./figures/africanum_all_depths.pdf')
plotLayers(glmThresholded[[min(indices):max(indices)]], 
           land = land, landCol = "black", title = "GLM Thresholded by depth for Poroderma africanum")
dev.off()

pdf('./figures/africanum_threshold_by_depth.pdf')
for (i in seq_along(names(glmThresholded))) {
plot(glmThresholded[[i]], main = paste0("Poroderma africanum thresholded glm at ", layerNames[i], " m"))
plot(land, add = T, col = "black")
plot(africanum_sp, add = T, col = "red")
}
dev.off()

# generate model with phosphate and aou excluded
glmModel_africanum2 <- glm(formula = response ~ Temperature * Oxygen * Nitrate * Salinity, 
                          family = binomial(link = "logit"),  data = datForMod)
summary(glmModel_africanum2)
glmModel_africanum2$coefficients
save(glmModel_africanum2, file = './results/glmModel_africanum2.Rdata')

layerNames <- as.numeric(gsub("[X]", "", names(temperature_3d_africanum)))
index <- seq(from = match(min(datForMod$depth), layerNames), 
             to = match(max(datForMod$depth), layerNames), by = 1)
depthPred <- NULL
for(j in index){
  depthPreds <- stack(temperature_3d_africanum[[j]], oxygen_3d_africanum[[j]],
                      nitrate_3d_africanum[[j]], salinity_3d_africanum[[j]])
  crs(depthPreds) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  names(depthPreds) <- c("Temperature", "Oxygen", "Nitrate", "Salinity")
  depthPred[[j]] <- mask(predict(depthPreds, glmModel_africanum2), africanum_accarea)
  depthPred[[j]] <- crop(depthPred[[j]], africanum_accarea)
  names(depthPred[[j]]) <- layerNames[[j]]
}
glmPred <- stack(depthPred[!unlist(lapply(depthPred, FUN = function(x) is.null(x)))])

glmThreshold <- quantile(xyzSample(datForMod[datForMod$response == 1,
                                             c("depth", "decimalLongitude", "decimalLatitude")], 
                                   brick(glmPred)), .1, na.rm = T)[[1]] # MS90
glmThresholded <- glmPred > glmThreshold
glmThresholded <- reclassify(glmThresholded, c(NA, NA, 0), include.lowest = T)

layerNames <- as.numeric(gsub("[X]", "", names(glmThresholded)))
datForMod$index <- unlist(lapply(datForMod$depth, FUN = function(x) which.min(abs(layerNames - x))))
indices <- unique(datForMod$index)
pdf('./figures/africanum_all_depths2.pdf')
plotLayers(glmThresholded[[min(indices):max(indices)]], 
           land = land, landCol = "black", title = "GLM Thresholded by depth for Poroderma africanum")
dev.off()

pdf('./figures/africanum_threshold_by_depth2.pdf')
for (i in seq_along(names(glmThresholded))) {
  plot(glmThresholded[[i]], main = paste0("Poroderma africanum thresholded glm at ", layerNames[i], " m"))
  plot(land, add = T, col = "black")
  plot(africanum_sp, add = T, col = "red")
}
dev.off()

# preliminary model for Poroderma pantherinum
# I'm going to remove all depth layers beyond a 100 m buffer beyond the deepest point at which pantherinum is found
# in order to avoid sampling at depths which are inaccessible to it

max(pantherinum_final$depth)
names(aou_3d)
# cutting off depth sampling for pantherinum at 475 m
aou_3d_pantherinum <- aou_3d[[1:36]]
nitrate_3d_pantherinum <- nitrate_3d[[1:36]]
oxygen_3d_pantherinum <- oxygen_3d[[1:36]]
phosphate_3d_pantherinum <- phosphate_3d[[1:36]]
salinity_3d_pantherinum <- salinity_3d[[1:36]]
temperature_3d_pantherinum <- temperature_3d[[1:36]]

# Gets the layer index for each occurrence by matching to depth
layerNames <- as.numeric(gsub("[X]", "", names(temperature_3d_pantherinum)))
pantherinum_final$index <- unlist(lapply(pantherinum_final$depth, FUN = function(x) 
  which.min(abs(layerNames - x))))
indices <- unique(pantherinum_final$index)

# downsampling occurrences
pantherinum_down <- data.frame()
for(i in indices){
  tempPoints <- pantherinum_final[pantherinum_final$index==i,]
  tempPoints <- downsample(tempPoints, temperature_3d_pantherinum[[1]])
  tempPoints$depth <- rep(layerNames[[i]], times = nrow(tempPoints))
  pantherinum_down <- rbind(pantherinum_down, tempPoints)
}

save(pantherinum_down, file = './data/pantherinum_down.Rdata')

print(paste0("Original number of points: ", nrow(pantherinum_final), "; number of downsampled occs: ", 
             nrow(pantherinum_down)))
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]
pointCompMap(occs1 = pantherinum_final, occs2 = pantherinum_down, 
             occs1Name = "Original", occs2Name = "Downsampled", 
             spName = "Poroderma pantherinum", 
             land = land)

aouVals <- xyzSample(occs = pantherinum_down, aou_3d_pantherinum)
tempVals <- xyzSample(occs = pantherinum_down, temperature_3d_pantherinum)
nitVals <- xyzSample(occs = pantherinum_down, nitrate_3d_pantherinum)
phosVals <- xyzSample(occs = pantherinum_down, phosphate_3d_pantherinum)
salVals <- xyzSample(occs = pantherinum_down, salinity_3d_pantherinum)
oxyVals <- xyzSample(occs = pantherinum_down, oxygen_3d_pantherinum)
vals <- cbind(pantherinum_down, oxyVals, tempVals, aouVals, nitVals, phosVals, salVals)
colnames(vals) <- c("decimalLongitude", "decimalLatitude", "depth", "Oxygen", "Temperature", "AOU", 
                    "Nitrate", "Phosphate", "Salinity")
vals <- vals[complete.cases(vals),]
row.names(vals) <- NULL
occsWdata <- vals

# Add response as a column
occsWdata$response <- rep(1, times = nrow(occsWdata))

# Background
backgroundVals <- mSampling3D(occs = pantherinum_down, 
                              envBrick = temperature_3d_pantherinum, 
                              mShp = pantherinum_accarea, 
                              depthLimit = "all")
aouVals <- xyzSample(occs = backgroundVals, aou_3d_pantherinum)
tempVals <- xyzSample(occs = backgroundVals, temperature_3d_pantherinum)
nitVals <- xyzSample(occs = backgroundVals, nitrate_3d_pantherinum)
phosVals <- xyzSample(occs = backgroundVals, phosphate_3d_pantherinum)
salVals <- xyzSample(occs = backgroundVals, salinity_3d_pantherinum)
oxyVals <- xyzSample(occs = backgroundVals, oxygen_3d_pantherinum)
vals <- cbind(backgroundVals, oxyVals, tempVals, aouVals, nitVals, phosVals, salVals)
colnames(vals) <- c("decimalLongitude", "decimalLatitude", "depth", "Oxygen", "Temperature", "AOU", 
                    "Nitrate", "Phosphate", "Salinity")
vals <- vals[complete.cases(vals),]
row.names(vals) <- NULL
backgroundWdata <- vals
# Add response as a column
backgroundWdata$response <- rep(0, times = nrow(backgroundWdata))

# Sample background points weighted by distance from centroid of occurrence environments
suitableCentroid <- apply(occsWdata[,c("Oxygen", "Temperature", "AOU", "Nitrate", "Phosphate", "Salinity")], 
                          MARGIN = 2, FUN = mean)
backgroundWdata$distance <- apply(backgroundWdata[,c("Oxygen", "Temperature", "AOU", 
                                                     "Nitrate", "Phosphate", "Salinity")], MARGIN = 1, 
                                  FUN = function(x) dist(rbind(suitableCentroid, x)))
backgroundWdata$sampleWeight <- (backgroundWdata$distance - 
                                   min(backgroundWdata$distance))/(max(backgroundWdata$distance)-
                                                                     min(backgroundWdata$distance))
sampleForAbsence <- sample(x = rownames(backgroundWdata), 
                           size = nrow(occsWdata) * 5, 
                           prob = backgroundWdata$sampleWeight)
backgroundWdata <- backgroundWdata[match(sampleForAbsence, 
                                         rownames(backgroundWdata)),]

# Unite datasets
datForMod <- rbind(occsWdata, backgroundWdata[,colnames(occsWdata)])
rm(suitableCentroid, sampleForAbsence, backgroundWdata, occsWdata)

# generate models
glmModel_pantherinum <- glm(formula = response ~ Temperature *  AOU * Oxygen * Nitrate * Phosphate * Salinity, 
                          family = binomial(link = "logit"),  data = datForMod)
summary(glmModel_pantherinum)
save(glmModel_pantherinum, file = './results/glmModel_pantherinum.Rdata')

layerNames <- as.numeric(gsub("[X]", "", names(temperature_3d_pantherinum)))
index <- seq(from = match(min(datForMod$depth), layerNames), 
             to = match(max(datForMod$depth), layerNames), by = 1)
depthPred <- NULL
for(j in index){
  depthPreds <- stack(temperature_3d_pantherinum[[j]], aou_3d_pantherinum[[j]], oxygen_3d_pantherinum[[j]],
                      nitrate_3d_pantherinum[[j]], phosphate_3d_pantherinum[[j]], salinity_3d_pantherinum[[j]])
  crs(depthPreds) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  names(depthPreds) <- c("Temperature", "AOU", "Oxygen", "Nitrate", "Phosphate", "Salinity")
  depthPred[[j]] <- mask(predict(depthPreds, glmModel_pantherinum), pantherinum_accarea)
  depthPred[[j]] <- crop(depthPred[[j]], pantherinum_accarea)
  names(depthPred[[j]]) <- layerNames[[j]]
}
glmPred <- stack(depthPred[!unlist(lapply(depthPred, FUN = function(x) is.null(x)))])

glmThreshold <- quantile(xyzSample(datForMod[datForMod$response == 1,
                                             c("depth", "decimalLongitude", "decimalLatitude")], 
                                   brick(glmPred)), .1, na.rm = T)[[1]] # MS90
glmThresholded <- glmPred > glmThreshold
glmThresholded <- reclassify(glmThresholded, c(NA, NA, 0), include.lowest = T)

layerNames <- as.numeric(gsub("[X]", "", names(glmThresholded)))
datForMod$index <- unlist(lapply(datForMod$depth, FUN = function(x) which.min(abs(layerNames - x))))
indices <- unique(datForMod$index)
pdf('./figures/pantherinum_all_depths.pdf')
plotLayers(glmThresholded[[min(indices):max(indices)]], 
           land = land, landCol = "black", title = "GLM Thresholded by depth for Poroderma pantherinum")
dev.off()

pantherinum_sp <- pantherinum_final
coordinates(pantherinum_sp) <- ~longitude+latitude
proj4string(pantherinum_sp) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

pdf('./figures/pantherinum_threshold_by_depth.pdf')
for (i in seq_along(names(glmThresholded))) {
  plot(glmThresholded[[i]], main = paste0("Poroderma pantherinum thresholded glm at ", layerNames[i], " m"))
  plot(land, add = T, col = "black")
  plot(pantherinum_sp, add = T, col = "red")
}
dev.off()

# generate models with aou and phosphate excluded
glmModel_pantherinum2 <- glm(formula = response ~ Temperature * Oxygen * Nitrate * Salinity, 
                            family = binomial(link = "logit"),  data = datForMod)
summary(glmModel_pantherinum2)
save(glmModel_pantherinum2, file = './results/glmModel_pantherinum2.Rdata')

layerNames <- as.numeric(gsub("[X]", "", names(temperature_3d_pantherinum)))
index <- seq(from = match(min(datForMod$depth), layerNames), 
             to = match(max(datForMod$depth), layerNames), by = 1)
depthPred <- NULL
for(j in index){
  depthPreds <- stack(temperature_3d_pantherinum[[j]], oxygen_3d_pantherinum[[j]],
                      nitrate_3d_pantherinum[[j]], salinity_3d_pantherinum[[j]])
  crs(depthPreds) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  names(depthPreds) <- c("Temperature", "Oxygen", "Nitrate", "Salinity")
  depthPred[[j]] <- mask(predict(depthPreds, glmModel_pantherinum2), pantherinum_accarea)
  depthPred[[j]] <- crop(depthPred[[j]], pantherinum_accarea)
  names(depthPred[[j]]) <- layerNames[[j]]
}
glmPred <- stack(depthPred[!unlist(lapply(depthPred, FUN = function(x) is.null(x)))])

glmThreshold <- quantile(xyzSample(datForMod[datForMod$response == 1,
                                             c("depth", "decimalLongitude", "decimalLatitude")], 
                                   brick(glmPred)), .1, na.rm = T)[[1]] # MS90
glmThresholded <- glmPred > glmThreshold
glmThresholded <- reclassify(glmThresholded, c(NA, NA, 0), include.lowest = T)

layerNames <- as.numeric(gsub("[X]", "", names(glmThresholded)))
datForMod$index <- unlist(lapply(datForMod$depth, FUN = function(x) which.min(abs(layerNames - x))))
indices <- unique(datForMod$index)
pdf('./figures/pantherinum_all_depths2.pdf')
plotLayers(glmThresholded[[min(indices):max(indices)]], 
           land = land, landCol = "black", title = "GLM Thresholded by depth for Poroderma pantherinum")
dev.off()

pantherinum_sp <- pantherinum_final
coordinates(pantherinum_sp) <- ~longitude+latitude
proj4string(pantherinum_sp) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

pdf('./figures/pantherinum_threshold_by_depth2.pdf')
for (i in seq_along(names(glmThresholded))) {
  plot(glmThresholded[[i]], main = paste0("Poroderma pantherinum thresholded glm at ", layerNames[i], " m"))
  plot(land, add = T, col = "black")
  plot(pantherinum_sp, add = T, col = "red")
}
dev.off()