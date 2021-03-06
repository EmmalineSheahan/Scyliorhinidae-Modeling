
library(voluModel) # Because of course
library(ggplot2) # For fancy plotting
library(rgdal) # For vector stuff. Will eventually be replaced with sf.
library(raster) # For raster stuff. Will eventually be replaced with terra.
library(viridisLite) # For high-contrast plotting palettes
# Load data
load(system.file("extdata/oxygenSmooth.RData", 
                 package='voluModel'))
td <- tempdir()
unzip(system.file("extdata/woa18_decav_t00mn01_cropped.zip", 
                  package = "voluModel"),
      exdir = paste0(td, "/temperature"), junkpaths = T)
temperature <- readOGR(dsn = paste0(td, "/temperature"), 
                       layer ="woa18_decav_t00mn01_cropped")
unlink(paste0(td, "/temperature"), recursive = T)
temperature@data[temperature@data == -999.999] <- NA
occs <- read.csv(system.file("extdata/Steindachneria_argentea.csv", 
                             package='voluModel'))
# Creating a RasterBrick
temperature <- rasterFromXYZ(cbind(temperature@coords,
                                   temperature@data))
# Get names of depths
envtNames <- gsub("[d,M]", "", names(temperature))
envtNames[[1]] <- "0"
names(temperature) <- envtNames
# Oxygen processing
names(oxygenSmooth) <- names(temperature)
# Clean points ----
occsClean <- occs[complete.cases(occs$depth),]
occsClean <- occsClean[occsClean$depth > 0.0,]
occsClean <- occsClean[occsClean$depth < 2000.0,]
occurrences <- occsClean[,c("decimalLatitude", "decimalLongitude", "depth")] 
# Preliminary cleaning
occurrences <- dplyr::distinct(occurrences)
occurrences <- occurrences[complete.cases(occurrences),]
# Gets the layer index for each occurrence by matching to depth
layerNames <- as.numeric(gsub("[X]", "", names(temperature)))
occurrences$index <- unlist(lapply(occurrences$depth, 
                                   FUN = function(x) which.min(abs(layerNames - x))))
indices <- unique(occurrences$index)
downsampledOccs <- data.frame()
for(i in indices){
  tempPoints <- occurrences[occurrences$index==i,]
  tempPoints <- downsample(tempPoints, temperature[[1]])
  tempPoints$depth <- rep(layerNames[[i]], times = nrow(tempPoints))
  downsampledOccs <- rbind(downsampledOccs, tempPoints)
}
occsWdata <- downsampledOccs[,c("decimalLatitude", "decimalLongitude", "depth")] 
# Extract data ----
occsWdata$temperature <- xyzSample(occs = occsWdata, temperature)
occsWdata$AOU <- xyzSample(occs = occsWdata, oxygenSmooth)
occsWdata <- occsWdata[complete.cases(occsWdata),]
# Study region
studyRegion <- marineBackground(occsWdata)
# Get limits
tempLims <- quantile(occsWdata$temperature,c(0, 1))
aouLims <- quantile(occsWdata$AOU,c(0, 1))
# Reclassify environmental bricks to presence/absence
temperaturePresence <- reclassify(temperature, 
                                  rcl = c(-Inf,tempLims[[1]],0,
                                          tempLims[[1]], tempLims[[2]], 1,
                                          tempLims[[2]], Inf, 0))
AOUpresence <- reclassify(oxygenSmooth, 
                          rcl = c(-Inf, aouLims[[1]],0,
                                  aouLims[[1]], aouLims[[2]], 1,
                                  aouLims[[2]], Inf, 0))
# Put it all together
envelopeModel3D <- temperaturePresence * AOUpresence
envelopeModel3D <- mask(crop(envelopeModel3D, studyRegion), 
                        mask = studyRegion)
names(envelopeModel3D) <- names(temperature)
rm(AOUpresence, downsampledOccs, occsClean, occurrences, temperaturePresence, 
   tempPoints, aouLims, envtNames, i, indices, layerNames, td, tempLims)

land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]

# plotting points
pointMap(occs = occs, land = land, landCol = "black", spName = "Steindachneria argentea", 
         ptSize = 2, ptCol = "orange")

pointCompMap(occs1 = occs, occs1Col = "red", occs1Name = "Raw", 
             occs2 = occsWdata, occs2Col = "orange", occs2Name = "Clean",
             spName = "Steindachneria argentea", agreeCol = "purple",
             land = land, landCol = "black", ptSize = 2)

# Plotting rasters  
oneRasterPlot(rast = temperature[[1]],
              land = land, title = "Sea Surface Temperature, WOA 2018",
              landCol = "black", option = "mako")

rasterComp(rast1 = envelopeModel3D[[1]], rast1Name = "Surface",
           rast2 = envelopeModel3D[[10]], rast2Name = "45m", 
           land = land, landCol = "black", 
           title = "Suitability of Habitat for Luminous Hake\nAt Two Different Depths")

# 3D visualization
layerNames <- as.numeric(gsub("[X]", "", names(envelopeModel3D)))
occsWdata$index <- unlist(lapply(occsWdata$depth, FUN = function(x) which.min(abs(layerNames - x))))
indices <- unique(occsWdata$index)
plotLayers(envelopeModel3D[[min(indices):max(indices)]], 
           title = "Envelope Model of Luminous Hake,\n 20 to 700m",
           land = land, landCol = "black")
