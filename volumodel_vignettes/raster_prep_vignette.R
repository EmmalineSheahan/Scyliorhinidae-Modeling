library(voluModel) # Since this is the package this vignette is about.
library(rgdal) # For vector stuff. Will eventually be replaced with sf.
library(raster) # For raster stuff. Will eventually be replaced with terra.
library(tibble) # For data organization
library(ggplot2) # For supplementary visualization
library(fields) # For raster interpolation
library(latticeExtra) # Some fancy plotting

load(system.file("extdata/oxygenSmooth.RData", 
                 package='voluModel'))

td <- tempdir()
unzip(system.file("extdata/woa18_decav_t00mn01_cropped.zip", 
                  package = "voluModel"),
      exdir = paste0(td, "/temperature"), junkpaths = T)
temperature <- readOGR(dsn = paste0(td, "/temperature"), 
                       layer ="woa18_decav_t00mn01_cropped")
unlink(paste0(td, "/temperature"), recursive = T)
# Specifying "no data" value
temperature@data[temperature@data == -999.999] <- NA
# Looking at the dataset
head(temperature@data)
# Plotting the dataset
layout(matrix(c(1, 2), ncol=2, byrow=TRUE), widths=c(4, 1))
land <- rnaturalearth::ne_countries(scale = "medium", 
                                    returnclass = "sf")[1]
ext <- extent(temperature@coords)
plot(temperature, main = "Distribution of voluModel Subset\nof WOA Temperature 2018",
     pch = 20, col = "red", xlim = ext[1:2], ylim = ext[3:4], cex = .6)
plot(land, col = "black", add = T)


# What does the WOA depth structure look like?
depths <- colnames(temperature@data)
depths <- as.numeric(gsub(depths[-1], pattern = "[d,M]", replacement = ""))
plot(0, xlim = c(0,1), ylim = c(0-max(depths), 0), axes=FALSE, type = "n", xlab = "", ylab = "Depth Intervals (m)")
axis(2, at = 0-depths, labels = depths)

#Creating a bottom raster from the point shapefile
temperatureBottom <- bottomRaster(temperature)
# Creating a 3D temperature RasterBrick from the point shapefile
temperature <- rasterFromXYZ(cbind(temperature@coords,
                                   temperature@data))
# Get names of depths
envtNames <- gsub("[d,M]", "", names(temperature))
envtNames[[1]] <- "0"
names(temperature) <- envtNames
# How do these files look?
p1 <- oneRasterPlot(temperature[[1]], land = land, landCol = "black", 
                    title= "Temperature (C)")
# R won't let you name the layers numbers only--it pastes an X at the beginning
p2 <- oneRasterPlot(temperature[["X100"]], land = land, landCol = "black", 
                    title= "Temperature (C)")
p3 <- oneRasterPlot(temperatureBottom,land = land, landCol = "black", 
                    title = "Temperature (C)")



temp <- c("Surface" = p1, "100m Depth" = p2, "Bottom" = p3)
update(temp, strip  = strip.custom(strip.levels = TRUE,
                                   horizontal = TRUE,
                                   bg = "black",
                                   fg = "white",
                                   par.strip.text = list(col = "white", cex = 1.2, font = 2)))

# AOU interpolation
td <- tempdir()
unzip(system.file("extdata/woa18_all_A00mn01_cropped.zip", 
                  package = "voluModel"),
      exdir = paste0(td, "/oxygen"), junkpaths = T)
# do something with the files
oxygen <- readOGR(dsn = paste0(td, "/oxygen"), 
                  layer = "woa18_all_A00mn01_cropped")
unlink(paste0(td, "/oxygen"), recursive = T)
oxygen@data[oxygen@data == -999.999] <- NA
plot(oxygen, main = "Distribution of voluModel subset of WOA AOU 2018",
     pch = 20, col = "red", xlim = ext[1:2], ylim = ext[3:4], cex = .6)
plot(land, col = "black", add = T)

# Creating a RasterBrick
oxygen <- rasterFromXYZ(cbind(oxygen@coords, oxygen@data))
for (i in 1:nlayers(oxygen)){
  oxygen[[i]] <- interpolateRaster(oxygen[[i]], lon.lat = T, fast = T, aRange = 5) #Thin plate spline interpolation
  oxygen[[i]] <- crop(mask(x = oxygen[[i]], mask = temperature[[i]]), temperature[[i]])
}
# Change names to match temperature
names(oxygen) <- envtNames
# How do these files look?
p1 <- oneRasterPlot(oxygen[[1]], land = land, landCol = "black", 
                    title= "Apparent Oxygen Utilization (µmol/kg), interpolated")
p2 <- oneRasterPlot(oxygen[["X100"]], land = land, landCol = "black", 
                    title= "Apparent Oxygen Utilization (µmol/kg), interpolated")
temp <- c("Surface" = p1, "100m Depth" = p2)
update(temp, strip  = strip.custom(strip.levels = TRUE,
                                   horizontal = TRUE,
                                   bg = "black",
                                   fg = "white",
                                   par.strip.text = list(col = "white", cex = 1.2, font = 2)))

oxygenSmooth <- oxygen
for (i in 1:nlayers(oxygen)){
  oxygenSmooth[[i]] <- smoothRaster(oxygenSmooth[[i]], lon.lat = T) #Thin plate spline interpolation
  oxygenSmooth[[i]] <- crop(mask(x = oxygenSmooth[[i]], mask = temperature[[i]]), temperature[[i]])
}
# Change names to match temperature
names(oxygenSmooth) <- names(temperature)
p1 <- oneRasterPlot(oxygenSmooth[[1]], land = land, landCol = "black", 
                    title= "Apparent Oxygen Utilization (µmol/kg), interpolated and smoothed")
p2 <- oneRasterPlot(oxygenSmooth[["X100"]], land = land, landCol = "black", 
                    title= "Apparent Oxygen Utilization (µmol/kg), interpolated and smoothed")
temp <- c("Surface" = p1, "100m Depth" = p2)
update(temp, strip  = strip.custom(strip.levels = TRUE,
                                   horizontal = TRUE,
                                   bg = "black",
                                   fg = "white",
                                   par.strip.text = list(col = "white", cex = 1.2, font = 2)))

p1 <- oneRasterPlot(oxygenSmooth[[1]], land = land, landCol = "black", 
                    title= "Apparent Oxygen Utilization (µmol/kg), interpolated and smoothed")
p2 <- oneRasterPlot(oxygenSmooth[["X100"]], land = land, landCol = "black", 
                    title= "Apparent Oxygen Utilization (µmol/kg), interpolated and smoothed")
temp <- c("Surface" = p1, "100m Depth" = p2)
update(temp, strip  = strip.custom(strip.levels = TRUE,
                                   horizontal = TRUE,
                                   bg = "black",
                                   fg = "white",
                                   par.strip.text = list(col = "white", cex = 1.2, font = 2)))
