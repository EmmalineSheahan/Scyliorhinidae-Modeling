library(utils)
library(raster)
library(voluModel)
library(dplyr)
library(rgdal)
library(sp)
library(sf)
library(rnaturalearth)

untar('./data/environment/woa18_all_A00mn01_shape.tar.gz', exdir = './data/environment/')

# investigate AOU data set from world ocean atlas
aou_raw <- shapefile('./data/environment/woa18_all_A00mn01.shp')

land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]

View(aou_raw@data)
aou_raw@data[aou_raw@data == -999.999] <- NA

# What does the WOA depth structure look like?
depths <- colnames(aou_raw@data)
depths <- as.numeric(gsub(depths[-1], pattern = "[d,M]", replacement = ""))
plot(0, xlim = c(0,1), ylim = c(0-max(depths), 0), axes=FALSE, type = "n", xlab = "", ylab = "Depth Intervals (m)")
axis(2, at = 0-depths, labels = depths)

#Creating a bottom raster from the point shapefile
aouBottom <- bottomRaster(aou_raw)
plot(aouBottom)
plot(land, add = T)

# Creating a 3D temperature RasterBrick from the point shapefile
aou <- rasterFromXYZ(cbind(aou_raw@coords,
                                   aou_raw@data))
# Get names of depths
envtNames <- gsub("[d,M]", "", names(aou_raw))
envtNames[[1]] <- "0"
names(aou) <- envtNames
# How do these files look?
p1 <- oneRasterPlot(aou[[1]], land = land, landCol = "black", 
                    title= "Apparent oxygen utilization")
# R won't let you name the layers numbers only--it pastes an X at the beginning
p2 <- oneRasterPlot(aou[["X100"]], land = land, landCol = "black", 
                    title= "apparent oxygen utilization")
p3 <- oneRasterPlot(aouBottom,land = land, landCol = "black", 
                    title = "apparent oxygen utilization")

# I'm going to create alpha hull buffers around the two species to 
# see what extent I need to cut the environmental rasters to
load('./data/africanum_final.Rdata')
load('./data/pantherinum_final.Rdata')

africanum_accarea <- marineBackground(africanum_final)
save(africanum_accarea, file = './data/africanum_accarea.Rdata')
pantherinum_accarea <- marineBackground(pantherinum_final)
save(pantherinum_accarea, file = './data/pantherinum_accarea.Rdata')

plot(africanum_accarea, col = "red")
plot(pantherinum_accarea, col = "green", add = T)
plot(land, add = T, col = "black")

extent(africanum_accarea)
extent(pantherinum_accarea)

ext_cutx <- c(8, 40)
ext_cuty <- c(-45, -20)
ext_cut <- rbind(ext_cutx, ext_cuty)
ext_cut <- extent(ext_cut)

aou_raw_cut <- crop(aou_raw, ext_cut)
plot(aou_raw_cut, col = "red")
plot(land, add = T)

# build function to clean and prepare rasters
clean_env <- function(wanted_file) {
untar(paste0('./data/environment/woa18_', wanted_file, '00mn01_shape.tar.gz'), exdir = './data/environment/')
temp_spdf <- shapefile(paste0('./data/environment/woa18_', wanted_file, '00mn01.shp'))
temp_spdf@data[temp_spdf@data == -999.999] <- NA
envtNames <- gsub("[d,M]", "", names(temp_spdf))
envtNames[[1]] <- "0"
names(temp_spdf) <- envtNames
temp_spdf <- crop(temp_spdf, ext_cut)
plot(temp_spdf, col = "red")
plot(land, add = T)
return(temp_spdf)
}

# prepare apparent oxygen utilization
aou_raw <- clean_env('all_A')
aou_Bottom <- bottomRaster(aou_raw)
aou_3d <- rasterFromXYZ(cbind(aou_raw@coords, aou_raw@data))
aou_Surface <- aou_3d[[1]]
save(aou_3d, file = './data/rasters/aou_3d.Rdata')
save(aou_Bottom, file = './data/rasters/aou_Bottom.Rdata')
save(aou_Surface, file = './data/rasters/aou_Surface.Rdata')

# prepare nitrate
nitrate_raw <- clean_env('all_n')
nitrate_Bottom <- bottomRaster(nitrate_raw)
nitrate_3d <- rasterFromXYZ(cbind(nitrate_raw@coords, nitrate_raw@data))
nitrate_Surface <- nitrate_3d[[1]]
save(nitrate_3d, file = './data/rasters/nitrate_3d.Rdata')
save(nitrate_Bottom, file = './data/rasters/nitrate_Bottom.Rdata')
save(nitrate_Surface, file = './data/rasters/nitrate_Surface.Rdata')

# prepare dissolved oxygen
oxygen_raw <- clean_env('all_o')
oxygen_Bottom <- bottomRaster(oxygen_raw)
oxygen_3d <- rasterFromXYZ(cbind(oxygen_raw@coords, oxygen_raw@data))
oxygen_Surface <- oxygen_3d[[1]]
save(oxygen_3d, file = './data/rasters/oxygen_3d.Rdata')
save(oxygen_Bottom, file = './data/rasters/oxygen_Bottom.Rdata')
save(oxygen_Surface, file = './data/rasters/oxygen_Surface.Rdata')

# prepare phosphate
phosphate_raw <- clean_env('all_p')
phosphate_Bottom <- bottomRaster(phosphate_raw)
phosphate_3d <- rasterFromXYZ(cbind(phosphate_raw@coords, phosphate_raw@data))
phosphate_Surface <- phosphate_3d[[1]]
save(phosphate_3d, file = './data/rasters/phosphate_3d.Rdata')
save(phosphate_Bottom, file = './data/rasters/phosphate_Bottom.Rdata')
save(phosphate_Surface, file = './data/rasters/phosphate_Surface.Rdata')

# prepare salinity
salinity_raw <- clean_env('decav_s')
salinity_Bottom <- bottomRaster(salinity_raw)
salinity_3d <- rasterFromXYZ(cbind(salinity_raw@coords, salinity_raw@data))
salinity_Surface <- salinity_3d[[1]]
save(salinity_3d, file = './data/rasters/salinity_3d.Rdata')
save(salinity_Bottom, file = './data/rasters/salinity_Bottom.Rdata')
save(salinity_Surface, file = './data/rasters/salinity_Surface.Rdata')

# prepare temperature
temperature_raw <- clean_env('decav_t')
temperature_Bottom <- bottomRaster(temperature_raw)
temperature_3d <- rasterFromXYZ(cbind(temperature_raw@coords, temperature_raw@data))
temperature_Surface <- temperature_3d[[1]]
save(temperature_3d, file = './data/rasters/temperature_3d.Rdata')
save(temperature_Bottom, file = './data/rasters/temperature_Bottom.Rdata')
save(temperature_Surface, file = './data/rasters/temperature_Surface.Rdata')
