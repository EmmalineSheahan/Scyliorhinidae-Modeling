# preparing rasters and occurrences to run 2D model in maxent GUI

library(voluModel)
library(raster)
library(dplyr)
library(rnaturalearth)
library(terra)

load('./data/africanum_down.Rdata')
load('./data/pantherinum_down.Rdata')
load('./data/africanum_accarea.Rdata')
load('./data/pantherinum_accarea.Rdata')
load('./data/rasters/aou_3d.Rdata')
load('./data/rasters/nitrate_3d.Rdata')
load('./data/rasters/oxygen_3d.Rdata')
load('./data/rasters/phosphate_3d.Rdata')
load('./data/rasters/salinity_3d.Rdata')
load('./data/rasters/temperature_3d.Rdata')

land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]

# new accessible area
load('./data/africanum_final.Rdata')
load('./data/pantherinum_final.Rdata')

pDist <- terra::distance(cbind(africanum_final$longitude, africanum_final$latitude), lonlat = TRUE)
buff1 <- mean(quantile(pDist, c(.1, .9), na.rm = TRUE))/2

# standard buffer of marineBackground on africanum is 227753.8

pDist <- terra::distance(cbind(pantherinum_final$longitude, pantherinum_final$latitude), lonlat = TRUE)
buff2 <- mean(quantile(pDist, c(.1, .9), na.rm = TRUE))/2

# standard buffer of marineBackground on pantherinum is 329091.1

# creating 2 new accessible areas for each, one with a 20% larger buffer and one with a 40% larger buffer
africanum_accarea2 <- marineBackground(africanum_final, buff = (buff1+(buff1*0.2)))
plot(africanum_accarea2)
plot(land, add = T)
save(africanum_accarea2, file = './data/africanum_accarea2.Rdata')

africanum_accarea3 <- marineBackground(africanum_final, buff = (buff1+(buff1*0.4)))
save(africanum_accarea3, file = './data/africanum_accarea3.Rdata')

africanum_accarea4 <- marineBackground(africanum_final, buff = (buff1+(buff1*0.6)))
save(africanum_accarea4, file = './data/africanum_accarea4.Rdata')

africanum_accarea5 <- marineBackground(africanum_final, buff = (buff1+(buff1*0.8)))
save(africanum_accarea5, file = './data/africanum_accarea5.Rdata')

pantherinum_accarea2 <- marineBackground(pantherinum_final, buff = (buff2+(buff2*0.2)))
save(pantherinum_accarea2, file = './data/pantherinum_accarea2.Rdata')

pantherinum_accarea3 <- marineBackground(pantherinum_final, buff = (buff2+(buff2*0.4)))
save(pantherinum_accarea3, file = './data/pantherinum_accarea3.Rdata')

pantherinum_accarea4 <- marineBackground(pantherinum_final, buff = (buff2+(buff2*0.6)))
save(pantherinum_accarea4, file = './data/pantherinum_accarea4.Rdata')

pantherinum_accarea5 <- marineBackground(pantherinum_final, buff = (buff2+(buff2*0.8)))
save(pantherinum_accarea5, file = './data/pantherinum_accarea5.Rdata')

# finish_raster function
finish_raster <- function(wanted_brick, wanted_mask, file_name) {
  ras_mean <- calc(wanted_brick, fun = mean, na.rm = T)
  ras_mean <- mask(ras_mean, mask = wanted_mask)
  plot(ras_mean)
  plot(land, add = T, col = "black")
  writeRaster(x = ras_mean, filename = file_name, format = "ascii", overwrite = T) 
}

# Mean rasters for africanum
aou_3d_africanum <- aou_3d[[1:38]]
nitrate_3d_africanum <- nitrate_3d[[1:38]]
oxygen_3d_africanum <- oxygen_3d[[1:38]]
phosphate_3d_africanum <- phosphate_3d[[1:38]]
salinity_3d_africanum <- salinity_3d[[1:38]]
temperature_3d_africanum <- temperature_3d[[1:38]]

aou_africanum_mean <- finish_raster(aou_3d_africanum, africanum_accarea, 
                                    './data/rasters/for_maxent/africanum/aou')
nitrate_africanum_mean <- finish_raster(nitrate_3d_africanum, africanum_accarea,
                                        './data/rasters/for_maxent/africanum/nitrate')
oxygen_africanum_mean <- finish_raster(oxygen_3d_africanum, africanum_accarea,
                                        './data/rasters/for_maxent/africanum/oxygen')
phosphate_africanum_mean <- finish_raster(phosphate_3d_africanum, africanum_accarea,
                                        './data/rasters/for_maxent/africanum/phosphate')
salinity_africanum_mean <- finish_raster(salinity_3d_africanum, africanum_accarea,
                                        './data/rasters/for_maxent/africanum/salinity')
temperature_africanum_mean <- finish_raster(temperature_3d_africanum, africanum_accarea,
                                        './data/rasters/for_maxent/africanum/temperature')

# Mean rasters for africanum 20% buffer
aou_africanum_mean2 <- finish_raster(aou_3d_africanum, africanum_accarea2, 
                                    './data/rasters/for_maxent/africanum2/aou')
nitrate_africanum_mean2 <- finish_raster(nitrate_3d_africanum, africanum_accarea2,
                                        './data/rasters/for_maxent/africanum2/nitrate')
oxygen_africanum_mean2 <- finish_raster(oxygen_3d_africanum, africanum_accarea2,
                                       './data/rasters/for_maxent/africanum2/oxygen')
phosphate_africanum_mean2 <- finish_raster(phosphate_3d_africanum, africanum_accarea2,
                                          './data/rasters/for_maxent/africanum2/phosphate')
salinity_africanum_mean2 <- finish_raster(salinity_3d_africanum, africanum_accarea2,
                                         './data/rasters/for_maxent/africanum2/salinity')
temperature_africanum_mean2 <- finish_raster(temperature_3d_africanum, africanum_accarea2,
                                            './data/rasters/for_maxent/africanum2/temperature')

# Mean rasters for africanum 40% buffer
aou_africanum_mean3 <- finish_raster(aou_3d_africanum, africanum_accarea3, 
                                     './data/rasters/for_maxent/africanum3/aou')
nitrate_africanum_mean3 <- finish_raster(nitrate_3d_africanum, africanum_accarea3,
                                         './data/rasters/for_maxent/africanum3/nitrate')
oxygen_africanum_mean3 <- finish_raster(oxygen_3d_africanum, africanum_accarea3,
                                        './data/rasters/for_maxent/africanum3/oxygen')
phosphate_africanum_mean3 <- finish_raster(phosphate_3d_africanum, africanum_accarea3,
                                           './data/rasters/for_maxent/africanum3/phosphate')
salinity_africanum_mean3 <- finish_raster(salinity_3d_africanum, africanum_accarea3,
                                          './data/rasters/for_maxent/africanum3/salinity')
temperature_africanum_mean3 <- finish_raster(temperature_3d_africanum, africanum_accarea3,
                                             './data/rasters/for_maxent/africanum3/temperature')

# Mean rasters for africanum 60% buffer
aou_africanum_mean4 <- finish_raster(aou_3d_africanum, africanum_accarea4, 
                                     './data/rasters/for_maxent/africanum4/aou')
nitrate_africanum_mean4 <- finish_raster(nitrate_3d_africanum, africanum_accarea4,
                                         './data/rasters/for_maxent/africanum4/nitrate')
oxygen_africanum_mean4 <- finish_raster(oxygen_3d_africanum, africanum_accarea4,
                                        './data/rasters/for_maxent/africanum4/oxygen')
phosphate_africanum_mean4 <- finish_raster(phosphate_3d_africanum, africanum_accarea4,
                                           './data/rasters/for_maxent/africanum4/phosphate')
salinity_africanum_mean4 <- finish_raster(salinity_3d_africanum, africanum_accarea4,
                                          './data/rasters/for_maxent/africanum4/salinity')
temperature_africanum_mean4 <- finish_raster(temperature_3d_africanum, africanum_accarea4,
                                             './data/rasters/for_maxent/africanum4/temperature')

# Mean rasters for africanum 80% buffer
aou_africanum_mean5 <- finish_raster(aou_3d_africanum, africanum_accarea5, 
                                     './data/rasters/for_maxent/africanum5/aou')
nitrate_africanum_mean5 <- finish_raster(nitrate_3d_africanum, africanum_accarea5,
                                         './data/rasters/for_maxent/africanum5/nitrate')
oxygen_africanum_mean5 <- finish_raster(oxygen_3d_africanum, africanum_accarea5,
                                        './data/rasters/for_maxent/africanum5/oxygen')
phosphate_africanum_mean5 <- finish_raster(phosphate_3d_africanum, africanum_accarea5,
                                           './data/rasters/for_maxent/africanum5/phosphate')
salinity_africanum_mean5 <- finish_raster(salinity_3d_africanum, africanum_accarea5,
                                          './data/rasters/for_maxent/africanum5/salinity')
temperature_africanum_mean5 <- finish_raster(temperature_3d_africanum, africanum_accarea5,
                                             './data/rasters/for_maxent/africanum5/temperature')


# occurrences for africanum
africanum_maxent <- africanum_down[,1:2]
nrow(africanum_maxent)
species_name <- rep("Poroderma africanum", times = 73)
africanum_maxent <- cbind(species_name, africanum_maxent)
colnames(africanum_maxent) <- c("species name", "longitude", "latitude")
write.csv(africanum_maxent, file = './data/rasters/for_maxent/africanum/africanum_maxent.csv')

# Mean rasters for pantherinum
aou_3d_pantherinum <- aou_3d[[1:36]]
nitrate_3d_pantherinum <- nitrate_3d[[1:36]]
oxygen_3d_pantherinum <- oxygen_3d[[1:36]]
phosphate_3d_pantherinum <- phosphate_3d[[1:36]]
salinity_3d_pantherinum <- salinity_3d[[1:36]]
temperature_3d_pantherinum <- temperature_3d[[1:36]]

aou_pantherinum_mean <- finish_raster(aou_3d_pantherinum, pantherinum_accarea, 
                                    './data/rasters/for_maxent/pantherinum/aou')
nitrate_pantherinum_mean <- finish_raster(nitrate_3d_pantherinum, pantherinum_accarea,
                                        './data/rasters/for_maxent/pantherinum/nitrate')
oxygen_pantherinum_mean <- finish_raster(oxygen_3d_pantherinum, pantherinum_accarea,
                                       './data/rasters/for_maxent/pantherinum/oxygen')
phosphate_pantherinum_mean <- finish_raster(phosphate_3d_pantherinum, pantherinum_accarea,
                                          './data/rasters/for_maxent/pantherinum/phosphate')
salinity_pantherinum_mean <- finish_raster(salinity_3d_pantherinum, pantherinum_accarea,
                                         './data/rasters/for_maxent/pantherinum/salinity')
temperature_pantherinum_mean <- finish_raster(temperature_3d_pantherinum, pantherinum_accarea,
                                            './data/rasters/for_maxent/pantherinum/temperature')

# mean rasters for pantherinum for 20% buffer
aou_pantherinum_mean2 <- finish_raster(aou_3d_pantherinum, pantherinum_accarea2, 
                                      './data/rasters/for_maxent/pantherinum2/aou')
nitrate_pantherinum_mean2 <- finish_raster(nitrate_3d_pantherinum, pantherinum_accarea2,
                                          './data/rasters/for_maxent/pantherinum2/nitrate')
oxygen_pantherinum_mean2 <- finish_raster(oxygen_3d_pantherinum, pantherinum_accarea2,
                                         './data/rasters/for_maxent/pantherinum2/oxygen')
phosphate_pantherinum_mean2 <- finish_raster(phosphate_3d_pantherinum, pantherinum_accarea2,
                                            './data/rasters/for_maxent/pantherinum2/phosphate')
salinity_pantherinum_mean2 <- finish_raster(salinity_3d_pantherinum, pantherinum_accarea2,
                                           './data/rasters/for_maxent/pantherinum2/salinity')
temperature_pantherinum_mean2 <- finish_raster(temperature_3d_pantherinum, pantherinum_accarea2,
                                              './data/rasters/for_maxent/pantherinum2/temperature')

# mean rasters for pantherinum for 40% buffer
aou_pantherinum_mean3 <- finish_raster(aou_3d_pantherinum, pantherinum_accarea3, 
                                       './data/rasters/for_maxent/pantherinum3/aou')
nitrate_pantherinum_mean3 <- finish_raster(nitrate_3d_pantherinum, pantherinum_accarea3,
                                           './data/rasters/for_maxent/pantherinum3/nitrate')
oxygen_pantherinum_mean3 <- finish_raster(oxygen_3d_pantherinum, pantherinum_accarea3,
                                          './data/rasters/for_maxent/pantherinum3/oxygen')
phosphate_pantherinum_mean3 <- finish_raster(phosphate_3d_pantherinum, pantherinum_accarea3,
                                             './data/rasters/for_maxent/pantherinum3/phosphate')
salinity_pantherinum_mean3 <- finish_raster(salinity_3d_pantherinum, pantherinum_accarea3,
                                            './data/rasters/for_maxent/pantherinum3/salinity')
temperature_pantherinum_mean3 <- finish_raster(temperature_3d_pantherinum, pantherinum_accarea3,
                                               './data/rasters/for_maxent/pantherinum3/temperature')
# mean rasters for pantherinum for 60% buffer
aou_pantherinum_mean4 <- finish_raster(aou_3d_pantherinum, pantherinum_accarea4, 
                                       './data/rasters/for_maxent/pantherinum4/aou')
nitrate_pantherinum_mean4 <- finish_raster(nitrate_3d_pantherinum, pantherinum_accarea4,
                                           './data/rasters/for_maxent/pantherinum4/nitrate')
oxygen_pantherinum_mean4 <- finish_raster(oxygen_3d_pantherinum, pantherinum_accarea4,
                                          './data/rasters/for_maxent/pantherinum4/oxygen')
phosphate_pantherinum_mean4 <- finish_raster(phosphate_3d_pantherinum, pantherinum_accarea4,
                                             './data/rasters/for_maxent/pantherinum4/phosphate')
salinity_pantherinum_mean4 <- finish_raster(salinity_3d_pantherinum, pantherinum_accarea3,
                                            './data/rasters/for_maxent/pantherinum4/salinity')
temperature_pantherinum_mean4 <- finish_raster(temperature_3d_pantherinum, pantherinum_accarea4,
                                               './data/rasters/for_maxent/pantherinum4/temperature')

# mean rasters for pantherinum for 80% buffer
aou_pantherinum_mean5 <- finish_raster(aou_3d_pantherinum, pantherinum_accarea5, 
                                       './data/rasters/for_maxent/pantherinum5/aou')
nitrate_pantherinum_mean5 <- finish_raster(nitrate_3d_pantherinum, pantherinum_accarea5,
                                           './data/rasters/for_maxent/pantherinum5/nitrate')
oxygen_pantherinum_mean5 <- finish_raster(oxygen_3d_pantherinum, pantherinum_accarea5,
                                          './data/rasters/for_maxent/pantherinum5/oxygen')
phosphate_pantherinum_mean5 <- finish_raster(phosphate_3d_pantherinum, pantherinum_accarea5,
                                             './data/rasters/for_maxent/pantherinum5/phosphate')
salinity_pantherinum_mean5 <- finish_raster(salinity_3d_pantherinum, pantherinum_accarea5,
                                            './data/rasters/for_maxent/pantherinum5/salinity')
temperature_pantherinum_mean5 <- finish_raster(temperature_3d_pantherinum, pantherinum_accarea5,
                                               './data/rasters/for_maxent/pantherinum5/temperature')

# occurrences for pantherinum
pantherinum_maxent <- pantherinum_down[,1:2]
nrow(pantherinum_maxent)
species_name <- rep("Poroderma pantherinum", times = 46)
pantherinum_maxent <- cbind(species_name, pantherinum_maxent)
colnames(pantherinum_maxent) <- c("species name", "longitude", "latitude")
write.csv(pantherinum_maxent, file = './data/rasters/for_maxent/pantherinum/pantherinum_maxent.csv')


