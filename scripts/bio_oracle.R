# preparing bio-oracle layers for maxent
library(sdmpredictors)
library(leaflet)
library(raster)
library(rnaturalearth)

land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]

list_datasets()
list_layers(dataset = "Bio-ORACLE")

# pull wanted layers
bio_temp <- load_layers("BO_sstmean")
bio_salinity <- load_layers("BO_salinity")
bio_nitrate <- load_layers("BO_nitrate")
bio_oxy <- load_layers("BO_dissox")

# reading in WOA data and accesible areas
load('./data/rasters/nitrate_3d.Rdata')
load('./data/rasters/oxygen_3d.Rdata')
load('./data/rasters/salinity_3d.Rdata')
load('./data/rasters/temperature_3d.Rdata')
load('./data/africanum_accarea.Rdata')
load('./data/africanum_accarea.Rdata')
load('./data/africanum_accarea2.Rdata')
load('./data/africanum_accarea3.Rdata')
load('./data/africanum_accarea4.Rdata')
load('./data/africanum_accarea5.Rdata')
load('./data/pantherinum_accarea.Rdata')
load('./data/pantherinum_accarea2.Rdata')
load('./data/pantherinum_accarea3.Rdata')
load('./data/pantherinum_accarea4.Rdata')
load('./data/pantherinum_accarea5.Rdata')

bio_temp <- crop(bio_temp, nitrate_3d[[1]])
save(bio_temp, file = './data/environment/BioOracle/bio_temp.Rdata')
bio_salinity <- crop(bio_salinity, nitrate_3d[[1]])
save(bio_salinity, file = './data/environment/BioOracle/bio_salinity.Rdata')
bio_nitrate <- crop(bio_nitrate, nitrate_3d[[1]])
save(bio_nitrate, file = './data/environment/BioOracle/bio_nitrate.Rdata')
bio_oxy <- crop(bio_oxy, nitrate_3d[[1]])
save(bio_oxy, file = './data/environment/BioOracle/bio_oxy.Rdata')

plot(bio_temp)
plot(land, col = "black", add = T)

plot(temperature_3d[[1]])
plot(land, col = "black", add = T)

plot(bio_salinity)
plot(land, col = "black", add = T)

plot(salinity_3d[[1]])
plot(land, col = "black", add = T)

plot(bio_nitrate)
plot(land, col = "black", add = T)

plot(bio_oxy)
plot(land, col = "black", add = T)

# examining bio-oracle bottom layers
bio_temp_bottom <- raster('./data/environment/BioOracle/Present.Benthic.Max.Depth.Temperature.Mean.asc')
bio_temp_bottom <- crop(bio_temp_bottom, temperature_3d[[1]])

plot(temperature_3d[[38]])
plot(land, col = "black", add = T)

plot(bio_temp_bottom)
plot(land, col = "black", add = T)

# find range of values at each layer