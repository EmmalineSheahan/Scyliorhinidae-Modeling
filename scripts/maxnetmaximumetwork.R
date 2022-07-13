# 3D modeling with maxnet on Poroderma africanum and Poroderma pantherinum

library(raster)
library(dplyr)
library(rgdal)
library(maxnet)
library(rnaturalearth)
library(dismo)
library(ENMeval)
library(voluModel)

# create land for plotting
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")[1]

# load in occurrnces
load('./data/africanum_down.Rdata')
load('./data/pantherinum_down.Rdata')

# load in environmental predictors
load('./data/rasters/nitrate_3d.Rdata')
load('./data/rasters/oxygen_3d.Rdata')
load('./data/rasters/salinity_3d.Rdata')
load('./data/rasters/temperature_3d.Rdata')

# load in accessible area
load('./data/africanum_accarea5.Rdata')
load('./data/pantherinum_accarea5.Rdata')

# separate environmental predictors for each species depending on depth limit
# rasters for africanum
nitrate_3d_africanum <- nitrate_3d[[1:38]]
oxygen_3d_africanum <- oxygen_3d[[1:38]]
salinity_3d_africanum <- salinity_3d[[1:38]]
temperature_3d_africanum <- temperature_3d[[1:38]]

# rasters for pantherinum
nitrate_3d_pantherinum <- nitrate_3d[[1:36]]
oxygen_3d_pantherinum <- oxygen_3d[[1:36]]
salinity_3d_pantherinum <- salinity_3d[[1:36]]
temperature_3d_pantherinum <- temperature_3d[[1:36]]

# crop environmental layers to accessible area
# crop_stack function to crop all layers in stack
# enviro = environmenatl layer brick
# accarea = acessible area polygon
crop_stack <- function(enviro, accarea) {
stack_list <- vector("list", length = dim(enviro)[3])
   for (i in 1:dim(enviro)[3]) {
        enviro[[i]] <- interpolateRaster(enviro[[i]])
        temp <- mask(enviro[[i]], accarea)
        stack_list[[i]] <- temp
   }
t <- stack(stack_list)
return(t)
}

africanum_nitrate <- crop_stack(nitrate_3d_africanum, africanum_accarea5)
africanum_oxygen <- crop_stack(oxygen_3d_africanum, africanum_accarea5)
africanum_salinity <- crop_stack(salinity_3d_africanum, africanum_accarea5)
africanum_temperature <- crop_stack(temperature_3d_africanum, africanum_accarea5)

pantherinum_nitrate <- crop_stack(nitrate_3d_pantherinum, pantherinum_accarea5)
pantherinum_oxygen <- crop_stack(oxygen_3d_pantherinum, pantherinum_accarea5)
pantherinum_salinity <- crop_stack(salinity_3d_pantherinum, pantherinum_accarea5)
pantherinum_temperature <- crop_stack(temperature_3d_pantherinum, pantherinum_accarea5)

# transform occurrence into list of spatial points by depth slice
depth_slices_a <- gsub(names(nitrate_3d_africanum), pattern = 'X', replacement = '')
sp_list_africanum <- vector("list", length = length(depth_slices_a))
for (i in seq_along(depth_slices_a)) {
test <- africanum_down %>% filter(depth == depth_slices_a[i])
if (nrow(test) == 0) {
  test <- NA 
  } else {
coordinates(test) <- ~longitude+latitude
}
sp_list_africanum[[i]] <- test
}

depth_slices_p <- gsub(names(nitrate_3d_pantherinum), pattern = 'X', replacement = '')
sp_list_pantherinum <- vector("list", length = length(depth_slices_p))
for (i in seq_along(depth_slices_p)) {
  test <- pantherinum_down %>% filter(depth == depth_slices_p[i])
  if (nrow(test) == 0) {
    test <- NA 
  } else {
    coordinates(test) <- ~longitude+latitude
  }
  sp_list_pantherinum[[i]] <- test
}

# create background points at each depth slice for each species
plot(africanum_temperature[[38]])
bg_lista <- vector("list", length = 38)
for (i in 1:38) {
bg <- dismo::randomPoints(africanum_temperature[[i]], n = 800) %>% as.data.frame()
coordinates(bg) <- ~x+y
bg_lista[[i]] <- bg
}

bg_lista[[1]]

plot(pantherinum_temperature[[36]])
bg_listp <- vector("list", length = 36)
for (i in 1:36) {
  bg <- dismo::randomPoints(pantherinum_temperature[[i]], n = 800) %>% as.data.frame()
  coordinates(bg) <- ~x+y
  bg_listp[[i]] <- bg
}

# Extract the environmental variables at the occurrence points for each depth value
# brick_extract function outputs vector of environmental values at each coordinate at each depth value
# wanted_brick = raster brick of environmental variables
# wanted_sp = list of coordinates by depth slice
brick_extract <- function(wanted_brick, wanted_sp) {
test_list <- vector("list", length = dim(wanted_brick)[3])
for (i in 1:dim(wanted_brick)[3]) {
test <- extract(wanted_brick[[i]], wanted_sp[[i]])
if (is.null(test) == T) {
  test <- NA
}
test_list[[i]] <- test
}
t <- unlist(test_list)
return(t)
}

# extracting for occs for africanum
temp_af <- brick_extract(africanum_temperature, sp_list_africanum)
nit_af <- brick_extract(africanum_nitrate, sp_list_africanum)
ox_af <- brick_extract(africanum_oxygen, sp_list_africanum)
sal_af <- brick_extract(africanum_salinity, sp_list_africanum)

present <- rep(1, times = length(temp_af))
occs_env_vals <- data.frame(present, temp_af, nit_af, ox_af, sal_af)
colnames(occs_env_vals) <- c("Present", "Temperature", "Nitrate", "Oxygen", "Salinity")

# extracting for background for africanum
temp_afbg <- brick_extract(africanum_temperature, bg_lista)
nit_afbg <- brick_extract(africanum_nitrate, bg_lista)
ox_afbg <- brick_extract(africanum_oxygen, bg_lista)
sal_afbg <- brick_extract(africanum_salinity, bg_lista)

absent <- rep(0, times = length(temp_afbg))
bg_env_vals <- data.frame(absent, temp_afbg, nit_afbg, ox_afbg, sal_afbg)
colnames(bg_env_vals) <- c("Present", "Temperature", "Nitrate", "Oxygen", "Salinity")

env_vals <- rbind(occs_env_vals, bg_env_vals)
env_vals <- env_vals %>% filter(!(is.na(Temperature)))

africanum_env <- env_vals[,2:5]
africanum_pres <- env_vals$Present

# attempting Maxent model for africanum
africanum_maxent <- maxent(x = africanum_env, p = africanum_pres)

# stack predictors by depth slice for predict function
stack_list <- vector("list", length = dim(africanum_temperature)[3])
for (i in 1:dim(africanum_temperature)[3]) {
  st <- stack(africanum_temperature[[i]], africanum_nitrate[[i]], africanum_oxygen[[i]], africanum_salinity[[i]])
  names(st) <- c("Temperature", "Nitrate", "Oxygen", "Salinity")
  stack_list[[i]] <- st
}

surface_test <- predict(africanum_maxent, stack_list[[1]])
plot(surface_test)  
plot(land, add = T, col = "black")

pdf('./figures/africanum_first_maxent.pdf')
for (i in 1:length(stack_list)) {
  pred <- predict(africanum_maxent, stack_list[[i]])
  plot(pred, main = depth_slices_a[[i]])
  plot(land, add = T, col = "black")
  #plot(sp_list_africanum[[i]], col = "red", add = T)
}
dev.off()

# attempt k-fold validation
present_only <- env_vals[1:59,]
africanum_k <- kfold(present_only$Present, k = 5)

africanum_model_list <- vector("list", length = 5)
africanum_eval_list <- vector("list", length = 5)

for (i in (1:5)) {
test <- present_only[africanum_k == i,]
train <- present_only[africanum_k != i,]
test_full <- rbind(test, env_vals[60:5379,])
train_full <- rbind(train, env_vals[60:5379,])
train_model <- maxent(x = train_full[,2:5], p = train_full$Present)
africanum_model_list[[i]] <- train_model
test_full_present <- test_full %>% filter(Present == 1)
test_full_absent <- test_full %>% filter(Present == 0)
ev <- evaluate(p = test_full_present[,2:5], a = test_full_absent[,2:5], model = train_model)
africanum_eval_list[[i]] <- ev
}

save(africanum_model_list, file = './results/africanum_model_list.Rdata')
save(africanum_eval_list, file = './results/africanum_eval_list.Rdata')
save(env_vals, file = './data/env_vals.Rdata')

# Now models for pantherinum
# extracting for occs for pantherinum
temp_pa <- brick_extract(pantherinum_temperature, sp_list_pantherinum)
nit_pa <- brick_extract(pantherinum_nitrate, sp_list_pantherinum)
ox_pa <- brick_extract(pantherinum_oxygen, sp_list_pantherinum)
sal_pa <- brick_extract(pantherinum_salinity, sp_list_pantherinum)

present <- rep(1, times = length(temp_pa))
occs_env_vals <- data.frame(present, temp_pa, nit_pa, ox_pa, sal_pa)
colnames(occs_env_vals) <- c("Present", "Temperature", "Nitrate", "Oxygen", "Salinity")

# extracting for background for pantherinum
temp_pabg <- brick_extract(pantherinum_temperature, bg_listp)
nit_pabg <- brick_extract(pantherinum_nitrate, bg_listp)
ox_pabg <- brick_extract(pantherinum_oxygen, bg_listp)
sal_pabg <- brick_extract(pantherinum_salinity, bg_listp)

absent <- rep(0, times = length(temp_pabg))
bg_env_vals <- data.frame(absent, temp_pabg, nit_pabg, ox_pabg, sal_pabg)
colnames(bg_env_vals) <- c("Present", "Temperature", "Nitrate", "Oxygen", "Salinity")

env_vals_p <- rbind(occs_env_vals, bg_env_vals)
env_vals_p <- env_vals_p %>% filter(!(is.na(Temperature)))

pantherinum_env <- env_vals_p[,2:5]
pantherinum_pres <- env_vals_p$Present

# attempting Maxent model for pantherinum
pantherinum_maxent <- maxent(x = pantherinum_env, p = pantherinum_pres)

# stack predictors by depth slice for predict function
stack_list_p <- vector("list", length = dim(pantherinum_temperature)[3])
for (i in 1:dim(pantherinum_temperature)[3]) {
  st <- stack(pantherinum_temperature[[i]], pantherinum_nitrate[[i]], 
              pantherinum_oxygen[[i]], pantherinum_salinity[[i]])
  names(st) <- c("Temperature", "Nitrate", "Oxygen", "Salinity")
  stack_list_p[[i]] <- st
}

surface_test <- predict(pantherinum_maxent, stack_list_p[[1]])
plot(surface_test)  
plot(land, add = T, col = "black")

pdf('./figures/pantherinum_first_maxent.pdf')
for (i in 1:length(stack_list_p)) {
  pred <- predict(pantherinum_maxent, stack_list_p[[i]])
  plot(pred, main = depth_slices_p[[i]])
  plot(land, add = T, col = "black")
  #plot(sp_list_pantherinum[[i]], col = "red", add = T)
}
dev.off()

# attempt k-fold validation
present_only <- env_vals_p[1:37,]
pantherinum_k <- kfold(present_only$Present, k = 5)

pantherinum_model_list <- vector("list", length = 5)
pantherinum_eval_list <- vector("list", length = 5)

for (i in (1:5)) {
  test <- present_only[pantherinum_k == i,]
  train <- present_only[pantherinum_k != i,]
  test_full <- rbind(test, env_vals_p[38:7057,])
  train_full <- rbind(train, env_vals_p[38:7057,])
  train_model <- maxent(x = train_full[,2:5], p = train_full$Present)
  pantherinum_model_list[[i]] <- train_model
  test_full_present <- test_full %>% filter(Present == 1)
  test_full_absent <- test_full %>% filter(Present == 0)
  ev <- evaluate(p = test_full_present[,2:5], a = test_full_absent[,2:5], model = train_model)
  pantherinum_eval_list[[i]] <- ev
}

save(pantherinum_model_list, file = './results/pantherinum_model_list.Rdata')
save(pantherinum_eval_list, file = './results/pantherinum_eval_list.Rdata')
save(env_vals_p, file = './data/env_evals_p.Rdata')
