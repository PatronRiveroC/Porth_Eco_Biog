
# ------------------------------------------------------------------------------------------------ #

### Title: Bias layer and background points ####
### Author: Patron-Rivero, C. ####
### Date: 11/14/2024 ###
### Project: "Ecological and Biogeographical Drivers of Speciation in Neotropical Hognose Pit Vipers, Porthidium (Squamata, Viperidae)" ###

# ------------------------------------------------------------------------------------------------ #

# Packages #

# ------------------------------------------------------------------------------------------------ #

library(terra)
library(dplyr)
library(adehabitatHR)
library(flexsdm)

# ------------------------------------------------------------------------------------------------ #

# Clean data #

# ------------------------------------------------------------------------------------------------ #

P_arc <- rast("E:/Porth_Ecol/env/M/P_arc_bio1.tif")
P_dun <- rast("E:/Porth_Ecol/env/M/P_dun_bio1.tif")
P_hes <- rast("E:/Porth_Ecol/env/M/P_hes_bio1.tif")
P_lan <- rast("E:/Porth_Ecol/env/M/P_lan_bio1.tif")
P_nas <- rast("E:/Porth_Ecol/env/M/P_nas_bio1.tif")
P_oph <- rast("E:/Porth_Ecol/env/M/P_oph_bio1.tif")
P_por <- rast("E:/Porth_Ecol/env/M/P_por_bio1.tif")
P_yuc <- rast("E:/Porth_Ecol/env/M/P_yuc_bio1.tif")

setwd("E:/Porth_Ecol/occ/4_bias")
press <- list.files(pattern = ".csv")
data_list <- list()

for (i in 1:length(press)) {
	
	setwd("E:/Porth_Ecol/occ/4_bias")
		
		if (grepl("\t", readLines(press[[i]], n = 1))) {
			
			data <- read.csv(press[[i]], header = TRUE, sep = "\t", row.names = NULL)
		
		} else if (grepl(",", readLines(press[[i]], n = 1))) {
			
			data <- read.csv(press[[i]], header = TRUE, sep = ",", row.names = NULL)
		
		}
	
		if ("species" %in% colnames(data) && "decimalLongitude" %in% colnames(data) && "decimalLatitude" %in% colnames(data)) {
			
			data <- data[, c("species", "decimalLongitude", "decimalLatitude")]
	
		} else if ("scientificname" %in% colnames(data) && "decimallongitude" %in% colnames(data) && "decimallatitude" %in% colnames(data)) {
			
			data <- data[, c("scientificname", "decimallongitude", "decimallatitude")]
		}
	
	data_1 <- na.omit(data)
	colnames(data_1) <- c("Spp", "Long", "Lat")
	data_2 <- data_1[!duplicated(data_1[, c("Lat", "Long")]), ]
	env <- raster::extract(r_merge, data_2[, c("Long", "Lat")])
	f_occ_env <- data.frame(data_2, env)
	occ_env <- na.omit(f_occ_env)
	data_clean <- occ_env[, c("Spp", "Long", "Lat")]
	data_list[[i]] <- data_clean

}

combined_data <- do.call(rbind, data_list)
env <- terra::extract(r_merge, combined_data[, c("Long", "Lat")])
f_occ_env <- data.frame(combined_data, env)
occ_env <- na.omit(f_occ_env)
data_clean <- occ_env[, c("Spp", "Long", "Lat")]
write.csv(data_clean, "E:/Porth_Ecol/occ/4_bias/bias.csv", row.names = F)

# ------------------------------------------------------------------------------------------------ #

# Inputs #

# ------------------------------------------------------------------------------------------------ #

den <- data_clean[, c("Long", "Lat")]
den <- SpatialPoints(den)
r <- rasterToPoints(r_merge)
r1 <- SpatialPoints(r[, 1:2])
r2 <- SpatialPixels(r1)

# ------------------------------------------------------------------------------------------------ #

# Kernel Utilization Distribution probability #

# ------------------------------------------------------------------------------------------------ #

kde <- kernelUD(den, h = "href", grid = r2)
k <- raster(kde)
projection(k) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
k1 <- crop(k, r_merge)
projection(r_merge) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
r_merge2 <- crop(r_merge, k1)
k3 <- mask(k1, r_merge2)
setwd("E:/Porth_Ecol/env/bias")
writeRaster(k3, filename = "bias.tif", format = "GTiff", overwrite = TRUE)

# ------------------------------------------------------------------------------------------------ #

# Background points #

# ------------------------------------------------------------------------------------------------ #

bg <- function(r_layer, output_file) {
	
	r_layer <- crop(r_layer, k3)
	k_mask <- mask(crop(k3, r_layer), r_layer)
  	bg_data <- sample_background(data = data_clean, x = "Long", y = "Lat", n = 10000, method = "biased", rlayer = r_layer, rbias = k_mask)
	write.csv(bg_data, output_file, row.names = FALSE)

}

bg(P_arc, "E:/Porth_Ecol/env/bias/P_arc.csv")
bg(P_dun, "E:/Porth_Ecol/env/bias/P_dun.csv")
bg(P_hes, "E:/Porth_Ecol/env/bias/P_hes.csv")
bg(P_lan, "E:/Porth_Ecol/env/bias/P_lan.csv")
bg(P_nas, "E:/Porth_Ecol/env/bias/P_nas.csv")
bg(P_oph, "E:/Porth_Ecol/env/bias/P_oph.csv")
bg(P_por, "E:/Porth_Ecol/env/bias/P_por.csv")
bg(P_yuc, "E:/Porth_Ecol/env/bias/P_yuc.csv")

# ------------------------------------------------------------------------------------------------ #

### EndNotRun
