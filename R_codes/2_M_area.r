
# ------------------------------------------------------------------------------------------------ #

### Title: M hypothesis and bias layer ####
### Author: Patron-Rivero, C. ####
### Date: 11/14/2024 ###
### Project: "Ecological and Biogeographical Drivers of Speciation in Neotropical Hognose Pit Vipers, Porthidium (Squamata, Viperidae)" ###

# ------------------------------------------------------------------------------------------------ #

# Packages #

# ------------------------------------------------------------------------------------------------ #

library(rgdal)
library(raster)
library(terra)
library(ellipsenm)

# ------------------------------------------------------------------------------------------------ #

# Inputs #

# ------------------------------------------------------------------------------------------------ #

setwd("E:/Porth_Ecol/")
jjm <- readOGR("./ecoreg/Ecoregions_Morrone.shp")
pred_p <- list.files("E:/Porth_Ecol/env", pattern = ".tif", full.names = TRUE)
stack <- terra::stack(pred_p)
spp <- c("P_arc", "P_dun", "P_hes", "P_lan", "P_nas", "P_oph", "P_por", "P_yuc")
press <- list.files("E:/Porth_Ecol/occ/joint", pattern = ".csv")
press <- press[-1]

# ------------------------------------------------------------------------------------------------ #

# Area M #

# ------------------------------------------------------------------------------------------------ #


for(i in 1:length(press)){
	setwd("E:/Porth_Ecol/occ/joint")
	occ <- read.delim(press[[i]], header = T, sep = ",")
	occ1 <- occ[c(2, 3)]
	occ2 <- SpatialPoints(occ1)
	CRS.new <- CRS ("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
	proj4string(occ2) <- CRS.new
	jjm_subset <- jjm[occ2, ]
	polys1 <- aggregate(jjm_subset, dissolve = TRUE)
	raster::shapefile(polys1, paste0("E:/Porth_Ecol/env/", press[[i]], ".shp"))
	Ms <- mask(crop(stack, jjm_subset), jjm_subset)
	setwd("E:/Porth_Ecol/env/M")
	writeRaster(Ms, filename = paste(spp[[i]], "_", names(stack), sep = ""), bylayer = TRUE, format = "GTiff")
	writeRaster(Ms, filename = paste(spp[[i]], "_", names(stack), sep = ""), bylayer = TRUE, format = "ascii")
	setwd("E:/Porth_Ecol/occ/joint")
}

# ------------------------------------------------------------------------------------------------ #

### EndNotRun
