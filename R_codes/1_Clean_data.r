
# ------------------------------------------------------------------------------------------------ #

### Title: Clean databases ####
### Author: Patron-Rivero, C. ####
### Date: 11/14/2024 ###
### Project: "Ecological and Biogeographical Drivers of Speciation in Neotropical Hognose Pit Vipers, Porthidium (Squamata, Viperidae)" ###

# ------------------------------------------------------------------------------------------------ #

# Packages #

# ------------------------------------------------------------------------------------------------ #

library(terra)
library(dplyr)
library(digitize)

# ------------------------------------------------------------------------------------------------ #

# Digitize data #

# ------------------------------------------------------------------------------------------------ #

setwd("E:/Porth_Ecol/occ/digitize")

dig <- list(c("P_arc1.jpg", "P_arc2.jpg", "P_dun1.jpg", "P_dun2.jpg", "P_dun3.jpg", "P_hes1.jpg", "P_lan1.jpg",
				"P_lan2.jpg", "P_lan3.jpg", "P_lan4.jpg", "P_lan5.jpg", "P_lan6.jpg", "P_lan7.jpg", "P_nas1.jpg",
				"P_nas2.jpg", "P_nas3.jpg", "P_nas4.jpg", "P_nas5.jpg", "P_nas6.jpg", "P_nas7.jpg", "P_nas8.jpg",
				"P_nas9.jpg", "P_oph1.jpg", "P_oph2.jpg", "P_oph3.jpg", "P_vol1.jpg")
dig_all <- list()

for(i in 1:length(dig) {

	dig_all[i] <- digitize(dig[i])
	
}	

combined_data <- do.call(rbind, dig_all)
write.csv(combined_data,"E:/Porth_Ecol/occ/raw/digitize.csv", row.names = FALSE)

# ------------------------------------------------------------------------------------------------ #

# Inputs #

# ------------------------------------------------------------------------------------------------ #

pred <- list.files("E:/Porth_Ecol/env", pattern = ".tif", full.names = TRUE)
stack <- terra::stack(pred)
setwd("E:/Porth_Ecol/occ/raw")
press <- list.files(pattern = ".csv")
data_list <- list()

# ------------------------------------------------------------------------------------------------ #

# Cleaning duplicated and obvious errors #

# ------------------------------------------------------------------------------------------------ #

for (i in 1:length(press)) {
	
	setwd("E:/Porth_Ecol/occ/raw")
		
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
	env <- raster::extract(stack, data_2[, c("Long", "Lat")])
	f_occ_env <- data.frame(data_2, env)
	occ_env <- na.omit(f_occ_env)
	data_clean <- occ_env[, c("Spp", "Long", "Lat")]
	data_list[[i]] <- data_clean

}

combined_data <- do.call(rbind, data_list)
data_filtered <- anti_join(combined_data, test, by = c("Long", "Lat"))
write.csv(data_filtered, "E:/Porth_Ecol/occ/joint/clean.csv", row.names = F)

# ------------------------------------------------------------------------------------------------ #

# Cleaning statistical outiers #

# ------------------------------------------------------------------------------------------------ #

data <- split(data_filtered, data_filtered$Spp)
spp <- c("P_arc", "P_dun", "P_hes", "P_lan", "P_nas", "P_oph", "P_por", "P_yuc")

for(i in 1:length(data)){
	
	occ_1 <- data[[i]]
	dups <- rep(0, nrow(occ_1))
	occ_e <- raster::extract(stack, occ_1[, 2:3], df = TRUE)
	full <- occ_1
		
		for (j in 2:13) {
			
			outliers <- outliers(occ_e$ID, occ_1$Spp, dups, occ_e[, j])
			full <- cbind(full, outliers)
			
		}

	column_indices <- 4:27
	final <- full
	
		for (k in column_indices) {
		
			final <- final[!(final[, k] == 1), ]
		
		}

	setwd("E:/Porth_Ecol/occ/joint")
	write.csv(final[, 1:3], spp[[i]], row.names = F)

}

# ------------------------------------------------------------------------------------------------ #

### EndNotRun
