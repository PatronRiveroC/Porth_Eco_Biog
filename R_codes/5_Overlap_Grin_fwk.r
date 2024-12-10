
# ------------------------------------------------------------------------------------------------ #

### Title: Overlap test by Grinnellian framework ####
### Author: Patrón-Rivero, C. ####
### Date: 11/14/2024 ###
### Project: "Ecological and Biogeographical Drivers of Speciation in Neotropical Hognose Pit Vipers, Porthidium (Squamata, Viperidae)" ###

# ------------------------------------------------------------------------------------------------ #

# Packages #

# ------------------------------------------------------------------------------------------------ #

library(terra)
library(raster)
library(ellipsenm)
library(rgl)
library(ntbox)

# ------------------------------------------------------------------------------------------------ #

# Inputs #

# ------------------------------------------------------------------------------------------------ #

var_dir <- "E:/Porth_Ecol/env/PCA"
spp_dir <- "E:/Porth_Ecol/Occ/Joint"
M_dir <- "E:E:/Porth_Ecol/Ms"
dir_save <- "E:/Porth_Ecol/Grinn"
spp <- list.files(spp_dir, pattern = ".csv")
spp_com <- combn(length(spp), 2)
var <- list.files(var_dir, pattern = ".tif", full.names = TRUE)
env <- terra::rast(c(var))
crs(env) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
sp_data_list <- list()
niche_spp <- list()
names <- vector()
cov <- list()
cen <- list()
nv <- vector()
sim_data <- data.frame()
overlap_df <- data.frame()

# ------------------------------------------------------------------------------------------------ #

# Functions #

# ------------------------------------------------------------------------------------------------ #

calc_species_elipsoid <- function(sp_data, Ms) {
	
	niche <- overlap_object(sp_data, species = "Species", longitude = "Long", latitude = "Lat", method = "covmat", level = 95, variables = Ms)
	
}

# ------------------------------------------------------------------------------------------------ #

# Prep data #

# ------------------------------------------------------------------------------------------------ #

for (i in 1:length(spp)) {
		
	sp <- read.csv(paste0(spp_dir, "/", spp[[i]]))
	Ms <- list.files(M_dir, pattern = ".shp")
	Ms_1 <- terra::vect(paste0(M_dir, "/", Ms[[i]]))
	crs(Ms_1) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
	Ms_2f <- terra::crop(env, Ms_1)
	Ms_2f <- terra::mask(Ms_2f, Ms_1)
	M <- raster::stack(Ms_2f)
	sp_data_list[[i]] <- list(sp_data = sp[, 1:3], Ms = M)
	setTxtProgressBar(txtProgressBar(min = 0, max = length(spp), style = 3, width = 40, char = "="), i)
}
	
for (i in 1:length(sp_data_list)) {
	
	niche_spp[[i]] <- calc_species_elipsoid(sp_data_list[[i]]$sp_data, sp_data_list[[i]]$Ms)
	setTxtProgressBar(txtProgressBar(min = 0, max = length(sp_data_list), style = 3, width = 40, char = "="), i)

}
	
for (i in 1:length(spp)) {
		
	names[[i]] <- substr(spp[[i]], 1, nchar(spp[[i]])-4)
	setTxtProgressBar(txtProgressBar(min = 0, max = length(spp), style = 3, width = 40, char = "="), i)
	
}

names_com <- combn(names, 2)
	
for (i in 1:length(niche_spp)) {
	
	data <- raster::extract(sp_data_list[[i]]$Ms, sp_data_list[[i]]$sp_data[, 2:3])
	data <- cbind(sp_data_list[[i]]$sp_data, data)
	data <- na.omit(data)
	data <- as.data.frame(data)
	covar_centroid <- cov_center(data, mve = TRUE, level = 0.95, vars = c(4, 5, 6, 7, 8))
	cov_n1 <- covar_centroid$covariance
	cov <- rbind(cov, cov_n1)
	cen_n1 <-  covar_centroid$centroid
	cen[[i]]  <- as.vector(c(as.numeric(cen_n1[[1]]), as.numeric(cen_n1[[2]]), as.numeric(cen_n1[[3]]), as.numeric(cen_n1[[4]]),
	as.numeric(cen_n1[[5]])))
	nv[[i]] <- covar_centroid$niche_volume
	setTxtProgressBar(txtProgressBar(min = 0, max = length(niche_spp), style = 3, width = 40, char = "="), i)	

}
	
covm <- lapply(seq(1, nrow(cov), by = 5), function(i) {cov[i:min(i + 4, nrow(cov)), ]})

# ------------------------------------------------------------------------------------------------ #

# Final analysis #

# ------------------------------------------------------------------------------------------------ #

for (i in 1:ncol(spp_com)) {
		
	overlap_st <- ellipsoid_overlap(niche_spp[[spp_com[1, i]]], niche_spp[[spp_com[2, i]]], overlap_type = "all", significance_test = TRUE, replicates = 1000)
	full <- overlap_st@significance_results$full_random$Niche_1_vs_2[, 3]
	union <- overlap_st@significance_results$union_random$Niche_1_vs_2[, 3]
	sim <- rep(paste0("test", i), 1000)
	sim_data <- rbind(sim_data, data.frame(full = full, union = union, simulation = sim))
	test <- paste0(names_com[1, i],"→", names_com[2, i])
	full_o <- overlap_st@full_overlap[, 3]
	union_o <- overlap_st@union_overlap[, 3]
	full_p <- overlap_st@full_overlap[, 4]
	union_p <- overlap_st@union_overlap[, 4]
	overlap_df <- rbind(overlap_df, data.frame(test = test, spp1 = names_com[1, i], spp2 = names_com[2, i], full = full_o, union = union_o, p_full = full_p, p_union = union_p))
	write.csv(overlap_df, paste0(dir_save, "/overlap_df.csv"), row.names = FALSE)
	write.csv(sim_data, paste0(dir_save, "/sim_data.csv"), row.names = FALSE)
	setTxtProgressBar(txtProgressBar(min = 0, max = ncol(spp_com), style = 3, width = 40, char = "="), i)
	
}

# ------------------------------------------------------------------------------------------------ #

### EndNotRun

