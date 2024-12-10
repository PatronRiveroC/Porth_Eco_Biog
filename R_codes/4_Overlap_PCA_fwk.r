
# ------------------------------------------------------------------------------------------------ #

### Title: Overlap test by PCA-environment framework ####
### Author: Patrón-Rivero, C. ####
### Date: 11/14/2024 ###
### Project: "Ecological and Biogeographical Drivers of Speciation in Neotropical Hognose Pit Vipers, Porthidium (Squamata, Viperidae)" ###

# ------------------------------------------------------------------------------------------------ #

# Packages #

# ------------------------------------------------------------------------------------------------ #

library(sp)
library(raster)
library(dismo)
library(ENMTools)
library(ecospat)
library(ellipsenm)
library(gtools)
library(rgdal)
library(ggplot2)
library(patchwork)
library(sf)
library(terra)

# ------------------------------------------------------------------------------------------------ #

# Function plots #

# ------------------------------------------------------------------------------------------------ #

histo_D <- function(dat, variable, line, title) {
			
			ggplot(dat, aes(x = variable, fill = variable)) +
			geom_histogram(fill = "red", alpha = 0.2, bins = 30) +
			geom_vline(xintercept = line, color = "black", linetype = "dashed") +
			ggtitle(title) +
			labs(x = NULL, y = NULL) +
			theme_bw() +
			theme(plot.title = element_text(hjust = 0, margin = margin(b = 10)))

}

histo_I <- function(dat, variable, line, title) {
			
			ggplot(dat, aes(x = variable, fill = variable)) +
			geom_histogram(fill = "blue", alpha = 0.2, bins = 30) +
			geom_vline(xintercept = line, color = "black", linetype = "dashed") +
			ggtitle(title) +
			labs(x = NULL, y = NULL) +
			theme_bw() +
			theme(plot.title = element_text(hjust = 0, margin = margin(b = 10)))

}

# ------------------------------------------------------------------------------------------------ #

# Function for similarity test #

# ------------------------------------------------------------------------------------------------ #

sim_spp_test <- function(spp_dir, var_dir, M_dir, n_back, dir_save) {
	
	spp <- list.files(spp_dir, pattern = ".csv")
	sp_data_list <- list()
	overlap <- list()
	
	var <- list.files(var_dir, pattern = ".tif", full.names = TRUE)
	env <- terra::rast(c(var))
	crs(env) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
		
	for (i in 1:length(spp)) {
		
		sp <- read.csv(paste0(spp_dir, "/", spp[[i]]))
		Ms <- list.files(M_dir, pattern = ".shp")
		Ms_1 <- terra::vect(paste0(M_dir, "/", Ms[[i]]))
		crs(Ms_1) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
		Ms_2f <- terra::crop(env, Ms_1)
		Ms_2f <- terra::mask(Ms_2f, Ms_1)
		M <- raster::stack(Ms_2f)
		Ms_3 <- rasterToPoints(M[[1]])
			
			if(nrow(Ms_3) < n_back) {
				
				Ms_comp <- Ms_3
				
				} else {
					
					Ms_comp <- Ms_3[sample(nrow(Ms_3), n_back), ]
				
				}
		
		sp1 <- raster::extract(M, sp[, 2:3])
		sp2 <- cbind(sp[, 1:3], sp1)
		sp3 <- na.omit(sp2)
		sp4 <- st_as_sf(sp3, coords = c("Long", "Lat"), crs = 4326)
		sp5 <- vect(sp4)
		Ms_comp <- as.data.frame(Ms_comp)
		Ms_comp <- st_as_sf(Ms_comp, coords = c("x", "y"), crs = 4326)
		Mss <- vect(Ms_comp)
		sp_data_list[[i]] <- list(sp_data = sp5, Ms = M, bg = Mss)
	
	}
	
	names <- vector()
	
	for (i in 1:length(spp)) {
		
		names[[i]] <- substr(spp[[i]], 1, nchar(spp[[i]])-4)
	
	}
	
	spp_com <- cbind(combn(length(spp), 2), combn(length(spp), 2, rev))
	names_com <- cbind(combn(names, 2), combn(names, 2, rev))
	
	for(i in 1:ncol(names_com)) {
		
		focal_sp <- enmtools.species()
		comp_sp <- enmtools.species()
		focal_sp$species.name <- names_com[1, i]
		comp_sp$species.name <- names_com[2, i]
		focal_sp$presence.points <- sp_data_list[[spp_com[1, i]]]$sp_data[, 2:3]
		comp_sp$presence.points <- sp_data_list[[spp_com[2, i]]]$sp_data[, 2:3]
		focal_sp$background.points <- sp_data_list[[spp_com[1, i]]]$bg[, 1:2]
		comp_sp$background.points <- sp_data_list[[spp_com[2, i]]]$bg[, 1:2]
		test <- paste0(names_com[1, i],"→", names_com[2, i])
		focal_m <- enmtools.ecospat.bg(focal_sp, comp_sp, env[[1:2]], nreps = 1000, bg.source = "points")
			
		D <- focal_m[[15]]$obs$D
		pD <- focal_m[[15]]$p.D
		simD <- as.data.frame(focal_m[[15]]$sim$D)
		colnames(simD)[1] <- "simulate"
		p_d <- as.numeric(focal_m[[15]]$obs$D)
		p_D <- round(focal_m[[15]]$p.D, 4)
		tD <- paste0(names_com[1, i], " vs ", names_com[2, i], " p = ", p_D)
		hD <- histo_D(dat = simD, variable = simD$simulate, line = p_d, title = tD)
		
		I <- focal_m[[15]]$obs$I
		pI <- focal_m[[15]]$p.I
		simI <- as.data.frame(focal_m[[15]]$sim$I)
		colnames(simI)[1] <- "simulate"
		p_i <- as.numeric(focal_m[[15]]$obs$I)
		p_I <- round(focal_m[[15]]$p.I, 4)
		tI <- paste0(names_com[1, i], " vs ", names_com[2, i], " p = ", p_I)
		hI <- histo_I(dat = simI, variable = simI$simulate, line = p_i, title = tI)
		
		qw <- hD / hI
		p1 <- qw + plot_annotation(theme = theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold")), 
				tag_levels = 'A', tag_suffix = ')')
		
		setwd(dir_save)
		ggsave(file = paste0(names_com[1, i], "_", names_com[2, i], ".jpeg"), plot = p1, width = 20, height = 20, dpi = 300, 
			   units = "cm", device = "jpeg")
		
		overlap[[i]] <- list(test = test, D = D, pD = pD, simD = simD, I = I, pI = pI, simI = simI)
		setTxtProgressBar(txtProgressBar(min = 0, max = ncol(names_com), style = 3, width = 40, char = "="), i)

	}					

	overlap_df <- data.frame(test = character(0), D = numeric(0), pD = numeric(0), I = numeric(0), pI = numeric(0))

	for (i in 1:length(overlap)) {
		
		test_column <- overlap[[i]]$test
		D_column <- overlap[[i]]$D
		pD_column <- overlap[[i]]$pD
		I_column <- overlap[[i]]$I
		pI_column <- overlap[[i]]$pI
		overlap_df <- rbind(overlap_df, data.frame(test = test_column, D = D_column, pD = pD_column, I = I_column, pI = pI_column))
	
	}

	sim_data <- data.frame(simD = numeric(0), simI = numeric(0), simulation = character(0))

	for (i in 1:length(overlap)) {
		
		simD_column <- overlap[[i]]$simD
		simI_column <- overlap[[i]]$simI
		simulation_factor <- paste("test", i, sep = "")
		sim_data <- rbind(sim_data, data.frame(simD = simD_column, simI = simI_column, simulation = simulation_factor))
	
	}
	
	colnames(sim_data) <- c("simD", "simI", "simulation")
	write.csv(overlap_df, "overlap_df.csv", row.names = FALSE)
	write.csv(sim_data, "sim_data.csv", row.names = FALSE)
	
}

# ------------------------------------------------------------------------------------------------ #

# Final analysis #

# ------------------------------------------------------------------------------------------------ #

set.seed(20)
spp_dir <- "E:/Porth_Ecol/occ/joint"
var_dir <- "E:/Porth_Ecol/env/PCA"
M_dir <- "E:/Porth_Ecol/env/M"
n_back <- 10000
dir_save <- "E:/Porth_Ecol/PCA_fwk"

sim_spp_test(spp_dir = spp_dir, var_dir = var_dir, M_dir = M_dir, n_back = n_back, dir_save = dir_save)

# ------------------------------------------------------------------------------------------------ #

### EndNotRun