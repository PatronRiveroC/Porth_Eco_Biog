
# ------------------------------------------------------------------------------------------------ #

### Title: Biogeographic stochastic models ####
### Author: Patrón-Rivero, C. ####
### Date: 11/14/2024 ###
### Project: "Ecological and Biogeographical Drivers of Speciation in Neotropical Hognose Pit Vipers, Porthidium (Squamata, Viperidae)" ###

# ------------------------------------------------------------------------------------------------ #

# Packages #

# ------------------------------------------------------------------------------------------------ #

library(terra)
library(ape)
library(optimx)
library(GenSA)
library(rexpokit)
library(cladoRcpp)
library(snow)
library(parallel)
library(BioGeoBEARS)

# ------------------------------------------------------------------------------------------------ #

# Preprocess #

# ------------------------------------------------------------------------------------------------ #

bio_morrone <- rast("E:/Porth_Ecol/bsm/ecor.asc")
setwd("E:/1_Porth_Ecol/1_csv/1_spp")
csv_files <- list.files(pattern = ".csv")
spp <- list()

for (file in csv_files) {
		
	data <- read.csv(file, header = TRUE)
	data <- extract(bio_morrone, data[2:3])
	data <- na.omit(data)
	data <- unique(data$geograf)
	spp[[file]] <- data
		
}

all_eco <- unique(unlist(spp))

# Mexican transition zone / Mesoamerican dominion (Tehuantepec left) = M
# 46 = Mexican transition zone = T = Sierra Madre "O"ccidental province = O
# 5 = Mesoamerican dominion = M = "B"alsas Basin province = B

# Mesoamerican dominion / Mexican transition zone (Tehuantepec rigth - Motagua–Polochic left) = T
# 12 = Mexican transition zone = T = Chiapas Highlands province = H
# 35 = Mesoamerican dominion = Pacific "L"owlands province = L
# 54 = Mesoamerican dominion = Veracruzan province = V

# Mesoamerican dominion (Yucatan Peninsula) = Y
# 57 = "Y"ucatan Peninsula province = Y

# Mesoamerican dominion (Motagua–Polochic rigth - Nicaraguan Depression left) = N
# 32 = Mos"Q"uito province = Q

# Pacific dominion (Nicaraguan Depression rigth - Western Andes Range left) = P
# 7 = C"A"uca province = A
# 13 = Chocó-Darién province = D
# 21 = Gua"J"ira province = J
# 22 = Guatuso-"T"alamanca province = T
# 30 = Magdalena province = M
# 42 = Puntarenas-"C"hiriquí province = C
# 45 = "S"abana province = S 
# 53 = Vene"Z"uelan province = Z
# 55 = "W"estern Ecuador province = W

# South American transition zone = S
# 38 = "P"áramo province = P

# ------------------------------------------------------------------------------------------------ #

# Inputs #

# ------------------------------------------------------------------------------------------------ #

wd <- np("E:/Porth_Ecol/bsm")
setwd(wd)

trfn <- np(paste(addslash("E:/Porth_Ecol/bsm"), "P_tree.newick", sep = ""))
tr = read.tree(trfn)
geogfn <- np(paste(addslash("E:/Porth_Ecol/bsm"), "Porth_geog.data", sep = ""))
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
max_range_size <- 6

areas = getareas_from_tipranges_object(tipranges)
states_list_0based = rcpp_areas_list_to_states_list(areas = areas, maxareas = max_range_size, include_null_range = TRUE)

ranges_list = NULL

for (i in 1:length(states_list_0based)) {    
		
	if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) ) {
        
		tmprange <- "_"
        
	} else {
        
		tmprange <- paste(areas[states_list_0based[[i]]+1], collapse="")
	}
			
	ranges_list <- c(ranges_list, tmprange)
    
}

# ------------------------------------------------------------------------------------------------ #

# Biogeographical assumptions #

# ------------------------------------------------------------------------------------------------ #

# 3.1 - 3.5: Isthmus of Tehuantepec
# 7.7: the Motagua–Polochic Faults
# 10: back-arc formation of the Nicaraguan Depression
# 15: oldest than the tree age

# STRATUM 1 (youngest)
disallowed1 <- c()
keepTF1 <- ranges_list %in% disallowed1 == FALSE
ranges_list_NEW1 <- ranges_list[keepTF1]

# STRATUM 2
ranges_to_lose <- ranges_list_NEW1[grep(pattern = "M", x = ranges_list)]
keepTF2 <- ranges_list_NEW1 %in% ranges_to_lose == FALSE
ranges_list_NEW2 <- ranges_list_NEW1[keepTF2]

# STRATUM 3
ranges_to_lose <- ranges_list_NEW2[grep(pattern = "T", x = ranges_list_NEW2)]
keepTF3 <- ranges_list_NEW2 %in% ranges_to_lose == FALSE
ranges_list_NEW3 <- ranges_list_NEW2[keepTF3]

# STRATUM 4 (oldest)
ranges_to_lose <- ranges_list_NEW3[grep(pattern = "N", x = ranges_list_NEW3)]
keepTF4 <- ranges_list_NEW3 %in% ranges_to_lose == FALSE
ranges_list_NEW4 <- ranges_list_NEW3[keepTF4]

states_list_0based_NEW1 <- states_list_0based[keepTF1]
states_list_0based_NEW2 <- states_list_0based_NEW1[keepTF2]
states_list_0based_NEW3 <- states_list_0based_NEW2[keepTF3]
states_list_0based_NEW4 <- states_list_0based_NEW3[keepTF4]
states_list_0based_NEW5 <- states_list_0based_NEW4

lists_of_states_lists_0based <- list()
lists_of_states_lists_0based[[1]] <- states_list_0based_NEW1
lists_of_states_lists_0based[[2]] <- states_list_0based_NEW2
lists_of_states_lists_0based[[3]] <- states_list_0based_NEW3
lists_of_states_lists_0based[[4]] <- states_list_0based_NEW4
lists_of_states_lists_0based[[5]] <- states_list_0based_NEW5

# ------------------------------------------------------------------------------------------------ #

# Run DEC model #

# ------------------------------------------------------------------------------------------------ #

BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn <- trfn
BioGeoBEARS_run_object$geogfn <- geogfn
BioGeoBEARS_run_object$max_range_size <- max_range_size
BioGeoBEARS_run_object$min_branchlength <- 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range <- TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
BioGeoBEARS_run_object$timesfn <- paste0(getwd(), "/timeperiods.txt")
BioGeoBEARS_run_object$dispersal_multipliers_fn <- paste0(getwd(), "/dispersal_multipliers.txt")
BioGeoBEARS_run_object$on_NaN_error <- -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup <- TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx <- TRUE    # if FALSE, use optim() instead of optimx();
BioGeoBEARS_run_object$num_cores_to_use <- 1
BioGeoBEARS_run_object$force_sparse <- FALSE    # force_sparse <- TRUE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object <- section_the_tree(inputs = BioGeoBEARS_run_object, make_master_table = TRUE, plot_pieces = FALSE, 
										   fossils_older_than = 0.001, cut_fossils = FALSE)
BioGeoBEARS_run_object$return_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_ancprobs <- TRUE    # get ancestral states from optim run
BioGeoBEARS_run_object$lists_of_states_lists_0based <- lists_of_states_lists_0based
BioGeoBEARS_run_object <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow <- TRUE
resfn <- "P_strat_DEC.Rdata"

if (runslow) {
    
	res <- bears_optim_run(BioGeoBEARS_run_object)
    save(res, file = resfn)
    resDEC <- res
    
	} else {
	
		load(resfn)
		resDEC <- res
	
    }

# ------------------------------------------------------------------------------------------------ #

# Run DEC + J model #

# ------------------------------------------------------------------------------------------------ #

BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn <- trfn
BioGeoBEARS_run_object$geogfn <- geogfn
BioGeoBEARS_run_object$max_range_size <- max_range_size
BioGeoBEARS_run_object$min_branchlength <- 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range <- TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
BioGeoBEARS_run_object$timesfn <- paste0(getwd(), "/timeperiods.txt")
BioGeoBEARS_run_object$dispersal_multipliers_fn <- paste0(getwd(), "/dispersal_multipliers.txt")
BioGeoBEARS_run_object$on_NaN_error <- -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup <- TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx <- TRUE    # if FALSE, use optim() instead of optimx();
BioGeoBEARS_run_object$num_cores_to_use <- 1
BioGeoBEARS_run_object$force_sparse <- FALSE    # force_sparse <- TRUE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object <- section_the_tree(inputs = BioGeoBEARS_run_object, make_master_table = TRUE, plot_pieces = FALSE,
										   fossils_older_than = 0.001, cut_fossils = FALSE)
BioGeoBEARS_run_object$return_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_ancprobs <- TRUE    # get ancestral states from optim run
dstart <- resDEC$outputs@params_table["d","est"]
estart <- resDEC$outputs@params_table["e","est"]
jstart = 0.0001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] <- dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] <- dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] <- estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] <- estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] <- "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] <- jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] <- jstart
BioGeoBEARS_run_object$lists_of_states_lists_0based <- lists_of_states_lists_0based
BioGeoBEARS_run_object <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn <- "P_strat_DEC_J.Rdata"
runslow <- TRUE

if (runslow) {
    
	res <- bears_optim_run(BioGeoBEARS_run_object)
    save(res, file = resfn)
    resDECj <- res
    
	} else {
    
    load(resfn)
    resDECj <- res
    
	}

# ------------------------------------------------------------------------------------------------ #

# Run DIVALIKE model #

# ------------------------------------------------------------------------------------------------ #

BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn <- trfn
BioGeoBEARS_run_object$geogfn <- geogfn
BioGeoBEARS_run_object$max_range_size <- max_range_size
BioGeoBEARS_run_object$min_branchlength <- 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range <- TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
BioGeoBEARS_run_object$timesfn <- paste0(getwd(), "/timeperiods.txt")
BioGeoBEARS_run_object$dispersal_multipliers_fn <- paste0(getwd(), "/dispersal_multipliers.txt")
BioGeoBEARS_run_object$on_NaN_error <- -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup <- TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx <- TRUE    # if FALSE, use optim() instead of optimx();
BioGeoBEARS_run_object$num_cores_to_use <- 1
BioGeoBEARS_run_object$force_sparse <- FALSE    # force_sparse <- TRUE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object <- section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE, fossils_older_than=0.001, cut_fossils=FALSE)
BioGeoBEARS_run_object$return_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_ancprobs <- TRUE    # get ancestral states from optim run
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] <- 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] <- 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] <- "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] <- "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] <- "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] <- "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] <- 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] <- 0.5
BioGeoBEARS_run_object$lists_of_states_lists_0based <- lists_of_states_lists_0based
BioGeoBEARS_run_object <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn <- "P_strat_DIVA.Rdata"

if (runslow) {
    res <- bears_optim_run(BioGeoBEARS_run_object)
    save(res, file = resfn)
    resDIVALIKE <- res
    
	} else {
    
    load(resfn)
    resDIVALIKE <- res
	
    }

# ------------------------------------------------------------------------------------------------ #

# Run DIVALIKE + J model #

# ------------------------------------------------------------------------------------------------ #

BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn <- trfn
BioGeoBEARS_run_object$geogfn <- geogfn
BioGeoBEARS_run_object$max_range_size <- max_range_size
BioGeoBEARS_run_object$min_branchlength <- 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range <- TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
BioGeoBEARS_run_object$timesfn <- paste0(getwd(), "/timeperiods.txt")
BioGeoBEARS_run_object$dispersal_multipliers_fn <- paste0(getwd(), "/dispersal_multipliers.txt")
BioGeoBEARS_run_object$on_NaN_error <- -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup <- TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx <- TRUE    # if FALSE, use optim() instead of optimx();
BioGeoBEARS_run_object$num_cores_to_use <- 1
BioGeoBEARS_run_object$force_sparse <- FALSE    # force_sparse <- TRUE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object <- section_the_tree(inputs = BioGeoBEARS_run_object, make_master_table = TRUE, plot_pieces = FALSE, 
										   fossils_older_than = 0.001, cut_fossils = FALSE)
BioGeoBEARS_run_object$return_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_ancprobs <- TRUE    # get ancestral states from optim run
dstart <- resDIVALIKE$outputs@params_table["d","est"]
estart <- resDIVALIKE$outputs@params_table["e","est"]
jstart <- 0.0001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] <- dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] <- dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] <- estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] <- estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] <- 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] <- 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] <- "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] <- "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] <- "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] <- "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] <- 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] <- 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] <- "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] <- jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] <- jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] <- 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] <- 1.99999
BioGeoBEARS_run_object$lists_of_states_lists_0based <- lists_of_states_lists_0based
BioGeoBEARS_run_object <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn <- "P_strat_DIVA_J.Rdata"
runslow <- TRUE

if (runslow) {
    
	res <- bears_optim_run(BioGeoBEARS_run_object)
    save(res, file = resfn)
    resDIVALIKEj <- res
    
	} else {
    
    load(resfn)
    resDIVALIKEj <- res
	
    }

# ------------------------------------------------------------------------------------------------ #

# Run BAYAREALIKE model #

# ------------------------------------------------------------------------------------------------ #

BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn <- trfn
BioGeoBEARS_run_object$geogfn <- geogfn
BioGeoBEARS_run_object$max_range_size <- max_range_size
BioGeoBEARS_run_object$min_branchlength <- 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range <- TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
BioGeoBEARS_run_object$timesfn <- paste0(getwd(), "/timeperiods.txt")
BioGeoBEARS_run_object$dispersal_multipliers_fn <- paste0(getwd(), "/dispersal_multipliers.txt")
BioGeoBEARS_run_object$on_NaN_error <- -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup <- TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx <- TRUE    # if FALSE, use optim() instead of optimx();
BioGeoBEARS_run_object$num_cores_to_use <- 1
BioGeoBEARS_run_object$force_sparse <- FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object <- section_the_tree(inputs = BioGeoBEARS_run_object, make_master_table = TRUE, plot_pieces = FALSE, 
										   fossils_older_than = 0.001, cut_fossils = FALSE)
BioGeoBEARS_run_object$return_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_ancprobs <- TRUE    # get ancestral states from optim run
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] <- 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] <- 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] <- 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] <- 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] <- "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] <- "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] <- "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] <- 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] <- 0.9999
BioGeoBEARS_run_object$lists_of_states_lists_0based <- lists_of_states_lists_0based
BioGeoBEARS_run_object <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow <- TRUE
resfn <- "P_strat_BAY.Rdata"

if (runslow) {
    
	res <- bears_optim_run(BioGeoBEARS_run_object)
    save(res, file = resfn)
    resBAYAREALIKE <- res
    
	} else {
    
	load(resfn)
    resBAYAREALIKE <- res
	
    }

# ------------------------------------------------------------------------------------------------ #

# Run BAYAREALIKE + J model #

# ------------------------------------------------------------------------------------------------ #

BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn <- trfn
BioGeoBEARS_run_object$geogfn <- geogfn
BioGeoBEARS_run_object$max_range_size <- max_range_size
BioGeoBEARS_run_object$min_branchlength <- 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range <- TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
BioGeoBEARS_run_object$timesfn <- paste0(getwd(), "/timeperiods.txt")
BioGeoBEARS_run_object$dispersal_multipliers_fn <- paste0(getwd(), "/dispersal_multipliers.txt")
BioGeoBEARS_run_object$on_NaN_error <- -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup <- TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx <- "GenSA"
BioGeoBEARS_run_object$num_cores_to_use <- 1
BioGeoBEARS_run_object$force_sparse <- FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object <- section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE, fossils_older_than=0.001, cut_fossils=FALSE)
BioGeoBEARS_run_object$return_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_ancprobs <- TRUE    # get ancestral states from optim run
dstart <- resBAYAREALIKE$outputs@params_table["d","est"]
estart <- resBAYAREALIKE$outputs@params_table["e","est"]
jstart <- 0.0001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] <- dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] <- dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] <- estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] <- estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] <- 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] <- 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] <- 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] <- 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] <- "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] <- jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] <- jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] <- 0.99999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] <- "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] <- "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] <- "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] <- 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] <- 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] <- 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] <- 4.9999999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] <- 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] <- 4.9999999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] <- 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] <- 0.99999
BioGeoBEARS_run_object$lists_of_states_lists_0based <  lists_of_states_lists_0based
BioGeoBEARS_run_object <- fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object = BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn <- "P_strat_BAY_J.Rdata"
runslow <- TRUE

if (runslow) {
    
	res <- bears_optim_run(BioGeoBEARS_run_object)
    save(res, file = resfn)
	resBAYAREALIKEj <- res
    
	} else {
    
    load(resfn)
    resBAYAREALIKEj <- res
	
    }

# ------------------------------------------------------------------------------------------------ #

# Summary of statistics to compare #

# ------------------------------------------------------------------------------------------------ #

restable <- NULL
teststable <- NULL
LnL_2 <- get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_1 <- get_LnL_from_BioGeoBEARS_results_object(resDECj)

numparams1 <- 3
numparams2 <- 2
stats <- AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

res2 <- extract_params_from_BioGeoBEARS_results_object(results_object = resDEC, returnwhat = "table", addl_params = c("j"), paramsstr_digits = 4)
res1 <- extract_params_from_BioGeoBEARS_results_object(results_object = resDECj, returnwhat = "table", addl_params = c("j"), paramsstr_digits = 4)

rbind(res2, res1)
tmp_tests <- conditional_format_table(stats)

restable <- rbind(restable, res2, res1)
teststable <- rbind(teststable, tmp_tests)

LnL_2 <- get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
LnL_1 <- get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)

numparams1 <- 3
numparams2 <- 2
stats <- AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

res2 <- extract_params_from_BioGeoBEARS_results_object(results_object = resDIVALIKE, returnwhat = "table", addl_params = c("j"), paramsstr_digits = 4)
res1 <- extract_params_from_BioGeoBEARS_results_object(results_object = resDIVALIKEj, returnwhat = "table", addl_params = c("j"), paramsstr_digits = 4)
rbind(res2, res1)
tmp_tests <- conditional_format_table(stats)

restable <- rbind(restable, res2, res1)
teststable <- rbind(teststable, tmp_tests)

LnL_2 <- get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
LnL_1 <- get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)
numparams1 <- 3
numparams2 <- 2
stats <- AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
res2 <- extract_params_from_BioGeoBEARS_results_object(results_object = resBAYAREALIKE, returnwhat = "table", addl_params = c("j"), paramsstr_digits = 4)
res1 <- extract_params_from_BioGeoBEARS_results_object(results_object = resBAYAREALIKEj, returnwhat = "table", addl_params = c("j"), paramsstr_digits = 4)
rbind(res2, res1)
tmp_tests <- conditional_format_table(stats)

restable <- rbind(restable, res2, res1)
teststable <- rbind(teststable, tmp_tests)

# ------------------------------------------------------------------------------------------------ #

# Essamble results #

# ------------------------------------------------------------------------------------------------ #

teststable$alt <- c("DEC+J", "DIVALIKE+J", "BAYAREALIKE+J")
teststable$null <- c("DEC", "DIVALIKE", "BAYAREALIKE")
row.names(restable) <- c("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J")
restable <- put_jcol_after_ecol(restable)

# ------------------------------------------------------------------------------------------------ #

# Model weigths #

# ------------------------------------------------------------------------------------------------ #

restable2 <- restable
AICtable <- calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
restable <- cbind(restable, AICtable)
restable_AIC_rellike <- AkaikeWeights_on_summary_table(restable = restable, colname_to_use = "AIC")
restable_AIC_rellike <- put_jcol_after_ecol(restable_AIC_rellike)
samplesize <- length(tr$tip.label)
AICtable <- calc_AICc_column(LnL_vals = restable$LnL, nparam_vals = restable$numparams, samplesize = samplesize)
restable2 <- cbind(restable2, AICtable)
restable_AICc_rellike <- AkaikeWeights_on_summary_table(restable = restable2, colname_to_use = "AICc")
restable_AICc_rellike <- put_jcol_after_ecol(restable_AICc_rellike)

write.table(restable_AICc_rellike, file = "restable.csv", quote = FALSE, sep = ",")
write.table(unlist_df(teststable), file = "teststable.csv", quote = FALSE, sep = ",")

# ------------------------------------------------------------------------------------------------ #

### EndNotRun
