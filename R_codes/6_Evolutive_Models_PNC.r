
# ------------------------------------------------------------------------------------------------ #

### Title: Univariate and Multivariate Evolutive Models ####
### Author: Patrón-Rivero, C. & J D Sánchez-R. ####
### Date: 16/14/2024 ###
### Project: "Ecological and Biogeographical Drivers of Speciation in Neotropical Hognose Pit Vipers, Porthidium (Squamata, Viperidae)" ###

# ------------------------------------------------------------------------------------------------ #

# Packages #

# ------------------------------------------------------------------------------------------------ #

library(maps)
library(ape)
library(phytools)
library(geiger)
library(readr)
library(mvMORPH)

# ------------------------------------------------------------------------------------------------ #

# Univariate Inputs #

# ------------------------------------------------------------------------------------------------ #

Tree <- read.tree("E:/1_Porth_Ecol/1_bsm/P_tree.newick")
Tree$tip.label
is.ultrametric(Tree)
TreeUlt <- force.ultrametric(Tree, method = c("nnls","extend"))
is.ultrametric(TreeUlt)
Cha <- read_delim("E:/1_Porth_Ecol/1_Paper/SM/SM_traits.csv")
rnames <- Cha[[1]] 
Dat <- data.matrix(Cha[2:17], rownames.force = TRUE)
rownames(Dat) <- unlist(Cha[1])             

# ------------------------------------------------------------------------------------------------ #

# Univariate function #

# ------------------------------------------------------------------------------------------------ #

fitBestModel <- function(trait = c("NP1", "NP2", "NP3", "NB", "bio1", "bio2", "bio3", "bio4", "bio7", "bio10", "bio11", "bio13",
								"bio14", "bio15", "bio16", "bio17")) {
  
	trait <- match.arg(trait, c("NP1", "NP2", "NP3", "NB", "bio1", "bio2", "bio3", "bio4", "bio7", "bio10", "bio11", "bio13",
								"bio14", "bio15", "bio16", "bio17"))
	models <- c("BM", "OU", "EB")
	summaries <- c("diffusion", "Ornstein-Uhlenbeck", "early burst")
    aic.se <- numeric(length(models))
	lnl.se <- numeric(length(models))
  
		for(m in 1:length(models)) {
    
			cat("\n\n\n\n\t*** ", paste(toupper(summaries[m]),": fitting ", sep = ""), models[m], " with SE *** \n", sep = "")
			tmp <- fitContinuous(TreeUlt, Dat[, trait], SE = NA, model = models[m], niter = 1000, ncores = 2)
			print(tmp)
			aic.se[m] <- tmp$opt$aicc
			lnl.se[m]<- tmp$opt$lnL
		
		}
  
	aic <- numeric(length(models))
	lnl <- numeric(length(models))
  
		for(m in 1:length(models)){
			
			cat("\n\n\n\n\t*** ", paste(toupper(summaries[m]),": fitting ", sep = ""), models[m], " *** \n", sep = "")
			tmp <- fitContinuous(TreeUlt,Dat[, trait], SE = 0, model = models[m], niter = 1000, ncores = 2)
			print(tmp)
			aic[m] <- tmp$opt$aicc
			lnl[m] <- tmp$opt$lnL
		
		}

	names(aic.se) <-names (lnl.se) <- names(aic) <- names(lnl) <- models
	delta_aic <- function(x) x-x[which(x==min(x))]
    daic<- delta_aic(aic)
	cat("\n\n\n\t\t\t\t*** MODEL COMPARISON: ", trait," *** \n", sep = "")
	cat("\tdelta-AIC values for models assuming no measurement error\t\t\t\t zero indicates the best model\n\n")
	print(daic, digits = 2)
	daic.se <- delta_aic(aic.se)
	cat("\n\n\n\n\t\t\t\t*** MODEL COMPARISON: ", trait," ***\n", sep = "")
	cat("\t\t   delta-AIC values for models estimating SE\t\t\t\t zero indicates the best model\n\n")
	print(daic.se, digits = 2)
	cat("\n\n\n")
	res_aicc <- rbind(aic, aic.se, daic, daic.se)
	rownames(res_aicc) <- c("AICc","AICc_SE","dAICc", "dAICc_SE")
	return(res_aicc)

}

# ------------------------------------------------------------------------------------------------ #

# Univariate results #

# ------------------------------------------------------------------------------------------------ #

res <- list()

for(i in 1:ncol(Dat)) {
	
	res[[i]] <- fitBestModel(colnames(Dat)[i])

}

# ------------------------------------------------------------------------------------------------ #

# Multivariate Inputs #

# ------------------------------------------------------------------------------------------------ #

NP <- data.matrix(Dat[1:3], rownames.force = TRUE)
NB <- data.matrix(Dat[, 4], rownames.force = TRUE)
TEM <- data.matrix(Dat[5:11], rownames.force = TRUE)
PREC <- data.matrix(Dat[12:16], rownames.force = TRUE)
Dat <- as.data.frame(Dat)

se <- function(data) {
  
  	for (col in colnames(data)) {
    
		se <- sd(data[[col]]) / sqrt(length(data[[col]]))
		data[[col]] <- se
	
	}
  
  return(data)

}

Dat_se <- se(Dat)

NP_se <- data.matrix(Dat_se[1:3], rownames.force = TRUE)
NB_se <- data.matrix(Dat_se[, 4], rownames.force = TRUE)
TEM_se <- data.matrix(Dat_se[5:11], rownames.force = TRUE)
PREC_se <- data.matrix(Dat_se[12:16], rownames.force = TRUE)

t <- list(NP, NB, TEM, PREC)
t_se <- list(NP_se, NB_se, TEM_se, PREC_se)

bm <- vector()
ou <- vector()
eb <- vector()

# ------------------------------------------------------------------------------------------------ #

# Multivariate Results #

# ------------------------------------------------------------------------------------------------ #

for(i in 1:length(t)) {

	bm_fit <- mvBM(TreeUlt, t[[i]], error = t_se[[i]], model = "BM1", echo = TRUE, control = list(maxit = 1000))
	ou_fit <- mvOU(TreeUlt, t[[i]], error = t_se[[i]], model = "OU1", echo = TRUE, control = list(maxit = 1000))
	eb_fit <- mvEB(TreeUlt, t[[i]], error = t_se[[i]], echo = TRUE, control = list(maxit = 1000))
	bm[i] <- bm_fit[3]
	ou[i] <- ou_fit[3]
	eb[i] <- eb_fit[3]

}

names <- c("NP", "NB", "TEM", "PREC")
AICc <- as.data.frame(cbind(names, as.numeric(bm), as.numeric(ou), as.numeric(eb)))

DAICc <- list()

for(i in 1:nrow(AICc)) {

	ind_min <- min(as.numeric(AICc[i, 2:4]))
	DAICc[[i]] <- as.numeric(AICc[i, 2:4]) - ind_min

}

D <- do.call(rbind, DAICc)

setwd("E:/1_Porth_Ecol/1_Paper/SM")
write.csv(AICc, "./AICc.csv", row.names = FALSE)
write.csv(D, "./DAICc.csv", row.names = FALSE)

# ------------------------------------------------------------------------------------------------ #

### EndNotRun
