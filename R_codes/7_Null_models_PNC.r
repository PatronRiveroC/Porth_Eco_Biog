
# ------------------------------------------------------------------------------------------------ #

### Title: Null models Spp ####
### Author: Patrón-Rivero, C. ####
### Date: 11/14/2024 ###
### Project: "Ecological and Biogeographical Drivers of Speciation in Neotropical Hognose Pit Vipers, Porthidium (Squamata, Viperidae)" ###

# ------------------------------------------------------------------------------------------------ #

# Packages #

# ------------------------------------------------------------------------------------------------ #

library(dplyr)
library(ggplot2)
library(patchwork)
library(car)

# ------------------------------------------------------------------------------------------------ #

# Prepare data #

# ------------------------------------------------------------------------------------------------ #

DI <- read.csv("E:/1_Porth_Ecol/1_Paper/data/all.csv")
FU <- read.csv("E:/1_Porth_Ecol/1_Paper/data/overlap_df.csv")
combinations <- c("P_arc", "P_lan", "P_por", "P_nas", "P_yuc")
target_spp2 <- c("P_hes", "P_dun", "P_oph")

process_data <- function(data, combinations, target_spp2) {
  
  df1 <- subset(data, Spp1 %in% combinations & Spp2 %in% target_spp2)
  df2 <- subset(data, Spp2 %in% combinations & Spp1 %in% target_spp2)

  df_nor <- rbind(df1, df2)
  df_nor$pair <- "Not related"
  
  df_rel <- anti_join(data, df_nor)
  df_rel <- na.omit(df_rel)
  df_rel$pair <- "Related"
  
  result <- rbind(df_nor, df_rel)
  return(result)
  
}

DI <- process_data(DI, combinations, target_spp2)
FU <- process_data(FU, combinations, target_spp2)

DI_s <- split(DI, DI$pair)
DI_nor <- DI_s[[1]]
DI_rel <- DI_s[[2]]

FU_s <- split(FU, FU$pair)
FU_nor <- FU_s[[1]]
FU_rel <- FU_s[[2]]

mean_data <- function(data, columns) {
	
	mean <- data %>%
    group_by(Spp1) %>%
    summarise(across(all_of(columns), mean))
	return(mean)
	
}

col <- c("D", "I")
DI_mean_n <- mean_data(DI_nor, col)
DI_mean_r <- mean_data(DI_rel, col)

col <- c("full", "union")
FU_mean_n <- mean_data(FU_nor, col)
FU_mean_r <- mean_data(FU_rel, col)

# ------------------------------------------------------------------------------------------------ #

# Random distribution of data #

# ------------------------------------------------------------------------------------------------ #

nor_D <- replicate(1000, {
  
  mm <- sample(DI_nor$D, size = 8, replace = TRUE)
  mean(mm)
  
})

nor_I <- replicate(1000, {
  
  mm <- sample(DI_nor$I, size = 8, replace = TRUE)
  mean(mm)
  
})

rel_D <- replicate(1000, {
  
  mm <- sample(DI_rel$D, size = 8, replace = TRUE)
  mean(mm)
  
})

rel_I <- replicate(1000, {
  
  mm <- sample(DI_rel$I, size = 8, replace = TRUE)
  mean(mm)
  
})

nor_F <- replicate(1000, {
  
  mm <- sample(FU_nor$full, size = 8, replace = TRUE)
  mean(mm)
  
})

nor_U <- replicate(1000, {
  
  mm <- sample(FU_nor$union, size = 8, replace = TRUE)
  mean(mm)
  
})

rel_F <- replicate(1000, {
  
  mm <- sample(FU_rel$full, size = 8, replace = TRUE)
  mean(mm)
  
})

rel_U <- replicate(1000, {
  
  mm <- sample(FU_rel$union, size = 8, replace = TRUE)
  mean(mm)
  
})

D_sim <- data.frame(pair = rep(c("Unrelated", "Related"), each = 1000), D = c(nor_D, rel_D))

I_sim <- data.frame(pair = rep(c("Unrelated", "Related"), each = 1000), I = c(nor_I, rel_I))

F_sim <- data.frame(pair = rep(c("Unrelated", "Related"), each = 1000), F = c(nor_F, rel_F))

U_sim <- data.frame(pair = rep(c("Unrelated", "Related"), each = 1000), U = c(nor_U, rel_U))

# ------------------------------------------------------------------------------------------------ #

# Stat analysis #

# ------------------------------------------------------------------------------------------------ #

mean_rel_D <- mean(DI_rel$D)
mean_rel_I <- mean(DI_rel$I)
mean_rel_F <- mean(FU_rel$full)
mean_rel_U <- mean(FU_rel$union)

# Normality test #

shapiro.test(DI_rel$D) # p.value < 0.05
shapiro.test(DI_rel$I) # p.value < 0.05
shapiro.test(FU_rel$full) # p.value < 0.05
shapiro.test(FU_rel$union) # p.value < 0.05

shapiro.test(DI_nor$D) # p.value < 0.05
shapiro.test(DI_nor$I) # p.value < 0.05
shapiro.test(FU_nor$full) # p.value < 0.05
shapiro.test(FU_nor$union) # p.value < 0.05

# Homogeneity test #

fligner.test(x = list(DI_rel$D, DI_nor$D)) # p.value > 0.05
fligner.test(x = list(DI_rel$I, DI_nor$I)) # p.value > 0.05
fligner.test(x = list(FU_rel$full, FU_nor$full)) # p.value > 0.05
fligner.test(x = list(FU_rel$union, FU_nor$union)) # p.value > 0.05

# Wilcoxon-Mann-Whitney Test #

# Index #

wilcox.test(x = DI$D, y = DI$I, conf.int = 0.95)
wilcox.test(x = FU$full, y = FU$union, conf.int = 0.95)

# Related vs not Related #

wilcox.test(x = DI_rel$D, y = DI_nor$D, conf.int = 0.95) # p.value > 0.05
wilcox.test(x = DI_rel$I, y = DI_nor$I, conf.int = 0.95) # p.value > 0.05
wilcox.test(x = FU_rel$full, y = FU_nor$full, conf.int = 0.95) # p.value > 0.05
wilcox.test(x = FU_rel$union, y = FU_nor$union, conf.int = 0.95) # p.value > 0.05

# ------------------------------------------------------------------------------------------------ #

# Histogram plot #

# ------------------------------------------------------------------------------------------------ #

create_histogram <- function(data, x, bins, segment_x, show_x_text = TRUE, show_y_text = TRUE) {

	p <- ggplot(data, aes(x = x)) +
		geom_histogram(fill = "#74c476", lwd = 0.2, color = "black", bins = bins, alpha = 0.5) +
		geom_segment(aes(x = segment_x, y = Inf, xend = segment_x, yend = 0.5),
					 color = "#ef3b2c", size = 0.4, arrow = arrow(length = unit(0.25, "cm"))) +
		labs(title = NULL, x = "Unrelated pairs niche overlap", y = "Frequency") +
		theme_bw() +
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  
	
	if (!show_x_text) {
	
		p <- p + labs(x = NULL)
	
	}
  
	if (!show_y_text) {
    
		p <- p + labs(y = NULL)
  
	}
  
	return(p)
}

nor_D <- as.data.frame(nor_D)
nor_I <- as.data.frame(nor_I)
nor_F <- as.data.frame(nor_F)
nor_U <- as.data.frame(nor_U)

D <- create_histogram(nor_D, nor_D$nor_D, 30, mean_rel_D, show_x_text = FALSE, show_y_text = TRUE)
I <- create_histogram(nor_I, nor_I$nor_I, 30, mean_rel_I, show_x_text = FALSE, show_y_text = FALSE)
F <- create_histogram(nor_F, nor_F$nor_F, 30, mean_rel_F, show_x_text = TRUE, show_y_text = TRUE)
U <- create_histogram(nor_U, nor_U$nor_U, 30, mean_rel_U, show_x_text = TRUE, show_y_text = FALSE)

p <- (D + I) / (F + U)
p1 <- p + plot_annotation(theme = theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold")), 
						  tag_levels = 'A', tag_suffix = ')')

setwd("E:/1_Porth_Ecol/1_Paper/Figures/Final")
ggsave(file = "Fig2.tiff", plot = p1, width = 15, height = 15, dpi = 600, units = "cm", device = "tiff")
ggsave(file = "Fig2.jpg", plot = p1, width = 15, height = 15, dpi = 600, units = "cm", device = "jpg")
ggsave(file = "Fig2.pdf", plot = p1, width = 15, height = 15, dpi = 600, units = "cm", device = "pdf")

# ------------------------------------------------------------------------------------------------ #

### EndNotRun
