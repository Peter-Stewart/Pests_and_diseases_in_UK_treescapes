# Packages ####
library(sf)
library(inlabru)
library(terra)
library(ggplot2)
library(MASS)
library(dplyr)
library(viridis)
library(cowplot)

# Standardize function from 'rethinking' package
standardize <- function(x){
  x <- scale(x)
  z <- as.numeric(x)
  attr(z, "scaled:center") <- attr(x, "scaled:center")
  attr(z, "scaled:scale") <- attr(x, "scaled:scale")
  return(z)
}

# Load data ####
setwd("C:/temp/pests_analysis/INLA_inputs")
GB_bound <- get(load("gb_mainland_boundary.Rdata"))
meshGB5 <- get(load("meshGB5.Rdata"))
covariate_stack <- rast("covariate_stack.tif")
correct <- covariate_stack$broadleaf_area_Ha
covariate_stack[is.na(covariate_stack)] <- 0
covariate_stack <- mask(covariate_stack, 
                        correct,
                        maskvalues = NA,
                        updatevalue = NA)

# Load saved model predictions ####
species <- c("Acer pseudoplatanus",
             "Betula pendula",
             "Fagus sylvatica",
             "Fraxinus excelsior",
             "Picea abies",
             "Picea sitchensis",
             "Pinus sylvestris",
             "Quercus robur",
             "Sorbus aucuparia")

setwd("D:/INLA_outputs/mesh5")
file_list <- list.files(pattern = "predictions.Rdata")
predictions_list <- list()
for(i in 1:length(file_list)){
  predictions_list[[i]] <- get(load(file_list[i]))
}
predictions_list <- predictions_list[-10]
names(predictions_list) <- species


# Load predicted host abundance layers ####
species_short <- c("Aps", # Acer pseudoplatanus
                   "Bpe", # Betula pendula
                   "Fsy", # Fagus sylvatica
                   "Fex", # Fraxinus excelsior
                   "Qro") # Quercus robur

files_short <- paste0(species_short, ".grd")
setwd("C:/temp/pests_analysis/Predicted_abundance_maps_rasters")
host_abundance <- rast(files_short)

# Project and clip predicted host abundance to match model predictions ####
dm <- c(106,194)
pred.df5 <- fm_pixels(meshGB5, 
                      dims = dm,
                      format = "terra",
                      mask = GB_bound)

host_abundance <- project(host_abundance, correct)
host_abundance[is.na(host_abundance)] <- 0
host_abundance <- mask(host_abundance, 
                       correct,
                       maskvalues = NA,
                       updatevalue = NA)


host_abundance_proj <- project(host_abundance, pred.df5)

predictions_list$`Acer pseudoplatanus`$abundance <- 
  eval_spatial(host_abundance_proj,
               predictions_list$`Acer pseudoplatanus`$geometry)

predictions_list$`Betula pendula`$abundance <- 
  eval_spatial(host_abundance_proj,
               predictions_list$`Betula pendula`$geometry)

predictions_list$`Fagus sylvatica`$abundance <- 
  eval_spatial(host_abundance_proj,
               predictions_list$`Fagus sylvatica`$geometry)

predictions_list$`Fraxinus excelsior`$abundance <- 
  eval_spatial(host_abundance_proj,
               predictions_list$`Fraxinus excelsior`$geometry)

predictions_list$`Quercus robur`$abundance <- 
  eval_spatial(host_abundance_proj,
               predictions_list$`Quercus robur`$geometry)


# Remove cells with NA abundance ####
predictions_list$`Acer pseudoplatanus` <- predictions_list$`Acer pseudoplatanus` %>% filter(!is.na(abundance))
predictions_list$`Betula pendula` <- predictions_list$`Betula pendula` %>% filter(!is.na(abundance))
predictions_list$`Fagus sylvatica` <- predictions_list$`Fagus sylvatica` %>% filter(!is.na(abundance))
predictions_list$`Fraxinus excelsior` <- predictions_list$`Fraxinus excelsior` %>% filter(!is.na(abundance))
predictions_list$`Quercus robur` <- predictions_list$`Quercus robur` %>% filter(!is.na(abundance))


# Use robust regression to predict P+D intensity from abundance, then take residuals ####
predictions_list$`Acer pseudoplatanus`$resid <- resid(rlm(standardize(mean) ~ standardize(abundance), data = predictions_list$`Acer pseudoplatanus`))
predictions_list$`Betula pendula`$resid <- resid(rlm(standardize(mean) ~ standardize(abundance), data = predictions_list$`Betula pendula`))
predictions_list$`Fagus sylvatica`$resid <- resid(rlm(standardize(mean) ~ standardize(abundance), data = predictions_list$`Fagus sylvatica`))
predictions_list$`Fraxinus excelsior`$resid <- resid(rlm(standardize(mean) ~ standardize(abundance), data = predictions_list$`Fraxinus excelsior`))
predictions_list$`Quercus robur`$resid <- resid(rlm(standardize(mean) ~ standardize(abundance), data = predictions_list$`Quercus robur`))


# Plot residuals ####
sp_to_plot <- species[c(1:4, 8)]

colpal <- c("purple4", "darkblue", "blue", "cyan3", "white", "orange", "red","red3", "darkred")
breaks <- c(-6 ,-3, -1.5, -0.5, 0.5, 1.5, 3, 6)

p_list <- list()
for(i in 1:length(sp_to_plot)){
  tmpdf <- predictions_list[[sp_to_plot[i]]]
  tmpdf <- tmpdf %>% select(resid, geometry)
  p_temp <- ggplot() +
    gg(tmpdf, geom = "tile") +
    geom_sf(data = GB_bound, fill = NA, lwd = 0.8, col = "black") +
   
    #ggtitle(paste(sp_to_plot[i])) +

    xlab("") + ylab("") +
    theme(legend.position = "none") +
    
    scale_fill_gradientn(colors = colpal, 
                         #values = scales::rescale(breaks),
                         limits = range(breaks))
    
  p_list[[i]] <- p_temp
  
}
names(p_list) <- sp_to_plot

# Get legend ####
tmpdf <- predictions_list$`Fraxinus excelsior`
tmpdf <- tmpdf %>% select(resid, geometry)
p_temp <- ggplot() +
  gg(tmpdf, geom = "tile") +
  geom_sf(data = GB_bound, fill = NA, lwd = 0.8, col = "black") +
  
  #ggtitle(paste(sp_to_plot[i])) +
  
  xlab("") + ylab("") +
  #theme(legend.position = "none") +
  
  scale_fill_gradientn(colors = colpal, 
                       values = scales::rescale(breaks),
                       limits = range(breaks),
                       name = "Residuals") 

resid_legend <- cowplot::get_legend(p_temp)

# Save plot ####
p_list[[6]] <- resid_legend

grid1 <- plot_grid(plotlist = p_list, nrow = 1, labels = c("A","B","C","D","E",""), rel_widths = c(1,1,1,1,1,0.4))

setwd("C:/temp/pests_analysis/figures_revisions")
ggsave2(filename = "intensity_vs_abundance.tiff",
        plot = grid1,
        bg = "white",
        dpi = 600,
        width = 165,
        height = 55,
        units = "mm",
        scale = 1.7)
