# Packages ####
library(sf)
library(inlabru)
library(terra)
library(ggplot2)
library(MASS)
library(dplyr)
library(viridis)
library(cowplot)
library(tidyr)
library(exactextractr)

# Load NFES Data ####
setwd("C:/temp/pests_analysis/National_Forest_Estate_Subcompartments_GB_2016")

NFES <- st_read("National_Forest_Estate_Subcompartments_GB_2016.shp")

NFES_Wales <- st_read("NRW_PRODUCTISED_SCDB_LLE.shp")

species <- c("Acer pseudoplatanus",
             "Betula pendula",
             "Fagus sylvatica",
             "Fraxinus excelsior",
             "Picea abies",
             "Picea sitchensis",
             "Pinus sylvestris",
             "Quercus robur",
             "Sorbus aucuparia")

species_NFES <- c("Sycamore",
                  "Birch (downy/silver)",
                  "Silver birch",
                  "Beech",
                  "Ash",
                  "Norway spruce",
                  "Sitka spruce",
                  "Scots pine",
                  "Oak (robur/petraea)",
                  "Pedunculate/common oak",
                  "Rowan")

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

dm <- c(106,194)
pred.df5 <- fm_pixels(meshGB5, 
                      dims = dm,
                      format = "terra",
                      mask = GB_bound)
# Process NFES data ####
NFES_sub <- NFES %>% filter(PRISPECIES %in% species_NFES |
                              SECSPECIES %in% species_NFES |
                              TERSPECIES %in% species_NFES) %>% select(OBJECTID,
                                                                       PRISPECIES,
                                                                       SECSPECIES,
                                                                       TERSPECIES,
                                                                       PRIPCTAREA,
                                                                       SECPCTAREA,
                                                                       TERPCTAREA,
                                                                       geometry) 

NFES_Wales_sub <- NFES_Wales %>% filter(prispecies %in% species_NFES |
                                          secspecies %in% species_NFES |
                                          terspecies %in% species_NFES) %>% select(objectid,
                                                                                   prispecies,
                                                                                   secspecies,
                                                                                   terspecies,
                                                                                   pripctarea,
                                                                                   secpctarea,
                                                                                   terpctarea,
                                                                                   geometry)

colnames(NFES_Wales_sub) <- colnames(NFES_sub)

NFES_sub <- rbind(NFES_sub, NFES_Wales_sub)

NFES_sub_long <- NFES_sub %>% 
  pivot_longer(
    cols = c(PRISPECIES, SECSPECIES, TERSPECIES),
    names_to = "level",
    values_to = "SPECIES"
  ) %>%
  filter(SPECIES %in% species_NFES)

NFES_sub_long$AREA <- NA
pb = txtProgressBar(min = 0, max = nrow(NFES_sub_long), initial = 0, style = 3) # Set up progress bar
for(i in 1:nrow(NFES_sub_long)){
  if(NFES_sub_long$level[i] == "PRISPECIES"){
    NFES_sub_long$AREA[i] <- NFES_sub_long$PRIPCTAREA[i]
  }else if(NFES_sub_long$level[i] == "SECSPECIES"){
    NFES_sub_long$AREA[i] <- NFES_sub_long$SECPCTAREA[i]
  }else if(NFES_sub_long$level[i] == "TERSPECIES"){
    NFES_sub_long$AREA[i] <- NFES_sub_long$TERPCTAREA[i]
  }
  setTxtProgressBar(pb, i)
  close(pb)
}

# Recode species variable ####
NFES_sub_long$SPECIES[NFES_sub_long$SPECIES == "Ash"] <- "Fraxinus excelsior"
NFES_sub_long$SPECIES[NFES_sub_long$SPECIES == "Beech"] <- "Fagus sylvatica"
NFES_sub_long$SPECIES[NFES_sub_long$SPECIES == "Birch (downy/silver)"] <- "Betula pendula"
NFES_sub_long$SPECIES[NFES_sub_long$SPECIES == "Norway spruce"] <- "Picea abies"
NFES_sub_long$SPECIES[NFES_sub_long$SPECIES == "Oak (robur/petraea)"] <- "Quercus robur"
NFES_sub_long$SPECIES[NFES_sub_long$SPECIES == "Pedunculate/common oak"] <- "Quercus robur"
NFES_sub_long$SPECIES[NFES_sub_long$SPECIES == "Rowan"] <- "Sorbus aucuparia"
NFES_sub_long$SPECIES[NFES_sub_long$SPECIES == "Scots pine"] <- "Pinus sylvestris"
NFES_sub_long$SPECIES[NFES_sub_long$SPECIES == "Silver birch"] <- "Betula pendula"
NFES_sub_long$SPECIES[NFES_sub_long$SPECIES == "Sitka spruce"] <- "Picea sitchensis"
NFES_sub_long$SPECIES[NFES_sub_long$SPECIES == "Sycamore"] <- "Acer pseudoplatanus"

# Calculate species cover for each cell ####
setwd("C:/temp/pests_analysis/National_Forest_Estate_Subcompartments_GB_2016")
for(sp in 1:length(species)){
  
  NFES_sp <- NFES_sub_long %>% filter(SPECIES == species[sp]) %>% select(AREA, geometry) %>% st_transform(crs(pred.df5))
  result <- exact_extract(pred.df5, NFES_sp, include_xy = T)
  
  area_values <- NFES_sp %>% st_drop_geometry() %>% select(AREA)
  
  to_remove <- NA
  for(i in seq_along(result)){
    if(nrow(result[[i]]) == 0){
      to_remove <- c(to_remove, i)
    }
  }
  to_remove <- to_remove[-1]
  
  if(length(to_remove) > 0){
    result <- result[-to_remove]
    area_values <- area_values[-to_remove,]
  }
  
  for(i in seq_along(result)){
    val <- area_values$AREA[i]
    result[[i]] <- cbind(result[[i]], val)
  }
  
  result <- do.call(rbind.data.frame, result)
  
  result %>% group_by(x, y) %>% 
    mutate(pixel_cov_fraction = sum(coverage_fraction),
           frac_value = val * coverage_fraction/pixel_cov_fraction) %>% 
    summarise(final_value = sum(frac_value, na.rm = T)) %>% 
    ungroup() -> for_raster
  
  coords <- as.matrix(for_raster[,c('x','y')])
  
  r <- rasterize(x = coords, y = pred.df5,
                 values = for_raster$final_value,
                 background = NA)
  r <- crop(r, pred.df5, mask = TRUE)
  
  writeRaster(r, filename = paste0(species[sp], "_NFES_cover.tif"))
}

# Load cover rasters ####
setwd("C:/temp/pests_analysis/National_Forest_Estate_Subcompartments_GB_2016")
file_list <- list.files(pattern = ".tif")
NFES_cover <- rast(file_list)
names(NFES_cover) <- gsub("_NFES_cover.tif", "", file_list)


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



predictions_list$`Acer pseudoplatanus`$abundance <- 
  eval_spatial(NFES_cover$`Acer pseudoplatanus`,
               predictions_list$`Acer pseudoplatanus`$geometry)

predictions_list$`Betula pendula`$abundance <- 
  eval_spatial(NFES_cover$`Betula pendula`,
               predictions_list$`Betula pendula`$geometry)

predictions_list$`Fagus sylvatica`$abundance <- 
  eval_spatial(NFES_cover$`Fagus sylvatica`,
               predictions_list$`Fagus sylvatica`$geometry)

predictions_list$`Fraxinus excelsior`$abundance <- 
  eval_spatial(NFES_cover$`Fraxinus excelsior`,
               predictions_list$`Fraxinus excelsior`$geometry)

predictions_list$`Picea abies`$abundance <-
  eval_spatial(NFES_cover$`Picea abies`,
              predictions_list$`Picea abies`$geometry)

predictions_list$`Picea sitchensis`$abundance <-
  eval_spatial(NFES_cover$`Picea sitchensis`,
               predictions_list$`Picea sitchensis`$geometry)

predictions_list$`Pinus sylvestris`$abundance <- 
  eval_spatial(NFES_cover$`Pinus sylvestris`,
               predictions_list$`Pinus sylvestris`$geometry)

predictions_list$`Quercus robur`$abundance <- 
  eval_spatial(NFES_cover$`Quercus robur`,
               predictions_list$`Quercus robur`$geometry)

predictions_list$`Sorbus aucuparia`$abundance <- 
  eval_spatial(NFES_cover$`Sorbus aucuparia`,
               predictions_list$`Sorbus aucuparia`$geometry)

# Rescale intensity to 0-1 ####
for(i in 1:length(species)){
  predictions_list[[i]]$mean <- (predictions_list[[i]]$mean - min(predictions_list[[i]]$mean)) / (max(predictions_list[[i]]$mean) - min(predictions_list[[i]]$mean))
}


# Remove cells with NA abundance ####
predictions_list$`Acer pseudoplatanus` <- predictions_list$`Acer pseudoplatanus` %>% filter(!is.na(abundance))
predictions_list$`Betula pendula` <- predictions_list$`Betula pendula` %>% filter(!is.na(abundance))
predictions_list$`Fagus sylvatica` <- predictions_list$`Fagus sylvatica` %>% filter(!is.na(abundance))
predictions_list$`Fraxinus excelsior` <- predictions_list$`Fraxinus excelsior` %>% filter(!is.na(abundance))
predictions_list$`Picea abies` <- predictions_list$`Picea abies` %>% filter(!is.na(abundance))
predictions_list$`Picea sitchensis` <- predictions_list$`Picea sitchensis` %>% filter(!is.na(abundance))
predictions_list$`Pinus sylvestris` <- predictions_list$`Pinus sylvestris` %>% filter(!is.na(abundance))
predictions_list$`Quercus robur` <- predictions_list$`Quercus robur` %>% filter(!is.na(abundance))
predictions_list$`Sorbus aucuparia` <- predictions_list$`Sorbus aucuparia` %>% filter(!is.na(abundance))


# Use robust regression to predict P+D intensity from abundance, then take residuals ####
predictions_list$`Acer pseudoplatanus`$resid <- resid(rlm(mean ~ abundance, data = predictions_list$`Acer pseudoplatanus`))
predictions_list$`Betula pendula`$resid <- resid(rlm(mean ~ abundance, data = predictions_list$`Betula pendula`))
predictions_list$`Fagus sylvatica`$resid <- resid(rlm(mean ~ abundance, data = predictions_list$`Fagus sylvatica`))
predictions_list$`Fraxinus excelsior`$resid <- resid(rlm(mean ~ abundance, data = predictions_list$`Fraxinus excelsior`))
predictions_list$`Picea abies`$resid <- resid(rlm(mean ~ abundance, data = predictions_list$`Picea abies`))
predictions_list$`Picea sitchensis`$resid <- resid(rlm(mean ~ abundance, data = predictions_list$`Picea sitchensis`))
predictions_list$`Pinus sylvestris`$resid <- resid(rlm(mean ~ abundance, data = predictions_list$`Pinus sylvestris`))
predictions_list$`Quercus robur`$resid <- resid(rlm(mean ~ abundance, data = predictions_list$`Quercus robur`))
predictions_list$`Sorbus aucuparia`$resid <- resid(rlm(mean ~ abundance, data = predictions_list$`Sorbus aucuparia`))


# Plot residuals ####
sp_to_plot <- species

colpal <- c("purple4", "darkblue", "blue", "cyan3", "white", "orange", "red","red3", "darkred")
breaks <- c(-0.8 ,-0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8)


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
p_list[[10]] <- resid_legend

grid1 <- plot_grid(plotlist = p_list[c(1,2,3,4,8,9,5,6,7,10)], nrow = 2, labels = c("A","B","C","D","E","F","G","H","I",""))

setwd("C:/temp/pests_analysis/figures_revisions")
ggsave2(filename = "intensity_vs_abundance_NFES.tiff",
        plot = grid1,
        bg = "white",
        dpi = 600,
        width = 165,
        height = 110,
        units = "mm",
        scale = 1.7)






