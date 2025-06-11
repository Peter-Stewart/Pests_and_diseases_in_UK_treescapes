# Packages ####
library(PointedSDMs)
library(terra)
library(ggplot2)
library(dplyr)
library(viridis)
library(cowplot)

# Function to automatically run ISDM multiple times with different input options
# and save all outputs to a single list
run_multi_ISDM <- function(df_list, 
                           species, 
                           mesh,
                           covs, 
                           covs_select = NULL,
                           priors_fixed = c(0, 1),
                           priors_pc = c(5, 0.01, 5000, 0.01),
                           predict_data = NULL,
                           save_predictions = FALSE,
                           only_fixed_effects = FALSE
){
  
  # for sp in species
  # for m in meshes
  # if covs != NULL then for cv in covs_select
  # for pr in priors?
  
  # Create list to store results
  results_list <- list()
  
  # Set up progress bar
  progress_bar = txtProgressBar(min=0, max=length(species), style = 3, char="=")
  
  # Loop over
  for(sp in 1:length(species)){
    
    # Create second list to store results for each species so results_list ends up nested
    sp_results_list <- list()
    
    # Subset df_list to given species
    sp_data <- df_list[species[sp]][[1]]
    
    # If no data for a species then that species' results list is NA
    if(length(sp_data) != 2){
      sp_results_list <- NA 
      results_list[[sp]] <- sp_results_list
      next
    }
    
    # Loop over list of covariate options
    for(cv in 1:length(covs_select)){
      cov_selected <- terra::subset(covs, covs_select[[cv]])
      
      # Create R6 object
      dat <- startISDM(sp_data,
                      Coordinates = c("X", "Y"), 
                      Projection = "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs",
                      Mesh = mesh,
                      responsePA = "Present",
                      trialsPA = "Surveyed",
                      spatialCovariates = cov_selected)
      
      # Specify type and priors for spatial field
      dat$specifySpatial(sharedSpatial = TRUE,
                         PC = TRUE,
                         prior.sigma = c(priors_pc[1], priors_pc[2]),
                         prior.range = c(priors_pc[3], priors_pc[4]))
      
      # Assign weakly regularising priors to intercept terms
      dat$priorsFixed(Effect = "intercept", 
                      datasetName = "THDAS_obs", 
                      mean.linear = priors_fixed[1],
                      prec.linear = priors_fixed[2])
      dat$priorsFixed(Effect = "intercept", 
                      datasetName = "APHA_obs",
                      mean.linear = priors_fixed[1],
                      prec.linear = priors_fixed[2])
      
      # Assign weakly regularising priors to effect sizes
      for(i in 1:length(names(cov_selected))){
        dat$priorsFixed(Effect = names(cov_selected)[i],
                        mean.linear = priors_fixed[1],
                        prec.linear = priors_fixed[2])
      }
      
      # Add bias field for the presence-only data
      dat$addBias("THDAS_obs")
      
      # Fit model
      mod <- fitISDM(data = dat, 
                     options = list(control.inla = list(int.strategy = "eb")))
      
      # Optionally generate and save predictions on predict_data
      if(save_predictions == TRUE){
        predictions <- predict(mod, predict_data, 
                               spatial = FALSE, 
                               biasfield = FALSE, 
                               predictor = TRUE,
                               fun = 'linear', n.samples = 1000)
        predictions_df <- as.data.frame(predictions$predictions)
        predictions_df_sf <- st_as_sf(predictions_df)
        save(predictions_df_sf, file = paste0(species[sp],"_","option",cv,"_predictions.Rdata"))
        save(mod, file = paste0(species[sp],"_","option",cv,"_model.Rdata"))
      }
      
      # Append model results to species results list
      if(only_fixed_effects == TRUE){
        sp_results_list[[cv]] <- mod$summary.fixed
      }else{
        sp_results_list[[cv]] <- mod
      }
      
      # Name of list element is covariate option number
      names(sp_results_list)[cv] <- paste("Option", cv)
    }
    
    # Append species results to main results list
    results_list[[sp]] <- sp_results_list
    
    # Update progress bar
    setTxtProgressBar(progress_bar, value = sp)
    
  }
  # Names for list elements are species names
  names(results_list) <- species
  
  # Return the results list
  close(progress_bar) 
  return(results_list)
}

# Function to plot effect size estimates for multiple models
effects_plot <- function(..., matrix_input = TRUE, effect_name = NULL, pchs = 17:23, v_off = 0.2, title = NULL){
  mods <- list(...)
  
  mat_list <- list()
  
  # Extract fixed effects estimates from each model
  for(m in 1:length(mods)){
    mat <- matrix(NA, nrow = length(mods[[m]]), ncol = 3)
    
    if(matrix_input == TRUE){
      for(i in 1:length(mods[[m]])){
        mat[i,1] <- mods[[m]][[i]][effect_name,"mean"]
        mat[i,2] <- mods[[m]][[i]][effect_name,"0.025quant"]
        mat[i,3] <- mods[[m]][[i]][effect_name,"0.975quant"]
      }
    }else{
      for(i in 1:length(mods[[m]])){
        mat[i,1] <- mods[[m]][[i]]$summary.fixed[effect_name,"mean"]
        mat[i,2] <- mods[[m]][[i]]$summary.fixed[effect_name,"0.025quant"]
        mat[i,3] <- mods[[m]][[i]]$summary.fixed[effect_name,"0.975quant"]
      }
    }
    mat_list[[m]] <- mat
  }
  
  # Create empty plot
  par(mar = c(5.1, 7, 4.1, 2.1))
  plot(NULL, 
       ylim = c(0, nrow(mat_list[[1]])+1), 
       xlim = c(min(unlist(mat_list))-0.025, 
                max(unlist(mat_list))+0.025),
       xlab = "Estimate", ylab = "",
       yaxt = "n",
       main = title)
  abline(v = 0, lty = 2)
  labs <- vector("character", length=nrow(mat_list[[1]]))
  for(i in 1:nrow(mat_list[[1]])){
    labs[i] <- paste0("m",i)
  }
  axis(side = 2, 
       at = seq(1,nrow(mat_list[[1]]),1),
       labels = labs,
       tick = TRUE,
       las = 2)
  
  for(m in 1:length(mods)){
    off <- -v_off*(1-m)
    for(i in 1:nrow(mat_list[[m]])){
      points(x = mat_list[[m]][i,1], y = i+off, pch = pchs[m], col = "black")
      lines(x = c(mat_list[[m]][i,2], mat_list[[m]][i,3]),
            y = c(i+off,i+off),
            lwd = 1)
    }
  }
}


# Load data ####
setwd("C:/temp/pests_analysis/INLA_inputs")

pest_data <- get(load("pest_data_ISDM_input_list.Rdata"))
meshGB1 <- get(load("meshGB1.Rdata"))
meshGB5 <- get(load("meshGB5.Rdata"))


covariate_stack <- rast("covariate_stack.tif")

GB_bound <- get(load("gb_mainland_boundary.Rdata"))
setwd("C:/temp/pests_analysis/NUTS_Level_1_January_2018_GCB_in_the_United_Kingdom_2022_4034932022524175836")
gb <- st_read("NUTS_Level_1_January_2018_GCB_in_the_United_Kingdom.shp")
gb <- st_transform(gb, "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs")
gb_main <- gb %>% filter(nuts118nm != "Northern Ireland") 
gb_all <- st_union(gb_main) # Combine the regions


# Ensure consistent NA values across raster layers ####
# separate out a layer which is correct
correct <- covariate_stack$broadleaf_area_Ha

covariate_stack[is.na(covariate_stack)] <- 0

covariate_stack <- mask(covariate_stack, 
                        correct,
                        maskvalues = NA,
                        updatevalue = NA)

# Standardise raster layers ####
# 0-centred standardise 
covariate_stack_sd <- terra::scale(covariate_stack)

# Load covariate lists ####
setwd("C:/temp/pests_analysis/ISDM_flists")
file_list <- list.files(pattern = "flist")
for(i in 1:length(file_list)){
  obj_name <- gsub(".Rdata","",file_list[i])
  temp_obj <- get(load(file_list[i]))
  assign(obj_name, temp_obj)
}

# Fit models - N.B. that this was done on JASMIN for each species separately ####
species <- c("Acer pseudoplatanus",
             "Betula pendula",
             "Fagus sylvatica",
             "Fraxinus excelsior",
             "Picea abies",
             "Picea sitchensis",
             "Pinus sylvestris",
             "Quercus robur",
             "Sorbus aucuparia")


# Recreation ####
setwd("F:/INLA_outputs/ISDM_drivers")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = species,
                      mesh = meshGB1,
                      covs = covariate_stack_sd,
                      covs_select = R_flist, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      save_predictions = FALSE,
                      only_fixed_effects = FALSE)
mod.time<-proc.time()-ptm
mod.time[3]/60
save(mod, file = "R_mods1.Rdata")


# Urbanisation ####
setwd("F:/INLA_outputs/ISDM_drivers")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = species,
                      mesh = meshGB1,
                      covs = covariate_stack_sd,
                      covs_select = Ur_flist, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      save_predictions = FALSE,
                      only_fixed_effects = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60
save(mod, file = "Ur_mods1.Rdata")

# Human pop. density ####
setwd("F:/INLA_outputs/ISDM_drivers")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = species,
                      mesh = meshGB1,
                      covs = covariate_stack_sd,
                      covs_select = H_flist, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      save_predictions = FALSE,
                      only_fixed_effects = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60
save(mod, file = "H_mods1.Rdata")

# Afforestation ####
setwd("F:/INLA_outputs/ISDM_drivers")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = species,
                      mesh = meshGB1,
                      covs = covariate_stack_sd,
                      covs_select = Af_flist, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      save_predictions = FALSE,
                      only_fixed_effects = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60
save(mod, file = "Af_mods1.Rdata")

setwd("F:/INLA_outputs/ISDM_drivers")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = species,
                      mesh = meshGB5,
                      covs = covariate_stack_sd,
                      covs_select = Af_flist, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      save_predictions = FALSE,
                      only_fixed_effects = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60
save(mod, file = "Af_mods5.Rdata")

# Deforestation ####
setwd("F:/INLA_outputs/ISDM_drivers")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = species,
                      mesh = meshGB1,
                      covs = covariate_stack_sd,
                      covs_select = Df_flist, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      save_predictions = FALSE,
                      only_fixed_effects = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60
save(mod, file = "Df_mods1.Rdata")

# Dist. to BCP ####
setwd("F:/INLA_outputs/ISDM_drivers")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = species,
                      mesh = meshGB1,
                      covs = covariate_stack_sd,
                      covs_select = d_BCP_flist, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      save_predictions = FALSE,
                      only_fixed_effects = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60
save(mod, file = "d_BCP_mods1.Rdata")

# Conifer area ####
setwd("F:/INLA_outputs/ISDM_drivers")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = species,
                      mesh = meshGB1,
                      covs = covariate_stack_sd,
                      covs_select = C_a_flist, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      save_predictions = FALSE,
                      only_fixed_effects = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60
save(mod, file = "C_a_mods1.Rdata")

# Woodland connectivity ####
setwd("F:/INLA_outputs/ISDM_drivers")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = species,
                      mesh = meshGB1,
                      covs = covariate_stack_sd,
                      covs_select = W_c_flist, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      save_predictions = FALSE,
                      only_fixed_effects = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60
save(mod, file = "W_c_mods1.Rdata")

# Plot effects - mesh 1 ####
species <- c("Acer pseudoplatanus",
             "Betula pendula",
             "Fagus sylvatica",
             "Fraxinus excelsior",
             "Picea abies",
             "Picea sitchensis",
             "Pinus sylvestris",
             "Quercus robur",
             "Sorbus aucuparia")

setwd("D:/INLA_outputs/ISDM_drivers/species_results/mesh1")
file_list <- list.files(pattern = "_mods")
for(i in 1:length(file_list)){
  obj_name <- gsub(".Rdata","",file_list[i])
  temp_obj <- get(load(file_list[i]))
  assign(obj_name, temp_obj)
}

R_mods1 <- unlist(unlist(list(mget(ls(pattern = "R_mods1"))), recursive = FALSE, use.names = FALSE), recursive = FALSE, use.names = TRUE)
R2_mods1 <- unlist(unlist(list(mget(ls(pattern = "R2_mods1"))), recursive = FALSE, use.names = FALSE), recursive = FALSE, use.names = TRUE)
Ur_mods1 <- unlist(unlist(list(mget(ls(pattern = "Ur_mods1"))), recursive = FALSE, use.names = FALSE), recursive = FALSE, use.names = TRUE)
H_mods1 <- unlist(unlist(list(mget(ls(pattern = "H_mods1"))), recursive = FALSE, use.names = FALSE), recursive = FALSE, use.names = TRUE)
Af_mods1 <- unlist(unlist(list(mget(ls(pattern = "Af_mods1"))), recursive = FALSE, use.names = FALSE), recursive = FALSE, use.names = TRUE)
Df_mods1 <- unlist(unlist(list(mget(ls(pattern = "Df_mods1"))), recursive = FALSE, use.names = FALSE), recursive = FALSE, use.names = TRUE)
d_BCP_mods1 <- unlist(unlist(list(mget(ls(pattern = "d_BCP_mods1"))), recursive = FALSE, use.names = FALSE), recursive = FALSE, use.names = TRUE)
C_a_mods1 <- unlist(unlist(list(mget(ls(pattern = "C_a_mods1"))), recursive = FALSE, use.names = FALSE), recursive = FALSE, use.names = TRUE)
W_c_mods1 <- unlist(unlist(list(mget(ls(pattern = "W_c_mods1"))), recursive = FALSE, use.names = FALSE), recursive = FALSE, use.names = TRUE)

species <- c("Acer pseudoplatanus",
             "Betula pendula",
             "Fagus sylvatica",
             "Fraxinus excelsior",
             "Quercus robur",
             "Sorbus aucuparia",
             "Picea abies",
             "Picea sitchensis",
             "Pinus sylvestris")

# Recreation
mod <- R_mods1
par(mfrow = c(2,5))
for(sp in 1:length(species)){
  effects_plot(mod[[species[sp]]], 
               effect_name = "recreation_weekly_1km",
               pchs = c(15,16,17,18,1),
               title = species[sp])
}

# Urbanisation
mod <- Ur_mods1
par(mfrow = c(2,5))
for(sp in 1:length(species)){
  effects_plot(mod[[species[sp]]], 
               effect_name = "urb_suburb_areaHa",
               pchs = c(15,16,17,18,1),
               title = species[sp])
}

# Human pop. density
mod <- H_mods1
par(mfrow = c(2,5))
for(sp in 1:length(species)){
  effects_plot(mod[[species[sp]]], 
               effect_name = "humanPop",
               pchs = c(15,16,17,18,1),
               title = species[sp])
}

# Afforestation 
mod <- Af_mods1
par(mfrow = c(2,5))
for(sp in 1:length(species)){
  effects_plot(mod[[species[sp]]], 
               effect_name = "LCC_afforestation_1km",
               pchs = c(15,16,17,18,1),
               title = species[sp])
}

# Deforestation 
mod <- Df_mods1
par(mfrow = c(2,5))
for(sp in 1:length(species)){
  effects_plot(mod[[species[sp]]], 
               effect_name = "LCC_deforestation_1km",
               pchs = c(15,16,17,18,1),
               title = species[sp])
}

# Dist. to BCP
mod <- d_BCP_mods1
par(mfrow = c(2,5))
for(sp in 1:length(species)){
  effects_plot(mod[[species[sp]]], 
               effect_name = "dist_BCP",
               pchs = c(15,16,17,18,1),
               title = species[sp])
}

# Conifer area
mod <- C_a_mods1
par(mfrow = c(2,5))
for(sp in 1:length(species)){
  effects_plot(mod[[species[sp]]], 
               effect_name = "conifer_area_Ha",
               pchs = c(15,16,17,18,1),
               title = species[sp])
}

# Woodland connectivity
mod <- W_c_mods1
par(mfrow = c(2,5))
for(sp in 1:length(species)){
  effects_plot(mod[[species[sp]]], 
               effect_name = "woodland_connectivity_sigma3",
               pchs = c(15,16,17,18,1),
               title = species[sp])
}

# Plot effects from maximal DAG for all key drivers
effect_list <- c("recreation_weekly_1km",
                 "recreation_yearly_1km",
                 "urb_suburb_areaHa",
                 "humanPop",
                 "LCC_afforestation_1km",
                 "LCC_deforestation_1km",
                 "dist_BCP",
                 "conifer_area_Ha",
                 "woodland_connectivity_sigma3")

effect_labels <- c("Recreation (w)",
                   "Recreation (y)",
                   "Urban area",
                   "Human pop.",
                   "Afforestation",
                   "Deforestation",
                   "Dist. to BCP",
                   "Conifer area",
                   "Connectivity")

mod_list <- c("R_mods1",
              "R2_mods1",
              "Ur_mods1",
              "H_mods1",
              "Af_mods1",
              "Df_mods1",
              "d_BCP_mods1",
              "C_a_mods1",
              "W_c_mods1")

species_labels <- c("A. pseudoplatanus",
                    "B. pendula",
                    "F. sylvatica",
                    "F. excelsior",
                    "Q. robur",
                    "S. aucuparia",
                    "P. abies",
                    "P. sitchensis",
                    "P. sylvestris")

# Plot with separate panel for each effect, showing all species per effect
setwd("C:/temp/pests_analysis/figures_new")
tiff("drivers_mesh1.tiff", width = 165, height = 160, units = 'mm', res = 600)
#par(pr)

layout(matrix(c(1:9), nrow = 3, ncol = 3, byrow = TRUE),
       widths = c(30,19,19))

for(m in 1:length(effect_list)){
  # Create matrix of effects for each species
  mat <- matrix(NA, nrow = length(species), ncol = 3)
  
  for(sp in 1:length(species)){
    mod <- get(mod_list[m], envir = .GlobalEnv)
    mod <- mod[[species[sp]]]
    
    mat[sp,1] <- mod[[length(mod)]][effect_list[m], "mean"]
    mat[sp,2] <- mod[[length(mod)]][effect_list[m], "0.025quant"]
    mat[sp,3] <- mod[[length(mod)]][effect_list[m], "0.975quant"]
  }
  
  if(m %in% c(1,4,7)){
    par(mar = c(4.1, 8, 2.1, 0.5))
  }else{
    par(mar = c(4.1, 1, 2.1, 0.5))
  }
  
  # Create empty plot
  plot(NULL, 
       ylim = c(0, nrow(mat)+1), 
       xlim = c(min(unlist(mat))-0.025, 
                max(unlist(mat))+0.025),
       xlab = "", ylab = "",
       yaxt = "n",
       main = effect_labels[m])
  abline(v = 0, lty = 2)
  title(LETTERS[m], adj = 0)
  
  if(m %in% 7:9){
    title(xlab="Estimate", line=2.6, cex.lab=1.2)
  }
  
  if(m %in% c(1,4,7)){
    axis(side = 2, 
         at = seq(nrow(mat),1,-1),
         labels = parse(text = paste0("italic('", species_labels, "')")),
         tick = TRUE,
         las = 2,
         cex.axis = 1)
  }else{
    axis(side = 2,
         at = seq(nrow(mat),1,-1),
         labels = rep("",9),
         tick = TRUE,
         las = 2)
  }
  
  for(i in 1:nrow(mat)){
    points(x = mat[i,1], y = 10-i, pch = 18, col = "black")
    lines(x = c(mat[i,2], mat[i,3]),
          y = rep(10-i, 2),
          lwd = 1)
  }
}
dev.off()


# Plot effects - mesh 5 ####
species <- c("Acer pseudoplatanus",
             "Betula pendula",
             "Fagus sylvatica",
             "Fraxinus excelsior",
             "Picea abies",
             "Picea sitchensis",
             "Pinus sylvestris",
             "Quercus robur",
             "Sorbus aucuparia")

setwd("D:/INLA_outputs/ISDM_drivers/species_results/mesh5")
file_list <- list.files(pattern = "_mods")
for(i in 1:length(file_list)){
  obj_name <- gsub(".Rdata","",file_list[i])
  temp_obj <- get(load(file_list[i]))
  assign(obj_name, temp_obj)
}

R_mods5 <- unlist(unlist(list(mget(ls(pattern = "R_mods5"))), recursive = FALSE, use.names = FALSE), recursive = FALSE, use.names = TRUE)
R2_mods5 <- unlist(unlist(list(mget(ls(pattern = "R2_mods5"))), recursive = FALSE, use.names = FALSE), recursive = FALSE, use.names = TRUE)
Ur_mods5 <- unlist(unlist(list(mget(ls(pattern = "Ur_mods5"))), recursive = FALSE, use.names = FALSE), recursive = FALSE, use.names = TRUE)
H_mods5 <- unlist(unlist(list(mget(ls(pattern = "H_mods5"))), recursive = FALSE, use.names = FALSE), recursive = FALSE, use.names = TRUE)
Af_mods5 <- unlist(unlist(list(mget(ls(pattern = "Af_mods5"))), recursive = FALSE, use.names = FALSE), recursive = FALSE, use.names = TRUE)
Df_mods5 <- unlist(unlist(list(mget(ls(pattern = "Df_mods5"))), recursive = FALSE, use.names = FALSE), recursive = FALSE, use.names = TRUE)
d_BCP_mods5 <- unlist(unlist(list(mget(ls(pattern = "d_BCP_mods5"))), recursive = FALSE, use.names = FALSE), recursive = FALSE, use.names = TRUE)
C_a_mods5 <- unlist(unlist(list(mget(ls(pattern = "C_a_mods5"))), recursive = FALSE, use.names = FALSE), recursive = FALSE, use.names = TRUE)
W_c_mods5 <- unlist(unlist(list(mget(ls(pattern = "W_c_mods5"))), recursive = FALSE, use.names = FALSE), recursive = FALSE, use.names = TRUE)

species <- c("Acer pseudoplatanus",
             "Betula pendula",
             "Fagus sylvatica",
             "Fraxinus excelsior",
             "Quercus robur",
             "Sorbus aucuparia",
             "Picea abies",
             "Picea sitchensis",
             "Pinus sylvestris")

# Recreation (weekly)
mod <- R_mods5
par(mfrow = c(2,5))
for(sp in 1:length(species)){
  effects_plot(mod[[species[sp]]], 
               effect_name = "recreation_weekly_1km",
               pchs = c(15,16,17,18,1),
               title = species[sp])
}

# Recreation (yearly)
mod <- R2_mods5
par(mfrow = c(2,5))
for(sp in 1:length(species)){
  effects_plot(mod[[species[sp]]], 
               effect_name = "recreation_yearly_1km",
               pchs = c(15,16,17,18,1),
               title = species[sp])
}

# Urbanisation
mod <- Ur_mods5
par(mfrow = c(2,5))
for(sp in 1:length(species)){
  effects_plot(mod[[species[sp]]], 
               effect_name = "urb_suburb_areaHa",
               pchs = c(15,16,17,18,1),
               title = species[sp])
}

# Human pop. density
mod <- H_mods5
par(mfrow = c(2,5))
for(sp in 1:length(species)){
  effects_plot(mod[[species[sp]]], 
               effect_name = "humanPop",
               pchs = c(15,16,17,18,1),
               title = species[sp])
}

# Afforestation 
mod <- Af_mods5
par(mfrow = c(2,5))
for(sp in 1:length(species)){
  effects_plot(mod[[species[sp]]], 
               effect_name = "LCC_afforestation_1km",
               pchs = c(15,16,17,18,1),
               title = species[sp])
}

# Deforestation 
mod <- Df_mods5
par(mfrow = c(2,5))
for(sp in 1:length(species)){
  effects_plot(mod[[species[sp]]], 
               effect_name = "LCC_deforestation_1km",
               pchs = c(15,16,17,18,1),
               title = species[sp])
}

# Dist. to BCP
mod <- d_BCP_mods5
par(mfrow = c(2,5))
for(sp in 1:length(species)){
  effects_plot(mod[[species[sp]]], 
               effect_name = "dist_BCP",
               pchs = c(15,16,17,18,1),
               title = species[sp])
}

# Conifer area
mod <- C_a_mods5
par(mfrow = c(2,5))
for(sp in 1:length(species)){
  effects_plot(mod[[species[sp]]], 
               effect_name = "conifer_area_Ha",
               pchs = c(15,16,17,18,1),
               title = species[sp])
}

# Woodland connectivity
mod <- W_c_mods5
par(mfrow = c(2,5))
for(sp in 1:length(species)){
  effects_plot(mod[[species[sp]]], 
               effect_name = "woodland_connectivity_sigma3",
               pchs = c(15,16,17,18,1),
               title = species[sp])
}

# Plot effects from maximal DAG for all key drivers
effect_list <- c("recreation_weekly_1km",
                 "recreation_yearly_1km",
                 "urb_suburb_areaHa",
                 "humanPop",
                 "LCC_afforestation_1km",
                 "LCC_deforestation_1km",
                 "dist_BCP",
                 "conifer_area_Ha",
                 "woodland_connectivity_sigma3")

effect_labels <- c("Recreation (w)",
                   "Recreation (y)",
                   "Urban area",
                   "Human pop.",
                   "Afforestation",
                   "Deforestation",
                   "Dist. to BCP",
                   "Conifer area",
                   "Connectivity")

mod_list <- c("R_mods5",
              "R2_mods5",
              "Ur_mods5",
              "H_mods5",
              "Af_mods5",
              "Df_mods5",
              "d_BCP_mods5",
              "C_a_mods5",
              "W_c_mods5")

species_labels <- c("A. pseudoplatanus",
                    "B. pendula",
                    "F. sylvatica",
                    "F. excelsior",
                    "Q. robur",
                    "S. aucuparia",
                    "P. abies",
                    "P. sitchensis",
                    "P. sylvestris")

# Plot with separate panel for each effect, showing all species per effect
setwd("C:/temp/pests_analysis/figures_new")
tiff("drivers_mesh5.tiff", width = 165, height = 160, units = 'mm', res = 600)
#par(pr)

layout(matrix(c(1:9), nrow = 3, ncol = 3, byrow = TRUE),
       widths = c(30,19,19))

for(m in 1:length(effect_list)){
  # Create matrix of effects for each species
  mat <- matrix(NA, nrow = length(species), ncol = 3)
  
  for(sp in 1:length(species)){
    mod <- get(mod_list[m], envir = .GlobalEnv)
    mod <- mod[[species[sp]]]
    
    mat[sp,1] <- mod[[length(mod)]][effect_list[m], "mean"]
    mat[sp,2] <- mod[[length(mod)]][effect_list[m], "0.025quant"]
    mat[sp,3] <- mod[[length(mod)]][effect_list[m], "0.975quant"]
  }
  
  if(m %in% c(1,4,7)){
    par(mar = c(4.1, 8, 2.1, 0.5))
  }else{
    par(mar = c(4.1, 1, 2.1, 0.5))
  }
  
  # Create empty plot
  plot(NULL, 
       ylim = c(0, nrow(mat)+1), 
       xlim = c(min(unlist(mat))-0.025, 
                max(unlist(mat))+0.025),
       xlab = "", ylab = "",
       yaxt = "n",
       main = effect_labels[m])
  abline(v = 0, lty = 2)
  title(LETTERS[m], adj = 0)
  
  if(m %in% 7:9){
    title(xlab="Estimate", line=2.6, cex.lab=1.2)
  }
  
  if(m %in% c(1,4,7)){
    axis(side = 2, 
         at = seq(nrow(mat),1,-1),
         labels = parse(text = paste0("italic('", species_labels, "')")),
         tick = TRUE,
         las = 2,
         cex.axis = 1)
  }else{
    axis(side = 2,
         at = seq(nrow(mat),1,-1),
         labels = rep("",9),
         tick = TRUE,
         las = 2)
  }
  
  for(i in 1:nrow(mat)){
    points(x = mat[i,1], y = 10-i, pch = 18, col = "black")
    lines(x = c(mat[i,2], mat[i,3]),
          y = rep(10-i, 2),
          lwd = 1)
  }
}
dev.off()
