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
      dat <- intModel(sp_data,
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


# Load data ####
setwd("C:/temp/pests_analysis/INLA_inputs")

pest_data_all <- get(load("pest_data_all_ISDM_input_list.Rdata"))
pest_data <- get(load("pest_data_ISDM_input_list.Rdata"))
meshGB1 <- get(load("meshGB1.Rdata"))
meshGB2 <- get(load("meshGB2.Rdata"))
meshGB3 <- get(load("meshGB3.Rdata"))
meshGB4 <- get(load("meshGB4.Rdata"))
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

# 0-1 standardise
#min_max_vals <- minmax(covariate_stack)
#covariate_stack_sd <- (covariate_stack - min_max_vals[1,]) / (min_max_vals[2,] - min_max_vals[1,])  

# Create grid to predict over ####
# dm <- c(53,97)
dm <- c(106,194)

pred.df1 <- fm_pixels(meshGB1, 
                      dims = dm,
                      format = "sf",
                      mask = GB_bound)

pred.df2 <- fm_pixels(meshGB2, 
                      dims = dm,
                      format = "sf",
                      mask = GB_bound)

pred.df3 <- fm_pixels(meshGB3, 
                      dims = dm,
                      format = "sf",
                      mask = GB_bound)

pred.df4 <- fm_pixels(meshGB4, 
                      dims = dm,
                      format = "sf",
                      mask = GB_bound)

pred.df5 <- fm_pixels(meshGB5, 
                      dims = dm,
                      format = "sf",
                      mask = GB_bound)

# Load WAIC tables ####
setwd("C:/temp/pests_analysis/WAIC_tables")
file_list <- list.files()
WAIC_list <- list()
for(i in 1:length(file_list)){
  WAIC_list[[i]] <- get(load(file_list[i]))
}
names_list <- gsub("WAIC_table_","",file_list)
names_list <- gsub(".Rdata","",names_list)
names(WAIC_list) <- names_list

# Fit models ####
#bru_options_set(bru_options_default())

# All newleaf hosts
best <- WAIC_list$newleaf %>% filter(delta_WAIC == 0)
cov_sets <- list(strsplit(best$covariates, ", ")[[1]])

# Mesh 1
setwd("F:/INLA_outputs/mesh1")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data_all_list,
                      species = "all",
                      mesh = meshGB1,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df1,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 2
setwd("F:/INLA_outputs/mesh2")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data_all_list,
                      species = "all",
                      mesh = meshGB2,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df2,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 3
setwd("F:/INLA_outputs/mesh3")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data_all_list,
                      species = "all",
                      mesh = meshGB3,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df3,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60


# Acer pseudoplatanus ####
best <- WAIC_list$Acer_pseudoplatanus %>% filter(delta_WAIC == 0)
cov_sets <- list(strsplit(best$covariates, ", ")[[1]])

# Mesh 1
setwd("F:/INLA_outputs/mesh1")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Acer pseudoplatanus",
                      mesh = meshGB1,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df1,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 2
setwd("F:/INLA_outputs/mesh2")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Acer pseudoplatanus",
                      mesh = meshGB2,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df2,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 3
setwd("F:/INLA_outputs/mesh3")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Acer pseudoplatanus",
                      mesh = meshGB3,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df3,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 4
setwd("F:/INLA_outputs/mesh4")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Acer pseudoplatanus",
                      mesh = meshGB4,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df4,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 5
setwd("F:/INLA_outputs/mesh5")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Acer pseudoplatanus",
                      mesh = meshGB5,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df5,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Betula pendula ####
best <- WAIC_list$Betula_pendula %>% filter(delta_WAIC == 0)
cov_sets <- list(strsplit(best$covariates, ", ")[[1]])

# Mesh 1
setwd("F:/INLA_outputs/mesh1")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Betula pendula",
                      mesh = meshGB1,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df1,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 2
setwd("F:/INLA_outputs/mesh2")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Betula pendula",
                      mesh = meshGB2,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df2,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 3
setwd("F:/INLA_outputs/mesh3")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Betula pendula",
                      mesh = meshGB3,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df3,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 4
setwd("F:/INLA_outputs/mesh4")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Betula pendula",
                      mesh = meshGB4,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df4,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 5
setwd("F:/INLA_outputs/mesh5")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Betula pendula",
                      mesh = meshGB5,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df5,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Fagus sylvatica ####
best <- WAIC_list$Fagus_sylvatica %>% filter(delta_WAIC == 0)
cov_sets <- list(strsplit(best$covariates, ", ")[[1]])

# Mesh 1
setwd("F:/INLA_outputs/mesh1")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Fagus sylvatica",
                      mesh = meshGB1,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df1,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 2
setwd("F:/INLA_outputs/mesh2")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Fagus sylvatica",
                      mesh = meshGB2,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df2,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 3
setwd("F:/INLA_outputs/mesh3")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Fagus sylvatica",
                      mesh = meshGB3,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df3,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 4
setwd("F:/INLA_outputs/mesh4")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Fagus sylvatica",
                      mesh = meshGB4,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df4,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 5
setwd("F:/INLA_outputs/mesh5")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Fagus sylvatica",
                      mesh = meshGB5,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df5,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Fraxinus excelsior ####
best <- WAIC_list$Fraxinus_excelsior %>% filter(delta_WAIC == 0)
cov_sets <- list(strsplit(best$covariates, ", ")[[1]])

# Mesh 1
setwd("F:/INLA_outputs/mesh1")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Fraxinus excelsior",
                      mesh = meshGB1,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df1,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 2
setwd("F:/INLA_outputs/mesh2")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Fraxinus excelsior",
                      mesh = meshGB2,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df2,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 3
setwd("F:/INLA_outputs/mesh3")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Fraxinus excelsior",
                      mesh = meshGB3,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df3,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 4
setwd("F:/INLA_outputs/mesh4")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Fraxinus excelsior",
                      mesh = meshGB4,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df4,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 5
setwd("F:/INLA_outputs/mesh5")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Fraxinus excelsior",
                      mesh = meshGB5,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df5,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Picea abies ####
best <- WAIC_list$Picea_abies %>% filter(delta_WAIC == 0)
cov_sets <- list(strsplit(best$covariates, ", ")[[1]])

# Mesh 1
setwd("F:/INLA_outputs/mesh1")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Picea abies",
                      mesh = meshGB1,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df1,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 2
setwd("F:/INLA_outputs/mesh2")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Picea abies",
                      mesh = meshGB2,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df2,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 3
setwd("F:/INLA_outputs/mesh3")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Picea abies",
                      mesh = meshGB3,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df3,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 4
setwd("F:/INLA_outputs/mesh4")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Picea abies",
                      mesh = meshGB4,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df4,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 5
setwd("F:/INLA_outputs/mesh5")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Picea abies",
                      mesh = meshGB5,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df5,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Picea sitchensis ####
best <- WAIC_list$Picea_sitchensis %>% filter(delta_WAIC == 0)
cov_sets <- list(strsplit(best$covariates, ", ")[[1]])

# Mesh 1
setwd("F:/INLA_outputs/mesh1")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Picea sitchensis",
                      mesh = meshGB1,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df1,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 2
setwd("F:/INLA_outputs/mesh2")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Picea sitchensis",
                      mesh = meshGB2,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df2,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 3
setwd("F:/INLA_outputs/mesh3")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Picea sitchensis",
                      mesh = meshGB3,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df3,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 4
setwd("F:/INLA_outputs/mesh4")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Picea sitchensis",
                      mesh = meshGB4,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df4,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 5
setwd("F:/INLA_outputs/mesh5")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Picea sitchensis",
                      mesh = meshGB5,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df5,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Pinus sylvestris ####
best <- WAIC_list$Pinus_sylvestris %>% filter(delta_WAIC == 0)
cov_sets <- list(strsplit(best$covariates, ", ")[[1]])

# Mesh 1
setwd("F:/INLA_outputs/mesh1")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Pinus sylvestris",
                      mesh = meshGB1,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df1,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 2
setwd("F:/INLA_outputs/mesh2")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Pinus sylvestris",
                      mesh = meshGB2,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df2,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 3
setwd("F:/INLA_outputs/mesh3")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Pinus sylvestris",
                      mesh = meshGB3,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df3,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 4
setwd("F:/INLA_outputs/mesh4")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Pinus sylvestris",
                      mesh = meshGB4,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df4,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 5
setwd("F:/INLA_outputs/mesh5")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Pinus sylvestris",
                      mesh = meshGB5,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df5,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Quercus robur ####
best <- WAIC_list$Quercus_robur %>% filter(delta_WAIC == 0)
cov_sets <- list(strsplit(best$covariates, ", ")[[1]])

# Mesh 1
setwd("F:/INLA_outputs/mesh1")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Quercus robur",
                      mesh = meshGB1,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df1,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 2
setwd("F:/INLA_outputs/mesh2")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Quercus robur",
                      mesh = meshGB2,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df2,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 3
setwd("F:/INLA_outputs/mesh3")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Quercus robur",
                      mesh = meshGB3,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df3,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 4
setwd("F:/INLA_outputs/mesh4")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Quercus robur",
                      mesh = meshGB4,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df4,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 5
setwd("F:/INLA_outputs/mesh5")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Quercus robur",
                      mesh = meshGB5,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df5,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60


# Sorbus aucuparia ####
best <- WAIC_list$Sorbus_aucuparia %>% filter(delta_WAIC == 0)
cov_sets <- list(strsplit(best$covariates, ", ")[[1]])


# Mesh 1
setwd("F:/INLA_outputs/mesh1")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Sorbus aucuparia",
                      mesh = meshGB1,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df1,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 2
setwd("F:/INLA_outputs/mesh2")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Sorbus aucuparia",
                      mesh = meshGB2,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df2,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 3
setwd("F:/INLA_outputs/mesh3")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Sorbus aucuparia",
                      mesh = meshGB3,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df3,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 4
setwd("F:/INLA_outputs/mesh4")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Sorbus aucuparia",
                      mesh = meshGB4,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df4,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Mesh 5
setwd("F:/INLA_outputs/mesh5")
ptm<-proc.time()
mod <- run_multi_ISDM(df_list = pest_data,
                      species = "Sorbus aucuparia",
                      mesh = meshGB5,
                      covs = covariate_stack_sd,
                      covs_select = cov_sets, 
                      priors_pc = c(1, 0.05, 10, 0.05),
                      predict_data = pred.df5,
                      save_predictions = TRUE)
mod.time<-proc.time()-ptm
mod.time[3]/60

# Plot maps ####
# Mesh 1 ####
# Load saved model predictions
species <- c("Acer pseudoplatanus",
             "Betula pendula",
             "Fagus sylvatica",
             "Fraxinus excelsior",
             "Picea abies",
             "Picea sitchensis",
             "Pinus sylvestris",
             "Quercus robur",
             "Sorbus aucuparia")

setwd("D:/INLA_outputs/mesh1")
file_list <- list.files(pattern = "predictions.Rdata")
predictions_list <- list()
for(i in 1:length(file_list)){
  predictions_list[[i]] <- get(load(file_list[i]))
}
names(predictions_list) <- species

#mylegend <- get(load("C:/temp/pests_analysis/mylegend.Rdata"))

# Plot prediction mean
p_list <- list()
for(sp in 1:length(species)){
  predictions_df_sf <- predictions_list[[species[sp]]]
  
  pred.mask <- predictions_df_sf
  
  pred.mask$mean <- (pred.mask$mean - min(pred.mask$mean)) / (max(pred.mask$mean) - min(pred.mask$mean))
  
  #pred.mask <- st_filter(predictions_df_sf, gb_all)
  pred.mask <- pred.mask %>% select(geometry, mean)
  #pred.mask <- pred.mask %>% select(geometry, sd)
  #pred.mask <- pred.mask %>% select(geometry, q0.025)
  #pred.mask <- pred.mask %>% select(geometry, q0.975)
  
  # Remove islands as islands were not modeled
  pred.mask.clip <- st_filter(pred.mask, GB_bound)
  
  p_temp <- ggplot() +
    gg(pred.mask.clip, geom = "tile") +
    geom_sf(data = GB_bound, fill = NA, lwd = 0.8, col = "black") +
    #geom_sf(data = df_sf) +
    #ggtitle(expression(paste(italic("Fraxinus excelsior")))) +
    #ggtitle("Predictions", subtitle = "(Linear predictor scale)") +
    #ggtitle(paste(species[sp])) +
    #ggtitle(paste(species[sp]), 
    #        subtitle = paste0("PB = ", nrow(pest_data[[species[sp]]]$THDAS_obs),", ",
    #                          "PA = ", nrow(pest_data[[species[sp]]]$APHA_obs))) +
    #ggtitle("All Newleaf hosts",
    #        subtitle = "PB = 1724, PA = 1995") +
    xlab("") + ylab("") +
    theme(legend.position = "none") +
    #labs(fill = "Mean") +
    #theme(legend.key.size = unit(1, "cm")) +
    scale_fill_viridis(option = "inferno", direction = 1)
  
  p_list[[sp]] <- p_temp
  print(sp)
}
names(p_list) <- species

p_list[[10]] = mylegend

grid1 <- plot_grid(plotlist = p_list[c(1,2,3,4,8,9,5,6,7,10)],
                   nrow = 2,
                   labels = c("A","B","C","D","E","F","G","H","I"),
                   hjust = -1)

setwd("C:/temp/pests_analysis/figures_new")
ggsave2(filename = "predictions_mesh1.tiff",
        plot = grid1,
        bg = "white",
        dpi = 600,
        width = 165,
        height = 110,
        units = "mm",
        scale = 1.7)


# Plot prediction uncertainty
p_list <- list()
for(sp in 1:length(species)){
  predictions_df_sf <- predictions_list[[species[sp]]]
  
  # Adjust alpha values to communicate prediction uncertainty
  #predictions_df_sf$alpha <- 1 - (abs(predictions_df_sf$sd)/abs(predictions_df_sf$mean))
  #predictions_df_sf$alpha[predictions_df_sf$alpha < 0] <- 0
  
  pred.mask <- predictions_df_sf
  
  pred.mask$sd <- (pred.mask$sd - min(pred.mask$sd)) / (max(pred.mask$mean) - min(pred.mask$mean))
  
  #pred.mask <- st_filter(predictions_df_sf, gb_all)
  #pred.mask <- pred.mask %>% select(geometry, alpha)
  pred.mask <- pred.mask %>% select(geometry, sd)
  #pred.mask <- pred.mask %>% select(geometry, q0.025)
  #pred.mask <- pred.mask %>% select(geometry, q0.975)
  
  # Remove islands as islands were not modeled
  pred.mask.clip <- st_filter(pred.mask, GB_bound)
  
  p_temp <- ggplot() +
    gg(pred.mask.clip, geom = "tile") +
    geom_sf(data = GB_bound, fill = NA, lwd = 0.8, col = "black") +
    #geom_sf(data = df_sf) +
    #ggtitle(expression(paste(italic("Fraxinus excelsior")))) +
    #ggtitle("Predictions", subtitle = "(Linear predictor scale)") +
    #ggtitle(paste(species[sp])) +
    #ggtitle(paste(species[sp]), 
    #        subtitle = paste0("PB = ", nrow(pest_data[[species[sp]]]$THDAS_obs),", ",
    #                          "PA = ", nrow(pest_data[[species[sp]]]$APHA_obs))) +
    #ggtitle("All Newleaf hosts",
    #        subtitle = "PB = 1724, PA = 1995") +
    xlab("") + ylab("") +
    theme(legend.position = "none") +
    #labs(fill = "SD") +
    #theme(legend.key.size = unit(1, "cm")) +
    
    scale_fill_viridis(option = "viridis", direction = 1)
  
  p_list[[sp]] <- p_temp
  print(sp)
}
names(p_list) <- species

p_list[[10]] = mylegend2

grid2 <- plot_grid(plotlist = p_list[c(1,2,3,4,8,9,5,6,7,10)],
                   nrow = 2,
                   labels = c("A","B","C","D","E","F","G","H","I"),
                   hjust = -1)

setwd("C:/temp/pests_analysis/figures_new")
ggsave2(filename = "uncertainty_mesh1.tiff",
        plot = grid2,
        bg = "white",
        dpi = 600,
        width = 165,
        height = 110,
        units = "mm",
        scale = 1.7)


# Mesh 2 ####
# Load saved model predictions
species <- c("Acer pseudoplatanus",
             "Betula pendula",
             "Fagus sylvatica",
             "Fraxinus excelsior",
             "Picea abies",
             "Picea sitchensis",
             "Pinus sylvestris",
             "Quercus robur",
             "Sorbus aucuparia")

setwd("D:/INLA_outputs/mesh2")
file_list <- list.files(pattern = "predictions.Rdata")
predictions_list <- list()
for(i in 1:length(file_list)){
  predictions_list[[i]] <- get(load(file_list[i]))
}
names(predictions_list) <- species

# Plot prediction mean
p_list <- list()
for(sp in 1:length(species)){
  predictions_df_sf <- predictions_list[[species[sp]]]
  
  pred.mask <- predictions_df_sf
  
  pred.mask$mean <- (pred.mask$mean - min(pred.mask$mean)) / (max(pred.mask$mean) - min(pred.mask$mean))
  
  #pred.mask <- st_filter(predictions_df_sf, gb_all)
  pred.mask <- pred.mask %>% select(geometry, mean)
  #pred.mask <- pred.mask %>% select(geometry, sd)
  #pred.mask <- pred.mask %>% select(geometry, q0.025)
  #pred.mask <- pred.mask %>% select(geometry, q0.975)
  
  # Remove islands as islands were not modeled
  pred.mask.clip <- st_filter(pred.mask, GB_bound)
  
  p_temp <- ggplot() +
    gg(pred.mask.clip, geom = "tile") +
    geom_sf(data = GB_bound, fill = NA, lwd = 0.8, col = "black") +
    #geom_sf(data = df_sf) +
    #ggtitle(expression(paste(italic("Fraxinus excelsior")))) +
    #ggtitle("Predictions", subtitle = "(Linear predictor scale)") +
    #ggtitle(paste(species[sp])) +
    #ggtitle(paste(species[sp]), 
    #        subtitle = paste0("PB = ", nrow(pest_data[[species[sp]]]$THDAS_obs),", ",
    #                          "PA = ", nrow(pest_data[[species[sp]]]$APHA_obs))) +
    #ggtitle("All Newleaf hosts",
    #        subtitle = "PB = 1724, PA = 1995") +
    xlab("") + ylab("") +
    theme(legend.position = "none") +
    #theme(legend.key.size = unit(0.7, "cm")) +
    scale_fill_viridis(option = "inferno", direction = 1)
  
  p_list[[sp]] <- p_temp
  print(sp)
}
names(p_list) <- species

p_list[[10]] = mylegend

grid1 <- plot_grid(plotlist = p_list[c(1,2,3,4,8,9,5,6,7,10)],
                   nrow = 2,
                   labels = c("A","B","C","D","E","F","G","H","I"),
                   hjust = -1)

setwd("C:/temp/pests_analysis/figures_new")
ggsave2(filename = "predictions_mesh2.tiff",
        plot = grid1,
        bg = "white",
        dpi = 600,
        width = 165,
        height = 110,
        units = "mm",
        scale = 1.7)

# Plot prediction uncertainty
p_list <- list()
for(sp in 1:length(species)){
  predictions_df_sf <- predictions_list[[species[sp]]]
  
  # Adjust alpha values to communicate prediction uncertainty
  #predictions_df_sf$alpha <- 1 - (abs(predictions_df_sf$sd)/abs(predictions_df_sf$mean))
  #predictions_df_sf$alpha[predictions_df_sf$alpha < 0] <- 0
  
  pred.mask <- predictions_df_sf
  
  pred.mask$sd <- (pred.mask$sd - min(pred.mask$sd)) / (max(pred.mask$mean) - min(pred.mask$mean))
  
  #pred.mask <- st_filter(predictions_df_sf, gb_all)
  #pred.mask <- pred.mask %>% select(geometry, alpha)
  pred.mask <- pred.mask %>% select(geometry, sd)
  #pred.mask <- pred.mask %>% select(geometry, q0.025)
  #pred.mask <- pred.mask %>% select(geometry, q0.975)
  
  # Remove islands as islands were not modeled
  pred.mask.clip <- st_filter(pred.mask, GB_bound)
  
  p_temp <- ggplot() +
    gg(pred.mask.clip, geom = "tile") +
    geom_sf(data = GB_bound, fill = NA, lwd = 0.8, col = "black") +
    #geom_sf(data = df_sf) +
    #ggtitle(expression(paste(italic("Fraxinus excelsior")))) +
    #ggtitle("Predictions", subtitle = "(Linear predictor scale)") +
    #ggtitle(paste(species[sp])) +
    #ggtitle(paste(species[sp]), 
    #        subtitle = paste0("PB = ", nrow(pest_data[[species[sp]]]$THDAS_obs),", ",
    #                          "PA = ", nrow(pest_data[[species[sp]]]$APHA_obs))) +
    #ggtitle("All Newleaf hosts",
    #        subtitle = "PB = 1724, PA = 1995") +
    xlab("") + ylab("") +
    theme(legend.position = "none") +
    #theme(legend.key.size = unit(0.7, "cm")) +
    scale_fill_viridis(option = "viridis", direction = 1)
  
  p_list[[sp]] <- p_temp
  print(sp)
}
names(p_list) <- species

p_list[[10]] = mylegend2

grid2 <- plot_grid(plotlist = p_list[c(1,2,3,4,8,9,5,6,7,10)],
                   nrow = 2,
                   labels = c("A","B","C","D","E","F","G","H","I"),
                   hjust = -1)

setwd("C:/temp/pests_analysis/figures_new")
ggsave2(filename = "uncertainty_mesh2.tiff",
        plot = grid2,
        bg = "white",
        dpi = 600,
        width = 165,
        height = 110,
        units = "mm",
        scale = 1.7)

# Mesh 3 ####
# Load saved model predictions
species <- c("Acer pseudoplatanus",
             "Betula pendula",
             "Fagus sylvatica",
             "Fraxinus excelsior",
             "Picea abies",
             "Picea sitchensis",
             "Pinus sylvestris",
             "Quercus robur",
             "Sorbus aucuparia")

setwd("D:/INLA_outputs/mesh3")
file_list <- list.files(pattern = "predictions.Rdata")
predictions_list <- list()
for(i in 1:length(file_list)){
  predictions_list[[i]] <- get(load(file_list[i]))
}
names(predictions_list) <- species

# Plot prediction mean
p_list <- list()
for(sp in 1:length(species)){
  predictions_df_sf <- predictions_list[[species[sp]]]
  
  pred.mask <- predictions_df_sf
  
  pred.mask$mean <- (pred.mask$mean - min(pred.mask$mean)) / (max(pred.mask$mean) - min(pred.mask$mean))
  
  #pred.mask <- st_filter(predictions_df_sf, gb_all)
  pred.mask <- pred.mask %>% select(geometry, mean)
  #pred.mask <- pred.mask %>% select(geometry, sd)
  #pred.mask <- pred.mask %>% select(geometry, q0.025)
  #pred.mask <- pred.mask %>% select(geometry, q0.975)
  
  # Remove islands as islands were not modeled
  pred.mask.clip <- st_filter(pred.mask, GB_bound)
  
  p_temp <- ggplot() +
    gg(pred.mask.clip, geom = "tile") +
    geom_sf(data = GB_bound, fill = NA, lwd = 0.8, col = "black") +
    #geom_sf(data = df_sf) +
    #ggtitle(expression(paste(italic("Fraxinus excelsior")))) +
    #ggtitle("Predictions", subtitle = "(Linear predictor scale)") +
    #ggtitle(paste(species[sp])) +
    #ggtitle(paste(species[sp]), 
    #        subtitle = paste0("PB = ", nrow(pest_data[[species[sp]]]$THDAS_obs),", ",
    #                          "PA = ", nrow(pest_data[[species[sp]]]$APHA_obs))) +
    #ggtitle("All Newleaf hosts",
    #        subtitle = "PB = 1724, PA = 1995") +
    xlab("") + ylab("") +
    theme(legend.position = "none") +
    #theme(legend.key.size = unit(0.7, "cm")) +
    scale_fill_viridis(option = "inferno", direction = 1)
  
  p_list[[sp]] <- p_temp
  print(sp)
}
names(p_list) <- species

p_list[[10]] = mylegend

grid1 <- plot_grid(plotlist = p_list[c(1,2,3,4,8,9,5,6,7,10)],
                   nrow = 2,
                   labels = c("A","B","C","D","E","F","G","H","I"),
                   hjust = -1)

setwd("C:/temp/pests_analysis/figures_new")
ggsave2(filename = "predictions_mesh3.tiff",
        plot = grid1,
        bg = "white",
        dpi = 600,
        width = 165,
        height = 110,
        units = "mm",
        scale = 1.7)


# Plot prediction uncertainty
p_list <- list()
for(sp in 1:length(species)){
  predictions_df_sf <- predictions_list[[species[sp]]]
  
  # Adjust alpha values to communicate prediction uncertainty
  #predictions_df_sf$alpha <- 1 - (abs(predictions_df_sf$sd)/abs(predictions_df_sf$mean))
  #predictions_df_sf$alpha[predictions_df_sf$alpha < 0] <- 0
  
  pred.mask <- predictions_df_sf
  
  pred.mask$sd <- (pred.mask$sd - min(pred.mask$sd)) / (max(pred.mask$mean) - min(pred.mask$mean))
  
  #pred.mask <- st_filter(predictions_df_sf, gb_all)
  #pred.mask <- pred.mask %>% select(geometry, alpha)
  pred.mask <- pred.mask %>% select(geometry, sd)
  #pred.mask <- pred.mask %>% select(geometry, q0.025)
  #pred.mask <- pred.mask %>% select(geometry, q0.975)
  
  # Remove islands as islands were not modeled
  pred.mask.clip <- st_filter(pred.mask, GB_bound)
  
  p_temp <- ggplot() +
    gg(pred.mask.clip, geom = "tile") +
    geom_sf(data = GB_bound, fill = NA, lwd = 0.8, col = "black") +
    #geom_sf(data = df_sf) +
    #ggtitle(expression(paste(italic("Fraxinus excelsior")))) +
    #ggtitle("Predictions", subtitle = "(Linear predictor scale)") +
    #ggtitle(paste(species[sp])) +
    #ggtitle(paste(species[sp]), 
    #        subtitle = paste0("PB = ", nrow(pest_data[[species[sp]]]$THDAS_obs),", ",
    #                          "PA = ", nrow(pest_data[[species[sp]]]$APHA_obs))) +
    #ggtitle("All Newleaf hosts",
    #        subtitle = "PB = 1724, PA = 1995") +
    xlab("") + ylab("") +
    theme(legend.position = "none") +
    #theme(legend.key.size = unit(0.7, "cm")) +
    scale_fill_viridis(option = "viridis", direction = 1)
  
  p_list[[sp]] <- p_temp
  print(sp)
}
names(p_list) <- species

p_list[[10]] = mylegend2

grid2 <- plot_grid(plotlist = p_list[c(1,2,3,4,8,9,5,6,7,10)],
                   nrow = 2,
                   labels = c("A","B","C","D","E","F","G","H","I"),
                   hjust = -1)

setwd("C:/temp/pests_analysis/figures_new")
ggsave2(filename = "uncertainty_mesh3.tiff",
        plot = grid2,
        bg = "white",
        dpi = 600,
        width = 165,
        height = 110,
        units = "mm",
        scale = 1.7)

# Mesh 4 ####
# Load saved model predictions
species <- c("Acer pseudoplatanus",
             "Betula pendula",
             "Fagus sylvatica",
             "Fraxinus excelsior",
             "Picea abies",
             "Picea sitchensis",
             "Pinus sylvestris",
             "Quercus robur",
             "Sorbus aucuparia")

setwd("D:/INLA_outputs/mesh4")
file_list <- list.files(pattern = "predictions.Rdata")
predictions_list <- list()
for(i in 1:length(file_list)){
  predictions_list[[i]] <- get(load(file_list[i]))
}
names(predictions_list) <- species

# Plot prediction mean
p_list <- list()
for(sp in 1:length(species)){
  predictions_df_sf <- predictions_list[[species[sp]]]
  
  pred.mask <- predictions_df_sf
  
  pred.mask$mean <- (pred.mask$mean - min(pred.mask$mean)) / (max(pred.mask$mean) - min(pred.mask$mean))
  
  #pred.mask <- st_filter(predictions_df_sf, gb_all)
  pred.mask <- pred.mask %>% select(geometry, mean)
  #pred.mask <- pred.mask %>% select(geometry, sd)
  #pred.mask <- pred.mask %>% select(geometry, q0.025)
  #pred.mask <- pred.mask %>% select(geometry, q0.975)
  
  # Remove islands as islands were not modeled
  pred.mask.clip <- st_filter(pred.mask, GB_bound)
  
  p_temp <- ggplot() +
    gg(pred.mask.clip, geom = "tile") +
    geom_sf(data = GB_bound, fill = NA, lwd = 0.8, col = "black") +
    #geom_sf(data = df_sf) +
    #ggtitle(expression(paste(italic("Fraxinus excelsior")))) +
    #ggtitle("Predictions", subtitle = "(Linear predictor scale)") +
    #ggtitle(paste(species[sp])) +
    #ggtitle(paste(species[sp]), 
    #        subtitle = paste0("PB = ", nrow(pest_data[[species[sp]]]$THDAS_obs),", ",
    #                          "PA = ", nrow(pest_data[[species[sp]]]$APHA_obs))) +
    #ggtitle("All Newleaf hosts",
    #        subtitle = "PB = 1724, PA = 1995") +
    xlab("") + ylab("") +
    theme(legend.position = "none") +
    #theme(legend.key.size = unit(0.7, "cm")) +
    scale_fill_viridis(option = "inferno", direction = 1)
  
  p_list[[sp]] <- p_temp
  print(sp)
}
names(p_list) <- species

p_list[[10]] = mylegend

grid1 <- plot_grid(plotlist = p_list[c(1,2,3,4,8,9,5,6,7,10)],
                   nrow = 2,
                   labels = c("A","B","C","D","E","F","G","H","I"),
                   hjust = -1)

setwd("C:/temp/pests_analysis/figures_new")
ggsave2(filename = "predictions_mesh4.tiff",
        plot = grid1,
        bg = "white",
        dpi = 600,
        width = 165,
        height = 110,
        units = "mm",
        scale = 1.7)


# Plot prediction uncertainty
p_list <- list()
for(sp in 1:length(species)){
  predictions_df_sf <- predictions_list[[species[sp]]]
  
  # Adjust alpha values to communicate prediction uncertainty
  #predictions_df_sf$alpha <- 1 - (abs(predictions_df_sf$sd)/abs(predictions_df_sf$mean))
  #predictions_df_sf$alpha[predictions_df_sf$alpha < 0] <- 0
  
  pred.mask <- predictions_df_sf
  
  pred.mask$sd <- (pred.mask$sd - min(pred.mask$sd)) / (max(pred.mask$mean) - min(pred.mask$mean))
  
  #pred.mask <- st_filter(predictions_df_sf, gb_all)
  #pred.mask <- pred.mask %>% select(geometry, alpha)
  pred.mask <- pred.mask %>% select(geometry, sd)
  #pred.mask <- pred.mask %>% select(geometry, q0.025)
  #pred.mask <- pred.mask %>% select(geometry, q0.975)
  
  # Remove islands as islands were not modeled
  pred.mask.clip <- st_filter(pred.mask, GB_bound)
  
  p_temp <- ggplot() +
    gg(pred.mask.clip, geom = "tile") +
    geom_sf(data = GB_bound, fill = NA, lwd = 0.8, col = "black") +
    #geom_sf(data = df_sf) +
    #ggtitle(expression(paste(italic("Fraxinus excelsior")))) +
    #ggtitle("Predictions", subtitle = "(Linear predictor scale)") +
    #ggtitle(paste(species[sp])) +
    #ggtitle(paste(species[sp]), 
    #        subtitle = paste0("PB = ", nrow(pest_data[[species[sp]]]$THDAS_obs),", ",
    #                          "PA = ", nrow(pest_data[[species[sp]]]$APHA_obs))) +
    #ggtitle("All Newleaf hosts",
    #        subtitle = "PB = 1724, PA = 1995") +
    xlab("") + ylab("") +
    theme(legend.position = "none") +
    #theme(legend.key.size = unit(0.7, "cm")) +
    scale_fill_viridis(option = "viridis", direction = 1)
  
  p_list[[sp]] <- p_temp
  print(sp)
}
names(p_list) <- species

p_list[[10]] = mylegend2

grid2 <- plot_grid(plotlist = p_list[c(1,2,3,4,8,9,5,6,7,10)],
                   nrow = 2,
                   labels = c("A","B","C","D","E","F","G","H","I"),
                   hjust = -1)

setwd("C:/temp/pests_analysis/figures_new")
ggsave2(filename = "uncertainty_mesh4.tiff",
        plot = grid2,
        bg = "white",
        dpi = 600,
        width = 165,
        height = 110,
        units = "mm",
        scale = 1.7)

# Mesh 5 ####
# Load saved model predictions
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
names(predictions_list) <- species

# Plot prediction mean
p_list <- list()
for(sp in 1:length(species)){
  predictions_df_sf <- predictions_list[[species[sp]]]
  
  pred.mask <- predictions_df_sf
  
  pred.mask$mean <- (pred.mask$mean - min(pred.mask$mean)) / (max(pred.mask$mean) - min(pred.mask$mean))
  
  
  #pred.mask <- st_filter(predictions_df_sf, gb_all)
  pred.mask <- pred.mask %>% select(geometry, mean)
  #pred.mask <- pred.mask %>% select(geometry, sd)
  #pred.mask <- pred.mask %>% select(geometry, q0.025)
  #pred.mask <- pred.mask %>% select(geometry, q0.975)
  
  # Remove islands as islands were not modeled
  pred.mask.clip <- st_filter(pred.mask, GB_bound)
  
  p_temp <- ggplot() +
    gg(pred.mask.clip, geom = "tile") +
    geom_sf(data = GB_bound, fill = NA, lwd = 0.8, col = "black") +
    #geom_sf(data = df_sf) +
    #ggtitle(expression(paste(italic("Fraxinus excelsior")))) +
    #ggtitle("Predictions", subtitle = "(Linear predictor scale)") +
    #ggtitle(paste(species[sp])) +
    #ggtitle(paste(species[sp]), 
    #        subtitle = paste0("PB = ", nrow(pest_data[[species[sp]]]$THDAS_obs),", ",
    #                          "PA = ", nrow(pest_data[[species[sp]]]$APHA_obs))) +
    #ggtitle("All Newleaf hosts",
    #        subtitle = "PB = 1724, PA = 1995") +
    xlab("") + ylab("") +
    theme(legend.position = "none") +
    #theme(legend.key.size = unit(0.7, "cm")) +
    scale_fill_viridis(option = "inferno", direction = 1)
  
  p_list[[sp]] <- p_temp
  print(sp)
}
names(p_list) <- species

p_list[[10]] = mylegend

grid1 <- plot_grid(plotlist = p_list[c(1,2,3,4,8,9,5,6,7,10)],
                   nrow = 2,
                   labels = c("A","B","C","D","E","F","G","H","I"),
                   hjust = -1)

setwd("C:/temp/pests_analysis/figures_new")
ggsave2(filename = "predictions_mesh5.tiff",
        plot = grid1,
        bg = "white",
        dpi = 600,
        width = 165,
        height = 110,
        units = "mm",
        scale = 1.7)


# Plot prediction uncertainty
p_list <- list()
for(sp in 1:length(species)){
  predictions_df_sf <- predictions_list[[species[sp]]]
  
  # Adjust alpha values to communicate prediction uncertainty
  #predictions_df_sf$alpha <- 1 - (abs(predictions_df_sf$sd)/abs(predictions_df_sf$mean))
  #predictions_df_sf$alpha[predictions_df_sf$alpha < 0] <- 0
  
  pred.mask <- predictions_df_sf
  
  pred.mask$sd <- (pred.mask$sd - min(pred.mask$sd)) / (max(pred.mask$mean) - min(pred.mask$mean))
  
  #pred.mask <- st_filter(predictions_df_sf, gb_all)
  #pred.mask <- pred.mask %>% select(geometry, alpha)
  pred.mask <- pred.mask %>% select(geometry, sd)
  #pred.mask <- pred.mask %>% select(geometry, q0.025)
  #pred.mask <- pred.mask %>% select(geometry, q0.975)
  
  # Remove islands as islands were not modeled
  pred.mask.clip <- st_filter(pred.mask, GB_bound)
  
  p_temp <- ggplot() +
    gg(pred.mask.clip, geom = "tile") +
    geom_sf(data = GB_bound, fill = NA, lwd = 0.8, col = "black") +
    #geom_sf(data = df_sf) +
    #ggtitle(expression(paste(italic("Fraxinus excelsior")))) +
    #ggtitle("Predictions", subtitle = "(Linear predictor scale)") +
    #ggtitle(paste(species[sp])) +
    #ggtitle(paste(species[sp]), 
    #        subtitle = paste0("PB = ", nrow(pest_data[[species[sp]]]$THDAS_obs),", ",
    #                          "PA = ", nrow(pest_data[[species[sp]]]$APHA_obs))) +
    #ggtitle("All Newleaf hosts",
    #        subtitle = "PB = 1724, PA = 1995") +
    xlab("") + ylab("") +
    theme(legend.position = "none") +
    #theme(legend.key.size = unit(0.7, "cm")) +
    scale_fill_viridis(option = "viridis", direction = 1)
  
  p_list[[sp]] <- p_temp
  print(sp)
}
names(p_list) <- species

p_list[[10]] = mylegend2

grid2 <- plot_grid(plotlist = p_list[c(1,2,3,4,8,9,5,6,7,10)],
                   nrow = 2,
                   labels = c("A","B","C","D","E","F","G","H","I"),
                   hjust = -1)

setwd("C:/temp/pests_analysis/figures_new")
ggsave2(filename = "uncertainty_mesh5.tiff",
        plot = grid2,
        bg = "white",
        dpi = 600,
        width = 165,
        height = 110,
        units = "mm",
        scale = 1.7)
