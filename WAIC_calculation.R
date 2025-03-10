# Packages ####
library(terra)
library(sf)
library(dplyr)
library(INLA)
library(inlabru)

# Custom functions ####
# Function to evaluate multiple spatial covariate values at observation coordinates
eval_spatial_multi <- function(data_sf, covariates){
  # Separate the geometry column to avoid overwriting it
  gm <- data_sf$geometry
  data_sf <- st_drop_geometry(data_sf)
  
  # Save original variable names for later
  names_original <- colnames(data_sf)
  
  # Loop over covariate layers in the raster
  for(i in 1:nlyr(covariates)){
    # Evaluate covariate values at observation locations
    cov_eval <- eval_spatial(covariates, gm, layer = names(covariates)[i])
    
    # For any missing values perform nearest neighbour imputation
    if (any(is.na(cov_eval))){
      warning(c("Imputed ", sum(is.na(cov_eval)), " missing covariate values for ", names(covariates)[i]))
      cov_eval <- bru_fill_missing(covariates, gm, cov_eval, layer = names(covariates)[i])
    }
    
    data_sf <- cbind(data_sf, cov_eval)
  }
  
  # Update variable names, re-attach geometry, and return dataset
  colnames(data_sf) <- c(names_original, names(covariates))
  data_sf <- st_set_geometry(data_sf, gm)
  return(data_sf)
}

# Function to run model for multiple species
run_multisp_bru <- function(comp, pest_data, covariates){
  # Create list to store results
  results_list <- list()
  
  # Loop over species
  for(sp in 1:length(names(pest_data))){
    
    # Extract APHA data for the species
    df <- pest_data[[sp]]
    df <- df$APHA_obs
    
    # If there are no data then skip to next iteration
    if(is.null(df) == TRUE) next
    
    # Convert to sf
    df_sf <- st_as_sf(df, coords = c("X","Y"))
    df_sf <- st_set_crs(df_sf, 27700) 
    
    # Evaluate covariate values at observation locations 
    df2 <- eval_spatial_multi(data_sf = df_sf, covariates = covariates)
    
    # Fit model
    mod <- bru(comp, family = "binomial", control.family=list(link="logit"), Ntrials = df_sf$Surveyed, data = df2)
    
    # Save estimates of fixed effects to results list
    results_list[[sp]] <- mod$summary.fixed
  }
}  


# Function to run multiple models for one species
run_multiform_bru <- function(comp_list, pest_data, covariates){
  
  df <- pest_data
  
  # Create list to store outputs
  out_list <- list()
  
  # Set up progress bar
  progress_bar = txtProgressBar(min=0, max=length(comp_list), style = 3, char="=")
  
  # Loop over formulae
  for(fo in 1:length(comp_list)){
    # Fit model
    mod <- bru(comp_list[[fo]], family = "binomial", control.family=list(link="logit"), Ntrials = df$Surveyed, data = df)
    
    # Save estimates of fixed effects to results list
    out_list[[fo]] <- mod
    
    # Update progress bar
    setTxtProgressBar(progress_bar, value = fo)
  }
  
  # Return output list
  close(progress_bar) 
  return(out_list)
} 

# Function to convert variable names into inlabru formula
names_to_formula <- function(sets, form_front, form_back){
  # Create empty list to store formulae
  form_list <- list()
  
  for(set in 1:length(sets)){
    # Create empty formula string
    form_str <- as.character(NULL)
    # Add variables in the adjustment set to the formula string
    for(i in 1:length(sets[[set]])){
      if(i ==1){
        form_str <- paste0(form_str, sets[[set]][i],"(", sets[[set]][i],', model = "linear", mean.linear = 0, prec.linear = 1)')
      }else{
        form_str <- paste0(form_str, " + ", sets[[set]][i],"(", sets[[set]][i],', model = "linear", mean.linear = 0, prec.linear = 1)')
      }
    }
    # Add front and back onto formula string
    form_full <- paste(form_front, form_str, form_back)
    form_full <- as.formula(form_full)
    
    # Append to the formula list
    form_list[[set]] <- form_full
  }
  
  # Return the formula list
  return(form_list)
}


# Load data ####
setwd("/data/notebooks/rstudio-peter/outputs/INLA_inputs")

all_pests_PA_all <- get(load("all_pests_PA_all.Rdata"))
all_pests_PA_newleaf <- get(load("all_pests_PA_newleaf.Rdata"))

pest_data <- get(load("pest_data_ISDM_input_list.Rdata"))
meshEW1 <- get(load("meshEW1.Rdata"))

covariate_stack <- rast("covariate_stack.tif")

# For species with enough APHA data, subset to enable separate model fitting
species <- c("Acer pseudoplatanus",
             "Betula pendula",
             "Fagus sylvatica",
             "Fraxinus excelsior",
             "Picea abies",
             "Picea sitchensis",
             "Pinus sylvestris",
             "Quercus robur",
             "Tsuga heterophylla")

pest_data_sp <- pest_data[species]

PA_list <- list()
for(i in 1:length(pest_data_sp)){
  PA_list[[i]] <- pest_data_sp[[i]]$APHA_obs
  PA_list[[i]] <- st_as_sf(PA_list[[i]], coords = c("X", "Y"))
  st_set_crs(PA_list[[i]], st_crs(all_pests_PA_all))
}
names(PA_list) <- names(pest_data_sp)

# Ensure consistent NA values across raster layers ####
# separate out a layer which is correct
correct <- covariate_stack$broadleaf_area_Ha

covariate_stack[is.na(covariate_stack)] <- 0

covariate_stack <- mask(covariate_stack, 
                        correct,
                        maskvalues = NA,
                        updatevalue = NA)


# Remove covariate layers which won't be used ####
covariate_stack <- terra::subset(covariate_stack, "footAccessProp1km", negate = TRUE)
covariate_stack <- terra::subset(covariate_stack, "gardens_GB", negate = TRUE)
covariate_stack <- terra::subset(covariate_stack, "GBR_rails_total", negate = TRUE)
covariate_stack <- terra::subset(covariate_stack, "humanPopInRange1km_1", negate = TRUE)
covariate_stack <- terra::subset(covariate_stack, "humanPopInRange1km_2", negate = TRUE)
covariate_stack <- terra::subset(covariate_stack, "MENEvisits1km_TickSolve_1", negate = TRUE)
covariate_stack <- terra::subset(covariate_stack, "MENEvisits1km_TickSolve_2", negate = TRUE)
covariate_stack <- terra::subset(covariate_stack, "OSgreenSpaceAccessPoints1km", negate = TRUE)
covariate_stack <- terra::subset(covariate_stack, "ShannonDiversityNatEnv", negate = TRUE)
covariate_stack <- terra::subset(covariate_stack, "urb_suburb_prop_cover", negate = TRUE)


# Standardise raster layers ####
# 0-centred standardise 
covariate_stack_sd <- terra::scale(covariate_stack)

# Generate model formulae for random subsets of the covariate stack ####
n_mods <- 50 # Number of models per covariate subset n
n <- 20:5 # Number of covariates in subset

sets_all <- list()
for(i in 1:length(n)){
  sets_temp <- list()
  for(j in 1:n_mods){
    sets_temp[[j]] <- sample(names(covariate_stack_sd), size = n[i], replace = FALSE)
  }
  sets_all[[i]] <- sets_temp
}

sets_all_unlisted <- unlist(sets_all, recursive = FALSE)
flist_all <- names_to_formula(sets_all_unlisted, 
                           form_front = "Present ~",
                           form_back = '+ spde(geometry, model = matern) + Intercept(1, model = "linear", mean.linear = 0, prec.linear = 1)')


# Evaluate spatial covariate values at the observation locations ####
all_pests_PA_all <- eval_spatial_multi(data_sf = all_pests_PA_all, covariates = covariate_stack_sd)
all_pests_PA_newleaf <- eval_spatial_multi(data_sf = all_pests_PA_newleaf, covariates = covariate_stack_sd)

PA_list_eval <- list()
for(i in 1:length(PA_list)){
  PA_list_eval[[i]] <- eval_spatial_multi(data_sf = PA_list[[i]], covariates = covariate_stack_sd)
}
names(PA_list_eval) <- names(PA_list)


# Fit models ####
# Mesh with max.edge = 30km
matern <- inla.spde2.pcmatern(meshEW1, 
                              prior.sigma = c(1, 0.01), # P(sigma > 1) = 0.01
                              prior.range = c(10*1000, 0.01)) # P(range < 5km) = 0.01

mods1 <- run_multiform_bru(comp_list = flist_all,
                          pest_data = all_pests_PA_all,
                          covariates = covariate_stack_sd)

mods2 <- run_multiform_bru(comp_list = flist_all,
                           pest_data = all_pests_PA_newleaf,
                           covariates = covariate_stack_sd)


Fs_mods <- run_multiform_bru(comp_list = flist_all,
                          pest_data = Fagus_sylvatica_APHA,
                          covariates = covariate_stack_sd)

Fe_mods <- run_multiform_bru(comp_list = flist_all,
                             pest_data = Fraxinus_excelsior_APHA,
                             covariates = covariate_stack_sd)

Qr_mods <- run_multiform_bru(comp_list = flist_all,
                             pest_data = Quercus_robur_APHA,
                             covariates = covariate_stack_sd)

# Produce table of WAIC values ####
WAIC_table1 <- matrix(NA, nrow = length(mods1), ncol = 2)
for(i in 1:nrow(WAIC_table1)){
  WAIC_table1[i,1] <- paste(sets_all_unlisted[[i]], collapse = ", ")
  WAIC_table1[i,2] <- mods1[[i]]$waic$waic
}
WAIC_table1 <- as.data.frame(WAIC_table1)
colnames(WAIC_table1) <- c("covariates","WAIC")
WAIC_table1$WAIC <- as.numeric(WAIC_table1$WAIC)
WAIC_table1$delta_WAIC <- WAIC_table1$WAIC - min(WAIC_table1$WAIC)

WAIC_table2 <- matrix(NA, nrow = length(mods2), ncol = 2)
for(i in 1:nrow(WAIC_table2)){
  WAIC_table2[i,1] <- paste(sets_all_unlisted[[i]], collapse = ", ")
  WAIC_table2[i,2] <- mods2[[i]]$waic$waic
}
WAIC_table2 <- as.data.frame(WAIC_table2)
colnames(WAIC_table2) <- c("covariates","WAIC")
WAIC_table2$WAIC <- as.numeric(WAIC_table2$WAIC)
WAIC_table2$delta_WAIC <- WAIC_table2$WAIC - min(WAIC_table2$WAIC)

WAIC_table_Fs <- matrix(NA, nrow = length(Fs_mods), ncol = 2)
for(i in 1:nrow(WAIC_table_Fs)){
  WAIC_table_Fs[i,1] <- paste(sets_all_unlisted[[i]], collapse = ", ")
  WAIC_table_Fs[i,2] <- Fs_mods[[i]]$waic$waic
}
WAIC_table_Fs <- as.data.frame(WAIC_table_Fs)
colnames(WAIC_table_Fs) <- c("covariates","WAIC")
WAIC_table_Fs$WAIC <- as.numeric(WAIC_table_Fs$WAIC)
WAIC_table_Fs$delta_WAIC <- WAIC_table_Fs$WAIC - min(WAIC_table_Fs$WAIC)

WAIC_table_Fe <- matrix(NA, nrow = length(Fe_mods), ncol = 2)
for(i in 1:nrow(WAIC_table_Fe)){
  WAIC_table_Fe[i,1] <- paste(sets_all_unlisted[[i]], collapse = ", ")
  WAIC_table_Fe[i,2] <- Fe_mods[[i]]$waic$waic
}
WAIC_table_Fe <- as.data.frame(WAIC_table_Fe)
colnames(WAIC_table_Fe) <- c("covariates","WAIC")
WAIC_table_Fe$WAIC <- as.numeric(WAIC_table_Fe$WAIC)
WAIC_table_Fe$delta_WAIC <- WAIC_table_Fe$WAIC - min(WAIC_table_Fe$WAIC)

WAIC_table_Qr <- matrix(NA, nrow = length(Qr_mods), ncol = 2)
for(i in 1:nrow(WAIC_table_Qr)){
  WAIC_table_Qr[i,1] <- paste(sets_all_unlisted[[i]], collapse = ", ")
  WAIC_table_Qr[i,2] <- Qr_mods[[i]]$waic$waic
}
WAIC_table_Qr <- as.data.frame(WAIC_table_Qr)
colnames(WAIC_table_Qr) <- c("covariates","WAIC")
WAIC_table_Qr$WAIC <- as.numeric(WAIC_table_Qr$WAIC)
WAIC_table_Qr$delta_WAIC <- WAIC_table_Qr$WAIC - min(WAIC_table_Qr$WAIC)


setwd("/data/notebooks/rstudio-peter/outputs/WAIC_tables")
save(WAIC_table1, file = "WAIC_table_all.Rdata")
save(WAIC_table2, file = "WAIC_table_newleaf.Rdata")

save(WAIC_table_Fs, file = "WAIC_table_Fs.Rdata")
save(WAIC_table_Fe, file = "WAIC_table_Fe.Rdata")
save(WAIC_table_Qr, file = "WAIC_table_Qr.Rdata")
