# Packages ####
library(terra)
library(sf)
library(PointedSDMs)
library(dplyr)

# Custom functions ####
# Opposite of %in% 
'%notin%' <- function(x,y)!('%in%'(x,y))

# Load data ####
# Pest data - APHA
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/pest_data/APHA/APHA_Rdata")
APHA_df_sf <- get(load("APHA_df_sf.Rdata"))

# Pest data - THDAS
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/pest_data/THDAS/THDAS_Rdata")
THDAS_df <- get(load("THDAS_df.Rdata"))
THDAS_df <- THDAS_df %>% st_set_geometry(THDAS_df$geometry)
THDAS_df <- THDAS_df %>% filter(Year > 2000) %>% filter(taxon_flag == 0)
  
# Covariates - construct covariate raster stack
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/processed1km")
rast_list <- list.files(pattern = ".tif")

covariate_list <- list()
for(i in 1:length(rast_list)){
  covariate_list[[i]] <- rast(rast_list[i])
}
names(covariate_list) <- rast_list

covariate_stack <- rast(covariate_list)
names(covariate_stack) <- gsub(".tif", "", names(covariate_stack))
  
# GB boundary
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/baselayers")
GB_bound <- get(load("gb_mainland_boundary.Rdata")) # Clean UK boundary from Dambly et al.

# Convex hull around Enland/Wales boundary for the PA models
EW_islands_removed <- st_read("EW_boundary_islands_removed.shp")
EW_bound <- st_convex_hull(EW_islands_removed)


# Convert CRS units for all data to km ####
covariate_stack <- project(covariate_stack, "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs")
GB_bound <- st_transform(GB_bound, "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs")
EW_bound <- st_transform(EW_bound, "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs")
APHA_df_sf <- st_transform(APHA_df_sf, "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs")
THDAS_df <- st_transform(THDAS_df, "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs")

# Create meshes for INLA ####
# UK for ISDM
use_crs <- st_crs(GB_bound)
meshGB1 <- fm_mesh_2d(loc.domain = st_as_sfc(GB_bound),
                    boundary = GB_bound,
                    max.edge = c(30, 150),
                    cutoff = 1,
                    crs = use_crs)

meshGB2 <- fm_mesh_2d(loc.domain = st_as_sfc(GB_bound),
                      boundary = GB_bound,
                      max.edge = c(20, 150),
                      cutoff = 1,
                      crs = use_crs)

meshGB3 <- fm_mesh_2d(loc.domain = st_as_sfc(GB_bound),
                      boundary = GB_bound,
                      max.edge = c(10, 150),
                      cutoff = 1,
                      crs = use_crs)

meshGB4 <- fm_mesh_2d(loc.domain = st_as_sfc(GB_bound),
                      boundary = GB_bound,
                      max.edge = c(40, 150),
                      cutoff = 1,
                      crs = use_crs)

meshGB5 <- fm_mesh_2d(loc.domain = st_as_sfc(GB_bound),
                      boundary = GB_bound,
                      max.edge = c(5, 150),
                      cutoff = 1,
                      crs = use_crs)


# England/Wales for spatial GLM
use_crs <- st_crs(EW_bound)
meshEW1 <- fm_mesh_2d(loc.domain = st_as_sfc(EW_bound),
                    boundary = EW_bound,
                    max.edge = c(30, 150),
                    cutoff = 1,
                    crs = use_crs)

meshEW2 <- fm_mesh_2d(loc.domain = st_as_sfc(EW_bound),
                      boundary = EW_bound,
                      max.edge = c(20, 150),
                      cutoff = 1,
                      crs = use_crs)

meshEW3 <- fm_mesh_2d(loc.domain = st_as_sfc(EW_bound),
                      boundary = EW_bound,
                      max.edge = c(10, 150),
                      cutoff = 1,
                      crs = use_crs)

meshEW4 <- fm_mesh_2d(loc.domain = st_as_sfc(EW_bound),
                      boundary = EW_bound,
                      max.edge = c(5, 150),
                      cutoff = 1,
                      crs = use_crs)

# Subset THDAS data to observations within UK boundary
THDAS_df <- st_filter(THDAS_df, GB_bound)

# Subset APHA data to priority host trees which are in THDAS
THDAS_hosts <- as.character(unique(THDAS_df$HOST.BINOMIAL))

APHA_df_sf$Host_in_THDAS <- ifelse(as.character(APHA_df_sf$Host_binomial) %in% THDAS_hosts, "Y", "N")
APHA_df_sf$Host_in_THDAS <- as.factor(APHA_df_sf$Host_in_THDAS)

APHA_df_sf_newleaf <- APHA_df_sf %>% filter(Host_in_THDAS == "Y")
APHA_df_sf_all <- APHA_df_sf

# Subset to premises and plant types which are plants (not e.g. seeds) in wider environment
APHA_df_sf_newleaf <- APHA_df_sf_newleaf %>% 
  filter(Plant_group %in% c("Grown plants", "Growing crop") | Job_desc %in% c("Estblshd plntgs", "Wild Plant Survey")) %>%
  filter(Premise_type %in% c("Farm", "Garden", "Heathland", "Watercourse", "woodland", "Recreational Facility")) %>%
  filter(Pest_found %notin% "Not Tested")

APHA_df_sf_all <- APHA_df_sf_all %>% 
  filter(Plant_group %in% c("Grown plants", "Growing crop") | Job_desc %in% c("Estblshd plntgs", "Wild Plant Survey")) %>%
  filter(Premise_type %in% c("Farm", "Garden", "Heathland", "Watercourse", "woodland", "Recreational Facility")) %>%
  filter(Pest_found %notin% "Not Tested")

# Process pest data into correct format for IDM
# Store in list with separate entry for each priority host tree
dlist <- list()
for(sp in 1:length(THDAS_hosts)){
  # Subset observations for the relevant host tree
  THDAS_sp <- THDAS_df %>% filter(HOST.BINOMIAL == THDAS_hosts[sp])
  APHA_sp <- APHA_df_sf_newleaf %>% filter(Host_binomial == THDAS_hosts[sp])
  
  if(nrow(THDAS_sp) == 0 | nrow(APHA_sp) == 0) next
  
  # THDAS observations are presence-only - only need coordinates
  THDAS_obs <- as.data.frame(st_coordinates(THDAS_sp))
  
  # APHA observations are presence-absence - need coordinates and count of presences with total count
  # First obtain APHA coordinates
  APHA_obs <- APHA_sp
  APHA_coords <- st_coordinates(APHA_obs)
  APHA_obs <- cbind(APHA_coords, APHA_obs)
  APHA_obs <- APHA_obs %>% dplyr::select("X", "Y", "Pest_found", "Visit_date") %>% st_drop_geometry()
  
  # Then create indicator variable for pest presence/absence in APHA data
  APHA_obs$Present <- NA
  for(i in 1:nrow(APHA_obs)){
    if(is.na(APHA_obs$Pest_found[i]) == TRUE){
      APHA_obs$Present[i] <- 0
    }else if(APHA_obs$Pest_found[i] %in% c("-ve")){
      APHA_obs$Present[i] <- 0
    }else{
      APHA_obs$Present[i] <- 1
    }
  }
  APHA_obs$Surveyed <- 1
  
  # Transform indicator variable into count of presences and total counts for each location
  APHA_obs <- APHA_obs %>% select(-Pest_found) %>%
    group_by(X, Y, Visit_date, .groups = "keep") %>%
    summarise(across(c(Present, Surveyed), max)) %>%
    group_by(X, Y, .groups = "keep") %>%
    summarise(across(c(Present, Surveyed), sum))
    
  APHA_obs <- as.data.frame(APHA_obs)
  APHA_obs <- APHA_obs %>% select(X, Y, Present, Surveyed)
  
  # Combine output data into single list
  combined_list <- list(THDAS_obs, APHA_obs)
  names(combined_list) <- c("THDAS_obs", "APHA_obs")
  # Append this list to the main output list
  dlist[[sp]] <- combined_list
}
names(dlist) <- THDAS_hosts

# For total burden drivers models, also group pests of all hosts together
# All hosts
all_pests_PA_all <- APHA_df_sf_all %>% select(Pest_found, Visit_date)
all_pests_PA_all$Present <- NA
all_pests_PA_all$Surveyed <- 1
for(i in 1:nrow(all_pests_PA_all)){
  # NA observations = no pest observed -> Present = 0
  if(is.na(all_pests_PA_all$Pest_found[i]) == TRUE){
    all_pests_PA_all$Present[i] <- 0
  }else if(all_pests_PA_all$Pest_found[i] %in% c("-ve")){
    # Pest sample is -ve -> Present = 0
    all_pests_PA_all$Present[i] <- 0
  }else{
    # Anything else is named pest -> Present = 1
    all_pests_PA_all$Present[i] <- 1
  }
}
coords <- st_coordinates(all_pests_PA_all)
all_pests_PA_all <- cbind(coords, all_pests_PA_all)

all_pests_PA_all <- all_pests_PA_all %>% select(-Pest_found) %>%
  group_by(X, Y, Visit_date, .groups = "keep") %>%
  summarise(across(c(Present, Surveyed), max)) %>%
  group_by(X, Y, .groups = "keep") %>%
  summarise(across(c(Present, Surveyed), sum))

# All Newleaf priority hosts
all_pests_PA_newleaf <- APHA_df_sf_newleaf %>% select(Pest_found, Visit_date)
all_pests_PA_newleaf$Present <- NA
all_pests_PA_newleaf$Surveyed <- 1
for(i in 1:nrow(all_pests_PA_newleaf)){
  # NA observations = no pest observed -> Present = 0
  if(is.na(all_pests_PA_newleaf$Pest_found[i]) == TRUE){
    all_pests_PA_newleaf$Present[i] <- 0
  }else if(all_pests_PA_newleaf$Pest_found[i] %in% c("-ve")){
    # Pest sample is -ve -> Present = 0
    all_pests_PA_newleaf$Present[i] <- 0
  }else{
    # Anything else is named pest -> Present = 1
    all_pests_PA_newleaf$Present[i] <- 1
  }
}
coords <- st_coordinates(all_pests_PA_newleaf)
all_pests_PA_newleaf <- cbind(coords, all_pests_PA_newleaf)

all_pests_PA_newleaf <- all_pests_PA_newleaf %>% select(-Pest_found) %>%
  group_by(X, Y, Visit_date, .groups = "keep") %>%
  summarise(across(c(Present, Surveyed), max)) %>%
  group_by(X, Y, .groups = "keep") %>%
  summarise(across(c(Present, Surveyed), sum))

# Individual newleaf hosts (APHA data)
hosts_PA <- list()
for(sp in 1:length(THDAS_hosts)){
  sp_PA <- APHA_df_sf_newleaf %>% 
    filter(Host_binomial == THDAS_hosts[sp]) %>% 
    select(Pest_found, Visit_date)
  
  if(nrow(sp_PA) == 0){
    next
  }
  
  sp_PA$Present <- NA
  sp_PA$Surveyed <- 1
  
  for(i in 1:nrow(sp_PA)){
    # NA observations = no pest observed -> Present = 0
    if(is.na(sp_PA$Pest_found[i]) == TRUE){
      sp_PA$Present[i] <- 0
    }else if(sp_PA$Pest_found[i] %in% c("-ve")){
      # Pest sample is -ve -> Present = 0
      sp_PA$Present[i] <- 0
    }else{
      # Anything else is named pest -> Present = 1
      sp_PA$Present[i] <- 1
    }
  }
  coords <- st_coordinates(sp_PA)
  sp_PA <- cbind(coords, sp_PA)
  
  hosts_PA[[sp]] <- sp_PA
}
names(hosts_PA) <- THDAS_hosts


# THDAS and APHA data with all hosts grouped together
THDAS_all <- as.data.frame(st_coordinates(THDAS_df))
APHA_coords <- st_coordinates(APHA_df_sf_newleaf)
APHA_coords <- cbind(APHA_coords, APHA_df_sf_newleaf)
APHA_all <- APHA_coords %>% dplyr::select("X", "Y", "Pest_found", "Visit_date") %>% st_drop_geometry()
APHA_all$Present <- NA
for(i in 1:nrow(APHA_all)){
  if(is.na(APHA_all$Pest_found[i]) == TRUE){
    APHA_all$Present[i] <- 0
  }else if(APHA_all$Pest_found[i] %in% c("-ve")){
    APHA_all$Present[i] <- 0
  }else{
    APHA_all$Present[i] <- 1
  }
}
APHA_all$Surveyed <- 1

APHA_all <- APHA_all %>% select(-Pest_found) %>%
  group_by(X, Y, Visit_date, .groups = "keep") %>%
  summarise(across(c(Present, Surveyed), max)) %>%
  group_by(X, Y, .groups = "keep") %>%
  summarise(across(c(Present, Surveyed), sum))

APHA_all <- as.data.frame(APHA_all)
APHA_all <- APHA_all %>% select(X, Y, Present, Surveyed)
pest_data_all <- list(THDAS_all, APHA_all)

names(pest_data_all) <- c("THDAS_obs", "APHA_obs")
pest_data_all_list <- list(pest_data_all)
names(pest_data_all_list) <- "all"

# Save inputs for PointedSDMs / inlabru ####
setwd("/data/notebooks/rstudio-peternewleafnew/outputs/INLA_inputs")

# Pest data
save(dlist, file = "pest_data_ISDM_input_list.Rdata")
save(pest_data_all_list, file = "pest_data_all_ISDM_input_list.Rdata")
save(all_pests_PA_all, file = "all_pests_PA_all.Rdata")
save(all_pests_PA_newleaf, file = "all_pests_PA_newleaf.Rdata")
save(hosts_PA, file = "hosts_PA.Rdata")
save(APHA_df_sf_all, file = "APHA_df_sf_all.Rdata")
save(APHA_df_sf_newleaf, file = "APHA_df_sf_newleaf.Rdata")

# GB meshes for ISDM
save(meshGB1, file = "meshGB1.Rdata")
save(meshGB2, file = "meshGB2.Rdata")
save(meshGB3, file = "meshGB3.Rdata")
save(meshGB4, file = "meshGB4.Rdata")
save(meshGB5, file = "meshGB5.Rdata")

# England/Wales meshes for drivers models
save(meshEW1, file = "meshEW1.Rdata")
save(meshEW2, file = "meshEW2.Rdata")
save(meshEW3, file = "meshEW3.Rdata")
save(meshEW4, file = "meshEW4.Rdata")

# Covariates
terra::writeRaster(covariate_stack, file = "covariate_stack.tif", overwrite = TRUE)

# GB/EW Boundaries
save(GB_bound, file = "gb_mainland_boundary.Rdata")
save(EW_bound, file = "EW_mainland_boundary.Rdata")





#################################################################
# Code for checking crs matches for all raster layers ####
crs_list <- list()
for(i in 1:length(covariate_list)){
  crs_list[[i]] <- terra::crs(covariate_list[[i]])
}
names(crs_list) <- rast_list

check_mat <- matrix(NA, length(crs_list), length(crs_list))
for(i in 1:nrow(check_mat)){
  for(j in 1:ncol(check_mat)){
    if(i == j){
      next
    }else{
      check_mat[i,j] <- ifelse(crs_list[[i]]==crs_list[[j]], TRUE, FALSE)
    }
  }
}
rownames(check_mat) <- names(crs_list)
colnames(check_mat) <- names(crs_list)
