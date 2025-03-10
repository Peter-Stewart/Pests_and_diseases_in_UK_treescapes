# Packages ####
library(readxl)
library(dplyr)
library(lubridate)
library(rnrfa)
library(sf)
library(stringr)
library(PostcodesioR)

# Custom fumctions ####
# Opposite of %in% 
'%notin%' <- function(x,y)!('%in%'(x,y))

# Get coordinates from postcodes
postcode_to_coords <- function(df){
  if("Postcode" %in% names(df) == FALSE){
    stop("Variable named 'Postcode' is required!")
  }
  
  # Create columns to store northing/easting obtained from postcode
  df$northing_post <- NA
  df$easting_post <- NA
  
  # If postcode is NA skip to next iteraction
  if(is.na(df$Postcode[i])) next
  
  # Loop over all observations
  for(i in 1:nrow(df)){
    if(i == 1){
      # If first observation then look up postcode
      dat <- postcode_lookup(paste(df$Postcode[i]))
      df$northing_post[i] <- dat$northings
      df$easting_post[i] <- dat$eastings
    }else{
      # If postcode same as previous obs then copy coordinates to save time
      if(paste(df$Postcode[i]) == paste(df$Postcode[i-1])){
        df$northing_post[i] <- df$northing_post[i-1]
        df$easting_post[i] <- df$easting_post[i-1]
      }else{
        # If postcode not same as previous obs then look it up
        dat <- postcode_lookup(paste(df$Postcode[i]))
        df$northing_post[i] <- dat$northings
        df$easting_post[i] <- dat$eastings
      }
    }
  }
  return(df)
}

# Calculate euclidean distance
euclid_dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

# Load data ####
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/pest_data/APHA")
file_list <- list.files(pattern = ".xlsx")

APHA_df1 <- read_excel(file_list[1], skip = 4)
APHA_df2 <- read_excel(file_list[2], skip = 4)

# GB baselayer for checking observations are within correct area
gb <- st_read("/data/notebooks/rstudio-phytophthora/spatial_data/baselayers/NUTS_Level_1_January_2018_GCB_in_the_United_Kingdom.shp")
# Remove Northern Ireland and Scotland (APHA data should only be from England and Wales)
gb_main <- gb %>% filter(nuts118cd != "UKN") %>% filter(nuts118cd != "UKM") 
gb_main_comb <- st_union(gb_main)

# OS 100km grid squares 
OS_grid <- st_read("/data/notebooks/rstudio-phytophthora/spatial_data/baselayers/OS_grid_squares_100km.shp")

# Remove last row of each dataframe as this contains summary stats
APHA_df1 <- APHA_df1[-nrow(APHA_df1),]
APHA_df2 <- APHA_df2[-nrow(APHA_df2),]

# Combine into one dataframe ####
APHA_df <- rbind(APHA_df1, APHA_df2)
APHA_df <- as.data.frame(APHA_df)
rm(APHA_df1, APHA_df2)

# Rename columns to remove spaces, special characters etc. ####
APHA_df <- APHA_df %>% rename(Month_int = `Month No.`) %>%
  rename(Premise_type = `Premise Type`) %>%
  rename(Postal_town = `Postal Town`) %>%
  rename(Grid_ref = `Grid Ref`) %>%
  rename(Risk_rating = `Risk Rating`) %>%
  rename(Area_code = `Area Code`) %>%
  rename(Registered_to_passport = `Registered to Passport Ticked on Client Record?`) %>%
  rename(Latest_PP_application_date = `Latest PP Application Date`) %>%
  rename(Visit_date = `Visit Date`) %>%
  rename(CRN_visit_date = `CRN + Visit Date`) %>%
  rename(Application_ID = `Application ID`) %>%
  rename(Work_category = `Work Category`) %>%
  rename(Job_code = `Job Code`) %>%
  rename(Job_desc = `Job Desc`) %>%
  rename(Origin_country = `Country of Origin`) %>%
  rename(Plant_group = `Plant Group`) %>%
  rename(Host_genus = `Host Genus`) %>%
  rename(Host_sp = `Host Species`) %>%
  rename(Inspection_ID = `Inspection ID`) %>%
  rename(Quantity_inspected = `Quantity Inspected`) %>%
  rename(Quantity_inspected_units = `Inspection Unit`) %>%
  rename(Sample_ID = `Sample Id`) %>%
  rename(Suspect_pest = `Suspect Pest`) %>%
  rename(Sample_comments = `Sample Comments`) %>%
  rename(Pest_found = `Pest Found`)

# Convert character variables to factors ####
APHA_df$Premise_type <- as.factor(APHA_df$Premise_type)
APHA_df$Postal_town <- as.factor(APHA_df$Postal_town)
APHA_df$Postcode <- as.factor(APHA_df$Postcode)
APHA_df$Risk_rating <- as.factor(APHA_df$Risk_rating)
APHA_df$Region <- as.factor(APHA_df$Region)
APHA_df$Area_code <- as.factor(APHA_df$Area_code)
APHA_df$Registered_to_passport <- as.factor(APHA_df$Registered_to_passport)
APHA_df$Application_ID <- as.factor(APHA_df$Application_ID)
APHA_df$Work_category <- as.factor(APHA_df$Work_category)
APHA_df$Job_code <- as.factor(APHA_df$Job_code)
APHA_df$Job_desc <- as.factor(APHA_df$Job_desc)
APHA_df$Origin_country <- as.factor(APHA_df$Origin_country)
APHA_df$Plant_group <- as.factor(APHA_df$Plant_group)

APHA_df$Inspection_ID <- as.factor(APHA_df$Inspection_ID)
APHA_df$Quantity_inspected_units <- as.factor(APHA_df$Quantity_inspected_units)
APHA_df$Sample_ID <- as.factor(APHA_df$Sample_ID)
APHA_df$Suspect_pest <- as.factor(APHA_df$Suspect_pest)
APHA_df$Sample_comments <- as.factor(APHA_df$Sample_comments)
APHA_df$Pest_found <- as.factor(APHA_df$Pest_found)


# Convert dates to lubridate format ####
APHA_df$Latest_PP_application_date <- lubridate::as_datetime(APHA_df$Latest_PP_application_date)
APHA_df$Visit_date <- lubridate::as_date(APHA_df$Visit_date)
  
# Convert other variables as required ####
APHA_df$Month_int <- as.integer(APHA_df$Month_int)

# Generate host binomials from genus and species names ####
APHA_df$Host_sp <- ifelse(is.na(APHA_df$Host_sp), "sp.", APHA_df$Host_sp)

APHA_df$Host_binomial <- paste(APHA_df$Host_genus, APHA_df$Host_sp)
APHA_df$Host_binomial <- as.factor(APHA_df$Host_binomial)
APHA_df$Host_genus <- as.factor(APHA_df$Host_genus)
APHA_df$Host_sp <- as.factor(APHA_df$Host_sp)

# Add row index for keeping track of observations ####
APHA_df$N <- 1:nrow(APHA_df)

# Data cleaning (not including removal of invalid coordinates) ####
# Some records contain e.g. "INPUT ERROR PLEASE IGNORE" in the comments - filter out
ignores1 <- grep("ERROR", APHA_df$Sample_comments, ignore.case = TRUE)
ignores2 <- grep("IGNORE", APHA_df$Sample_comments, ignore.case = TRUE)
ignores <- unique(c(ignores1,ignores2))
APHA_df <- APHA_df[-ignores,]

# Deal with invalid coordinates ####
# First convert grid references to easting / northing
coords <- osg_parse(APHA_df$Grid_ref)
APHA_df$northing <- coords$northing
APHA_df$easting <- coords$easting

# Remove points with no coordinates
APHA_df <- APHA_df %>% filter(!is.na(northing)) %>% filter(!is.na(easting))

# In some cases grid reference has been entered in the comments - extract from here
# Subset cases where grid reference is mentioned in comments
comment_coords <- APHA_df %>% filter(grepl("grid ref", Sample_comments, ignore.case = TRUE))
# Extract the grid reference from the comments
ref_str <- comment_coords$Sample_comments
new_ref <- str_match(ref_str, "[A-Z]{2}\\s*\\d{5}\\s*\\d{5}")
new_ref <- gsub(" ","", new_ref) # Remove spaces from grid references
new_ref <- ifelse(grepl("EF", new_ref) == TRUE, NA, new_ref) # Remove invalid grid square "EF"
# Replace invalid grid reference with new grid reference extracted from comments
for(i in 1:nrow(comment_coords)){
  if(is.na(new_ref[i])) next
  else{
    comment_coords$Grid_ref[i] <- new_ref[i]
  }
}
coords_new <- osg_parse(comment_coords$Grid_ref)
comment_coords$northing <- coords_new$northing
comment_coords$easting <- coords_new$easting
# Merge back in to main dataframe
APHA_df <- APHA_df %>% filter(as.factor(N) %notin% as.factor(comment_coords$N))
APHA_df <- rbind(APHA_df, comment_coords)

# Some other points have incorrect grid references which place them in the sea, or in Scotland!
# Flag these points for checking
APHA_df_sf_temp <- st_as_sf(APHA_df, coords = c("easting", "northing"))
APHA_df_sf_temp <- st_set_crs(APHA_df_sf_temp, 27700)
APHA_df_sf_temp2 <- st_filter(APHA_df_sf_temp, gb_main)
APHA_df_sf_temp <- APHA_df_sf_temp %>% filter(as.factor(N) %notin% as.factor(APHA_df_sf_temp2$N))
index_out_of_bounds <- APHA_df_sf_temp$N # Get indices for observations which are outside England/Wales
rm(APHA_df_sf_temp, APHA_df_sf_temp2) # Remove the temporary sf dataframes

APHA_df_out <- APHA_df %>% filter(as.factor(N) %in% as.factor(index_out_of_bounds))
APHA_df <- APHA_df %>% filter(as.factor(N) %notin% as.factor(index_out_of_bounds))

# Use postcode to generate coordinates
APHA_df_out <- postcode_to_coords(APHA_df_out)
APHA_df_out <- APHA_df_out %>% filter(!is.na(northing_post))

# Small number of points appear to actually be in Scotland
# remove as sampling effort not comparable
APHA_df_sf_temp <- st_as_sf(APHA_df_out, coords = c("easting_post", "northing_post"))
APHA_df_sf_temp <- st_set_crs(APHA_df_sf_temp, 27700)
APHA_df_sf_temp <- st_filter(APHA_df_sf_temp, gb_main)
APHA_df_out <- APHA_df_out %>% filter(as.factor(N) %in% as.factor(APHA_df_sf_temp$N))
rm(APHA_df_sf_temp)

# Replace northing/easting coordinate with postcode coordinate for APHA_df_out
APHA_df_out$northing <- APHA_df_out$northing_post
APHA_df_out$easting <- APHA_df_out$easting_post

# Merge back into main dataframe
APHA_df_out <- APHA_df_out %>% select(-northing_post, -easting_post)
APHA_df <- rbind(APHA_df, APHA_df_out)

# Check for other points where grid reference does not match postcode
# First get coordinates from postcodes for all observations
APHA_df <- postcode_to_coords(APHA_df)

# Then calculate euclidean distance between grid ref and postcode coordinates
APHA_df$coords_diff <- NA
for(i in 1:nrow(APHA_df)){
  APHA_df$coords_diff[i] <- euclid_dist(x1 = c(APHA_df$easting[i], APHA_df$northing[i]), 
                                        x2 = c(APHA_df$easting_post[i], APHA_df$northing_post[i]))
}

dist_threshold <- 50000 # Distance threshold (in m) for observation to be considered mismatch

# Subset observations which exceed distance threshold
APHA_to_check <- APHA_df %>% filter(coords_diff >= dist_threshold)

# Check for cases where grid ref. has grid cell code wrong
# Subset cases where line connecting the postcode/grid ref coordinates is vertical or horizontal
APHA_to_check$easting_diff <- abs(APHA_to_check$easting - APHA_to_check$easting_post)
APHA_to_check$northing_diff <- abs(APHA_to_check$northing - APHA_to_check$northing_post)

diff_threshold <- 5000 # Threshold (in m) for x or y distance
APHA_to_check2 <- APHA_to_check %>% filter(easting_diff <= diff_threshold | northing_diff <= diff_threshold)
# For these observations use OS square from postcode with numbers from grid ref
temp_sf <- st_as_sf(APHA_to_check2, coords = c("easting_post", "northing_post")) %>% st_set_crs(27700)
postcode_grids <- st_join(temp_sf, OS_grid, join = st_within) %>% select(N, tile_name)
rm(temp_sf)
postcode_grids <- st_drop_geometry(postcode_grids)
APHA_to_check2 <- merge(APHA_to_check2, postcode_grids, by = "N")
APHA_to_check2$Grid_ref_new <- paste0(APHA_to_check2$tile_name,
                                      substring(APHA_to_check2$Grid_ref, 3))
coords <- osg_parse(APHA_to_check2$Grid_ref_new)
APHA_to_check2$northing <- coords$northing
APHA_to_check2$easting <- coords$easting
APHA_to_check2 <- APHA_to_check2 %>% select(-Grid_ref_new)

# Check distances now
APHA_to_check2$coords_diff2 <- NA
for(i in 1:nrow(APHA_to_check2)){
  APHA_to_check2$coords_diff2[i] <- euclid_dist(x1 = c(APHA_to_check2$easting[i], APHA_to_check2$northing[i]), 
                                        x2 = c(APHA_to_check2$easting_post[i], APHA_to_check2$northing_post[i]))
}
# Most points are now close but there is a discontinuity with some points above 19km
# Merge these points into APHA_to_check3 for manual checking
APHA_to_check2b <- APHA_to_check2 %>% filter(coords_diff2 > 19000) %>% 
  select(-tile_name, -coords_diff2)

# Subset observations which are not likely due to grid cell error
APHA_to_check3 <- APHA_to_check %>% filter(N %notin% APHA_to_check2$N)
APHA_to_check3 <- rbind(APHA_to_check3, APHA_to_check2b)
APHA_to_check2 <- APHA_to_check2 %>% filter(N %notin% APHA_to_check3$N)

# Manual checking
# Create crosswalk spreadsheet with observation index (N) and action to be taken (action)
dat <- APHA_to_check3 %>% select(N, Postcode)
#write.csv(dat, file = "manual_check_crosswalk_blank.csv")

# Fill in crosswalk by manually checking the records
# as of 1st Dec have manually checked all observations for postcodes that appear >4 times
# Load filled crosswalk spreadsheet
crosswalk <- read.csv("manual_check_crosswalk_filled_v1.csv", 
                      stringsAsFactors = TRUE,
                      na.strings = c("","NA"))

crosswalk <- crosswalk %>% select(-Postcode)

# Merge crosswalk with data which were being checked
APHA_to_check3 <- merge(APHA_to_check3, crosswalk, by = "N")

# Correct coordinates based on action from crosswalk
# Remove observations which could not be / were not manually verified
APHA_to_check3 <- APHA_to_check3 %>% filter(action %notin% c("not_identifiable", "remove")) %>%
  filter(!is.na(action))

# Correct other observations
for(i in 1:nrow(APHA_to_check3)){
  # If grid ref is correct then do nothing
  if(APHA_to_check3$action[i] == "grid_ref_correct"){
    next
  } # If post code is correct then use postcode northing/easting
  else if(APHA_to_check3$action[i] == "post_code_correct"){
    APHA_to_check3$northing[i] <- APHA_to_check3$northing_post[i]
    APHA_to_check3$easting[i] <- APHA_to_check3$easting_post[i]
    
  } # If grid reference is manually provided then use that and recalculate northing/easting
  else if(grepl("[A-Z]{2}\\s*\\d{5}\\s*\\d{5}", APHA_to_check3$action[i]) == TRUE){
    APHA_to_check3$Grid_ref[i] <- paste(APHA_to_check3$action[i])
    coords <- osg_parse(APHA_to_check3$Grid_ref[i])
    APHA_to_check3$northing[i] <- coords$northing
    APHA_to_check3$easting[i] <- coords$easting
  }
}

# Merge all checked data back into main dataframe
APHA_df <- APHA_df %>% filter(N %notin% APHA_to_check$N)
APHA_to_check2 <- APHA_to_check2 %>% select(colnames(APHA_df))
APHA_to_check3 <- APHA_to_check3 %>% select(colnames(APHA_df))

APHA_df <- rbind(APHA_df, APHA_to_check2, APHA_to_check3)

# Convert to sf, and set CRS to British National Grid (EPSG 27700) ####
APHA_df_sf <- st_as_sf(APHA_df, coords = c("easting", "northing"))
APHA_df_sf <- st_set_crs(APHA_df_sf, 27700) 

# Save outputs ####
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/pest_data/APHA/APHA_Rdata")

save(APHA_df, file = "APHA_df.Rdata")
save(APHA_df_sf, file = "APHA_df_sf.Rdata")
