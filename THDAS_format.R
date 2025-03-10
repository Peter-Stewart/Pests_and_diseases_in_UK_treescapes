# Packages ####
library(lubridate)
library(dplyr)
library(sf)
library(stringr)

# Set directories ####
p1 <- "/data/notebooks/rstudio-phytophthora/spatial_data/pest_data/THDAS/THDAS_csv"
p2 <- "/data/notebooks/rstudio-phytophthora/spatial_data/pest_data/THDAS/THDAS_Rdata"

# List of improper taxonomic names to flag
improper_taxa <- c("Animals",
                   "Insects etc",
                   "Woodpecker",
                   #"Acute Oak Decline",
                   #"Oak dieback",
                   "sawfly"#,
                   #"Top dying",
                   #"Virus",
                   #"Complex"
                   )

# Get list of csv files ####
setwd(p1)
file_list <- list.files(p1)

# Format files and save as .Rdata ####
THDAS_list <- list() # Empty list to store all results
THDAS_names <- str_sub(file_list, end = -5) # Names of dataframes minus .csv suffix

# Loop over all csv files
for(i in 1:length(file_list)){
  
  df <- read.csv(file_list[i])
  
  # Convert character variables to factors
  df$ENQUIRY_ID <- as.factor(df$ENQUIRY_ID)
  df$COUNTRY <- as.factor(df$COUNTRY)
  df$ORIGIN <- as.factor(df$ORIGIN)
  df$HOST.GENUS <- as.factor(df$HOST.GENUS)
  df$HOST.SPECIES <- as.factor(df$HOST.SPECIES)
  df$HOST.VARIETY <- as.factor(df$HOST.VARIETY)
  df$COMMON_NAME <- as.factor(df$COMMON_NAME)
  df$FAMILY <- as.factor(df$FAMILY)
  df$AGENT.GENUS <- as.factor(df$AGENT.GENUS)
  df$AGENT.SPECIES <- as.factor(df$AGENT.SPECIES)
  df$AGENT.VARIETY <- as.factor(df$AGENT.VARIETY)
  df$PEST_PATH <- as.factor(df$PEST_PATH)
  
  # Convert date to lubridate format
  df$ENQUIRY_DATE <- as_datetime(df$ENQUIRY_DATE)
  df$ENQUIRY_DATE <- date(df$ENQUIRY_DATE)
  
  # Remove empty end column
  if("...18" %in% colnames(df)){
    df <- df %>% select(-...18)
  }
  if("...19" %in% colnames(df)){
    df <- df %>% select(-...19)
  }
  
  # Remove observations with missing or incorrect location data 
  df <- df %>% filter(!is.na(EASTING)) %>% filter(!is.na(NORTHING)) %>% 
    filter(NORTHING != 0)
  
  # Flag observations which do not have proper taxonomic names
  df$taxon_flag <- ifelse(df$AGENT.GENUS %in% improper_taxa, 1, 0)
  
  # Convert to sf, and set CRS to British National Grid (EPSG 27700)
  df <- st_as_sf(df, coords = c("EASTING", "NORTHING"))
  df <- st_set_crs(df, 27700) 
  
  # Save files as .Rdata
  save(df, file = paste0(p2, "/", THDAS_names[i], ".Rdata"))
  THDAS_list[[i]] <- df
  rm(df)
}

# Save list of all files
names(THDAS_list) <- THDAS_names
save(THDAS_list, file = paste0(p2,"/","THDAS_list.Rdata"))

# Save dataframe version of THDAS data
THDAS <- THDAS_list
THDAS_df <- as.data.frame(do.call(rbind, THDAS))
THDAS_df$HOST.BINOMIAL <- rep(names(THDAS), times = sapply(THDAS, nrow))

THDAS_df$AGENT.BINOMIAL <- NA
for(j in 1:nrow(THDAS_df)){
  if(is.na(THDAS_df$AGENT.SPECIES[j])){
    THDAS_df$AGENT.BINOMIAL[j] <- paste(THDAS_df$AGENT.GENUS[j], "sp.")
  }else{
    THDAS_df$AGENT.BINOMIAL[j] <- paste(THDAS_df$AGENT.GENUS[j], THDAS_df$AGENT.SPECIES[j])
  }
}
THDAS_df$AGENT.BINOMIAL <- as.factor(THDAS_df$AGENT.BINOMIAL)
THDAS_df$HOST.BINOMIAL <- as.factor(THDAS_df$HOST.BINOMIAL)

save(THDAS_df, file = paste0(p2,"/","THDAS_df.Rdata"))
