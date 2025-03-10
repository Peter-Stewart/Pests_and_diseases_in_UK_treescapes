setwd("/data/notebooks/rstudio-phytophthora/spatial_data/processed1km")
file_names <- list.files(pattern = "tif")
var_names <- gsub(".tif", "", file_names)

DAG_var_names <- names(R_d14a)

var_names_crosswalk <- as.data.frame(DAG_var_names)
var_names_crosswalk$rast_names <- NA

# Manually input
var_names_crosswalk$rast_names[1] <- "ancient_woodland_pctcover"
var_names_crosswalk$rast_names[2] <- "LCC_afforestation_1km"
var_names_crosswalk$rast_names[3] <- "broadleaf_area_Ha"
var_names_crosswalk$rast_names[4] <- "conifer_area_Ha"
var_names_crosswalk$rast_names[5] <- "LCC_deforestation_1km"
var_names_crosswalk$rast_names[6] <- "elevation_median_1km"
var_names_crosswalk$rast_names[7] <- "humanPop"
var_names_crosswalk$rast_names[8] <- "recreation_weekly_1km"
var_names_crosswalk$rast_names[9] <- "urb_suburb_areaHa"
var_names_crosswalk$rast_names[10] <- "vpd_total"
var_names_crosswalk$rast_names[11] <- "woodland_area_Ha"
var_names_crosswalk$rast_names[12] <- "woodland_connectivity_sigma3"
var_names_crosswalk$rast_names[13] <- "dist_BCP"
var_names_crosswalk$rast_names[14] <- "dist_garden_GB_m"
var_names_crosswalk$rast_names[15] <- NA
var_names_crosswalk$rast_names[16] <- NA
var_names_crosswalk$rast_names[17] <- NA

setwd("/data/notebooks/rstudio-peter/scripts/analysis")
write.csv(var_names_crosswalk, "var_names_crosswalk.csv")
