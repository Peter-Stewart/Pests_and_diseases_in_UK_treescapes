# Packages ####
library(terra)
library(sf)
library(dplyr)

# Create 1km grid squares for GB ####
# Use one of the layers from processed1km folder as a template
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/processed1km")
GB_grid <- rast("broadleaf_area_Ha.tif")

# Roads - total length ####
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/GBR_roads")
GBR_roads <- st_read("GBR_roads.shp")
GBR_roads_crs <- GBR_roads %>% st_transform(27700) # transform to British National Grid (EPSG 27700)
GBR_roads_vect <- vect(GBR_roads_crs) %>% as.lines()
GBR_roads_total <- rasterizeGeom(GBR_roads_vect, GB_grid, "length") # Calculate total length
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/processed1km")
writeRaster(GBR_roads_total, file = "GBR_roads_total.tif")  # Save results

# Railways - total length ####
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/GBR_rails")
GBR_rails <- st_read("GBR_rails.shp")
GBR_rails_crs <- GBR_rails %>% st_transform(27700) # transform to British National Grid (EPSG 27700)
GBR_rails_vect <- vect(GBR_rails_crs) %>% as.lines()
GBR_rails_total <- rasterizeGeom(GBR_rails_vect, GB_grid, "length") # Calculate total length
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/processed1km")
writeRaster(GBR_rails_total, file = "GBR_rails_total.tif") # Save results

# Parks and gardens - distance to nearest ####
# Also consider option to recalculate for older pest records which may pre-date newer gardens
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/Parks and Gardens")
gardens_England <- st_read("England/ParksAndGardens_04April2023.shp")
gardens_Scotland <- st_read("Scotland/Gardens_and_Designed_Landscapes.shp")
gardens_Wales <- st_read("Wales/cadw_rhpg_registeredareas.shp")

gardens_England_vect <- vect(gardens_England) %>% as.points()
gardens_Scotland_vect <- vect(gardens_Scotland) %>% as.points()
gardens_Wales_vect <- vect(gardens_Wales) %>% as.points()

# This works because max number of gardens in a 1km square is 1
gardens_England_rast <- rasterize(gardens_England_vect, GB_grid)
gardens_Scotland_rast <- rasterize(gardens_Scotland_vect, GB_grid)
gardens_Wales_rast <- rasterize(gardens_Wales_vect, GB_grid)

# Combine into one raster (as dist. to nearest garden for cell in England may be e.g. in Wales)
gardens_GB_rast <- terra::merge(gardens_England_rast, gardens_Scotland_rast)
gardens_GB_rast <- terra:: merge(gardens_GB_rast, gardens_Wales_rast)

# Calculate distances
dist_garden_GB_m <- terra::distance(gardens_GB_rast)

# Save results
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/processed1km")
writeRaster(gardens_GB_rast, file = "gardens_GB.tif")
writeRaster(dist_garden_GB_m, file = "dist_garden_GB_m.tif")

# Priority tree layers - convert crs only ####
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/priority_tree_layers_tif")
file_list <- list.files(getwd(), pattern = ".tif")
for(i in 1:length(file_list)){
  dat <- rast(file_list[i])
  crs(dat) <- crs(GB_grid)
  writeRaster(dat,
              file = paste0("EPSG_27700/",file_list[i]),
              overwrite = TRUE)
}

# Elevation - median ####
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/DEM")
DEM <- rast("DEM_OS36.tif")
DEM_1km <- terra::aggregate(DEM,
                            fact = 40,
                            fun = median)
DEM_1km <- terra::crop(DEM_1km, GB_grid)
DEM_1km <- extend(DEM_1km, GB_grid)
ext(DEM_1km) <- ext(GB_grid)
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/processed1km")
writeRaster(DEM_1km, file = "elevation_median_1km.tif")

# Human pop. and human pop. in range - transform CRS only ####
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/processed1km")
pop1 <- rast("humanPop.tif")
pop2 <- rast("humanPopInRange1km.tif")

pop1 <- project(pop1, GB_grid)
pop2 <- project(pop2, GB_grid)

writeRaster(pop1, file = "humanPop.tif", overwrite = TRUE)
writeRaster(pop2, file = "humanPopInRange1km.tif", overwrite = TRUE)


# Canopy height - transform CRS and crop ####
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/Canopy_height/ETH_GlobalCanopyHeight_2020")
file_list <- list.files(pattern = ".tif")
canopy_list <- list()
for(i in 1:length(file_list)){
  canopy_list[[i]] <- rast(file_list[i])
}
canopy_sprc <- sprc(canopy_list)
canopy <- merge(canopy_sprc)
canopy_GB <- project(canopy, GB_grid)
canopy_GB_crop <- terra::crop(canopy_GB, GB_grid, mask = TRUE)
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/processed1km")
writeRaster(canopy_GB_crop, file = "canopy_height.tif", overwrite = TRUE)


# Ancient woodlands - total area (ha) ####
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/Ancient_woodlands/England")
AW_England <- st_read("Ancient_Woodland___Natural_England.shp")
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/Ancient_woodlands/Wales")
AW_Wales <- st_read("NRW_ANCIENT_WOODLAND_INVENTORY_2021Polygon.shp")
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/Ancient_woodlands/Scotland")
AW_Scotland <- st_read("AWI_SCOTLAND.shp")

AW_England$country <- "England"
AW_Wales$country <- "Wales"
AW_Scotland$country <- "Scotland"

AW_England <- AW_England %>% select(NAME, THEMNAME, AREA, country, geometry)
AW_Wales <- AW_Wales %>% select(unique_id, catagory_n, area_hecta, country, geometry)
AW_Scotland <- AW_Scotland %>% select(WOOD_ID, ANTIQUITY, HECTARE, country, geometry)

cnames <- c("site", "type", "area_ha", "country", "geometry")
colnames(AW_England) <- cnames
colnames(AW_Wales) <- cnames
colnames(AW_Scotland) <- cnames

AW_all <- rbind(AW_England, AW_Wales, AW_Scotland)

AW_all_vect <- vect(AW_all) %>% as.polygons()

AW_total <- rasterize(AW_all_vect, GB_grid, cover = TRUE)
AW_total[is.na(AW_total)] <- 0
AW_total <- mask(AW_total, GB_grid, maskvalues = NA, updatevalue = NA)

setwd("/data/notebooks/rstudio-phytophthora/spatial_data/Ancient_woodlands")
save(AW_all, file = "Ancient_woodlands_all.Rdata")

setwd("/data/notebooks/rstudio-phytophthora/spatial_data/processed1km")
writeRaster(AW_total, file = "ancient_woodland_pctcover.tif")


# Flow accumulation - transform crs and crop ####
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/HydroSHEDS/Flow_accumulation_ACA")
flow_accumulation <- rast("hyd_eu_aca_15s.tif")
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/HydroSHEDS/HydroSHEDS_land_mask")
land_mask <- rast("hyd_eu_msk_15s.tif")

flow_accumulation_mask <- mask(flow_accumulation, 
                               land_mask, 
                               maskvalues = 2, 
                               updatevalue = NA)

flow_accumulation_GB <- project(flow_accumulation_mask, GB_grid)

setwd("/data/notebooks/rstudio-phytophthora/spatial_data/processed1km")
writeRaster(flow_accumulation_GB, file = "flow_accumulation_ha.tif")


# Rivers - total length ####
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/HydroSHEDS/HydroRIVERS")
rivers <- st_read("HydroRIVERS_v10_eu.shp")
rivers_crs <- rivers %>% st_transform(27700)
rivers_vect <- vect(rivers_crs) %>% as.lines()
rivers_total <- rasterizeGeom(rivers_vect, GB_grid, "length")
rivers_total_crop <- terra::crop(rivers_total, GB_grid, mask = TRUE)

setwd("/data/notebooks/rstudio-phytophthora/spatial_data/processed1km")
writeRaster(rivers_total_crop, file = "rivers_total.tif")

# Vascular plant alpha-diversity (Sabatini et al. 2022) - transform CRS and mask ####
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/Sabatini_vascular_plant_alpha_diversity")

alpha_diversity <- rast("w3_tile_joint_sr1000.tif")

alpha_diversity_crs <- project(alpha_diversity, GB_grid)
alpha_diversity_crs_crop <- terra::crop(alpha_diversity_crs, GB_grid, mask = TRUE)

setwd("/data/notebooks/rstudio-phytophthora/spatial_data/processed1km")
writeRaster(alpha_diversity_crs_crop, file = "vascular_plant_alpha_diversity.tif")


# Land cover change - aggregate to derive several measures ####
setwd("/data/notebooks/rstudio-phytophthora/assets/landcoverchange-1990-2015-25m-raster-gb-v1/data/07b6e5e9-b766-48e5-a28c-5b3e35abecc0")

LCC <- rast("LCC_GB_1990_to_2015.tif")

# Afforestation - layer 5 = 1
LCC_afforestation <- LCC$Layer_5
LCC_afforestation$afforestation <- terra::ifel(LCC_afforestation$Layer_5 == 1, 1, 0)
LCC_afforestation <- LCC_afforestation$afforestation

LCC_afforestation_1km <- terra::aggregate(LCC_afforestation,
                            fact = 40,
                            fun = mean)

LCC_afforestation_1km <- project(LCC_afforestation_1km, GB_grid)
LCC_afforestation_1km <- terra::crop(LCC_afforestation_1km, GB_grid, mask = TRUE)

# Deforestation - layer 4 = 1
LCC_deforestation <- LCC$Layer_4
LCC_deforestation$deforestation <- terra::ifel(LCC_deforestation$Layer_4 == 1, 1, 0)
LCC_deforestation <- LCC_deforestation$deforestation

LCC_deforestation_1km <- terra::aggregate(LCC_deforestation,
                                          fact = 40,
                                          fun = mean)
LCC_deforestation_1km <- project(LCC_deforestation_1km, GB_grid)
LCC_deforestation_1km <- terra::crop(LCC_deforestation_1km, GB_grid, mask = TRUE)

setwd("/data/notebooks/rstudio-phytophthora/spatial_data/processed1km")
writeRaster(LCC_afforestation_1km, file = "LCC_afforestation_1km.tif")
writeRaster(LCC_deforestation_1km, file = "LCC_deforestation_1km.tif")


# VPD - sum from 2000 to 2022 ####
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/monthly_VPD")
file_list <- list.files(pattern = ".tif$")
vpd <- rast(file_list[40:length(file_list)])
vpd_crs <- project(vpd, GB_grid)
vpd_total <- sum(vpd_crs)
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/processed1km")
writeRaster(vpd_total, file = "vpd_total.tif")


# Border Control Point (BCP) locations - distance to nearest ####
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/BCP_locations")
BCP <- st_read("Border Control Point (BCPs) Locations.kml")
BCP <- BCP %>% select(Name) %>% st_zm()
BCP_crs <- BCP %>% st_transform(27700)

BCP_vect <- vect(BCP_crs) %>% as.points()

BCP_rast <- rasterize(BCP_vect, GB_grid)
BCP_rast_mask <- crop(BCP_rast, GB_grid, mask = TRUE) # Need to mask to remove points in NI

dist_BCP <- terra::distance(BCP_rast_mask)
dist_BCP_mask <- crop(dist_BCP, GB_grid, mask = TRUE)

setwd("/data/notebooks/rstudio-phytophthora/spatial_data/processed1km")
writeRaster(dist_BCP_mask, file = "dist_BCP.tif")


# Woodland connectivity ####
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/processed1km")
woodland_area <- rast("woodland_area_Ha.tif") 
broadleaf_area <- rast("broadleaf_area_Ha.tif")
conifer_area <- rast("conifer_area_Ha.tif")

# Create weights matrix
m <- rast(ncols=180,nrows=180, xmin=0)
weights <- focalMat(m,
                    3, # Value of sigma
                    "Gauss")

# Calculate connectivity as area of woodland in adjacent cells weighted by the weights matrix
woodland_connectivity <- focal(woodland_area, 
                               w = weights, 
                               fun = "sum", 
                               na.policy = "omit",
                               na.rm = TRUE)

broadleaf_connectivity <- focal(broadleaf_area, 
                               w = weights, 
                               fun = "sum", 
                               na.policy = "omit",
                               na.rm = TRUE)

conifer_connectivity <- focal(conifer_area, 
                               w = weights, 
                               fun = "sum", 
                               na.policy = "omit",
                               na.rm = TRUE)

writeRaster(woodland_connectivity, file = "woodland_connectivity_sigma3.tif")
writeRaster(broadleaf_connectivity, file = "broadleaf_connectivity_sigma3.tif")
writeRaster(conifer_connectivity, file = "conifer_connectivity_sigma3.tif")


# Recreation ####
setwd("/data/notebooks/rstudio-phytophthora/spatial_data/Recreation maps_EMBARGO")
recreation_weekly <- rast("Weekly_250m_NonUrban_PathsOnly.tif")
recreation_yearly <- rast("Yearly_250m_NonUrban_PathsOnly.tif")

recreation_weekly_1km <- terra::aggregate(recreation_weekly,
                                          fact = 4,
                                          fun = sum)
recreation_yearly_1km <- terra::aggregate(recreation_yearly,
                                          fact = 4,
                                          fun = sum)
recreation_weekly_1km_crs <- project(recreation_weekly_1km, GB_grid)
recreation_yearly_1km_crs <- project(recreation_yearly_1km, GB_grid)

recreation_weekly_1km_crs[is.na(recreation_weekly_1km_crs)] <- 0
recreation_yearly_1km_crs[is.na(recreation_yearly_1km_crs)] <- 0

recreation_weekly_1km_crs_mask <- mask(recreation_weekly_1km_crs, GB_grid, maskvalues = NA, updatevalue = NA)
recreation_yearly_1km_crs_mask <- mask(recreation_yearly_1km_crs, GB_grid, maskvalues = NA, updatevalue = NA)

setwd("/data/notebooks/rstudio-phytophthora/spatial_data/processed1km")
writeRaster(recreation_weekly_1km_crs_mask, file = "recreation_weekly_1km.tif")
writeRaster(recreation_yearly_1km_crs_mask, file = "recreation_yearly_1km.tif")
            
