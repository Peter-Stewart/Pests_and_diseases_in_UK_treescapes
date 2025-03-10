# Packages ####
library(PointedSDMs)
library(terra)
library(ggplot2)
library(dplyr)
library(viridis)
library(cowplot)

# Load all predictions
species <- c("Acer pseudoplatanus",
             "Betula pendula",
             "Fagus sylvatica",
             "Fraxinus excelsior",
             "Picea abies",
             "Picea sitchensis",
             "Pinus sylvestris",
             "Quercus robur",
             "Sorbus aucuparia")

species_ord <- c("Acer pseudoplatanus",
                 "Betula pendula",
                 "Fagus sylvatica",
                 "Fraxinus excelsior",
                 "Quercus robur",
                 "Sorbus aucuparia",
                 "Picea abies",
                 "Picea sitchensis",
                 "Pinus sylvestris") 

setwd("D:/INLA_outputs/mesh1")
file_list <- list.files(pattern = "predictions.Rdata")
predictions_list1 <- list()
for(i in 1:length(file_list)){
  predictions_list1[[i]] <- get(load(file_list[i]))
}
names(predictions_list1) <- species

setwd("D:/INLA_outputs/mesh2")
file_list <- list.files(pattern = "predictions.Rdata")
predictions_list2 <- list()
for(i in 1:length(file_list)){
  predictions_list2[[i]] <- get(load(file_list[i]))
}
names(predictions_list2) <- species

setwd("D:/INLA_outputs/mesh3")
file_list <- list.files(pattern = "predictions.Rdata")
predictions_list3 <- list()
for(i in 1:length(file_list)){
  predictions_list3[[i]] <- get(load(file_list[i]))
}
names(predictions_list3) <- species

setwd("D:/INLA_outputs/mesh4")
file_list <- list.files(pattern = "predictions.Rdata")
predictions_list4 <- list()
for(i in 1:length(file_list)){
  predictions_list4[[i]] <- get(load(file_list[i]))
}
names(predictions_list4) <- species

setwd("D:/INLA_outputs/mesh5")
file_list <- list.files(pattern = "predictions.Rdata")
predictions_list5 <- list()
for(i in 1:length(file_list)){
  predictions_list5[[i]] <- get(load(file_list[i]))
}
names(predictions_list5) <- species


# Calculate correlation between predictions and 5km mesh predictions
out_table <- matrix(NA, nrow = length(species), ncol = 5)
out_table[,1] <- species_ord
for(sp in 1:length(species_ord)){
  out_table[sp,2] <- round(cor(predictions_list5[[species_ord[sp]]]$mean, predictions_list3[[species_ord[sp]]]$mean), 3)
  out_table[sp,3] <- round(cor(predictions_list5[[species_ord[sp]]]$mean, predictions_list2[[species_ord[sp]]]$mean), 3)
  out_table[sp,4] <- round(cor(predictions_list5[[species_ord[sp]]]$mean, predictions_list1[[species_ord[sp]]]$mean), 3)
  out_table[sp,5] <- round(cor(predictions_list5[[species_ord[sp]]]$mean, predictions_list4[[species_ord[sp]]]$mean), 3)
  
}
