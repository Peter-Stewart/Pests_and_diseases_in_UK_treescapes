# Packages ####
library(ggdag)
library(dagitty)
library(ggplot2)
library(tidyr)
library(cowplot)

# Function to calculate minimum adjustment sets for multiple DAGs
adjustmentSets_multi <- function(patt, effect = "total"){
  dag_names <- ls(pattern = patt, envir = .GlobalEnv)
  adj_sets_mat <- matrix(NA, nrow = length(dag_names), ncol = 2)
  for(i in 1:length(dag_names)){
    adj_sets_mat[i,1] <- dag_names[i]
    adj_sets_mat[i,2] <- paste(adjustmentSets(get(dag_names[i]), effect = effect), collapse = "; ")
  }
  adj_sets_df <- as.data.frame(adj_sets_mat)
  
  colnames(adj_sets_df) <- c("DAG", "adjustment_sets")
  
  return(adj_sets_df)
}

# Function to generate bru formulae automatically from minimal adjustment sets
create_covariate_list <- function(adj_sets, crosswalk, form_front, var_focal, form_back){
  # Get unique adjustment sets
  adj_sets_u <- unique(adj_sets$adjustment_sets)
  print(paste("Found",length(adj_sets_u),"unique adjustment sets."))
  
  # Create empty list to store formulae
  form_list <- list()
  
  for(set in 1:length(adj_sets_u)){
    # Create empty formula string
    form_str <- as.character(NULL)
    # Add variables in the adjustment set to the formula string
    for(i in 1:length(crosswalk$DAG_var_names)){
      if(grepl(crosswalk$DAG_var_names[i], adj_sets_u[set], fixed = TRUE) == TRUE){
        form_str <- paste0(form_str, " + ", crosswalk$rast_names[i],"(", crosswalk$rast_names[i],', model = "linear", mean.linear = 0, prec.linear = 1)')
      }
    }
    # Add front and back onto formula string
    form_front <- form_front
    form_focal <- paste0(var_focal, "(", var_focal, ', model = "linear", mean.linear = 0, prec.linear = 1)')
    form_back <- form_back
    
    form_full <- paste(form_front, form_focal, form_str, form_back)
    form_full <- as.formula(form_full)
    
    # Append to the formula list
    form_list[[set]] <- form_full
  }
  
  # Return the formula list
  return(form_list)
}

# Function to create run_multi_ISDM covariate list from adjustment sets
create_covariate_list <- function(adj_sets, crosswalk, var_focal){
  # Get unique adjustment sets
  adj_sets_u <- unique(adj_sets$adjustment_sets)
  print(paste("Found",length(adj_sets_u),"unique adjustment sets."))
  
  # Create empty list to store formulae
  form_list <- list()
  
  for(set in 1:length(adj_sets_u)){
    # Create empty formula string
    form_str <- as.character(NULL)
    # Add variables in the adjustment set to the formula string
    for(i in 1:length(crosswalk$DAG_var_names)){
      if(grepl(crosswalk$DAG_var_names[i], adj_sets_u[set], fixed = TRUE) == TRUE){
        if(length(form_str) == 0){
          form_str <- paste0(crosswalk$rast_names[i])
        }else{
          form_str <- paste0(form_str, ",", crosswalk$rast_names[i])
          
        }
      }
    }
    
    if(length(form_str) == 0){
      form_full <- list(var_focal)
    }else{
      form_full <- paste0(form_str,",", var_focal)
      form_full <- strsplit(form_full, split = ",")
    }
    
    # Append to the formula list
    form_list[[set]] <- form_full[[1]]
  }
  
  # Return the formula list
  return(form_list)
}


# Dataframe of node coordinates for DAGs ####
DAG_coords1 <- read.csv("C:/temp/pests_analysis/Model_flists/DAG_node_coordinates.csv")
DAG_coords1$y <- -1*DAG_coords1$y
DAG_coords1$y[DAG_coords1$name == "d_BCP"] <- 0.187

# DAG vs. variable names crosswalk
var_names_crosswalk <- read.csv("C:/temp/pests_analysis/Model_flists/var_names_crosswalk.csv")
var_names_crosswalk <- var_names_crosswalk[!is.na(var_names_crosswalk$rast_names),]

# Recreation ####
R_d01a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                 p_e ~ p_i, # Due to propagule pressure 
                 p_i ~ R,  # Recreation
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)


R_d02a <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d02b <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, p_e ~ H, # Human pop. density
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d03a <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d03a <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d04a <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 H ~ E, R ~ E, # Elevation
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d04b <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d04c <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)


R_d05a <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, # VPD
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d05b <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d06a <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, # Dist. to BCP
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d06b <- dagify(y ~ p_i + p_e,
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, # Dist. to BCP
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d06c <- dagify(y ~ p_i + p_e,
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d07a <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, # Woodland area
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d07b <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, # Woodland area
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d07c <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, # Woodland area
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)


R_d07c <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, # Woodland area
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d07c <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d08a <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, # Distance to park/garden
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d08b <- dagify(y ~ p_i + p_e,
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d09a <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, # Ancient woodland cover
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d09b <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d10a <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                 C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, # Conifer area
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d10b <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                 C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, # Conifer area
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d10c <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                 C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d11a <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                 C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                 B_a ~ E, W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d11b <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                 C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                 B_a ~ E, W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Ur, # Broadleaf area
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d11c <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                 C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                 B_a ~ E, W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Ur, R ~ B_a, # Broadleaf area
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d12a <- dagify(y ~ p_i + p_e,
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                 C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                 B_a ~ E, W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Ur, R ~ B_a, # Broadleaf area
                 W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d12b <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                 C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                 B_a ~ E, W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Ur, R ~ B_a, # Broadleaf area
                 W_c ~ W_a, p_i ~ W_c, R ~ W_c, # Woodland connectivity
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d13a <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                 C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                 B_a ~ E, W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Ur, R ~ B_a, # Broadleaf area
                 W_c ~ W_a, p_i ~ W_c, R ~ W_c, # Woodland connectivity
                 p_i ~ Af, W_a ~ Af, C_a ~ Af, B_a ~ Af, # Afforestation
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

R_d14a <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, 
                 p_i ~ R, # Recreation
                 R ~ H, p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                 H ~ E, R ~ E, Ur ~ E, # Elevation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                 C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                 B_a ~ E, W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Ur, R ~ B_a, # Broadleaf area
                 W_c ~ W_a, p_i ~ W_c, R ~ W_c, # Woodland connectivity
                 p_i ~ Af, W_a ~ Af, C_a ~ Af, B_a ~ Af, # Afforestation
                 p_i ~ Df, W_a ~ Df, C_a ~ Df, B_a ~ Df, A ~ Df, # Deforestation
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "R",
                 coords = DAG_coords1)

# Calculate minimum adjustment sets 
R_sets1 <- adjustmentSets_multi("R_d")
R_sets1$DAG_num <- 1:nrow(R_sets1)
R_sets1_sep <- separate_longer_delim(data = R_sets1, 
                                     cols = adjustment_sets, 
                                     delim = "; ")


# Create list of inlabru formulae from minimum adjustment sets
R_flist <- create_covariate_list(adj_sets = R_sets1_sep,
                                 crosswalk = var_names_crosswalk,
                                 var_focal = "recreation_weekly_1km")

setwd("C:/temp/pests_analysis/ISDM_flists")
save(R_flist, file = "R_flist.Rdata")

# Create tidy format DAGs
dag_names <- ls(pattern = "R_d", envir = .GlobalEnv)
for(i in 1:length(dag_names)){
  d <- get(dag_names[i])
  
  d_tidy <- d %>% tidy_dagitty()
  
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "p_i" | d_tidy$data$name == "p_e", "latent", "observed")
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "y", "outcome", d_tidy$data$colour)
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "R", "exposure", d_tidy$data$colour)
  
  assign(paste0(dag_names[i],"_tidy"), d_tidy)
}

# Plot DAGs
p_list <- list()
dag_names <- ls(pattern = "^R_d.+_tidy$", envir = .GlobalEnv)
for(i in 1:length(dag_names)){
  d_temp <- get(dag_names[i])
  p_temp <- d_temp %>% 
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges() +
    geom_dag_point(aes(colour = colour)) +
    geom_dag_text(size = 3.5) +
    theme_dag() +
    theme(legend.position="none") +
    scale_color_manual(values = c("#cc2055", "#8d8d8d", "black", "#14a0f0"))
  p_list[[i]] <- p_temp
}


# Save all DAGs to .tiff
setwd("F:/INLA_outputs/DAGs")
for(i in seq(1, length(p_list), by = 8)){
  grid <- plot_grid(plotlist = p_list[i:(i+7)],
                    ncol = 4,
                    labels = i:(i+7))
  
  ggsave(filename = paste0("R_DAGs_",i,"_to_",i+7,".tiff"),
         plot = grid,
         dpi = "print",
         width = 16, 
         height = 7.8, 
         units = "cm",
         bg = "white",
         scale = 2.5)
}


# Urbanisation ####
Ur_d01a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d02a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d03a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, # Recreation
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)


Ur_d03b <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)


Ur_d04a <- dagify(y ~ p_i + p_e, 
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, # Elevation
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d04b <- dagify(y ~ p_i + p_e, 
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)


Ur_d05a <- dagify(y ~ p_i + p_e, 
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, # VPD
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d05b <- dagify(y ~ p_i + p_e, 
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d06a <- dagify(y ~ p_i + p_e, 
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                  p_i ~ d_BCP, # Dist. to BCP
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d06b <- dagify(y ~ p_i + p_e,
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                  p_i ~ d_BCP, d_BCP ~ Ur, # Dist. to BCP
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d06c <- dagify(y ~ p_i + p_e,
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                  p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d07a <- dagify(y ~ p_i + p_e, 
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                  p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                  p_e ~ W_a, W_a ~ E, # Woodland area
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d07b <- dagify(y ~ p_i + p_e, 
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                  p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                  p_e ~ W_a, W_a ~ E, p_i ~ W_a, # Woodland area
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d07c <- dagify(y ~ p_i + p_e, 
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                  p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                  p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, # Woodland area
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)


Ur_d07c <- dagify(y ~ p_i + p_e, 
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                  p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                  p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, # Woodland area
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d07c <- dagify(y ~ p_i + p_e, 
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                  p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                  p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d08a <- dagify(y ~ p_i + p_e, 
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                  p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                  p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                  R ~ d_G, p_i ~ d_G, # Distance to park/garden
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d08b <- dagify(y ~ p_i + p_e,
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                  p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                  p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                  R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d09a <- dagify(y ~ p_i + p_e, 
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                  p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                  p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                  R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                  R ~ A, W_a ~ A, # Ancient woodland cover
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d09b <- dagify(y ~ p_i + p_e, 
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                  p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                  p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                  R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                  R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d10a <- dagify(y ~ p_i + p_e, 
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                  p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                  p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                  R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                  R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                  C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, # Conifer area
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d10b <- dagify(y ~ p_i + p_e, 
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                  p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                  p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                  R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                  R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                  C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, # Conifer area
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d10c <- dagify(y ~ p_i + p_e, 
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                  p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                  p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                  R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                  R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                  C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d11a <- dagify(y ~ p_i + p_e, 
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                  p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                  p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                  R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                  R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                  C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                  B_a ~ E, W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d11b <- dagify(y ~ p_i + p_e, 
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                  p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                  p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                  R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                  R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                  C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                  B_a ~ E, W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Ur, # Broadleaf area
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d11c <- dagify(y ~ p_i + p_e, 
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                  p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                  p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                  R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                  R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                  C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                  B_a ~ E, W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Ur, R ~ B_a, # Broadleaf area
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d12a <- dagify(y ~ p_i + p_e,
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                  p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                  p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                  R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                  R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                  C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                  B_a ~ E, W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Ur, R ~ B_a, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d12b <- dagify(y ~ p_i + p_e, 
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                  p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                  p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                  R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                  R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                  C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                  B_a ~ E, W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Ur, R ~ B_a, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, R ~ W_c, # Woodland connectivity
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d13a <- dagify(y ~ p_i + p_e, 
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                  p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                  p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                  R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                  R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                  C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                  B_a ~ E, W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Ur, R ~ B_a, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, R ~ W_c, # Woodland connectivity
                  p_i ~ Af, W_a ~ Af, C_a ~ Af, B_a ~ Af, # Afforestation
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

Ur_d14a <- dagify(y ~ p_i + p_e, 
                  p_e ~ p_i, 
                  p_i ~ Ur, p_e ~ Ur, # Urbanisation
                  Ur ~ H, p_i ~ H, # Human pop. density
                  R ~ Ur, p_i ~ R, R ~ H, # Recreation
                  H ~ E, R ~ E, Ur ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                  p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                  p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                  R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                  R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                  C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                  B_a ~ E, W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Ur, R ~ B_a, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, R ~ W_c, # Woodland connectivity
                  p_i ~ Af, W_a ~ Af, C_a ~ Af, B_a ~ Af, # Afforestation
                  p_i ~ Df, W_a ~ Df, C_a ~ Df, B_a ~ Df, A ~ Df, # Deforestation
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Ur",
                  coords = DAG_coords1)

# Calculate minimum adjustment sets 
Ur_sets1 <- adjustmentSets_multi("Ur_d")
Ur_sets1$DAG_num <- 1:nrow(Ur_sets1)
Ur_sets1_sep <- separate_longer_delim(data = Ur_sets1, 
                                      cols = adjustment_sets, 
                                      delim = "; ")


# Create list of inlabru formulae from minimum adjustment sets
Ur_flist <- create_covariate_list(adj_sets = Ur_sets1_sep, 
                                  crosswalk = var_names_crosswalk,
                                  var_focal = "urb_suburb_areaHa")


setwd("C:/temp/pests_analysis/ISDM_flists")
save(Ur_flist, file = "Ur_flist.Rdata")

# Create tidy format DAGs
dag_names <- ls(pattern = "Ur_d", envir = .GlobalEnv)
for(i in 1:length(dag_names)){
  d <- get(dag_names[i])
  
  d_tidy <- d %>% tidy_dagitty()
  
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "p_i" | d_tidy$data$name == "p_e", "latent", "observed")
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "y", "outcome", d_tidy$data$colour)
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "Ur", "exposure", d_tidy$data$colour)
  
  assign(paste0(dag_names[i],"_tidy"), d_tidy)
}

# Plot DAGs
p_list <- list()
dag_names <- ls(pattern = "^Ur_d.+_tidy$", envir = .GlobalEnv)
for(i in 1:length(dag_names)){
  d_temp <- get(dag_names[i])
  p_temp <- d_temp %>% 
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges() +
    geom_dag_point(aes(colour = colour)) +
    geom_dag_text(size = 3.5) +
    theme_dag() +
    theme(legend.position="none") +
    scale_color_manual(values = c("#cc2055", "#8d8d8d", "black", "#14a0f0"))
  p_list[[i]] <- p_temp
}


# Save all DAGs to .tiff
setwd("F:/INLA_outputs/DAGs")
for(i in seq(1, length(p_list), by = 8)){
  grid <- plot_grid(plotlist = p_list[i:(i+7)],
                    ncol = 4,
                    labels = i:(i+7))
  
  ggsave(filename = paste0("Ur_DAGs_",i,"_to_",i+7,".tiff"),
         plot = grid,
         dpi = "print",
         width = 16, 
         height = 7.8, 
         units = "cm",
         bg = "white",
         scale = 2.5)
}

# Human pop. density ####
H_d01a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, p_e ~ H, # Human pop. density
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d02a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d03a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d04a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d04b <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, # Recreation
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d05a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, # VPD
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d05b <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d06a <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, # Dist. to BCP
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d06b <- dagify(y ~ p_i + p_e,
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, # Dist. to BCP
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d06c <- dagify(y ~ p_i + p_e,
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d07a <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, # Woodland area
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d07b <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, # Woodland area
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d07c <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, # Woodland area
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)


H_d07c <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, # Woodland area
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d07c <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d08a <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, # Distance to park/garden
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d08b <- dagify(y ~ p_i + p_e,
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d09a <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, # Ancient woodland cover
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d09b <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d10a <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                 C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, # Conifer area
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d10b <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                 C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, # Conifer area
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d10c <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                 C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d11a <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                 C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                 B_a ~ E, W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d11b <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                 C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                 B_a ~ E, W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Ur, # Broadleaf area
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d11c <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                 C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                 B_a ~ E, W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Ur, R ~ B_a, # Broadleaf area
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d12a <- dagify(y ~ p_i + p_e,
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                 C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                 B_a ~ E, W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Ur, R ~ B_a, # Broadleaf area
                 W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d12b <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                 C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                 B_a ~ E, W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Ur, R ~ B_a, # Broadleaf area
                 W_c ~ W_a, p_i ~ W_c, R ~ W_c, # Woodland connectivity
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d13a <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                 C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                 B_a ~ E, W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Ur, R ~ B_a, # Broadleaf area
                 W_c ~ W_a, p_i ~ W_c, R ~ W_c, # Woodland connectivity
                 p_i ~ Af, W_a ~ Af, C_a ~ Af, B_a ~ Af, # Afforestation
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

H_d14a <- dagify(y ~ p_i + p_e, 
                 p_e ~ p_i, # Due to propagule pressure
                 p_i ~ H, # Human pop. density
                 Ur ~ H, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                 Ur ~ E, H ~ E, # Elevation
                 R ~ E, R ~ H, R ~ Ur, # Recreation
                 VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                 p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                 p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                 R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                 R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                 C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                 B_a ~ E, W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Ur, R ~ B_a, # Broadleaf area
                 W_c ~ W_a, p_i ~ W_c, R ~ W_c, # Woodland connectivity
                 p_i ~ Af, W_a ~ Af, C_a ~ Af, B_a ~ Af, # Afforestation
                 p_i ~ Df, W_a ~ Df, C_a ~ Df, B_a ~ Df, A ~ Df, # Deforestation
                 latent = c("p_i", "p_e"), 
                 outcome = "y", 
                 exposure = "H",
                 coords = DAG_coords1)

# Calculate minimum adjustment sets 
H_sets1 <- adjustmentSets_multi("H_d")
H_sets1$DAG_num <- 1:nrow(H_sets1)
H_sets1_sep <- separate_longer_delim(data = H_sets1, 
                                     cols = adjustment_sets, 
                                     delim = "; ")


# Create list of inlabru formulae from minimum adjustment sets
H_flist <- create_covariate_list(adj_sets = H_sets1_sep, 
                                 crosswalk = var_names_crosswalk,
                                 var_focal = "humanPop"
)


setwd("C:/temp/pests_analysis/ISDM_flists")
save(H_flist, file = "H_flist.Rdata")

# Create tidy format DAGs
dag_names <- ls(pattern = "H_d", envir = .GlobalEnv)
for(i in 1:length(dag_names)){
  d <- get(dag_names[i])
  
  d_tidy <- d %>% tidy_dagitty()
  
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "p_i" | d_tidy$data$name == "p_e", "latent", "observed")
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "y", "outcome", d_tidy$data$colour)
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "H", "exposure", d_tidy$data$colour)
  
  assign(paste0(dag_names[i],"_tidy"), d_tidy)
}

# Plot DAGs
p_list <- list()
dag_names <- ls(pattern = "^H_d.+_tidy$", envir = .GlobalEnv)
for(i in 1:length(dag_names)){
  d_temp <- get(dag_names[i])
  p_temp <- d_temp %>% 
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges() +
    geom_dag_point(aes(colour = colour)) +
    geom_dag_text(size = 3.5) +
    theme_dag() +
    theme(legend.position="none") +
    scale_color_manual(values = c("#cc2055", "#8d8d8d", "black", "#14a0f0"))
  p_list[[i]] <- p_temp
}


# Save all DAGs to .tiff
setwd("F:/INLA_outputs/DAGs")
for(i in seq(1, length(p_list), by = 8)){
  grid <- plot_grid(plotlist = p_list[i:(i+7)],
                    ncol = 4,
                    labels = i:(i+7))
  
  ggsave(filename = paste0("H_DAGs_",i,"_to_",i+7,".tiff"),
         plot = grid,
         dpi = "print",
         width = 16, 
         height = 7.8, 
         units = "cm",
         bg = "white",
         scale = 2.5)
}

# Afforestation ####
Af_d01a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Af, # Afforestation
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Af",
                  coords = DAG_coords1)

Af_d02a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Af, # Afforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Af, # Woodland area
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Af",
                  coords = DAG_coords1)

Af_d03a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Af, # Afforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Af, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Af, # Conifer area
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Af",
                  coords = DAG_coords1)

Af_d04a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Af, # Afforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Af, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Af, # Conifer area
                  W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Af, # Broadleaf area
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Af",
                  coords = DAG_coords1)

Af_d05a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Af, # Afforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Af, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Af, # Conifer area
                  W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Af, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Af",
                  coords = DAG_coords1)

Af_d06a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Af, # Afforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Af, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Af, # Conifer area
                  W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Af, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                  W_a ~ A, # Ancient woodland cover
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Af",
                  coords = DAG_coords1)

Af_d07a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Af, # Afforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Af, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Af, # Conifer area
                  W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Af, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                  W_a ~ A, # Ancient woodland cover
                  p_i ~ Df, W_a ~ Df, C_a ~ Df, B_a ~ Df, A ~ Df, # Deforestation
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Af",
                  coords = DAG_coords1)

Af_d08a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Af, # Afforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Af, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Af, # Conifer area
                  W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Af, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                  W_a ~ A, # Ancient woodland cover
                  p_i ~ Df, W_a ~ Df, C_a ~ Df, B_a ~ Df, A ~ Df, # Deforestation
                  W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Af",
                  coords = DAG_coords1)

Af_d09a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Af, # Afforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Af, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Af, # Conifer area
                  W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Af, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                  W_a ~ A, # Ancient woodland cover
                  p_i ~ Df, W_a ~ Df, C_a ~ Df, B_a ~ Df, A ~ Df, # Deforestation
                  W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, # VPD
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Af",
                  coords = DAG_coords1)

Af_d10a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Af, # Afforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Af, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Af, # Conifer area
                  W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Af, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                  W_a ~ A, # Ancient woodland cover
                  p_i ~ Df, W_a ~ Df, C_a ~ Df, B_a ~ Df, A ~ Df, # Deforestation
                  W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, # VPD
                  p_i ~ Ur, p_e ~ Ur, VPD ~ Ur, W_a ~ Ur, C_a ~ Ur, B_a ~ Ur, A ~ Ur, Ur ~ E, # Urbanisation
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Af",
                  coords = DAG_coords1)

Af_d11a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Af, # Afforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Af, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Af, # Conifer area
                  W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Af, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                  W_a ~ A, # Ancient woodland cover
                  p_i ~ Df, W_a ~ Df, C_a ~ Df, B_a ~ Df, A ~ Df, # Deforestation
                  W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, # VPD
                  p_i ~ Ur, p_e ~ Ur, VPD ~ Ur, W_a ~ Ur, C_a ~ Ur, B_a ~ Ur, A ~ Ur, Ur ~ E, # Urbanisation
                  p_i ~ H, Ur ~ H, H ~ E, # Human pop. density
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Af",
                  coords = DAG_coords1)

Af_d12a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Af, # Afforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Af, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Af, # Conifer area
                  W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Af, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                  W_a ~ A, # Ancient woodland cover
                  p_i ~ Df, W_a ~ Df, C_a ~ Df, B_a ~ Df, A ~ Df, # Deforestation
                  W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, # VPD
                  p_i ~ Ur, p_e ~ Ur, VPD ~ Ur, W_a ~ Ur, C_a ~ Ur, B_a ~ Ur, A ~ Ur, Ur ~ E, # Urbanisation
                  p_i ~ H, Ur ~ H, H ~ E, # Human pop. density
                  p_i ~ R, R ~ Ur, R ~ W_a, R ~ B_a, R ~ C_a, R ~ A, R ~ H, # Recreation
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Af",
                  coords = DAG_coords1)

Af_d13a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Af, # Afforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Af, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Af, # Conifer area
                  W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Af, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                  W_a ~ A, # Ancient woodland cover
                  p_i ~ Df, W_a ~ Df, C_a ~ Df, B_a ~ Df, A ~ Df, # Deforestation
                  W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, # VPD
                  p_i ~ Ur, p_e ~ Ur, VPD ~ Ur, W_a ~ Ur, C_a ~ Ur, B_a ~ Ur, A ~ Ur, Ur ~ E, # Urbanisation
                  p_i ~ H, Ur ~ H, H ~ E, # Human pop. density
                  p_i ~ R, R ~ Ur, R ~ W_a, R ~ B_a, R ~ C_a, R ~ A, R ~ H, # Recreation
                  R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Af",
                  coords = DAG_coords1)

Af_d14a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Af, # Afforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Af, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Af, # Conifer area
                  W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Af, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                  W_a ~ A, # Ancient woodland cover
                  p_i ~ Df, W_a ~ Df, C_a ~ Df, B_a ~ Df, A ~ Df, # Deforestation
                  W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, # VPD
                  p_i ~ Ur, p_e ~ Ur, VPD ~ Ur, W_a ~ Ur, C_a ~ Ur, B_a ~ Ur, A ~ Ur, Ur ~ E, # Urbanisation
                  p_i ~ H, Ur ~ H, H ~ E, # Human pop. density
                  p_i ~ R, R ~ Ur, R ~ W_a, R ~ B_a, R ~ C_a, R ~ A, R ~ H, # Recreation
                  R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                  p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Af",
                  coords = DAG_coords1)

# Calculate minimum adjustment sets 
Af_sets1 <- adjustmentSets_multi("Af_d")
Af_sets1$DAG_num <- 1:nrow(Af_sets1)
Af_sets1_sep <- separate_longer_delim(data = Af_sets1, 
                                      cols = adjustment_sets, 
                                      delim = "; ")


# Create list of inlabru formulae from minimum adjustment sets
Af_flist <- create_covariate_list(adj_sets = Af_sets1_sep, 
                                  crosswalk = var_names_crosswalk,
                                  var_focal = "LCC_afforestation_1km"
)


setwd("C:/temp/pests_analysis/ISDM_flists")
save(Af_flist, file = "Af_flist.Rdata")

# Create tidy format DAGs
dag_names <- ls(pattern = "Af_d", envir = .GlobalEnv)
for(i in 1:length(dag_names)){
  d <- get(dag_names[i])
  
  d_tidy <- d %>% tidy_dagitty()
  
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "p_i" | d_tidy$data$name == "p_e", "latent", "observed")
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "y", "outcome", d_tidy$data$colour)
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "Af", "exposure", d_tidy$data$colour)
  
  assign(paste0(dag_names[i],"_tidy"), d_tidy)
}

# Plot DAGs
p_list <- list()
dag_names <- ls(pattern = "^Af_d.+_tidy$", envir = .GlobalEnv)
for(i in 1:length(dag_names)){
  d_temp <- get(dag_names[i])
  p_temp <- d_temp %>% 
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges() +
    geom_dag_point(aes(colour = colour)) +
    geom_dag_text(size = 3.5) +
    theme_dag() +
    theme(legend.position="none") +
    scale_color_manual(values = c("#cc2055", "#8d8d8d", "black", "#14a0f0"))
  p_list[[i]] <- p_temp
}


# Save all DAGs to .tiff
setwd("F:/INLA_outputs/DAGs")
for(i in seq(1, length(p_list), by = 8)){
  grid <- plot_grid(plotlist = p_list[i:(i+7)],
                    ncol = 4,
                    labels = i:(i+7))
  
  ggsave(filename = paste0("Af_DAGs_",i,"_to_",i+7,".tiff"),
         plot = grid,
         dpi = "print",
         width = 16, 
         height = 7.8, 
         units = "cm",
         bg = "white",
         scale = 2.5)
}

# Deforestation ####
Df_d01a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Df, # Deforestation
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Df",
                  coords = DAG_coords1)

Df_d02a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Df, # Deforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Df, # Woodland area
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Df",
                  coords = DAG_coords1)

Df_d03a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Df, # Deforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Df, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Df, # Conifer area
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Df",
                  coords = DAG_coords1)

Df_d04a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Df, # Deforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Df, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Df, # Conifer area
                  W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Df, # Broadleaf area
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Df",
                  coords = DAG_coords1)

Df_d05a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Df, # Deforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Df, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Df, # Conifer area
                  W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Df, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Df",
                  coords = DAG_coords1)

Df_d06a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Df, # Deforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Df, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Df, # Conifer area
                  W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Df, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                  W_a ~ A, A ~ Df, # Ancient woodland cover
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Df",
                  coords = DAG_coords1)

Df_d07a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Df, # Deforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Df, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Df, # Conifer area
                  W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Df, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                  W_a ~ A, A ~ Df, # Ancient woodland cover
                  p_i ~ Af, W_a ~ Af, C_a ~ Af, B_a ~ Af, A ~ Af, # Afforestation
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Df",
                  coords = DAG_coords1)

Df_d08a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Df, # Deforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Df, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Df, # Conifer area
                  W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Df, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                  W_a ~ A, A ~ Df, # Ancient woodland cover
                  p_i ~ Af, W_a ~ Af, C_a ~ Af, B_a ~ Af, A ~ Af, # Afforestation
                  W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Df",
                  coords = DAG_coords1)

Df_d09a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Df, # Deforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Df, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Df, # Conifer area
                  W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Df, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                  W_a ~ A, A ~ Df, # Ancient woodland cover
                  p_i ~ Af, W_a ~ Af, C_a ~ Af, B_a ~ Af, A ~ Af, # Afforestation
                  W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, # VPD
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Df",
                  coords = DAG_coords1)

Df_d10a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Df, # Deforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Df, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Df, # Conifer area
                  W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Df, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                  W_a ~ A, A ~ Df, # Ancient woodland cover
                  p_i ~ Af, W_a ~ Af, C_a ~ Af, B_a ~ Af, A ~ Af, # Afforestation
                  W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, # VPD
                  p_i ~ Ur, p_e ~ Ur, VPD ~ Ur, W_a ~ Ur, C_a ~ Ur, B_a ~ Ur, A ~ Ur, Ur ~ E, # Urbanisation
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Df",
                  coords = DAG_coords1)

Df_d11a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Df, # Deforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Df, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Df, # Conifer area
                  W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Df, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                  W_a ~ A, A ~ Df, # Ancient woodland cover
                  p_i ~ Af, W_a ~ Af, C_a ~ Af, B_a ~ Af, A ~ Af, # Afforestation
                  W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, # VPD
                  p_i ~ Ur, p_e ~ Ur, VPD ~ Ur, W_a ~ Ur, C_a ~ Ur, B_a ~ Ur, A ~ Ur, Ur ~ E, # Urbanisation
                  p_i ~ H, Ur ~ H, H ~ E, # Human pop. density
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Df",
                  coords = DAG_coords1)

Df_d12a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Df, # Deforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Df, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Df, # Conifer area
                  W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Df, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                  W_a ~ A, A ~ Df, # Ancient woodland cover
                  p_i ~ Af, W_a ~ Af, C_a ~ Af, B_a ~ Af, A ~ Af, # Afforestation
                  W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, # VPD
                  p_i ~ Ur, p_e ~ Ur, VPD ~ Ur, W_a ~ Ur, C_a ~ Ur, B_a ~ Ur, A ~ Ur, Ur ~ E, # Urbanisation
                  p_i ~ H, Ur ~ H, H ~ E, # Human pop. density
                  p_i ~ R, R ~ Ur, R ~ W_a, R ~ B_a, R ~ C_a, R ~ A, R ~ H, # Recreation
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Df",
                  coords = DAG_coords1)

Df_d13a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Df, # Deforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Df, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Df, # Conifer area
                  W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Df, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                  W_a ~ A, A ~ Df, # Ancient woodland cover
                  p_i ~ Af, W_a ~ Af, C_a ~ Af, B_a ~ Af, A ~ Af, # Afforestation
                  W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, # VPD
                  p_i ~ Ur, p_e ~ Ur, VPD ~ Ur, W_a ~ Ur, C_a ~ Ur, B_a ~ Ur, A ~ Ur, Ur ~ E, # Urbanisation
                  p_i ~ H, Ur ~ H, H ~ E, # Human pop. density
                  p_i ~ R, R ~ Ur, R ~ W_a, R ~ B_a, R ~ C_a, R ~ A, R ~ H, # Recreation
                  R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Df",
                  coords = DAG_coords1)

Df_d14a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                  p_e ~ p_i, # Due to propagule pressure
                  p_i ~ Df, # Deforestation
                  p_e ~ W_a, p_i ~ W_a, W_a ~ Df, # Woodland area
                  W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Df, # Conifer area
                  W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Df, # Broadleaf area
                  W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                  W_a ~ A, A ~ Df, # Ancient woodland cover
                  p_i ~ Af, W_a ~ Af, C_a ~ Af, B_a ~ Af, A ~ Af, # Afforestation
                  W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                  VPD ~ E, p_e ~ VPD, # VPD
                  p_i ~ Ur, p_e ~ Ur, VPD ~ Ur, W_a ~ Ur, C_a ~ Ur, B_a ~ Ur, A ~ Ur, Ur ~ E, # Urbanisation
                  p_i ~ H, Ur ~ H, H ~ E, # Human pop. density
                  p_i ~ R, R ~ Ur, R ~ W_a, R ~ B_a, R ~ C_a, R ~ A, R ~ H, # Recreation
                  R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                  p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                  latent = c("p_i", "p_e"), 
                  outcome = "y", 
                  exposure = "Df",
                  coords = DAG_coords1)

# Calculate minimum adjustment sets 
Df_sets1 <- adjustmentSets_multi("Df_d")
Df_sets1$DAG_num <- 1:nrow(Df_sets1)
Df_sets1_sep <- separate_longer_delim(data = Df_sets1, 
                                      cols = adjustment_sets, 
                                      delim = "; ")


# Create list of inlabru formulae from minimum adjustment sets
Df_flist <- create_covariate_list(adj_sets = Df_sets1_sep, 
                                  crosswalk = var_names_crosswalk,
                                  var_focal = "LCC_deforestation_1km"
)


setwd("C:/temp/pests_analysis/ISDM_flists")
save(Df_flist, file = "Df_flist.Rdata")

# Create tidy format DAGs
dag_names <- ls(pattern = "Df_d", envir = .GlobalEnv)
for(i in 1:length(dag_names)){
  d <- get(dag_names[i])
  
  d_tidy <- d %>% tidy_dagitty()
  
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "p_i" | d_tidy$data$name == "p_e", "latent", "observed")
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "y", "outcome", d_tidy$data$colour)
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "Df", "exposure", d_tidy$data$colour)
  
  assign(paste0(dag_names[i],"_tidy"), d_tidy)
}

# Plot DAGs
p_list <- list()
dag_names <- ls(pattern = "^Df_d.+_tidy$", envir = .GlobalEnv)
for(i in 1:length(dag_names)){
  d_temp <- get(dag_names[i])
  p_temp <- d_temp %>% 
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges() +
    geom_dag_point(aes(colour = colour)) +
    geom_dag_text(size = 3.5) +
    theme_dag() +
    theme(legend.position="none") +
    scale_color_manual(values = c("#cc2055", "#8d8d8d", "black", "#14a0f0"))
  p_list[[i]] <- p_temp
}


# Save all DAGs to .tiff
setwd("F:/INLA_outputs/DAGs")
for(i in seq(1, length(p_list), by = 8)){
  grid <- plot_grid(plotlist = p_list[i:(i+7)],
                    ncol = 4,
                    labels = i:(i+7))
  
  ggsave(filename = paste0("Df_DAGs_",i,"_to_",i+7,".tiff"),
         plot = grid,
         dpi = "print",
         width = 16, 
         height = 7.8, 
         units = "cm",
         bg = "white",
         scale = 2.5)
}


# Dist. to BCP ####
d_BCP_d01a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                     p_e ~ p_i, # Due to propagule pressure
                     p_i ~ d_BCP, # Dist. to BCP
                     latent = c("p_i", "p_e"), 
                     outcome = "y", 
                     exposure = "d_BCP",
                     coords = DAG_coords1)

d_BCP_d02a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                     p_e ~ p_i, # Due to propagule pressure
                     p_i ~ d_BCP, # Dist. to BCP
                     d_BCP ~ Ur, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                     latent = c("p_i", "p_e"), 
                     outcome = "y", 
                     exposure = "d_BCP",
                     coords = DAG_coords1)

d_BCP_d03a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                     p_e ~ p_i, # Due to propagule pressure
                     p_i ~ d_BCP, # Dist. to BCP
                     d_BCP ~ Ur, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                     p_i ~ H, Ur ~ H, # Human pop. density
                     latent = c("p_i", "p_e"), 
                     outcome = "y", 
                     exposure = "d_BCP",
                     coords = DAG_coords1)

d_BCP_d04a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                     p_e ~ p_i, # Due to propagule pressure
                     p_i ~ d_BCP, # Dist. to BCP
                     d_BCP ~ Ur, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                     p_i ~ H, Ur ~ H, # Human pop. density
                     d_BCP ~ E, Ur ~ E, H ~ E, # Elevation
                     latent = c("p_i", "p_e"), 
                     outcome = "y", 
                     exposure = "d_BCP",
                     coords = DAG_coords1)

d_BCP_d05a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                     p_e ~ p_i, # Due to propagule pressure
                     p_i ~ d_BCP, # Dist. to BCP
                     d_BCP ~ Ur, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                     p_i ~ H, Ur ~ H, # Human pop. density
                     d_BCP ~ E, Ur ~ E, H ~ E, # Elevation
                     VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                     latent = c("p_i", "p_e"), 
                     outcome = "y", 
                     exposure = "d_BCP",
                     coords = DAG_coords1)

d_BCP_d06a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                     p_e ~ p_i, # Due to propagule pressure
                     p_i ~ d_BCP, # Dist. to BCP
                     d_BCP ~ Ur, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                     p_i ~ H, Ur ~ H, # Human pop. density
                     d_BCP ~ E, Ur ~ E, H ~ E, # Elevation
                     VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                     p_i ~ R, R ~ Ur, R ~ H, # Recreation
                     latent = c("p_i", "p_e"), 
                     outcome = "y", 
                     exposure = "d_BCP",
                     coords = DAG_coords1)

d_BCP_d07a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                     p_e ~ p_i, # Due to propagule pressure
                     p_i ~ d_BCP, # Dist. to BCP
                     d_BCP ~ Ur, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                     p_i ~ H, Ur ~ H, # Human pop. density
                     d_BCP ~ E, Ur ~ E, H ~ E, # Elevation
                     VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                     p_i ~ R, R ~ Ur, R ~ H, # Recreation
                     R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                     latent = c("p_i", "p_e"), 
                     outcome = "y", 
                     exposure = "d_BCP",
                     coords = DAG_coords1)

d_BCP_d08a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                     p_e ~ p_i, # Due to propagule pressure
                     p_i ~ d_BCP, # Dist. to BCP
                     d_BCP ~ Ur, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                     p_i ~ H, Ur ~ H, # Human pop. density
                     d_BCP ~ E, Ur ~ E, H ~ E, # Elevation
                     VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                     p_i ~ R, R ~ Ur, R ~ H, # Recreation
                     R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                     p_e ~ W_a, p_i ~ W_a, R ~ W_a, W_a ~ Ur, W_a ~ E, # Woodland area
                     latent = c("p_i", "p_e"), 
                     outcome = "y", 
                     exposure = "d_BCP",
                     coords = DAG_coords1)

d_BCP_d09a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                     p_e ~ p_i, # Due to propagule pressure
                     p_i ~ d_BCP, # Dist. to BCP
                     d_BCP ~ Ur, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                     p_i ~ H, Ur ~ H, # Human pop. density
                     d_BCP ~ E, Ur ~ E, H ~ E, # Elevation
                     VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                     p_i ~ R, R ~ Ur, R ~ H, # Recreation
                     R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                     p_e ~ W_a, p_i ~ W_a, R ~ W_a, W_a ~ Ur, W_a ~ E, # Woodland area
                     W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, R ~ C_a, C_a ~ Ur, # Conifer area
                     latent = c("p_i", "p_e"), 
                     outcome = "y", 
                     exposure = "d_BCP",
                     coords = DAG_coords1)

d_BCP_d10a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                     p_e ~ p_i, # Due to propagule pressure
                     p_i ~ d_BCP, # Dist. to BCP
                     d_BCP ~ Ur, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                     p_i ~ H, Ur ~ H, # Human pop. density
                     d_BCP ~ E, Ur ~ E, H ~ E, # Elevation
                     VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                     p_i ~ R, R ~ Ur, R ~ H,# Recreation
                     R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                     p_e ~ W_a, p_i ~ W_a, R ~ W_a, W_a ~ Ur, W_a ~ E, # Woodland area
                     W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, R ~ C_a, C_a ~ Ur, # Conifer area
                     W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, R ~ B_a, B_a ~ Ur, # Broadleaf area
                     latent = c("p_i", "p_e"), 
                     outcome = "y", 
                     exposure = "d_BCP",
                     coords = DAG_coords1)

d_BCP_d11a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                     p_e ~ p_i, # Due to propagule pressure
                     p_i ~ d_BCP, # Dist. to BCP
                     d_BCP ~ Ur, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                     p_i ~ H, Ur ~ H, # Human pop. density
                     d_BCP ~ E, Ur ~ E, H ~ E, # Elevation
                     VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                     p_i ~ R, R ~ Ur, R ~ H,# Recreation
                     R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                     p_e ~ W_a, p_i ~ W_a, R ~ W_a, W_a ~ Ur, W_a ~ E, # Woodland area
                     W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, R ~ C_a, C_a ~ Ur, # Conifer area
                     W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, R ~ B_a, B_a ~ Ur, # Broadleaf area
                     W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                     latent = c("p_i", "p_e"), 
                     outcome = "y", 
                     exposure = "d_BCP",
                     coords = DAG_coords1)

d_BCP_d12a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                     p_e ~ p_i, # Due to propagule pressure
                     p_i ~ d_BCP, # Dist. to BCP
                     d_BCP ~ Ur, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                     p_i ~ H, Ur ~ H, # Human pop. density
                     d_BCP ~ E, Ur ~ E, H ~ E, # Elevation
                     VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                     p_i ~ R, R ~ Ur, R ~ H,# Recreation
                     R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                     p_e ~ W_a, p_i ~ W_a, R ~ W_a, W_a ~ Ur, W_a ~ E, # Woodland area
                     W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, R ~ C_a, C_a ~ Ur, # Conifer area
                     W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, R ~ B_a, B_a ~ Ur, # Broadleaf area
                     W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                     W_a ~ A, R ~ A, # Ancient woodland cover
                     latent = c("p_i", "p_e"), 
                     outcome = "y", 
                     exposure = "d_BCP",
                     coords = DAG_coords1)

d_BCP_d13a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                     p_e ~ p_i, # Due to propagule pressure
                     p_i ~ d_BCP, # Dist. to BCP
                     d_BCP ~ Ur, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                     p_i ~ H, Ur ~ H, # Human pop. density
                     d_BCP ~ E, Ur ~ E, H ~ E, # Elevation
                     VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                     p_i ~ R, R ~ Ur, R ~ H,# Recreation
                     R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                     p_e ~ W_a, p_i ~ W_a, R ~ W_a, W_a ~ Ur, W_a ~ E, # Woodland area
                     W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, R ~ C_a, C_a ~ Ur, # Conifer area
                     W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, R ~ B_a, B_a ~ Ur, # Broadleaf area
                     W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                     W_a ~ A, R ~ A, # Ancient woodland cover
                     p_i ~ Af, W_a ~ Af, C_a ~ Af, B_a ~ Af, # Afforestation
                     latent = c("p_i", "p_e"), 
                     outcome = "y", 
                     exposure = "d_BCP",
                     coords = DAG_coords1)

d_BCP_d14a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                     p_e ~ p_i, # Due to propagule pressure
                     p_i ~ d_BCP, # Dist. to BCP
                     d_BCP ~ Ur, p_i ~ Ur, p_e ~ Ur, # Urbanisation
                     p_i ~ H, Ur ~ H, # Human pop. density
                     d_BCP ~ E, Ur ~ E, H ~ E, # Elevation
                     VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                     p_i ~ R, R ~ Ur, R ~ H,# Recreation
                     R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                     p_e ~ W_a, p_i ~ W_a, R ~ W_a, W_a ~ Ur, W_a ~ E, # Woodland area
                     W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, R ~ C_a, C_a ~ Ur, # Conifer area
                     W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, R ~ B_a, B_a ~ Ur, # Broadleaf area
                     W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                     W_a ~ A, R ~ A, # Ancient woodland cover
                     p_i ~ Af, W_a ~ Af, C_a ~ Af, B_a ~ Af, # Afforestation
                     p_i ~ Df, W_a ~ Df, C_a ~ Df, B_a ~ Df, A ~ Df, # Deforestation
                     latent = c("p_i", "p_e"), 
                     outcome = "y", 
                     exposure = "d_BCP",
                     coords = DAG_coords1)


# Calculate minimum adjustment sets 
d_BCP_sets1 <- adjustmentSets_multi("d_BCP_d")
d_BCP_sets1$DAG_num <- 1:nrow(d_BCP_sets1)
d_BCP_sets1_sep <- separate_longer_delim(data = d_BCP_sets1, 
                                         cols = adjustment_sets, 
                                         delim = "; ")


# Create list of inlabru formulae from minimum adjustment sets
d_BCP_flist <- create_covariate_list(adj_sets = d_BCP_sets1_sep, 
                                     crosswalk = var_names_crosswalk,
                                     var_focal = "dist_BCP"
)


setwd("C:/temp/pests_analysis/ISDM_flists")
save(d_BCP_flist, file = "d_BCP_flist.Rdata")

# Create tidy format DAGs
dag_names <- ls(pattern = "d_BCP_d", envir = .GlobalEnv)
for(i in 1:length(dag_names)){
  d <- get(dag_names[i])
  
  d_tidy <- d %>% tidy_dagitty()
  
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "p_i" | d_tidy$data$name == "p_e", "latent", "observed")
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "y", "outcome", d_tidy$data$colour)
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "d_BCP", "exposure", d_tidy$data$colour)
  
  assign(paste0(dag_names[i],"_tidy"), d_tidy)
}

# Plot DAGs
p_list <- list()
dag_names <- ls(pattern = "^d_BCP_d.+_tidy$", envir = .GlobalEnv)
for(i in 1:length(dag_names)){
  d_temp <- get(dag_names[i])
  p_temp <- d_temp %>% 
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges() +
    geom_dag_point(aes(colour = colour)) +
    geom_dag_text(size = 3.5) +
    theme_dag() +
    theme(legend.position="none") +
    scale_color_manual(values = c("#cc2055", "#8d8d8d", "black", "#14a0f0"))
  p_list[[i]] <- p_temp
}


# Save all DAGs to .tiff
setwd("F:/INLA_outputs/DAGs")
for(i in seq(1, length(p_list), by = 8)){
  grid <- plot_grid(plotlist = p_list[i:(i+7)],
                    ncol = 4,
                    labels = i:(i+7))
  
  ggsave(filename = paste0("d_BCP_DAGs_",i,"_to_",i+7,".tiff"),
         plot = grid,
         dpi = "print",
         width = 16, 
         height = 7.8, 
         units = "cm",
         bg = "white",
         scale = 2.5)
}


# Conifer area ####
C_a_d01a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ C_a, p_e ~ C_a, # Conifer area
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "C_a",
                   coords = DAG_coords1)

C_a_d02a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ C_a, p_e ~ C_a, # Conifer area
                   p_e ~ W_a, p_i ~ W_a, W_a ~ C_a, # Woodland area
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "C_a",
                   coords = DAG_coords1)

C_a_d03a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ C_a, p_e ~ C_a, # Conifer area
                   p_e ~ W_a, p_i ~ W_a, W_a ~ C_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "C_a",
                   coords = DAG_coords1)

C_a_d04a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ C_a, p_e ~ C_a, # Conifer area
                   p_e ~ W_a, p_i ~ W_a, W_a ~ C_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "C_a",
                   coords = DAG_coords1)

C_a_d05a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ C_a, p_e ~ C_a, # Conifer area
                   p_e ~ W_a, p_i ~ W_a, W_a ~ C_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                   p_i ~ Af, C_a ~ Af, B_a ~ Af, W_a ~ Af, # Afforestation
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "C_a",
                   coords = DAG_coords1)

C_a_d06a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ C_a, p_e ~ C_a, # Conifer area
                   p_e ~ W_a, p_i ~ W_a, W_a ~ C_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                   p_i ~ Af, C_a ~ Af, B_a ~ Af, W_a ~ Af, # Afforestation
                   p_i ~ Df, C_a ~ Df, B_a ~ Df, W_a ~ Df, # Deforestation
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "C_a",
                   coords = DAG_coords1)

C_a_d07a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ C_a, p_e ~ C_a, # Conifer area
                   p_e ~ W_a, p_i ~ W_a, W_a ~ C_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                   p_i ~ Af, C_a ~ Af, B_a ~ Af, W_a ~ Af, # Afforestation
                   p_i ~ Df, C_a ~ Df, B_a ~ Df, W_a ~ Df, # Deforestation
                   W_a ~ A, A ~ Df, # Ancient woodland cover
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "C_a",
                   coords = DAG_coords1)

C_a_d08a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ C_a, p_e ~ C_a, # Conifer area
                   p_e ~ W_a, p_i ~ W_a, W_a ~ C_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                   p_i ~ Af, C_a ~ Af, B_a ~ Af, W_a ~ Af, # Afforestation
                   p_i ~ Df, C_a ~ Df, B_a ~ Df, W_a ~ Df, # Deforestation
                   W_a ~ A, A ~ Df, # Ancient woodland cover
                   W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "C_a",
                   coords = DAG_coords1)


C_a_d09a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ C_a, p_e ~ C_a, # Conifer area
                   p_e ~ W_a, p_i ~ W_a, W_a ~ C_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                   p_i ~ Af, C_a ~ Af, B_a ~ Af, W_a ~ Af, # Afforestation
                   p_i ~ Df, C_a ~ Df, B_a ~ Df, W_a ~ Df, # Deforestation
                   W_a ~ A, A ~ Df, # Ancient woodland cover
                   W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                   VPD ~ E, p_e ~ VPD, # VPD
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "C_a",
                   coords = DAG_coords1)

C_a_d10a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ C_a, p_e ~ C_a, # Conifer area
                   p_e ~ W_a, p_i ~ W_a, W_a ~ C_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                   p_i ~ Af, C_a ~ Af, B_a ~ Af, W_a ~ Af, # Afforestation
                   p_i ~ Df, C_a ~ Df, B_a ~ Df, W_a ~ Df, # Deforestation
                   W_a ~ A, A ~ Df, # Ancient woodland cover
                   W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                   VPD ~ E, p_e ~ VPD, # VPD
                   p_i ~ Ur, p_e ~ Ur, VPD ~ Ur, W_a ~ Ur, C_a ~ Ur, B_a ~ Ur, A ~ Ur, Ur ~ E, # Urbanisation
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "C_a",
                   coords = DAG_coords1)

C_a_d11a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ C_a, p_e ~ C_a, # Conifer area
                   p_e ~ W_a, p_i ~ W_a, W_a ~ C_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                   p_i ~ Af, C_a ~ Af, B_a ~ Af, W_a ~ Af, # Afforestation
                   p_i ~ Df, C_a ~ Df, B_a ~ Df, W_a ~ Df, # Deforestation
                   W_a ~ A, A ~ Df, # Ancient woodland cover
                   W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                   VPD ~ E, p_e ~ VPD, # VPD
                   p_i ~ Ur, p_e ~ Ur, VPD ~ Ur, W_a ~ Ur, C_a ~ Ur, B_a ~ Ur, A ~ Ur, Ur ~ E, # Urbanisation
                   p_i ~ H, Ur ~ H, H ~ E, # Human pop. density
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "C_a",
                   coords = DAG_coords1)

C_a_d12a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ C_a, p_e ~ C_a, # Conifer area
                   p_e ~ W_a, p_i ~ W_a, W_a ~ C_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                   p_i ~ Af, C_a ~ Af, B_a ~ Af, W_a ~ Af, # Afforestation
                   p_i ~ Df, C_a ~ Df, B_a ~ Df, W_a ~ Df, # Deforestation
                   W_a ~ A, A ~ Df, # Ancient woodland cover
                   W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                   VPD ~ E, p_e ~ VPD, # VPD
                   p_i ~ Ur, p_e ~ Ur, VPD ~ Ur, W_a ~ Ur, C_a ~ Ur, B_a ~ Ur, A ~ Ur, Ur ~ E, # Urbanisation
                   p_i ~ H, Ur ~ H, H ~ E, # Human pop. density
                   p_i ~ R, R ~ Ur, R ~ W_a, R ~ B_a, R ~ C_a, R ~ A, R ~ H,# Recreation
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "C_a",
                   coords = DAG_coords1)

C_a_d13a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ C_a, p_e ~ C_a, # Conifer area
                   p_e ~ W_a, p_i ~ W_a, W_a ~ C_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                   p_i ~ Af, C_a ~ Af, B_a ~ Af, W_a ~ Af, # Afforestation
                   p_i ~ Df, C_a ~ Df, B_a ~ Df, W_a ~ Df, # Deforestation
                   W_a ~ A, A ~ Df, # Ancient woodland cover
                   W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                   VPD ~ E, p_e ~ VPD, # VPD
                   p_i ~ Ur, p_e ~ Ur, VPD ~ Ur, W_a ~ Ur, C_a ~ Ur, B_a ~ Ur, A ~ Ur, Ur ~ E, # Urbanisation
                   p_i ~ H, Ur ~ H, H ~ E, # Human pop. density
                   p_i ~ R, R ~ Ur, R ~ W_a, R ~ B_a, R ~ C_a, R ~ A, R ~ H,# Recreation
                   R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "C_a",
                   coords = DAG_coords1)

C_a_d14a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ C_a, p_e ~ C_a, # Conifer area
                   p_e ~ W_a, p_i ~ W_a, W_a ~ C_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   W_c ~ W_a, p_i ~ W_c, # Woodland connectivity
                   p_i ~ Af, C_a ~ Af, B_a ~ Af, W_a ~ Af, # Afforestation
                   p_i ~ Df, C_a ~ Df, B_a ~ Df, W_a ~ Df, # Deforestation
                   W_a ~ A, A ~ Df, # Ancient woodland cover
                   W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                   VPD ~ E, p_e ~ VPD, # VPD
                   p_i ~ Ur, p_e ~ Ur, VPD ~ Ur, W_a ~ Ur, C_a ~ Ur, B_a ~ Ur, A ~ Ur, Ur ~ E, # Urbanisation
                   p_i ~ H, Ur ~ H, H ~ E, # Human pop. density
                   p_i ~ R, R ~ Ur, R ~ W_a, R ~ B_a, R ~ C_a, R ~ A, R ~ H,# Recreation
                   R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                   p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "C_a",
                   coords = DAG_coords1)

# Calculate minimum adjustment sets 
C_a_sets1 <- adjustmentSets_multi("C_a_d")
C_a_sets1$DAG_num <- 1:nrow(C_a_sets1)
C_a_sets1_sep <- separate_longer_delim(data = C_a_sets1, 
                                       cols = adjustment_sets, 
                                       delim = "; ")


# Create list of inlabru formulae from minimum adjustment sets
C_a_flist <- create_covariate_list(adj_sets = C_a_sets1_sep, 
                                   crosswalk = var_names_crosswalk,
                                   var_focal = "conifer_area_Ha"
)


setwd("C:/temp/pests_analysis/ISDM_flists")
save(C_a_flist, file = "C_a_flist.Rdata")

# Create tidy format DAGs
dag_names <- ls(pattern = "C_a_d", envir = .GlobalEnv)
for(i in 1:length(dag_names)){
  d <- get(dag_names[i])
  
  d_tidy <- d %>% tidy_dagitty()
  
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "p_i" | d_tidy$data$name == "p_e", "latent", "observed")
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "y", "outcome", d_tidy$data$colour)
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "C_a", "exposure", d_tidy$data$colour)
  
  assign(paste0(dag_names[i],"_tidy"), d_tidy)
}

# Plot DAGs
p_list <- list()
dag_names <- ls(pattern = "^C_a_d.+_tidy$", envir = .GlobalEnv)
for(i in 1:length(dag_names)){
  d_temp <- get(dag_names[i])
  p_temp <- d_temp %>% 
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges() +
    geom_dag_point(aes(colour = colour)) +
    geom_dag_text(size = 3.5) +
    theme_dag() +
    theme(legend.position="none") +
    scale_color_manual(values = c("#cc2055", "#8d8d8d", "black", "#14a0f0"))
  p_list[[i]] <- p_temp
}


# Save all DAGs to .tiff
setwd("F:/INLA_outputs/DAGs")
for(i in seq(1, length(p_list), by = 8)){
  grid <- plot_grid(plotlist = p_list[i:(i+7)],
                    ncol = 4,
                    labels = i:(i+7))
  
  ggsave(filename = paste0("C_a_DAGs_",i,"_to_",i+7,".tiff"),
         plot = grid,
         dpi = "print",
         width = 16, 
         height = 7.8, 
         units = "cm",
         bg = "white",
         scale = 2.5)
}

# Woodland connectivity ####
W_c_d01a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ W_c, # Woodland connectivity
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "W_c",
                   coords = DAG_coords1)

W_c_d02a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ W_c, # Woodland connectivity
                   p_e ~ W_a, p_i ~ W_a, W_c ~ W_a, # Woodland area
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "W_c",
                   coords = DAG_coords1)

W_c_d03a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ W_c, # Woodland connectivity
                   p_e ~ W_a, p_i ~ W_a, W_c ~ W_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "W_c",
                   coords = DAG_coords1)

W_c_d04a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ W_c, # Woodland connectivity
                   p_e ~ W_a, p_i ~ W_a, W_c ~ W_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   p_i ~ C_a, p_e ~ C_a, W_a ~ C_a,# Conifer area
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "W_c",
                   coords = DAG_coords1)

W_c_d05a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ W_c, # Woodland connectivity
                   p_e ~ W_a, p_i ~ W_a, W_c ~ W_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   p_i ~ C_a, p_e ~ C_a, W_a ~ C_a,# Conifer area
                   p_i ~ Af, C_a ~ Af, B_a ~ Af, W_a ~ Af, # Afforestation
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "W_c",
                   coords = DAG_coords1)

W_c_d06a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ W_c, # Woodland connectivity
                   p_e ~ W_a, p_i ~ W_a, W_c ~ W_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   p_i ~ C_a, p_e ~ C_a, W_a ~ C_a,# Conifer area
                   p_i ~ Af, C_a ~ Af, B_a ~ Af, W_a ~ Af, # Afforestation
                   p_i ~ Df, C_a ~ Df, B_a ~ Df, W_a ~ Df, # Deforestation
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "W_c",
                   coords = DAG_coords1)

W_c_d07a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ W_c, # Woodland connectivity
                   p_e ~ W_a, p_i ~ W_a, W_c ~ W_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   p_i ~ C_a, p_e ~ C_a, W_a ~ C_a,# Conifer area
                   p_i ~ Af, C_a ~ Af, B_a ~ Af, W_a ~ Af, # Afforestation
                   p_i ~ Df, C_a ~ Df, B_a ~ Df, W_a ~ Df, # Deforestation
                   W_a ~ A, A ~ Df, # Ancient woodland cover
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "W_c",
                   coords = DAG_coords1)

W_c_d08a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ W_c, # Woodland connectivity
                   p_e ~ W_a, p_i ~ W_a, W_c ~ W_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   p_i ~ C_a, p_e ~ C_a, W_a ~ C_a,# Conifer area
                   p_i ~ Af, C_a ~ Af, B_a ~ Af, W_a ~ Af, # Afforestation
                   p_i ~ Df, C_a ~ Df, B_a ~ Df, W_a ~ Df, # Deforestation
                   W_a ~ A, A ~ Df, # Ancient woodland cover
                   W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "W_c",
                   coords = DAG_coords1)

W_c_d09a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ W_c, # Woodland connectivity
                   p_e ~ W_a, p_i ~ W_a, W_c ~ W_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   p_i ~ C_a, p_e ~ C_a, W_a ~ C_a,# Conifer area
                   p_i ~ Af, C_a ~ Af, B_a ~ Af, W_a ~ Af, # Afforestation
                   p_i ~ Df, C_a ~ Df, B_a ~ Df, W_a ~ Df, # Deforestation
                   W_a ~ A, A ~ Df, # Ancient woodland cover
                   W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                   VPD ~ E, p_e ~ VPD, # VPD
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "W_c",
                   coords = DAG_coords1)

W_c_d10a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ W_c, # Woodland connectivity
                   p_e ~ W_a, p_i ~ W_a, W_c ~ W_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   p_i ~ C_a, p_e ~ C_a, W_a ~ C_a,# Conifer area
                   p_i ~ Af, C_a ~ Af, B_a ~ Af, W_a ~ Af, # Afforestation
                   p_i ~ Df, C_a ~ Df, B_a ~ Df, W_a ~ Df, # Deforestation
                   W_a ~ A, A ~ Df, # Ancient woodland cover
                   W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                   VPD ~ E, p_e ~ VPD, # VPD
                   p_i ~ Ur, p_e ~ Ur, VPD ~ Ur, W_a ~ Ur, C_a ~ Ur, B_a ~ Ur, A ~ Ur, Ur ~ E, # Urbanisation
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "W_c",
                   coords = DAG_coords1)

W_c_d11a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ W_c, # Woodland connectivity
                   p_e ~ W_a, p_i ~ W_a, W_c ~ W_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   p_i ~ C_a, p_e ~ C_a, W_a ~ C_a,# Conifer area
                   p_i ~ Af, C_a ~ Af, B_a ~ Af, W_a ~ Af, # Afforestation
                   p_i ~ Df, C_a ~ Df, B_a ~ Df, W_a ~ Df, # Deforestation
                   W_a ~ A, A ~ Df, # Ancient woodland cover
                   W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                   VPD ~ E, p_e ~ VPD, # VPD
                   p_i ~ Ur, p_e ~ Ur, VPD ~ Ur, W_a ~ Ur, C_a ~ Ur, B_a ~ Ur, A ~ Ur, Ur ~ E, # Urbanisation
                   p_i ~ H, Ur ~ H, H ~ E, # Human pop. density
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "W_c",
                   coords = DAG_coords1)

W_c_d12a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ W_c, # Woodland connectivity
                   p_e ~ W_a, p_i ~ W_a, W_c ~ W_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   p_i ~ C_a, p_e ~ C_a, W_a ~ C_a,# Conifer area
                   p_i ~ Af, C_a ~ Af, B_a ~ Af, W_a ~ Af, # Afforestation
                   p_i ~ Df, C_a ~ Df, B_a ~ Df, W_a ~ Df, # Deforestation
                   W_a ~ A, A ~ Df, # Ancient woodland cover
                   W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                   VPD ~ E, p_e ~ VPD, # VPD
                   p_i ~ Ur, p_e ~ Ur, VPD ~ Ur, W_a ~ Ur, C_a ~ Ur, B_a ~ Ur, A ~ Ur, Ur ~ E, # Urbanisation
                   p_i ~ H, Ur ~ H, H ~ E, # Human pop. density
                   p_i ~ R, R ~ Ur, R ~ W_a, R ~ B_a, R ~ C_a, R ~ A, R ~ H,# Recreation
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "W_c",
                   coords = DAG_coords1)

W_c_d13a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ W_c, # Woodland connectivity
                   p_e ~ W_a, p_i ~ W_a, W_c ~ W_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   p_i ~ C_a, p_e ~ C_a, W_a ~ C_a,# Conifer area
                   p_i ~ Af, C_a ~ Af, B_a ~ Af, W_a ~ Af, # Afforestation
                   p_i ~ Df, C_a ~ Df, B_a ~ Df, W_a ~ Df, # Deforestation
                   W_a ~ A, A ~ Df, # Ancient woodland cover
                   W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                   VPD ~ E, p_e ~ VPD, # VPD
                   p_i ~ Ur, p_e ~ Ur, VPD ~ Ur, W_a ~ Ur, C_a ~ Ur, B_a ~ Ur, A ~ Ur, Ur ~ E, # Urbanisation
                   p_i ~ H, Ur ~ H, H ~ E, # Human pop. density
                   p_i ~ R, R ~ Ur, R ~ W_a, R ~ B_a, R ~ C_a, R ~ A, R ~ H,# Recreation
                   R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "W_c",
                   coords = DAG_coords1)

W_c_d14a <- dagify(y ~ p_i + p_e, # Pest presence (y) = f(introduction probability, establishment probability)
                   p_e ~ p_i, # Due to propagule pressure
                   p_i ~ W_c, # Woodland connectivity
                   p_e ~ W_a, p_i ~ W_a, W_c ~ W_a, # Woodland area
                   W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, # Broadleaf area
                   p_i ~ C_a, p_e ~ C_a, W_a ~ C_a,# Conifer area
                   p_i ~ Af, C_a ~ Af, B_a ~ Af, W_a ~ Af, # Afforestation
                   p_i ~ Df, C_a ~ Df, B_a ~ Df, W_a ~ Df, # Deforestation
                   W_a ~ A, A ~ Df, # Ancient woodland cover
                   W_a ~ E, C_a ~ E, B_a ~ E, # Elevation
                   VPD ~ E, p_e ~ VPD, # VPD
                   p_i ~ Ur, p_e ~ Ur, VPD ~ Ur, W_a ~ Ur, C_a ~ Ur, B_a ~ Ur, A ~ Ur, Ur ~ E, # Urbanisation
                   p_i ~ H, Ur ~ H, H ~ E, # Human pop. density
                   p_i ~ R, R ~ Ur, R ~ W_a, R ~ B_a, R ~ C_a, R ~ A, R ~ H,# Recreation
                   R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                   p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   exposure = "W_c",
                   coords = DAG_coords1)

# Calculate minimum adjustment sets 
W_c_sets1 <- adjustmentSets_multi("W_c_d")
W_c_sets1$DAG_num <- 1:nrow(W_c_sets1)
W_c_sets1_sep <- separate_longer_delim(data = W_c_sets1, 
                                       cols = adjustment_sets, 
                                       delim = "; ")


# Create list of inlabru formulae from minimum adjustment sets
W_c_flist <- create_covariate_list(adj_sets = W_c_sets1_sep, 
                                   crosswalk = var_names_crosswalk,
                                   var_focal = "woodland_connectivity_sigma3"
)

setwd("C:/temp/pests_analysis/ISDM_flists")
save(W_c_flist, file = "W_c_flist.Rdata")

# Create tidy format DAGs
dag_names <- ls(pattern = "W_c_d", envir = .GlobalEnv)
for(i in 1:length(dag_names)){
  d <- get(dag_names[i])
  
  d_tidy <- d %>% tidy_dagitty()
  
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "p_i" | d_tidy$data$name == "p_e", "latent", "observed")
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "y", "outcome", d_tidy$data$colour)
  d_tidy$data$colour <- ifelse(d_tidy$data$name == "W_c", "exposure", d_tidy$data$colour)
  
  assign(paste0(dag_names[i],"_tidy"), d_tidy)
}

# Plot DAGs
p_list <- list()
dag_names <- ls(pattern = "^W_c_d.+_tidy$", envir = .GlobalEnv)
for(i in 1:length(dag_names)){
  d_temp <- get(dag_names[i])
  p_temp <- d_temp %>% 
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges() +
    geom_dag_point(aes(colour = colour)) +
    geom_dag_text(size = 3.5) +
    theme_dag() +
    theme(legend.position="none") +
    scale_color_manual(values = c("#cc2055", "#8d8d8d", "black", "#14a0f0"))
  p_list[[i]] <- p_temp
}


# Save all DAGs to .tiff
setwd("F:/INLA_outputs/DAGs")
for(i in seq(1, length(p_list), by = 8)){
  grid <- plot_grid(plotlist = p_list[i:(i+7)],
                    ncol = 4,
                    labels = i:(i+7))
  
  ggsave(filename = paste0("W_c_DAGs_",i,"_to_",i+7,".tiff"),
         plot = grid,
         dpi = "print",
         width = 16, 
         height = 7.8, 
         units = "cm",
         bg = "white",
         scale = 2.5)
}


# Maximum DAG plot for main text ####
DAG_main <- dagify(y ~ p_i + p_e, 
                   p_e ~ p_i, 
                   p_i ~ R, # Recreation
                   R ~ H, p_i ~ H, # Human pop. density
                   Ur ~ H, p_i ~ Ur, p_e ~ Ur, R ~ Ur, # Urbanisation
                   H ~ E, R ~ E, Ur ~ E, # Elevation
                   VPD ~ E, p_e ~ VPD, VPD ~ Ur, # VPD
                   p_i ~ d_BCP, d_BCP ~ Ur, d_BCP ~ E, # Dist. to BCP
                   p_e ~ W_a, W_a ~ E, p_i ~ W_a, W_a ~ Ur, R ~ W_a, # Woodland area
                   R ~ d_G, p_i ~ d_G, p_e ~ d_G, # Distance to park/garden
                   R ~ A, W_a ~ A, A ~ Ur, # Ancient woodland cover
                   C_a ~ E, W_a ~ C_a, p_i ~ C_a, p_e ~ C_a, C_a ~ Ur, R ~ C_a, # Conifer area
                   B_a ~ E, W_a ~ B_a, p_i ~ B_a, p_e ~ B_a, B_a ~ Ur, R ~ B_a, # Broadleaf area
                   W_c ~ W_a, p_i ~ W_c, R ~ W_c, # Woodland connectivity
                   p_i ~ Af, W_a ~ Af, C_a ~ Af, B_a ~ Af, # Afforestation
                   p_i ~ Df, W_a ~ Df, C_a ~ Df, B_a ~ Df, A ~ Df, # Deforestation
                   latent = c("p_i", "p_e"), 
                   outcome = "y", 
                   coords = DAG_coords1)
DAG_main_tidy <- DAG_main %>% tidy_dagitty()

DAG_main_tidy$data$colour <- ifelse(DAG_main_tidy$data$name == "p_i" | DAG_main_tidy$data$name == "p_e", "latent", "observed")
DAG_main_tidy$data$colour <- ifelse(DAG_main_tidy$data$name == "y", "outcome", DAG_main_tidy$data$colour)

p_main <- DAG_main_tidy %>% 
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges() +
  geom_dag_node(aes(colour = colour)) +
  geom_dag_text(size = 3.2) +
  theme_dag() +
  theme(legend.position="none") +
  scale_colour_manual(values = c("#8d8d8d", "black", "#14a0f0"))

setwd("C:/temp/pests_analysis/figures_new")
ggsave(filename = "DAG_main.tiff",
       plot = p_main,
       dpi = 600,
       width = 165, 
       height = 110, 
       units = "mm",
       bg = "white",
       scale = 1)
