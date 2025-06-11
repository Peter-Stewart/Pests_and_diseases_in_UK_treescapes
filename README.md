# Pests_and_diseases_in_UK_treescapes
Code to reproduce the analyses in the manuscript _Patterns and drivers of pest and disease occurrence in UK treescapes._

## Description of files
- **APHA_format.R**: Clean and format the APHA dataset and save as Rdata file.
- **Create_var_names_crosswalk.R**: Create crosswalk between DAG nodes and variable names as they appear in the covariate data.
- **Drivers_DAGs.R**: Produce DAGs and calculate minimum adjustment sets. 
- **Drivers_analysis.R**: Fit ISDMs and obtain estimates of effects of focal drivers. 
- **INLA_dataprep.R**: Prepare datasets for fitting ISDMs.
- **Prediction_mesh_correlations.R**: Calculate correlations between spatial predictions from ISDMs using different resolution meshes for the spatial random effect.
- **Predictions_analysis.R**: Fit ISDMs and produce spatial predictions of pest and disease intensity.
- **Process_to_1km.R**: Produce spatial covariate layers at 1km resolution.
- **THDAS_format.R**: Clean and format the THDAS dataset and save as Rdata file.
- **WAIC_calculation.R**: Fit glms using random subsets of the covariates, and store the WAIC values for each model. 
