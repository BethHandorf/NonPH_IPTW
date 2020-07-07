# NonPH_IPTW
Code to run simulations for manuscript "Analysis of survival data with non-proportional hazards: A comparison of propensity score weighted methods"

Contents  
simulated_covariates.csv - This file contains simulated (synthetic) sets of covariates constructed to be similar to those found in the National Cancer Database.  
simulate_datasets_logTime.R - This R code generates simulated treatments and outcomes for patients.  Each simulated dataset is appended to a large file which will contain all simulated datasets for later processing.  This code is designed to be run in batch mode.  
run_all_methods_batch.R - This R code selects a single simulated dataset and obtains the true and estimated survival results.  It also runs a bootstrap to obtain variance estimates and get upper and lower confidence limits.  This code is designed to be run in batch mode.  
execute_all_methods.R - This R code obtains estimates of survival outcomes using each of the methods described in the manuscript.  

