

Code to run simulations for manuscript "Analysis of survival data with non-proportional hazards: A comparison of propensity score weighted methods" https://arxiv.org/abs/2009.00785


Notes: The user will need to set their working directory to a folder containing all scripts and datafiles.


Contents

1. Data files

simulated_covariates.csv - This file contains simulated (synthetic) sets of covariates constructed to be similar to those found in the National Cancer Database.

FNPH.1_SNPH.0.25_TE.-0.69_SS.5000_seedRange.1-5.csv - this file contains the simulated data used for the first 5 itterations of the base case.

Bias_results_graph_input.csv - this file contains the numerical results used to create Figure 1

Var_results_graph_input.csv - this file contains the numerical results used to create Figure A1


2. Example run of selected simulations 

Note: Because of the bootstrap step, each itteration of the simulation was highly time-intensive to run.  Therefore, all simulations were run in parallel.  We therefore provide a way for the user to run a few itterations with a subset of the simulated datasets.

run_all_methods_example.R - This R code uses one simulated dataset (of itterations 1-5) and obtains the true and estimated survival results. It also runs a bootstrap to obtain variance estimates and get upper and lower confidence limits. 

execute_all_methods.R - This R code obtains estimates of survival outcomes using each of the methods described in the manuscript.


3. Code which fully re-creates the simulation results presented in section 4

3A. Create datasets

Data was simulated using seeds 1-500.  Successive simulations were appended to a single file.  

simulate_datasets_logTime.R - This R code generates simulated treatments and outcomes for patients. The functional form of the non-proportionality is a log-time effect. Each simulated dataset is appended to a large file which will contain all simulated datasets for later processing. This code is designed to be run in batch mode.  Parameters can be modified by entering arguments into the command line call for this function.  Alternatively, this can be done manually in the code file.

Here are the base case values (scenario 1)
SNPH=0.25 (strength of the non-proportional hazards as a function of log time)
TE=-log(2) (treatment effect at time 0)
SS=5000 (sample size)

Scenario 3 uses the following inputs:
SNPH=0.125 
TE=-log(2) 
SS=5000 

Scenario 4 uses the following inputs:
SNPH=0.25 
TE=-0.41
SS=5000 

Scenario 5 uses the following inputs:
SNPH=0.25 
TE=-log(2) 
SS=1000 


simulate_datasets_PWC.R - This additional file generates simulated treatments and outcomes for patients. The functional form of the non-proportionality is a piece-wise effect (Scenario 2). Each simulated dataset is appended to a large file which will contain all simulated datasets for later processing. This code is designed to be run in batch mode.

3B. Run analysis methods

Each simulated set is run in parallel.  The input datafile is identified, as is the particular seed of interest.  These should be entered via the command line.  The script then looks up the rows for the respective seeds in the previously generated simulated datafiles.  

run_all_methods_batch.R - This R code selects a single simulated dataset and obtains the true and estimated survival results. It also runs a bootstrap to obtain variance estimates and get upper and lower confidence limits. This code is designed to be run in batch mode.  Results from this code are used to generate Table 1 and A2, and Figures 1 and A1.

execute_all_methods.R - This R code obtains estimates of survival outcomes using each of the methods described in the manuscript.

4. Code which recreates simulation figures

results_graphs_v2.R - this R code creates Figure 1 and Figure A1

5. Code which runs the clinical examples from section 5
Note: *Data to run this code is not available here due to restrictions from the NCDB data use agreement*

sarcoma.R - This file runs the analysis described in section 5.1 and creates Figure 2 and the data for Table 2. 

execute_all_sarcoma.R - This R code obtains estimates of survival outcomes using each of the methods described in the manuscript, with specific inputs and covariates for the sarcoma example.

kidney.R - This file runs the analysis described in section 5.1 and creates Figure 2 and the data for Table 2. 

execute_all_kidney.R - This R code obtains estimates of survival outcomes using each of the methods described in the manuscript, with specific inputs and covariates for the renal example.


The following software was used:

R version 3.5.1 (2018-07-02) 
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 6 (Core)

Matrix products: default
BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
 [1] khroma_1.2.0      wesanderson_0.3.6 gridExtra_2.3     reshape2_1.4.3    muhaz_1.2.6.1    
 [6] ggplot2_3.2.1     simsurv_0.2.3     pseudo_1.4.3      geepack_1.2-1     KMsurv_0.1-5     
[11] flexsurv_1.1.1    survMisc_0.5.5    survival_2.44-1   survRM2_1.0-2 

loaded via a namespace (and not attached):
 [1] deSolve_1.24       tinytex_0.16       zoo_1.8-4          tidyselect_0.2.5   xfun_0.9         
 [6] inum_1.0-1         purrr_0.3.3        splines_3.5.1      lattice_0.20-38    mstate_0.2.12     
[11] colorspace_1.4-1   vctrs_0.2.0       generics_0.0.2     rlang_0.4.0       pillar_1.4.2      
[16] withr_2.1.2        glue_1.3.1         RColorBrewer_1.1-2 plyr_1.8.5         lifecycle_0.1.0   
[21] stringr_1.4.0      munsell_0.5.0      gtable_0.3.0       mvtnorm_1.0-11      knitr_1.25        
[26] broom_0.5.6        Rcpp_1.0.7         xtable_1.8-4       scales_1.1.1       backports_1.1.8   
[31] km.ci_0.5-2        stringi_1.4.3      dplyr_0.8.3        quadprog_1.5-7     tools_3.5.1       
[36] magrittr_1.5       tibble_2.1.3       Formula_1.2-3      crayon_1.3.4       tidyr_1.0.0       
[41] pkgconfig_2.0.3    partykit_1.2-5     ellipsis_0.3.0     MASS_7.3-51.4      libcoin_1.0-5     
[46] Matrix_1.2-17      data.table_1.12.2  R6_2.4.0           rpart_4.1-15       nlme_3.1-141      
[51] compiler_3.5.1    
