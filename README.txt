Code to run simulations for manuscript "Analysis of survival data with non-proportional hazards: A comparison of propensity score weighted methods" https://arxiv.org/abs/2009.00785

Notes: The user will need to set their working directory to a folder containing all scripts and datafiles.

Contents
1. Simulation data input files

simulated_covariates.csv - This file contains simulated (synthetic) sets of covariates constructed to be similar to those found in the National Cancer Database.
This file has contents summarized in Data_dictionary_sim_data.xlsx, with additional information in the text of the manuscript

FNPH.1_SNPH.0.25_TE.-0.69_SS.5000_seedRange.1-5.csv - this file contains the simulated data used for the first 5 itterations of the base case.

###################################################


2. Example run of selected simulations

Note: Because of the bootstrap step, each itteration of the simulation was highly time-intensive to run. Therefore, all simulations were run in parallel. We therefore provide a way for the user to run a few itterations with a subset of the simulated datasets.

run_all_methods_example.R - This R code uses one simulated dataset (of itterations 1-5) and obtains the true and estimated survival results. 
                            It also runs a bootstrap to obtain variance estimates and get upper and lower confidence limits.
			    Please note the software version requirements at the end of this document.  This code may not run properly
			    with different versions as some required functions have changed over time.

execute_all_methods.R - This R code obtains estimates of survival outcomes using each of the methods described in the manuscript.

###################################################

3. Code which fully re-creates the simulation results presented in section 4

##
3A. Create datasets

Data was simulated using seeds 1-500. Successive simulations were appended to a single file.

simulate_datasets.R - This R code generates simulated treatments and outcomes for patients. Each simulated dataset is appended to a large file which will contain all simulated datasets for later processing. This code is designed to be run in batch mode. Parameters can be modified by entering arguments into the command line call for this function. Alternatively, this can be done manually in the code file.

Here are the base case values (scenario 1) SNPH=0.25 (strength of the non-proportional hazards as a function of log time) TE=-log(2) (treatment effect at time 0) SS=5000 (sample size)

Scenario 2 uses a piece-wise constant non-proportional hazards function

Scenario 3 uses the following inputs: SNPH=0.125 TE=-log(2) SS=5000

Scenario 4 uses the following inputs: SNPH=0.25 TE=-0.41 SS=5000

Scenario 5 uses the following inputs: SNPH=0.25 TE=-log(2) SS=1000

Scenario 6 uses a generalized gamma distribution with the following parameters: TE=-log(2) mu=2 (intercept) sigma=1.2 labmda=-0.5

Scenario 7 uses a generalized gamma distribution with the following parameters: TE=-log(2) mu=3.5 (intercept) sigma=1.2 labmda=2.5

##
3B. Run analysis methods

Each simulated set is run in parallel. The input datafile is identified, as is the particular seed of interest. These should be entered via the command line. The script then looks up the rows for the respective seeds in the previously generated simulated datafiles.

run_all_methods_batch.R - This R code selects a single simulated dataset and obtains the true and estimated survival results. It also runs a bootstrap to obtain variance estimates and get upper and lower confidence limits. 
                        This code is designed to be run in batch mode. Results from this code are used to generate Table 1 and A2, and Figures 1 and A1.

execute_all_methods.R - This R code obtains estimates of survival outcomes using each of the methods described in the manuscript.

###################################################

4. Intermediate results files

All disaggregated results of the simulation study are available in the folder "Intermediate_results"

These file names start with the estimand, followed by the scenario:

res.2y - 2 year survival
res.5y - 5 year survival
res.10y - 10 year survival
res.median - Median survival
res.rms - restricted mean survival

FNPH.1_SNPH.0.25_TE.0.69_SS.5000 - scenario 1
FNPH.2_SNPH.0.25_TE.0.69_SS.5000 - scenario 2
FNPH.1_SNPH.0.125_TE.0.69_SS.5000 - scenario 3
FNPH.1_SNPH.0.25_TE.0.41_SS.5000 - scenario 4
FNPH.1_SNPH.0.25_TE.0.69_SS.1000 - scenario 5
FNPH.3_SNPH.0.25_TE.0.69Q-0.5_SS.5000 - scenario 6
FNPH.3_SNPH.0.25_TE.0.69Q2.5_SS.5000 - scenario 7

So, for example, the file 
res.2y.FNPH.1_SNPH.0.25_TE.0.69_SS.5000.csv is the results for 2 year survival for scenario 1

A data dictionary explaining what each of the columns measures is available in the file:
Simulations_intermediate_results_data_dictionary.csv


###################################################

5. Code which recreates simulation figures and tables

create_MS_tables.R - this code summarizes the intermediate results files presented listed in (4)
		     and creates csv files.  Tables for confidence intervals are created directly.
		     The summarized results which create the figures are also saved.

Fig1_input_bias.csv - this file contains the numerical results used to create Figure 1

FigA1_input_var.csv - this file contains the numerical results used to create Figure A1

results_graphs_v2.R - this R code creates Figure 1 and Figure A1

###################################################


6. Code which runs the clinical examples from section 5

Note: The actual data use in this paper is not available here due to restrictions from the NCDB data use agreement

Section 5.1

sarcoma.R - This file runs the analysis described in section 5.1 and creates Figure 2 and the data for Table 2.

execute_all_sarcoma.R - This R code obtains estimates of survival outcomes using each of the methods described in the manuscript, with specific inputs and covariates for the sarcoma example.

sarcoma_sim_set.csv - This a synthetic set of data created to be similar to the actual sarcoma dataset so that the sarcoma code
                     can be run and tested.  THIS SET IS FULLY SIMULATED AND CONTAINS NO REAL PATIENT DATA.

This file has contents summarized in Data_dictionary_sim_data.xlsx

Notes for demonstration run: 
-This code will not reproduce the results from the mansucript as it uses a synthetic sample.  
 It is included for demonstration purposes only
-To make the run time of the demonstration reasonable, the number of bootstrap iterations 
 is set to 50.  In the paper we used 500.  This can be changed in line 166 in the file sarcoma.R

Section 5.2
kidney.R - This file runs the analysis described in section 5.1 and creates Figure 2 and the data for Table 2.

execute_all_kidney.R - This R code obtains estimates of survival outcomes using each of the methods described in the manuscript, with specific inputs and covariates for the renal example.

kidney_sim_set.csv - This a synthetic set of data created to be similar to the actual kidney dataset so that the kidney code
                     can be run and tested.  THIS SET IS FULLY SIMULATED AND CONTAINS NO REAL PATIENT DATA.

This file has contents summarized in Data_dictionary_sim_data.xlsx

Notes for demonstration run: 
-This code will not reproduce the results from the manuscript as it uses a synthetic sample.  
 It is included for demonstration purposes only
-The analysis presented in the code required a large amount of working memory for the pseudo
 observations approach.  to make this feasible to run on a desktop computer, this example
 analysis only uses 10% of the 32,000 samples.  The standard errors are therefore larger
 in the example run than they were in the original run.  The full sample can be used by
 commenting out line 375 and uncommenting line 377 in the execute_all_kidney.R file
-To make the run time of the demonstration reasonable, the number of bootstrap iterations 
 is set to 50.  In the paper we used 500.  This can be changed in line 180 in the file kidney.R

###################################################

7. Interim simulated datasets

Simulated survival times are available in the following files:
FNPH.1_SNPH.0.25_TE.0.69_SS.5000.csv - scenario 1 & 5
FNPH.2_SNPH.0.25_TE.0.69_SS.5000.csv - scenario 2
FNPH.1_SNPH.0.125_TE.0.69_SS.5000.csv - scenario 3
FNPH.1_SNPH.0.25_TE.0.41_SS.5000 - scenario 4
FNPH.3_SNPH.0.25_TE.0.69PSForm1Q-0.5_SS.5000.csv - scenario 6
FNPH.3_SNPH.0.25_TE.0.69PSForm1Q2.5_SS.5000.csv - scenario 7

These data files are too large to contain in the supplement but can be downloaded here:
https://drive.google.com/file/d/1nCpzlRJI0BrQhUDL6YA_T74Lk9wf8UZs/view?usp=sharing

Note: seeds 1-500 were used in the final manuscript

Note: scenario 5 uses the first 500 patients for each treatment arm simulated in scenario 1



###################################################

The following software was used (unix batch mode):

R version 3.5.1 (2018-07-02) Platform: x86_64-pc-linux-gnu (64-bit) Running under: CentOS Linux 6 (Core)

Matrix products: default BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so

locale: [1] LC_CTYPE=en_US.UTF-8 LC_NUMERIC=C [3] LC_TIME=en_US.UTF-8 LC_COLLATE=en_US.UTF-8 [5] LC_MONETARY=en_US.UTF-8 LC_MESSAGES=en_US.UTF-8 [7] LC_PAPER=en_US.UTF-8 LC_NAME=C [9] LC_ADDRESS=C LC_TELEPHONE=C [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages: [1] stats graphics grDevices utils datasets methods base

other attached packages: [1] khroma_1.2.0 wesanderson_0.3.6 gridExtra_2.3 reshape2_1.4.3 muhaz_1.2.6.1
[6] ggplot2_3.2.1 simsurv_0.2.3 pseudo_1.4.3 geepack_1.2-1 KMsurv_0.1-5
[11] flexsurv_1.1.1 survMisc_0.5.5 survival_2.44-1 survRM2_1.0-2

loaded via a namespace (and not attached): [1] deSolve_1.24 tinytex_0.16 zoo_1.8-4 tidyselect_0.2.5 xfun_0.9
[6] inum_1.0-1 purrr_0.3.3 splines_3.5.1 lattice_0.20-38 mstate_0.2.12
[11] colorspace_1.4-1 vctrs_0.2.0 generics_0.0.2 rlang_0.4.0 pillar_1.4.2
[16] withr_2.1.2 glue_1.3.1 RColorBrewer_1.1-2 plyr_1.8.5 lifecycle_0.1.0
[21] stringr_1.4.0 munsell_0.5.0 gtable_0.3.0 mvtnorm_1.0-11 knitr_1.25
[26] broom_0.5.6 Rcpp_1.0.7 xtable_1.8-4 scales_1.1.1 backports_1.1.8
[31] km.ci_0.5-2 stringi_1.4.3 dplyr_0.8.3 quadprog_1.5-7 tools_3.5.1
[36] magrittr_1.5 tibble_2.1.3 Formula_1.2-3 crayon_1.3.4 tidyr_1.0.0
[41] pkgconfig_2.0.3 partykit_1.2-5 ellipsis_0.3.0 MASS_7.3-51.4 libcoin_1.0-5
[46] Matrix_1.2-17 data.table_1.12.2 R6_2.4.0 rpart_4.1-15 nlme_3.1-141
[51] compiler_3.5.1


The following software was used for the example run (local on Windows PC)

R version 3.6.2 (2019-12-12)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 17134)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] pseudo_1.4.3   geepack_1.3-1  KMsurv_0.1-5   flexsurv_1.1.1 survRM2_1.0-3 
[6] survMisc_0.5.5 survival_3.2-3

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.7         pillar_1.4.4       compiler_3.6.2     RColorBrewer_1.1-2
 [5] tools_3.6.2        digest_0.6.25      nlme_3.1-142       evaluate_0.14     
 [9] lifecycle_1.0.0    tibble_3.0.1       gtable_0.3.0       lattice_0.20-38   
[13] pkgconfig_2.0.3    rlang_0.4.11       Matrix_1.4-1       yaml_2.2.1        
[17] mvtnorm_1.1-1      xfun_0.14          gridExtra_2.3      dplyr_1.0.0       
[21] knitr_1.28         generics_0.1.0     vctrs_0.3.8        grid_3.6.2        
[25] tidyselect_1.1.1   data.table_1.12.8  deSolve_1.28       mstate_0.2.12     
[29] glue_1.4.1         R6_2.5.0           km.ci_0.5-2        rmarkdown_2.3     
[33] tidyr_1.1.3        purrr_0.3.4        ggplot2_3.3.2      magrittr_2.0.1    
[37] backports_1.1.8    MASS_7.3-51.4      scales_1.1.1       ellipsis_0.3.1    
[41] htmltools_0.5.0    splines_3.6.2      xtable_1.8-4       colorspace_1.4-1  
[45] quadprog_1.5-8     muhaz_1.2.6.1      munsell_0.5.0      broom_0.5.6       
[49] crayon_1.4.1       zoo_1.8-8         