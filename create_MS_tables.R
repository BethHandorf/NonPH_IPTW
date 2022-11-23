

#This code takes a result file, summarizes it, and pulls out the summarized results presented in the manuscript


setwd("/")


dict<-read.csv("Intermediate_results/Simulations_intermediate_results_data_dictionary.csv", header=TRUE)

filenames<-list(c("Intermediate_results/res.2y.FNPH.1_SNPH.0.25_TE.0.69_SS.5000.csv", #base
             "Intermediate_results/res.2y.FNPH.2_SNPH.0.25_TE.0.69_SS.5000.csv", #LPW
             "Intermediate_results/res.2y.FNPH.1_SNPH.0.125_TE.0.69_SS.5000.csv", #Modest Non-ph
             "Intermediate_results/res.2y.FNPH.1_SNPH.0.25_TE.0.41_SS.5000.csv", #Modest Treatment effect
             "Intermediate_results/res.2y.FNPH.1_SNPH.0.25_TE.0.69_SS.1000.csv", #small sample size
             "Intermediate_results/res.2y.FNPH.3_SNPH.0.25_TE.0.69Q-0.5_SS.5000.csv", #GG lambda=-0.5
             "Intermediate_results/res.2y.FNPH.3_SNPH.0.25_TE.0.69Q2.5_SS.5000.csv"), # GG lambda=2.5
             
             c("Intermediate_results/res.5y.FNPH.1_SNPH.0.25_TE.0.69_SS.5000.csv", #base
             "Intermediate_results/res.5y.FNPH.2_SNPH.0.25_TE.0.69_SS.5000.csv", #LPW
             "Intermediate_results/res.5y.FNPH.1_SNPH.0.125_TE.0.69_SS.5000.csv", #Modest Non-ph
             "Intermediate_results/res.5y.FNPH.1_SNPH.0.25_TE.0.41_SS.5000.csv", #Modest Treatment effect
             "Intermediate_results/res.5y.FNPH.1_SNPH.0.25_TE.0.69_SS.1000.csv", #small sample size
             "Intermediate_results/res.5y.FNPH.3_SNPH.0.25_TE.0.69Q-0.5_SS.5000.csv", #GG lambda=-0.5
             "Intermediate_results/res.5y.FNPH.3_SNPH.0.25_TE.0.69Q2.5_SS.5000.csv"), # GG lambda=2.5
             
             c("Intermediate_results/res.10y.FNPH.1_SNPH.0.25_TE.0.69_SS.5000.csv", #base
             "Intermediate_results/res.10y.FNPH.2_SNPH.0.25_TE.0.69_SS.5000.csv", #LPW
             "Intermediate_results/res.10y.FNPH.1_SNPH.0.125_TE.0.69_SS.5000.csv", #Modest Non-ph
             "Intermediate_results/res.10y.FNPH.1_SNPH.0.25_TE.0.41_SS.5000.csv", #Modest Treatment effect
             "Intermediate_results/res.10y.FNPH.1_SNPH.0.25_TE.0.69_SS.1000.csv", #small sample size
             "Intermediate_results/res.10y.FNPH.3_SNPH.0.25_TE.0.69Q-0.5_SS.5000.csv", #GG lambda=-0.5
             "Intermediate_results/res.10y.FNPH.3_SNPH.0.25_TE.0.69Q2.5_SS.5000.csv"), # GG lambda=2.5
             
             c("Intermediate_results/res.median.FNPH.1_SNPH.0.25_TE.0.69_SS.5000.csv", #base
             "Intermediate_results/res.median.FNPH.2_SNPH.0.25_TE.0.69_SS.5000.csv", #LPW
             "Intermediate_results/res.median.FNPH.1_SNPH.0.125_TE.0.69_SS.5000.csv", #Modest Non-ph
             "Intermediate_results/res.median.FNPH.1_SNPH.0.25_TE.0.41_SS.5000.csv", #Modest Treatment effect
             "Intermediate_results/res.median.FNPH.1_SNPH.0.25_TE.0.69_SS.1000.csv", #small sample size
             "Intermediate_results/res.median.FNPH.3_SNPH.0.25_TE.0.69Q-0.5_SS.5000.csv", #GG lambda=-0.5
             "Intermediate_results/res.median.FNPH.3_SNPH.0.25_TE.0.69Q2.5_SS.5000.csv"), # GG lambda=2.5
             
             c("Intermediate_results/res.RMS.FNPH.1_SNPH.0.25_TE.0.69_SS.5000.csv", #base
             "Intermediate_results/res.RMS.FNPH.2_SNPH.0.25_TE.0.69_SS.5000.csv", #LPW
             "Intermediate_results/res.RMS.FNPH.1_SNPH.0.125_TE.0.69_SS.5000.csv", #Modest Non-ph
             "Intermediate_results/res.RMS.FNPH.1_SNPH.0.25_TE.0.41_SS.5000.csv", #Modest Treatment effect
             "Intermediate_results/res.RMS.FNPH.1_SNPH.0.25_TE.0.69_SS.1000.csv", #small sample size
             "Intermediate_results/res.RMS.FNPH.3_SNPH.0.25_TE.0.69Q-0.5_SS.5000.csv", #GG lambda=-0.5
             "Intermediate_results/res.RMS.FNPH.3_SNPH.0.25_TE.0.69Q2.5_SS.5000.csv") # GG lambda=2.5
             
             
             )


bias.res<-matrix(rep(NA,56*5),ncol=5)
var.res<-matrix(rep(NA,56*5),ncol=5)
cov.res<-matrix(rep(NA,49*5),ncol=5)

for ( i in 1:5) {
  for (j in 1:7) {

    readfilename<-filenames[[i]][j]  
    
    #Read in results and header file
    resultstab<-read.csv(readfilename, header=FALSE)
    
    #First 2 columns of results tab are labels and seed numbers 
    resultstab<-resultstab[,-c(1:2)]
    
    names(resultstab)<-dict$Varname
    
    #Summarize results
    res.sum<-apply(resultstab,2,mean, na.rm=TRUE)
    
    #################################################################
    #Select elements of these summaries used in the tables and figures
    
    #Bias (Figure 1)
    #Note: the order of results is the order in the dataset "bias_results_graph_input.csv".  
    #This data is re-ordered for the publication
    
    bias.tmp<-res.sum[c("AFT.GG_bias",
                        "AFT.WBL.LS_bias",
                        "CTV.LPW_bias",
                        "CTV.LT_bias",
                        "NCM_bias",
                        "NUA_bias",
                        "PO_bias",
                        "WKM_bias")]
    
    bias.res[(8*(j-1)+1):(8*j),i]<-bias.tmp
    
    
    
    #Variance (Figure A1)
    #Note: the order of results is the order in the dataset "var_results_graph_input.csv".  
    
    var.tmp<-res.sum[c("AFT.GG_var",
                       "AFT.WBL.LS_var",
                       "CTV.LPW_var",
                       "CTV.LT_var",
                       "NCM_var",
                       "NUA_var",
                       "PO_var",
                       "WKM_var")]
    
    var.res[(8*(j-1)+1):(8*j),i]<-var.tmp
    
    #Coverage (Table 1, Table A2)
    #Note: the order here is the order in which the results appear in the manuscript
  
    cov.tmp<-res.sum[c("NCM_cov",
                       "CTV.LPW_cov",
                       "CTV.LT_cov",
                       "AFT.GG_cov",
                       "AFT.WBL.LS_cov",
                       "PO_cov",
                       "WKM_cov")]
    
    cov.tmp<-round(cov.tmp,2)
    
    cov.res[(7*(j-1)+1):(7*j),i]<-cov.tmp
    
  }
}


#Bias

colnames(bias.res)<-c("2y","5y","10y","median","RMS")

bias.res<-cbind.data.frame(
  Scenario=c(rep("Base",8),
             rep("LPW",8),
             rep("Modest NPH",8),
             rep("Modest TE",8),
             rep("Small SS",8),
             rep("GG lambda -0.5",8),
             rep("GG lambda 2.5",8)),
  
  Method=rep(c("AFT.GG",
               "AFT.WBL.LS",
               "CTV.LT",
               "CTV.LPW",
               "NCM",
               "NUA",
               "PO",
               "WKM"),7), bias.res)

write.csv(bias.res,"Tables_graphs/Fig1_input_bias.csv", row.names=FALSE)


#Variance

colnames(var.res)<-c("2y","5y","10y","median","RMS")

var.res<-cbind.data.frame(
  Scenario=c(rep("Base",8),
             rep("LPW",8),
             rep("Modest NPH",8),
             rep("Modest TE",8),
             rep("Small SS",8),
             rep("GG lambda -0.5",8),
             rep("GG lambda 2.5",8)),
  
  Method=rep(c("AFT.GG",
  "AFT.WBL.LS",
  "CTV.LPW",
  "CTV.LT",
  "NCM",
  "NUA",
  "PO",
  "WKM"),7), var.res)

write.csv(var.res,"Tables_graphs/FigA1_input_var.csv", row.names=FALSE)


#Coverage
colnames(cov.res)<-c("2y","5y","10y","median","RMS")

cov.res<-cbind.data.frame(
  Scenario=c(rep("Base",7),
             rep("PWC",7),
             rep("Modest NPH",7),
             rep("Modest TE",7),
             rep("Small SS",7),
             rep("GG lambda -0.5",7),
             rep("GG lambda 2.5",7)),
  
  Method=rep(c("Cox",
               "CTV LT",
               "CTV PWC",
               "AFT GG",
               "AFT WBL LS",
               "Pseudo",
               "Wtd KM"),7), cov.res)

write.csv(cov.res[1:7,],"Tables_graphs/Table1_cov.csv", row.names=FALSE)
write.csv(cov.res[8:49,],"Tables_graphs/TableA2_cov.csv", row.names=FALSE)