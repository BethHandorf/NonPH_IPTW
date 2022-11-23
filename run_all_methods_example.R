# Author: Elizabeth Handorf
# Supplemental R code for 
# "Analysis of survival data with non-proportional hazards: 
# A comparison of propensity score weighted methods"
# 2020-07-07

#Run all methods and get bootstrap standard errors for one simulated dataset 
#Required inputs: "dfile" which contains example simulated datasets (with seeds 1-5)
#                 "simulated_covariates.csv" which contains the matrix of covariates
#                 "execute_all_methods.R" which contains the R code to run each survival method

library(survMisc)
library(survRM2)
library(flexsurv)
library(pseudo)

start.time<-Sys.time()

#Some code to extract the inputs from the command line, modified from
#http://yangfeng.wordpress.com/2009/09/03/including-arguments-in-r-cmd-batch-mode/
##First read in the arguments listed at the command line
#args=(commandArgs(TRUE))

##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
# if(length(args)==0){
#   print("No arguments supplied.")
# 
#   dfile<-"FNPH.2_SNPH.0.25_TE.0.69_SS.5000.csv"
#   seed=101
#   setwd("local input dir")
#   
#   resdir<-"local results dir"
#   
# }else{
# 
#   dfile <- args[1]
#   seed <- as.double(args[2])
#   setwd("remote input dir")
#   resdir<-"remote results dir"
# }

dfile<-"FNPH.1_SNPH.0.25_TE.-0.69_SS.5000_seedRange.1-5.csv"
seed=5   #For example run, 1-5 will work for seed
setwd("..")
   
resdir<-"./"

#Change the directory as needed
setwd("./")
resdir<-getwd()


print(dfile)
print(seed)
getwd()

source("execute_all_methods.R")

covariates<-read.csv("simulated_covariates.csv")
dim(covariates)
names(covariates)[1]<-"ID.cov"

#Save results

res.vec.names<-c(
  "seed",
  "TruePO_Z0",
  "TruePO_Z1",
  "TruePO_diff",
  "NUA_Z0",
  "NUA_Z1",
  "NUA_diff",
  "NUA_bias",
  "NUA_cov",
  "NUA_var",
  "NUA_sig",
  "NCM_Z0",
  "NCM_Z1",
  "NCM_diff",
  "NCM_bias",
  "NCM_cov",
  "NCM_var",
  "NCM_sig",
  "WKM_Z0",
  "WKM_Z1",
  "WKM_diff",
  "WKM_bias",
  "WKM_cov",
  "WKM_var",
  "WKM_sig",
  "CTV.LT_Z0",
  "CTV.LT_Z1",
  "CTV.LT_diff",
  "CTV.LT_bias",
  "CTV.LT_cov",
  "CTV.LT_var",
  "CTV.LT_sig",
  "CTV.LPW_Z0",
  "CTV.LPW_Z1",
  "CTV.LPW_diff",
  "CTV.LPW_bias",
  "CTV.LPW_cov",
  "CTV.LPW_var",
  "CTV.LPW_sig",
  "AFT.GG_Z0",
  "AFT.GG_Z1",
  "AFT.GG_diff",
  "AFT.GG_bias",
  "AFT.GG_cov",
  "AFT.GG_var",
  "AFT.GG_sig",
  "AFT.WBL.LS_Z0",
  "AFT.WBL.LS_Z1",
  "AFT.WBL.LS_diff",
  "AFT.WBL.LS_bias",
  "AFT.WBL.LS_cov",
  "AFT.WBL.LS_var",
  "AFT.WBL.LS_sig",
  "PO_Z0",
  "PO_Z1",
  "PO_diff",
  "PO_bias",
  "PO_cov",
  "PO_var",
  "PO_sig"
)

setwd(resdir)

#Create empty vectors to save results
res.med.all<-rep(NA,60)
res.rms.all<-rep(NA,60)
res.2y.all<-rep(NA,60)
res.5y.all<-rep(NA,60)
res.10y.all<-rep(NA,60)

res.med.all[1]<-seed
res.rms.all[1]<-seed
res.2y.all[1]<-seed
res.5y.all[1]<-seed
res.10y.all[1]<-seed


names(res.med.all)<-names(res.rms.all)<-names(res.2y.all)<-
  names(res.5y.all)<-names(res.10y.all)<-res.vec.names

#Only read in the sample we are using for this itteration
dataj<-read.csv(dfile, nrows=5000, skip=5000*(seed-1)+1, header=FALSE)

#Rename the data files as needed
if(dim(dataj)[2]==10){
  names(dataj)<-c("X","itter","ID","trt","timePO.0","timePO.1","censtime","eventtime","time","fail")
}
if(dim(dataj)[2]==9){
  names(dataj)<-c("itter","ID","trt","timePO.0","timePO.1","censtime","eventtime","time","fail")
}


dim(dataj)
print(dataj[1:10,])
print(dataj[4990:5000,])
dataj<-merge(dataj, covariates, by.x="ID", by.y="ID.cov")
dim(dataj)

#Truth
#Survival curve - potential outcomes of untreated
last.obs.time<-min(max(c(dataj$time[dataj$trt==0&dataj$fail==1]),c(dataj$time[dataj$trt==1&dataj$fail==1])))
S.0<-survfit(Surv(timePO.0,rep(1,length(timePO.0)))~1,data=dataj)

#True median, restricted mean survival, and survival at specific timepoints
S.0.median<-as.numeric(summary(S.0)$table['median'])
S.0.RMS<-as.numeric(summary(S.0,rmean=last.obs.time)$table['*rmean'])

S.0.2y<-summary(S.0,times=c(2))$surv
S.0.5y<-summary(S.0,times=c(5))$surv
S.0.10y<-summary(S.0,times=c(10))$surv


#Survival curve - potential outcomes of treated
S.1<-survfit(Surv(timePO.1,rep(1,length(timePO.1)))~1,data=dataj)

#True median, restricted mean survival, and survival at specific timepoints
S.1.median<-as.numeric(summary(S.1)$table['median'])
S.1.RMS<-as.numeric(summary(S.1,rmean=last.obs.time)$table['*rmean'])

S.1.2y<-summary(S.1,times=c(2))$surv
S.1.5y<-summary(S.1,times=c(5))$surv
S.1.10y<-summary(S.1,times=c(10))$surv


true.Potential.outcomes<-cbind(c(S.0.median,S.0.RMS,S.0.2y,S.0.5y,S.0.10y),
                               c(S.1.median,S.1.RMS,S.1.2y,S.1.5y,S.1.10y))
rownames(true.Potential.outcomes)<-c("median","RMS","2y","5y","10y")
true.Potential.outcomes<-cbind(true.Potential.outcomes, true.Potential.outcomes[,2]-true.Potential.outcomes[,1])

colnames(true.Potential.outcomes)<-c("untreated","treated","diff")


#list of all methods used to generate estimates

#NUA=Naove unadjusted	
#NCM=Naove cox model (weighted no TV covariate)	
#WKM=weighted KM	
#CVT_LT=Cox model with linear effect of log-time	
#CTV_LPW=Cox model with linear piecewise effect
#AFT_GG=AFT generalized gamma	
#AFT_WBL_LS=AFT weibull with location and scale	
#PO=Pseudo-observations

methods<-c("NUA","NCM","WKM","CTV_LT","CTV_LPW","AFT_GG","AFT_WBL_LS","PO")

#Get the estimates for the full sample

start.time<-Sys.time()
estimates<-runAllMethods(dataj,d=c(1:5000))
end.time<-Sys.time()
print(end.time-start.time)


#number of bootstrap replicates
M<-500

#Save matrices of bootstrap output
d.med.boot<-d.rms.boot<-d.2y.boot<-d.5y.boot<-d.10y.boot<-matrix(rep(NA,M*length(methods)),nrow=M)
colnames(d.med.boot)<-colnames(d.rms.boot)<-colnames(d.2y.boot)<-colnames(d.5y.boot)<-colnames(d.10y.boot)<-methods

#Run bootstrap resampling
set.seed(seed)
for (j in 1:M) {
  res.tmp<-runAllMethods(data=dataj,d=sample(seq(1:dim(dataj)[1]),replace=TRUE))
  
  d.med.boot[j,]<-res.tmp[[1]][,3]
  d.rms.boot[j,]<-res.tmp[[2]][,3]
  d.2y.boot[j,]<-res.tmp[[3]][,3]
  d.5y.boot[j,]<-res.tmp[[4]][,3]
  d.10y.boot[j,]<-res.tmp[[5]][,3]
}

#bootstrap variance estimate
var.d.med<-apply(d.med.boot,2,var)
var.d.rms<-apply(d.rms.boot,2,var)
var.d.2y<-apply(d.2y.boot,2,var)
var.d.5y<-apply(d.5y.boot,2,var)
var.d.10y<-apply(d.10y.boot,2,var)

d.med.CL<-apply(d.med.boot, 2, quantile, probs = c(0.025, 0.975))
d.rms.CL<-apply(d.rms.boot, 2, quantile, probs = c(0.025, 0.975))
d.2y.CL<-apply(d.2y.boot, 2, quantile, probs = c(0.025, 0.975))
d.5y.CL<-apply(d.5y.boot, 2, quantile, probs = c(0.025, 0.975))
d.10y.CL<-apply(d.10y.boot, 2, quantile, probs = c(0.025, 0.975))

#Bias and coverage

bias.d.med<-estimates[[1]][,3]-true.Potential.outcomes[1,3]
bias.d.rms<-estimates[[2]][,3]-true.Potential.outcomes[2,3]
bias.d.2y<-estimates[[3]][,3]-true.Potential.outcomes[3,3]
bias.d.5y<-estimates[[4]][,3]-true.Potential.outcomes[4,3]
bias.d.10y<-estimates[[5]][,3]-true.Potential.outcomes[5,3]

#rbind(bias.d.med, bias.d.rms, bias.d.2y, bias.d.5y, bias.d.10y)

cov.d.med<-d.med.CL[1,]<=true.Potential.outcomes[1,3]&d.med.CL[2,]>=true.Potential.outcomes[1,3]
cov.d.rms<-d.rms.CL[1,]<=true.Potential.outcomes[2,3]&d.rms.CL[2,]>=true.Potential.outcomes[2,3]
cov.d.2y<-d.2y.CL[1,]<=true.Potential.outcomes[3,3]&d.2y.CL[2,]>=true.Potential.outcomes[3,3]
cov.d.5y<-d.5y.CL[1,]<=true.Potential.outcomes[4,3]&d.5y.CL[2,]>=true.Potential.outcomes[4,3]
cov.d.10y<-d.10y.CL[1,]<=true.Potential.outcomes[5,3]&d.10y.CL[2,]>=true.Potential.outcomes[5,3]

#Statistical significance (<>0)
ssig.d.med<-!(d.med.CL[1,]<=0&d.med.CL[2,]>=0)
ssig.d.rms<-!(d.rms.CL[1,]<=0&d.rms.CL[2,]>=0)
ssig.d.2y<-!(d.2y.CL[1,]<=0&d.2y.CL[2,]>=0)
ssig.d.5y<-!(d.5y.CL[1,]<=0&d.5y.CL[2,]>=0)
ssig.d.10y<-!(d.10y.CL[1,]<=0&d.10y.CL[2,]>=0)

#reformat the output for multiple simulations
#output 1 vector for each of the 5 outcomes measures

res.med.all[1]<-seed

res.med.all[2:60]<-c(true.Potential.outcomes["median",],
                     estimates[[1]][1,],
                     bias.d.med[1],
                     cov.d.med[1],
                     var.d.med[1],
                     ssig.d.med[1],
                     estimates[[1]][2,],
                     bias.d.med[2],
                     cov.d.med[2],
                     var.d.med[2],
                     ssig.d.med[2],
                     estimates[[1]][3,],
                     bias.d.med[3],
                     cov.d.med[3],
                     var.d.med[3],
                     ssig.d.med[3],
                     estimates[[1]][4,],
                     bias.d.med[4],
                     cov.d.med[4],
                     var.d.med[4],
                     ssig.d.med[4],
                     estimates[[1]][5,],
                     bias.d.med[5],
                     cov.d.med[5],
                     var.d.med[5],
                     ssig.d.med[5],
                     estimates[[1]][6,],
                     bias.d.med[6],
                     cov.d.med[6],
                     var.d.med[6],
                     ssig.d.med[6],
                     estimates[[1]][7,],
                     bias.d.med[7],
                     cov.d.med[7],
                     var.d.med[7],
                     ssig.d.med[7],
                     estimates[[1]][8,],
                     bias.d.med[8],
                     cov.d.med[8],
                     var.d.med[8],
                     ssig.d.med[8])


res.rms.all[2:60]<-c(true.Potential.outcomes["RMS",],
                     estimates[[2]][1,],
                     bias.d.rms[1],
                     cov.d.rms[1],
                     var.d.rms[1],
                     ssig.d.rms[1],
                     estimates[[1]][2,],
                     bias.d.rms[2],
                     cov.d.rms[2],
                     var.d.rms[2],
                     ssig.d.rms[2],
                     estimates[[2]][3,],
                     bias.d.rms[3],
                     cov.d.rms[3],
                     var.d.rms[3],
                     ssig.d.rms[3],
                     estimates[[2]][4,],
                     bias.d.rms[4],
                     cov.d.rms[4],
                     var.d.rms[4],
                     ssig.d.rms[4],
                     estimates[[2]][5,],
                     bias.d.rms[5],
                     cov.d.rms[5],
                     var.d.rms[5],
                     ssig.d.rms[5],
                     estimates[[2]][6,],
                     bias.d.rms[6],
                     cov.d.rms[6],
                     var.d.rms[6],
                     ssig.d.rms[6],
                     estimates[[2]][7,],
                     bias.d.rms[7],
                     cov.d.rms[7],
                     var.d.rms[7],
                     ssig.d.rms[7],
                     estimates[[2]][8,],
                     bias.d.rms[8],
                     cov.d.rms[8],
                     var.d.rms[8],
                     ssig.d.rms[8])


res.2y.all[2:60]<-c(true.Potential.outcomes["2y",],
                    estimates[[3]][1,],
                    bias.d.2y[1],
                    cov.d.2y[1],
                    var.d.2y[1],
                    ssig.d.2y[1],
                    estimates[[3]][2,],
                    bias.d.2y[2],
                    cov.d.2y[2],
                    var.d.2y[2],
                    ssig.d.2y[2],
                    estimates[[3]][3,],
                    bias.d.2y[3],
                    cov.d.2y[3],
                    var.d.2y[3],
                    ssig.d.2y[3],
                    estimates[[3]][4,],
                    bias.d.2y[4],
                    cov.d.2y[4],
                    var.d.2y[4],
                    ssig.d.2y[4],
                    estimates[[3]][5,],
                    bias.d.2y[5],
                    cov.d.2y[5],
                    var.d.2y[5],
                    ssig.d.2y[5],
                    estimates[[3]][6,],
                    bias.d.2y[6],
                    cov.d.2y[6],
                    var.d.2y[6],
                    ssig.d.2y[6],
                    estimates[[3]][7,],
                    bias.d.2y[7],
                    cov.d.2y[7],
                    var.d.2y[7],
                    ssig.d.2y[7],
                    estimates[[3]][8,],
                    bias.d.2y[8],
                    cov.d.2y[8],
                    var.d.2y[8],
                    ssig.d.2y[8])


res.5y.all[2:60]<-c(true.Potential.outcomes["5y",],
                    estimates[[4]][1,],
                    bias.d.5y[1],
                    cov.d.5y[1],
                    var.d.5y[1],
                    ssig.d.5y[1],
                    estimates[[4]][2,],
                    bias.d.5y[2],
                    cov.d.5y[2],
                    var.d.5y[2],
                    ssig.d.5y[2],
                    estimates[[4]][3,],
                    bias.d.5y[3],
                    cov.d.5y[3],
                    var.d.5y[3],
                    ssig.d.5y[3],
                    estimates[[4]][4,],
                    bias.d.5y[4],
                    cov.d.5y[4],
                    var.d.5y[4],
                    ssig.d.5y[4],
                    estimates[[4]][5,],
                    bias.d.5y[5],
                    cov.d.5y[5],
                    var.d.5y[5],
                    ssig.d.5y[5],
                    estimates[[4]][6,],
                    bias.d.5y[6],
                    cov.d.5y[6],
                    var.d.5y[6],
                    ssig.d.5y[6],
                    estimates[[4]][7,],
                    bias.d.5y[7],
                    cov.d.5y[7],
                    var.d.5y[7],
                    ssig.d.5y[7],
                    estimates[[4]][8,],
                    bias.d.5y[8],
                    cov.d.5y[8],
                    var.d.5y[8],
                    ssig.d.5y[8])


res.10y.all[2:60]<-c(true.Potential.outcomes["10y",],
                     estimates[[5]][1,],
                     bias.d.10y[1],
                     cov.d.10y[1],
                     var.d.10y[1],
                     ssig.d.10y[1],
                     estimates[[5]][2,],
                     bias.d.10y[2],
                     cov.d.10y[2],
                     var.d.10y[2],
                     ssig.d.10y[2],
                     estimates[[5]][3,],
                     bias.d.10y[3],
                     cov.d.10y[3],
                     var.d.10y[3],
                     ssig.d.10y[3],
                     estimates[[5]][4,],
                     bias.d.10y[4],
                     cov.d.10y[4],
                     var.d.10y[4],
                     ssig.d.10y[4],
                     estimates[[5]][5,],
                     bias.d.10y[5],
                     cov.d.10y[5],
                     var.d.10y[5],
                     ssig.d.10y[5],
                     estimates[[5]][6,],
                     bias.d.10y[6],
                     cov.d.10y[6],
                     var.d.10y[6],
                     ssig.d.10y[6],
                     estimates[[5]][7,],
                     bias.d.10y[7],
                     cov.d.10y[7],
                     var.d.10y[7],
                     ssig.d.10y[7],
                     estimates[[5]][8,],
                     bias.d.10y[8],
                     cov.d.10y[8],
                     var.d.10y[8],
                     ssig.d.10y[8])


#Append to results file

write.table(data.frame(t(res.med.all)),paste("res.median.",dfile,sep=""), sep=",", col.names=FALSE, append=TRUE)
write.table(data.frame(t(res.rms.all)),paste("res.RMS.",dfile,sep=""), sep=",", col.names=FALSE, append=TRUE)
write.table(data.frame(t(res.2y.all)),paste("res.2y.",dfile,sep=""), sep=",", col.names=FALSE, append=TRUE)
write.table(data.frame(t(res.5y.all)),paste("res.5y.",dfile,sep=""), sep=",", col.names=FALSE, append=TRUE)
write.table(data.frame(t(res.10y.all)),paste("res.10y.",dfile,sep=""), sep=",", col.names=FALSE, append=TRUE)

end.time<-Sys.time()
print(end.time-start.time)

