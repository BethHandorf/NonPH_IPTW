# Author: Elizabeth Handorf
# Supplemental R code for 
# "Analysis of survival data with non-proportional hazards: 
# A comparison of propensity score weighted methods"
# 2020-07-07

###### Code to simulate treatments and outcomes for all different scenarios
###### Each scenario will be simulated m=500 times, saving the output in one large file 
###### with iteration numbers for later extracting sets of treatments and outcomes
###### For manuscript - seeds 1-500 used for the 500 simulations


#Required inputs: simulated_covariates.csv which contains the matrix of covariates 


library(simsurv)
library(survival)


setwd("directory")

#Some code to extrac the inputs from the command line, modified from
#http://yangfeng.wordpress.com/2009/09/03/including-arguments-in-r-cmd-batch-mode/
##First read in the arguments listed at the command line
args=(commandArgs(TRUE))

##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  FNPH=1  #1=log-time; Note - separate code for the piecewise linear form
  SNPH=0.25
  TE = -0.67
  SS<-2500*2
  startSeed<-1
  endSeed<-5

}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

# #Base case
# #FNPH.LT_SNPH.S_TE.L_SS.L
# FNPH="LT"  #Note - separate code for the piecewise linear form
# SNPH=0.25
# TE=-log(2)
# SS<-2500*2

#Read in file of simulated covariates

sim.covariates<-read.csv(
  "simulated_covariates.csv")
sim.covariates<-sim.covariates[1:SS,]



#Create an ID variable that will be used for linkage
names(sim.covariates)[1]<-"ID"


####################################################################################################


##########################
#Set up treatment model

betas.cov<-c(.03,.15,.35,0.15,-0.8,.4,0.3,-0.4,-.1,-.15,-.2,-.25,-.15,-.25,-.3,-.2,-.25,-.3,-.3,-.35)
design.matrix<-model.matrix(~age+as.factor(charlson)+male+low.stage+high.grade+histology+white+hispanic+
                              as.factor(facility)+as.factor(income)+as.factor(education)+as.factor(insurance),
                            data=sim.covariates)

#set intercept s.t. probability of treatment = some value
#for simplicity, use the average values of X here
p<-0.5

xavg = apply(design.matrix[,2:21],2,mean)

alpha = -1*(log(((1/p)-1)/exp(-1*xavg%*%betas.cov)))

#set betas for propensity score model 
betas.ps<-c(alpha,betas.cov)

cbind(colnames(design.matrix),betas.ps)
#Calculate probability of treatment given model parameters and simulated data
logitp<-design.matrix%*%betas.ps
expit<-function(x) { exp(x)/(1+exp(x))}

prob.trt<-expit(logitp)

#Covariate effects for survival
design.matrix.surv<-design.matrix[,-1]
#Potential outcomes:
design.matrix.surv.PO.0<-cbind(trt=0,design.matrix.surv)
design.matrix.surv.PO.1<-cbind(trt=1,design.matrix.surv)


#Survival treatment effect sizes
betas.surv<-c(TE,.04,.3,.8,.3,-.4,.2,0.2,-.1,-.2,0,0,0,-.1,-.1,-.2,-.3,-.2,-.15,-.1,-.2)
names(betas.surv)<-colnames(design.matrix.surv.PO.1)

cbind(colnames(design.matrix.surv.PO.1),betas.surv)

#Initialize data frame to store all simulated datasets
sim.dat<-data.frame()

TE_pos<-abs(TE)

filename<-paste("FNPH.", FNPH, "_SNPH.", SNPH, "_TE.", round(TE_pos,2), 
                "_SS.", SS,"_seedRange.", startSeed, "-", endSeed,  ".csv", sep="")


colname.flag<-TRUE
for (i in startSeed:endSeed) {

  #Simulate both potential survival outcomes with seed
  set.seed(i)
  Time.PO.0 <- simsurv(lambdas = .1, gammas = 0.8, betas =betas.surv,
                            x = as.data.frame(design.matrix.surv.PO.0),
                            tde = c(trt = SNPH), tdefunction = "log", interval = c(1e-15, 1000))

  set.seed(i)
  Time.PO.1 <- simsurv(lambdas = .1, gammas = 0.8, betas =betas.surv,
                            x = as.data.frame(design.matrix.surv.PO.1),
                            tde = c(trt = SNPH), tdefunction = "log")


  #simulate censoring time
  censtime<-runif(dim(sim.covariates)[1],min=8, max=15)

  #simulate actual treatment received
  trt=sapply(prob.trt,FUN=rbinom,n=1,size=1)

  sim.dat.tmp<-cbind.data.frame(itter=i, ID = Time.PO.0$id, trt, timePO.0 = Time.PO.0$eventtime,
                              timePO.1=Time.PO.1$eventtime,censtime)
  #eventtime = true survival time under treatment received
  #time = observed time
  #fail = failure indicator
  sim.dat.tmp$eventtime<- trt*sim.dat.tmp$timePO.1 + (1-trt)*sim.dat.tmp$timePO.0
  sim.dat.tmp$time<-apply(cbind(sim.dat.tmp$eventtime,sim.dat.tmp$censtime),1,min)
  sim.dat.tmp$fail<-sim.dat.tmp$censtime>sim.dat.tmp$eventtime

  #Write out results
  write.table(sim.dat.tmp,filename, append=TRUE,sep=",", col.names=colname.flag)
  sim.dat<-rbind.data.frame(sim.dat, sim.dat.tmp)
  colname.flag<-FALSE

}




