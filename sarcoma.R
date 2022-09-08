library(survival)
library(ggplot2)
#library(epicalc)
library(muhaz)
library(flexsurv)
library(pseudo)
library(twang)


setwd("./")

dat<-read.csv("R:\\NCDB\\sarcoma\\sarcoma 2012\\sarcomasurv.csv", header=T)

# survival plots
s.obj<-survfit(Surv(DX_LASTCONTACT_DEATH_MONTHS, PUF_VITAL_STATUS==0)~Chemo,data=dat)
survdiff(Surv(DX_LASTCONTACT_DEATH_MONTHS, PUF_VITAL_STATUS==0)~Chemo,data=dat)
s.obj

par(mfrow=c(1,2))
par(mar=c(4, 4, 4, 1) + 0.1)

plot(s.obj,col=c("firebrick4","navy"),xlab="Months",mark.time=F,lwd=2,
     main="Unadjusted ")
legend("topright",c("No chemo","Chemo"), lwd=2, lty=1, col=c("firebrick4","navy"))

##Proportional hazards assumption

plot(s.obj,col=c("firebrick4","navy"),xlab="Months",mark.time=F,lwd=2,
     main="Complementary log-log",  fun="cloglog")
legend("topleft",c("No chemo","Chemo"), lwd=2, lty=1, col=c("firebrick4","navy"))

died.ind<-as.numeric(dat$PUF_VITAL_STATUS==0)


#Assess the PH assumption

regression<-coxph(Surv(DX_LASTCONTACT_DEATH_MONTHS, PUF_VITAL_STATUS==0)~ Chemo, data=dat)

cox.zph(regression)


multiple.regression<-coxph(Surv(DX_LASTCONTACT_DEATH_MONTHS, PUF_VITAL_STATUS==0)~ Chemo +histgrp+
                             ageCat+SEX+RaceCat+ insurance+insurance+ income+ noHSD+ distance+ location+
                             charlson+  lowgrade+ sizeGrp+ YEAR_OF_DIAGNOSIS+ FACILITY_LOCATION_CD +Facility + site, data=dat)

cox.zph(multiple.regression)




#######Propensity score analysis for Chemo across all histologies

ps.mod<-glm(Chemo~ ageCat+histgrp+as.factor(SEX)+RaceCat+ insurance+income+ noHSD+ distance+ location+
              charlson+  as.factor(lowgrade)+ sizeGrp+ YEAR_OF_DIAGNOSIS+ as.factor(FACILITY_LOCATION_CD) +Facility + site, data=dat, family=binomial)
ps2<-predict(ps.mod,type="response")

dat$ps2<-ps2
p.chemo<-sum(dat$Chemo==1)/length(dat$Chemo)

dat$psweight2<-(dat$Chemo*p.chemo/ps2)+((1-dat$Chemo)*(1-p.chemo)/(1-ps2))

#Balance statistics 
datatmp<-dat
datatmp$ageCat<-as.factor(datatmp$ageCat)
datatmp$histgrp<-as.factor(datatmp$histgrp)
datatmp$SEX<-as.factor(datatmp$SEX)
datatmp$RaceCat<-as.factor(datatmp$RaceCat)
datatmp$insurance<-as.factor(datatmp$insurance)
datatmp$income<-as.factor(datatmp$income)
datatmp$noHSD<-as.factor(datatmp$noHSD)
datatmp$distance<-as.factor(datatmp$distance)
datatmp$location<-as.factor(datatmp$location)
datatmp$charlson<-as.factor(datatmp$charlson)
datatmp$lowgrade<-as.factor(datatmp$lowgrade)
datatmp$sizeGrp<-as.factor(datatmp$sizeGrp)
datatmp$YEAR_OF_DIAGNOSIS<-as.factor(datatmp$YEAR_OF_DIAGNOSIS)
datatmp$FACILITY_LOCATION_CD<-as.factor(datatmp$FACILITY_LOCATION_CD)
datatmp$Facility<-as.factor(datatmp$Facility)
datatmp$site<-as.factor(datatmp$site)



bal<-bal.stat(datatmp, vars=c("ageCat","histgrp","SEX","RaceCat","insurance","income","noHSD", "distance", "location",
                                     "charlson","lowgrade","sizeGrp", "YEAR_OF_DIAGNOSIS","FACILITY_LOCATION_CD","Facility", "site")
                   , treat.var="Chemo",
                   w.all=datatmp$psweight2,
                   #w.all=1, #check unweighted balance stats
                   sampw=1,
                   estimand="ATE",multinom=FALSE)

################### Create plots


#for IPTW paper
par(mfrow=c(1,1))
par(mar=c(5, 4, 4, 2) + 0.1)
s.obj<-survfit(Surv(DX_LASTCONTACT_DEATH_MONTHS, PUF_VITAL_STATUS==0)~Chemo,data=dat,weights=psweight2)
s.obj
plot(s.obj,col=c("firebrick4","navy"),xlab="Months",mark.time=F,lwd=2,
     main="All histologies", ylab="Survival")
legend("topright",c("No chemo (N=3883)","Chemo (N=1494)"), 
       lwd=2, lty=1, col=c("firebrick4","navy"))



#Graphs for pub

s.obj<-survfit(Surv(DX_LASTCONTACT_DEATH_MONTHS, PUF_VITAL_STATUS==0)~Chemo,data=dat)
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
plot(s.obj,lty=c(1,2),xlab="Months",mark.time=F,lwd=2,
     main="Unadjusted ", ylab="Survival")
legend("topright",c("No chemo","Chemo"), lwd=2, lty=c(1,2))

##Proportional hazards assumption

plot(s.obj,lty=c(1,2),xlab="log(Months)",ylab="log(-log(Survival))",mark.time=F,lwd=2,
     main="Complementary log-log",  fun="cloglog")
legend("topleft",c("No chemo","Chemo"), lwd=2, lty=c(1,2))

#Adjusted
s.obj<-survfit(Surv(DX_LASTCONTACT_DEATH_MONTHS, PUF_VITAL_STATUS==0)~Chemo,data=dat,weights=psweight2)
s.obj
plot(s.obj,lty=c(1,2),xlab="Months",mark.time=F,lwd=2,
     main="IPTW", ylab="Survival")
legend("topright",c("No chemo","Chemo"), 
       lwd=2, lty=c(1,2))


# Create dataset used for analysis by the IPTW method code
# file is modified to use relevant variables here
# remove extraneous variables and rename where necessary

dat.iptw.input<-subset(dat, select=c(Chemo, DX_LASTCONTACT_DEATH_MONTHS, PUF_VITAL_STATUS, ageCat,histgrp,SEX,RaceCat, insurance,income, noHSD, distance, location,
                                       charlson,  lowgrade, sizeGrp, YEAR_OF_DIAGNOSIS, FACILITY_LOCATION_CD ,Facility , site))
names(dat.iptw.input)[1]<-"trt"
names(dat.iptw.input)[3]<-"alive"
dat.iptw.input$fail<-as.numeric(!dat.iptw.input$alive)
dat.iptw.input$time<-dat.iptw.input$DX_LASTCONTACT_DEATH_MONTHS/12

setwd("./")

source("execute_all_sarcoma.R")

#maximum time
max(dat.iptw.input$time)

dat.iptw.input<-dat.iptw.input[dat.iptw.input$time >0,]
#Require time >3mo
#dat.iptw.input<-dat.iptw.input[dat.iptw.input$time >3/12,]
dat.iptw.input$ID<-c(1:dim(dat.iptw.input)[1])

start.time<-Sys.time()
estimates<-runAllMethods(dat.iptw.input,d=c(1:dim(dat.iptw.input)[1]))
end.time<-Sys.time()
print(end.time-start.time)

estimates

#Now create confidence limits

#number of bootstrap replicates
M<-500

#Save matrices of bootstrap output
d.med.boot<-d.rms.boot<-d.2y.boot<-d.5y.boot<-d.10y.boot<-matrix(rep(NA,M*8),nrow=M)
#colnames(d.med.boot)<-colnames(d.rms.boot)<-colnames(d.2y.boot)<-colnames(d.5y.boot)<-colnames(d.10y.boot)<-methods

set.seed(12456)

for (j in 1:M) {
  res.tmp<-runAllMethods(data=dat.iptw.input,d=sample(seq(1:dim(dat.iptw.input)[1]),replace=TRUE))
  
  d.med.boot[j,]<-res.tmp[[1]][,3]
  d.rms.boot[j,]<-res.tmp[[2]][,3]
  d.2y.boot[j,]<-res.tmp[[3]][,3]
  d.5y.boot[j,]<-res.tmp[[4]][,3]
  d.10y.boot[j,]<-res.tmp[[5]][,3]
}


d.med.CL<-apply(d.med.boot, 2, quantile, probs = c(0.025, 0.975))
d.rms.CL<-apply(d.rms.boot, 2, quantile, probs = c(0.025, 0.975))
d.2y.CL<-apply(d.2y.boot, 2, quantile, probs = c(0.025, 0.975))
d.5y.CL<-apply(d.5y.boot, 2, quantile, probs = c(0.025, 0.975))
#d.10y.CL<-apply(d.10y.boot, 2, quantile, probs = c(0.025, 0.975))
