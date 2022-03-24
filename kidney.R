library(survival)
library(ggplot2)
#library(epicalc)
library(muhaz)
library(flexsurv)
library(pseudo)


setwd("./")

dat<-read.csv("ps_data.csv",header=T)

#Just use T1a disease
dat<-dat[dat$Stg_grp=="T1a",]

dat$surgery<-as.numeric(dat$surgery=="PN")

#limit to 75+
dat<-dat[dat$AGE>50&dat$AGE<=60,]

#limit to 10 years of follow-up
dat$DX_LASTCONTACT_DEATH_MONTHS[dat$DX_LASTCONTACT_DEATH_MONTHS>120]<-120
dat$PUF_VITAL_STATUS[dat$DX_LASTCONTACT_DEATH_MONTHS==120]<-1


########################################

# survival plots
s.obj<-survfit(Surv(DX_LASTCONTACT_DEATH_MONTHS, PUF_VITAL_STATUS==0)~surgery,data=dat)
survdiff(Surv(DX_LASTCONTACT_DEATH_MONTHS, PUF_VITAL_STATUS==0)~surgery,data=dat)
s.obj

par(mfrow=c(1,2))
par(mar=c(4, 4, 4, 1) + 0.1)

plot(s.obj,col=c("firebrick4","navy"),xlab="Months",mark.time=F,lwd=2,
     main="Unadjusted ")
legend("bottomleft",c("Radical","Partial"), lwd=2, lty=1, col=c("firebrick4","navy"))


##Proportional hazards assumption

plot(s.obj,col=c("firebrick4","navy"),xlab="Months",mark.time=F,lwd=2,
     main="Complementary log-log",  fun="cloglog")
legend("topleft",c("Radical","Partial"), lwd=2, lty=1, col=c("firebrick4","navy"))

died.ind<-as.numeric(dat$PUF_VITAL_STATUS==0)


reg<-coxph(Surv(DX_LASTCONTACT_DEATH_MONTHS, PUF_VITAL_STATUS==0)~ surgery, data=dat)

cox.zph(reg)




#######Propensity score analysis for surgery 

ps.mod<-glm(surgery~ as.factor(FACILITY_TYPE_CD)+as.factor(FACILITY_LOCATION_CD)+
              as.factor(SEX)+as.factor(INSURANCE_STATUS)+as.factor(MED_INC_QUAR_12)+
              as.factor(NO_HSD_QUAR_12)+as.factor(CDCC_TOTAL)+
              as.factor(YEAR_OF_DIAGNOSIS)+
              hist_grp+raceCat+hispanic+urban+Grade_cat, data=dat, family=binomial)
ps2<-predict(ps.mod,type="response")

dat$ps2<-ps2
p.surgery<-sum(dat$surgery==1)/length(dat$surgery)

dat$psweight2<-(dat$surgery*p.surgery/ps2)+((1-dat$surgery)*(1-p.surgery)/(1-ps2))

#Trim dataset by weight <.1 and >.9

dat<-dat[dat$ps2>.1 & dat$ps2 <.9,]

write.csv(dat)

################### Create Figure 3

#Graphs for pub
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))

# survival plots
s.obj<-survfit(Surv(DX_LASTCONTACT_DEATH_MONTHS, PUF_VITAL_STATUS==0)~surgery,data=dat)
plot(s.obj,lty=c(1,2),mark.time=F,lwd=2,xlab="Months", ylab="Survival",
     main="Unadjusted ")
legend("bottomleft",c("Radical","Partial"), lwd=2, lty=c(1,2))


##Proportional hazards assumption

plot(s.obj,lty=c(1,2),xlab="log(Months)",ylab="log(-log(Survival))",mark.time=F,lwd=2,
     main="Complementary log-log",  fun="cloglog")
legend("topleft",c("Radical","Partial"), lwd=2, lty=c(1,2))

#Adjusted
s.obj<-survfit(Surv(DX_LASTCONTACT_DEATH_MONTHS, PUF_VITAL_STATUS==0)~surgery,
               data=dat,weights=psweight2)
s.obj
plot(s.obj,lty=c(1,2),xlab="Months",mark.time=F,lwd=2,
     main="IPTW", ylab="Survival")
legend("bottomleft",c("Radical ","Partial"), 
       lwd=2, lty=c(1,2))


# Create dataset used for analysis by the IPTW method code
# file is modified to use relevant variables here
# remove extraneous variables and rename where necessary

dat.iptw.input<-subset(dat, select=c(surgery, DX_LASTCONTACT_DEATH_MONTHS, PUF_VITAL_STATUS, hist_grp,SEX,raceCat, INSURANCE_STATUS, NO_HSD_QUAR_12,
                                     MED_INC_QUAR_12,FACILITY_TYPE_CD,FACILITY_LOCATION_CD, CDCC_TOTAL,
                                     YEAR_OF_DIAGNOSIS,hispanic,urban,Grade_cat ))
names(dat.iptw.input)[1]<-"trt"
names(dat.iptw.input)[3]<-"alive"
dat.iptw.input$fail<-as.numeric(!dat.iptw.input$alive)
dat.iptw.input$time<-dat.iptw.input$DX_LASTCONTACT_DEATH_MONTHS/12



source("execute_all_kidney.R")

#maximum time
max(dat.iptw.input$time)
dat.iptw.input<-dat.iptw.input[!is.na(dat.iptw.input$time),]
dat.iptw.input<-dat.iptw.input[dat.iptw.input$time >0,]
#Require time >3mo
#dat.iptw.input<-dat.iptw.input[dat.iptw.input$time >3/12,]
dat.iptw.input$ID<-c(1:dim(dat.iptw.input)[1])

start.time<-Sys.time()
estimates<-run_Cox_KM(dat.iptw.input,d=c(1:dim(dat.iptw.input)[1]))
end.time<-Sys.time()
print(end.time-start.time)

estimates

#Bootstrap to get confidence limits


#Now create confidence limits

#number of bootstrap replicates
M<-500

#Save matrices of bootstrap output
d.med.boot<-d.rms.boot<-d.2y.boot<-d.5y.boot<-d.10y.boot<-matrix(rep(NA,M*9),nrow=M)
#colnames(d.med.boot)<-colnames(d.rms.boot)<-colnames(d.2y.boot)<-colnames(d.5y.boot)<-colnames(d.10y.boot)<-methods


set.seed(12456)

for (j in 1:M) {
        res.tmp<-run_Cox_KM(data=dat.iptw.input,d=sample(seq(1:dim(dat.iptw.input)[1]),replace=TRUE))
        
       #d.med.boot[j,]<-res.tmp[[1]][,3]
        d.rms.boot[j,]<-res.tmp[[2]][,3]
        d.2y.boot[j,]<-res.tmp[[3]][,3]
        d.5y.boot[j,]<-res.tmp[[4]][,3]
        d.10y.boot[j,]<-res.tmp[[5]][,3]
}



# Confidence limits - read in file

#d.med.CL<-apply(d.med.boot, 2, quantile, probs = c(0.025, 0.975),na.rm=TRUE)
d.rms.CL<-apply(d.rms.boot, 2, quantile, probs = c(0.025, 0.975),na.rm=TRUE)
d.2y.CL<-apply(d.2y.boot, 2, quantile, probs = c(0.025, 0.975),na.rm=TRUE)
d.5y.CL<-apply(d.5y.boot, 2, quantile, probs = c(0.025, 0.975),na.rm=TRUE)
d.10y.CL<-apply(d.10y.boot, 2, quantile, probs = c(0.025, 0.975))

cbind(estimates[[2]],t(d.rms.CL))
cbind(estimates[[3]],t(d.2y.CL))
cbind(estimates[[4]],t(d.5y.CL))
