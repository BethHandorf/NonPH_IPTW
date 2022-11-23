library(survival)
library(ggplot2)
#library(epicalc)
library(muhaz)
library(flexsurv)
library(pseudo)
library(twang)

setwd("./")

dat<-read.csv("kidney_sim_set.csv",header=T)

#Just use T1a disease
#dat<-dat[dat$Stg_grp=="T1a",]

#dat$surgery<-as.numeric(dat$surgery=="PN")

#limit to 75+
#dat<-dat[dat$AGE>50&dat$AGE<=60,]

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

dat<-dat[!is.na(dat$DX_LASTCONTACT_DEATH_MONTHS),]

#dat<-dat[dat$DX_LASTCONTACT_DEATH_MONTHS>0,]
#died.ind<-as.numeric(dat$PUF_VITAL_STATUS==0)
#dat$time<-dat$DX_LASTCONTACT_DEATH_MONTHS/12
#aft.wbl.mod.L.S<-flexsurvreg(Surv(time,died.ind)~as.factor(surgery),anc = list(shape = ~ as.factor(surgery)),
#                             data=dat,  dist="weibull")


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

#Balance statistics
datatmp<-dat

datatmp$FACILITY_LOCATION_CD<-as.factor(datatmp$FACILITY_LOCATION_CD)
datatmp$FACILITY_TYPE_CD<-as.factor(datatmp$FACILITY_TYPE_CD)
datatmp$SEX<-as.factor(datatmp$SEX)
datatmp$INSURANCE_STATUS<-as.factor(datatmp$INSURANCE_STATUS)
datatmp$MED_INC_QUAR_12<-as.factor(datatmp$MED_INC_QUAR_12)
datatmp$NO_HSD_QUAR_12<-as.factor(datatmp$NO_HSD_QUAR_12)
datatmp$CDCC_TOTAL<-as.factor(datatmp$CDCC_TOTAL)
datatmp$YEAR_OF_DIAGNOSIS<-as.factor(datatmp$YEAR_OF_DIAGNOSIS)
datatmp$hist_grp<-as.factor(datatmp$hist_grp)
datatmp$raceCat<-as.factor(datatmp$raceCat)
datatmp$hispanic<-as.factor(datatmp$hispanic)
datatmp$urban<-as.factor(datatmp$urban)
datatmp$Grade_cat<-as.factor(datatmp$Grade_cat)



bal<-bal.stat(datatmp, vars=c("FACILITY_TYPE_CD","FACILITY_LOCATION_CD","SEX","INSURANCE_STATUS","MED_INC_QUAR_12",
"NO_HSD_QUAR_12","CDCC_TOTAL","YEAR_OF_DIAGNOSIS","hist_grp","raceCat","hispanic","urban","Grade_cat")
              , treat.var="surgery",
              w.all=datatmp$psweight2,
              #w.all=1, #check balance if unweighted
              sampw=1,
              estimand="ATE",multinom=FALSE)

bal

################### Create Figure 3

png("Tables_graphs/Figure_3_mock_results.png")

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


dev.off()

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
#In the paper this used 500.  Using 50 here to make the demo run in a reasonable time
M<-50

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

sum.rms<-cbind(estimates[[2]],t(d.rms.CL))
sum.2y<-cbind(estimates[[3]],t(d.2y.CL))
sum.5y<-cbind(estimates[[4]],t(d.5y.CL))

#Create Table 3
Table3<-cbind(sum.rms[c(2,4,5,6,7,3,9),3:5],
      sum.2y[c(2,4,5,6,7,3,9),3:5],
      sum.5y[c(2,4,5,6,7,3,9),3:5])

colnames(Table3)<-c("RMS Delta","RMS LCL","RMS UCL","2y Delta","2y LCL","2y UCL","5y Delta","5y LCL","5y UCL")
rownames(Table3)<-c("Cox","TV Cox: log-T", "TV Cox: PWC", "AFT Gen Gamma", "AFT Wbl TV shape", "Weighted K-M", "Pseudo-obs")

Table3<-as.data.frame(Table3)
write.csv(Table3, "Tables_graphs/Table3_mock_results.csv")
