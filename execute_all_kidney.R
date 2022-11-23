


run_Cox_KM<-function(data,d) {
  methods<-c("NUA","NCM","WKM","CTV_LT","CTV_LPW","AFT_GG","AFT_WBL_LS","PO","PO_adj_wt")
  
  dat<-data[d,]
  last.obs.time<-min(c(max(dat$time[dat$trt==0&dat$fail==1],na.rm=TRUE),
                       max(dat$time[dat$trt==1&dat$fail==1],na.rm=TRUE)))
  
  #Initalize tables of results of interest 
  res.median<-res.RMS<-res.2y<-res.5y<-res.10y<-matrix(rep(NA,2*length(methods)),ncol=2)  
  
  ############# Estimate propensity score weights #########
  
  ps.mod<-glm(trt~ as.factor(FACILITY_TYPE_CD)+as.factor(FACILITY_LOCATION_CD)+
                as.factor(SEX)+as.factor(INSURANCE_STATUS)+as.factor(MED_INC_QUAR_12)+
                as.factor(NO_HSD_QUAR_12)+as.factor(CDCC_TOTAL)+
                as.factor(YEAR_OF_DIAGNOSIS)+
                hist_grp+raceCat+hispanic+urban+Grade_cat, data=dat, family=binomial)
  ps2<-predict(ps.mod,type="response")
  
  dat$ps2<-ps2
  p.surgery<-sum(dat$trt==1)/length(dat$trt)
  
  dat$psweight2<-(dat$trt*p.surgery/ps2)+((1-dat$trt)*(1-p.surgery)/(1-ps2))
  
  
 
  
  #trim
  dat<-dat[dat$ps2>.1 & dat$ps2<0.9,]
  
  
  # ########## 1. Naive (unadjusted) results from observed data only
  # ########################
  # #start.time<-Sys.time()
  # 
  # i<-match("NUA",methods) #Index number of this method
  # #Unadjusted observed survival curves
  # S.obs<-survfit(Surv(time,fail)~trt, data=dat)
  # 
  # #medians
  # res.median[i,]<-as.numeric(summary(S.obs)$table[,'median'])
  # 
  # #RMS
  # res.RMS[i,]<-as.numeric(summary(S.obs,rmean=last.obs.time)$table[,'*rmean'])
  # 
  # #Survival
  # surv.prs<-matrix(as.numeric(summary(S.obs,times=c(2,5,10))$surv),nrow=3)
  # res.2y[i,]<-surv.prs[1,]
  # res.5y[i,]<-surv.prs[2,]
  # res.10y[i,]<-surv.prs[3,]
  
  #end.time<-Sys.time()
  #print(end.time-start.time)
  
  
  ############## 2. Naive cox model #########################
  #start.time<-Sys.time()
  
  i<-match("NCM",methods) #Index number of this method
  NCM.mod<-coxph(Surv(time,fail)~trt, weights=psweight2, data=dat, robust=TRUE)
  
  s.NCM.trt0<-survfit(NCM.mod,newdata= data.frame(trt=0))
  s.NCM.trt1<-survfit(NCM.mod,newdata= data.frame(trt=1))
  
  #medians
  res.median[i,1]<-as.numeric(summary(s.NCM.trt0)$table['median'])
  res.median[i,2]<-as.numeric(summary(s.NCM.trt1)$table['median'])
  #RMS
  res.RMS[i,1]<-as.numeric(summary(s.NCM.trt0,rmean=last.obs.time)$table['*rmean'])
  res.RMS[i,2]<-as.numeric(summary(s.NCM.trt1,rmean=last.obs.time)$table['*rmean'])
  
  #Survival
  surv0.prs<-matrix(as.numeric(summary(s.NCM.trt0,times=c(2,5,10))$surv),nrow=3)
  surv1.prs<-matrix(as.numeric(summary(s.NCM.trt1,times=c(2,5,10))$surv),nrow=3)
  
  res.2y[i,1]<-surv0.prs[1,]
  res.5y[i,1]<-surv0.prs[2,]
  res.10y[i,1]<-surv0.prs[3,]
  res.2y[i,2]<-surv1.prs[1,]
  res.5y[i,2]<-surv1.prs[2,]
  res.10y[i,2]<-surv1.prs[3,]
  
  #end.time<-Sys.time()
  #print(end.time-start.time)
  
  
  
  ############## 3. Weighted Kaplan-Meier curves
  #start.time<-Sys.time()
  
  
  i<-match("WKM",methods) #Index number of this method
  
  S.WKM<-survfit(Surv(time,fail)~trt, data=dat, weights=psweight2)
  
  #medians
  res.median[i,]<-as.numeric(summary(S.WKM)$table[,'median'])
  
  #RMS
  res.RMS[i,]<-as.numeric(summary(S.WKM,rmean=last.obs.time)$table[,'*rmean'])
  
  #Survival
  
  surv.prs<-matrix(as.numeric(summary(S.WKM,times=c(2,5,10))$surv),nrow=3)
  res.2y[i,]<-surv.prs[1,]
  res.5y[i,]<-surv.prs[2,]
  res.10y[i,]<-surv.prs[3,]
  
  #end.time<-Sys.time()
  #print(end.time-start.time)
  
  
  ################## 4. Cox model with TV effect of log-time
  #start.time<-Sys.time()
  
  i<-match("CTV_LT",methods) #Index number of this method
  #Split the dataset at the failure times
  #cut.points <- unique(dat$time[dat$fail == 1])
  
  cut.points<-seq(from=0,to=11,by=1/12)
  SURV2 <- survSplit(data = dat, cut = cut.points, end = "time",
                     start = "time0", event = "fail")
  SURV2$trtLT <-SURV2$trt*log(SURV2$time)
  
  #Remove very small time differences - these cause errors
  tdiff<-(SURV2$time-SURV2$time0)
  SURV2<-SURV2[tdiff>10^-5,]
  CTV.mod<-coxph(Surv(time0,time,fail)~trt+trtLT,data=SURV2, weights=psweight2)
  
  #Fitted curves from this cox model
  #Get a list of all the time intervals (from the pt with longest follow-up)
  last <- SURV2$ID[which.max(SURV2$time)]
  intervals <- SURV2[SURV2$ID == last, c("time0", "time", "fail")]
  
  #curve for control
  covs<-data.frame(trt = 0, intervals)
  covs$trtLT <- covs$trt * log(covs$time)
  
  s.CTV.trt0<-survfit(CTV.mod, newdata = covs, individual = TRUE)
  
  #redo for treated
  covs<-data.frame(trt = 1, intervals)
  covs$trtLT <- covs$trt * log(covs$time)
  
  s.CTV.trt1<-survfit(CTV.mod, newdata = covs, individual = TRUE)
  
  #medians
  res.median[i,1]<-as.numeric(summary(s.CTV.trt0)$table['median'])
  res.median[i,2]<-as.numeric(summary(s.CTV.trt1)$table['median'])
  #RMS
  res.RMS[i,1]<-as.numeric(summary(s.CTV.trt0,rmean=last.obs.time)$table['*rmean'])
  res.RMS[i,2]<-as.numeric(summary(s.CTV.trt1,rmean=last.obs.time)$table['*rmean'])
  #Survival
  surv0.prs<-matrix(as.numeric(summary(s.CTV.trt0,times=c(2,5,10))$surv),nrow=3)
  surv1.prs<-matrix(as.numeric(summary(s.CTV.trt1,times=c(2,5,10))$surv),nrow=3)
  
  res.2y[i,1]<-surv0.prs[1,]
  res.5y[i,1]<-surv0.prs[2,]
  res.10y[i,1]<-surv0.prs[3,]
  res.2y[i,2]<-surv1.prs[1,]
  res.5y[i,2]<-surv1.prs[2,]
  res.10y[i,2]<-surv1.prs[3,]
  
  #end.time<-Sys.time()
  #print(end.time-start.time)
  
  
  # ################## 5. Cox model with piecewise TV effect
  # #start.time<-Sys.time()
  # 
  i<-match("CTV_LPW",methods) #Index number of this method
  #Split the dataset at times 2 and 5
  cut.points <- c(2,5)
  SURV2 <- survSplit(data = dat, cut = cut.points, end = "time",
                     start = "time0", event = "fail")

  #Remove very small time differences - these cause errors
  tdiff<-(SURV2$time-SURV2$time0)
  SURV2<-SURV2[tdiff>10^-5,]
  SURV2$trt_2y <-SURV2$trt*as.numeric(SURV2$time0==2)
  SURV2$trt_5y <-SURV2$trt*as.numeric(SURV2$time0==5)
  CTV.mod.LPW<-coxph(Surv(time0,time,fail)~trt+trt_2y+trt_5y,data=SURV2, weights=psweight2)

  #Fitted curves from this cox model

  #Get a list of all the time intervals (from the pt with longest follow-up)
  last <- SURV2$ID[which.max(SURV2$time)]
  intervals <- SURV2[SURV2$ID == last, c("time0", "time", "fail")]

  #curve for control
  covs<-data.frame(trt = 0, intervals)
  covs$trt_2y <- covs$trt * as.numeric(covs$time0==2)
  covs$trt_5y <- covs$trt * as.numeric(covs$time0==5)
  s.CTV.trt0<-survfit(CTV.mod.LPW, newdata = covs, individual = TRUE)

  #redo for treated
  covs<-data.frame(trt = 1, intervals)
  covs$trt_2y <- covs$trt * as.numeric(covs$time0==2)
  covs$trt_5y <- covs$trt * as.numeric(covs$time0==5)
  s.CTV.trt1<-survfit(CTV.mod.LPW, newdata = covs, individual = TRUE)

  #medians
  res.median[i,1]<-as.numeric(summary(s.CTV.trt0)$table['median'])
  res.median[i,2]<-as.numeric(summary(s.CTV.trt1)$table['median'])
  #RMS
  res.RMS[i,1]<-as.numeric(summary(s.CTV.trt0,rmean=last.obs.time)$table['*rmean'])
  res.RMS[i,2]<-as.numeric(summary(s.CTV.trt1,rmean=last.obs.time)$table['*rmean'])
  #Survival
  surv0.prs<-matrix(as.numeric(summary(s.CTV.trt0,times=c(2,5,10))$surv),nrow=3)
  surv1.prs<-matrix(as.numeric(summary(s.CTV.trt1,times=c(2,5,10))$surv),nrow=3)

  res.2y[i,1]<-surv0.prs[1,]
  res.5y[i,1]<-surv0.prs[2,]
  res.10y[i,1]<-surv0.prs[3,]
  res.2y[i,2]<-surv1.prs[1,]
  res.5y[i,2]<-surv1.prs[2,]
  res.10y[i,2]<-surv1.prs[3,]
  # 
  # #end.time<-Sys.time()
  # #print(end.time-start.time)
  # 
  # 
  # ################## 6. parametric AFT model - Generalized Gamma

  #start.time<-Sys.time()
  i<-match("AFT_GG",methods) #Index number of this method
  aft.ggam.mod<-flexsurvreg(Surv(time,fail)~as.factor(trt), data=dat, weights = psweight2, dist="gengamma")

  #Median
  medians<-summary(aft.ggam.mod, fn = median.ggam, t = 1, B = 10000) #Helper function defined above
  res.median[i,1]<-medians[['as.factor(trt)=0']]$est
  res.median[i,2]<-medians[['as.factor(trt)=1']]$est

  #Mean restricted survival time
  rmst<-summary(aft.ggam.mod, fn = rmst_gengamma, t = last.obs.time, B = 10000)
  res.RMS[i,1]<-rmst[['as.factor(trt)=0']]$est
  res.RMS[i,2]<-rmst[['as.factor(trt)=1']]$est

  #predicted survival function
  sum.ggam<-summary(aft.ggam.mod, t=c(2,5,10))
  res.2y[i,1]<-sum.ggam[['as.factor(trt)=0']]$est[1]
  res.5y[i,1]<-sum.ggam[['as.factor(trt)=0']]$est[2]
  res.10y[i,1]<-sum.ggam[['as.factor(trt)=0']]$est[3]
  res.2y[i,2]<-sum.ggam[['as.factor(trt)=1']]$est[1]
  res.5y[i,2]<-sum.ggam[['as.factor(trt)=1']]$est[2]
  res.10y[i,2]<-sum.ggam[['as.factor(trt)=1']]$est[3]

  #end.time<-Sys.time()
  #print(end.time-start.time)


  ################## 7. parametric AFT model - Weibull allowing location and scale to vary
  #start.time<-Sys.time()

  i<-match("AFT_WBL_LS",methods) #Index number of this method

  #wt.tmp<-rep(1,3774)
  wt.tmp<-as.numeric(dat$psweight2)

  aft.wbl.mod.L.S<-flexsurvreg(Surv(time,fail)~as.factor(trt), anc = list(shape = ~ as.factor(trt)),
                               data=dat, weights= wt.tmp,
                               dist="weibull", inits=c(0.034, 3.75, -0.126,0.218 ))


  #aft.wbl.mod.L.S<-flexsurvreg(Surv(time,fail)~as.factor(dat$trt),anc = list(shape = ~ as.factor(dat$trt)),
  #                             data=dat, weights = ps.IPTW, dist="weibull")



  medians<-summary(aft.wbl.mod.L.S, fn = median.weibull, t = 1, B = 10000)
  res.median[i,1]<-medians[['as.factor(trt)=0']]$est
  res.median[i,2]<-medians[['as.factor(trt)=1']]$est

  #Mean restricted survival time
  rmst<-summary(aft.wbl.mod.L.S, fn = rmst_weibull, t = last.obs.time, B = 10000)
  res.RMS[i,1]<-rmst[['as.factor(trt)=0']]$est
  res.RMS[i,2]<-rmst[['as.factor(trt)=1']]$est

  sum.wbl<-summary(aft.wbl.mod.L.S, t=c(2,5,10))
  res.2y[i,1]<-sum.wbl[['as.factor(trt)=0']]$est[1]
  res.5y[i,1]<-sum.wbl[['as.factor(trt)=0']]$est[2]
  res.10y[i,1]<-sum.wbl[['as.factor(trt)=0']]$est[3]
  res.2y[i,2]<-sum.wbl[['as.factor(trt)=1']]$est[1]
  res.5y[i,2]<-sum.wbl[['as.factor(trt)=1']]$est[2]
  res.10y[i,2]<-sum.wbl[['as.factor(trt)=1']]$est[3]

  #end.time<-Sys.time()
  #print(end.time-start.time)


  ################## 8. Pseudo-observations method
  #start.time<-Sys.time()

  # i<-match("PO",methods) #Index number of this method
  # 
  # #Old correction for variance stabilized weights
  # #Centering the weights at 1 makes results more consistent
  # #weight.correction<-as.numeric(dat$trt==1)/mean(dat$ps.IPTW[dat$trt==1]) +
  # #  as.numeric(dat$trt==0)/mean(dat$ps.IPTW[dat$trt==0])
  # #  dat$ps.IPTW.PO<-dat$ps.IPTW*weight.correction
  # 
  # dat$ps.IPTW.PO<-dat$psweight2
  # 
  # #Find median survival - need to use the whole curve
  # 
  # #pseudo.obs<-pseudosurv(dat$time, dat$fail) #,tmax=cutoffs)
  # cutoffs<-seq(from=1/12,to=last.obs.time,by=1/12)
  # cutoffs.obs<-cutoffs[cutoffs>min(dat$time[dat$fail==1]) & cutoffs<max(dat$time[dat$fail==1])]
  # pseudo.obs<-pseudosurv(dat$time, dat$fail,tmax=cutoffs.obs)
  # 
  # 
  # #get the pseudo-observations
  # pseudo.obs.mat<-as.matrix(pseudo.obs$pseudo)
  # #Weight them by the propensity score weights
  # pseudo.obs.mat.wt<-diag(as.numeric(dat$ps.IPTW.PO))%*%pseudo.obs.mat
  # #Survival at each time is the average of the weighted pseudo observations
  # S_t_0<-colSums(pseudo.obs.mat.wt[dat$trt==0,])/dim(dat)[1]
  # S_t_1<-colSums(pseudo.obs.mat.wt[dat$trt==1,])/dim(dat)[1]
  # #Find the first time where the survival curve is <0.5
  # res.median[i,1]<-min(pseudo.obs$time[S_t_0<0.5])
  # res.median[i,2]<-min(pseudo.obs$time[S_t_1<0.5])
  # 
  # #Pseudo mean - built in function
  # pseudo.RMS = pseudomean(dat$time, dat$fail,tmax=last.obs.time)
  # 
  # res.RMS[i,1]<-sum(pseudo.RMS[dat$trt==0]*dat$ps.IPTW.PO[dat$trt==0])/dim(dat)[1]
  # res.RMS[i,2]<-sum(pseudo.RMS[dat$trt==1]*dat$ps.IPTW.PO[dat$trt==1])/dim(dat)[1]
  # 
  # 
  # #Survival at specific timepoints
  # #cutoffs <- c(1:14)
  # #pseudo.obs<-pseudosurv(dat$time, dat$fail,tmax=cutoffs)
  # 
  # pseudo.obs.df<-data.frame(pseudo.obs$pseudo)
  # 
  # #Directly estimate the probability of surviving to time t
  # #This is just the mean of the pseudo observations
  # #Applying propensity score weights
  # 
  # res.2y[i,1]<-sum(pseudo.obs.df$time.2[dat$trt==0]*dat$ps.IPTW.PO[dat$trt==0])/dim(dat)[1]
  # res.2y[i,2]<-sum(pseudo.obs.df$time.2[dat$trt==1]*dat$ps.IPTW.PO[dat$trt==1])/dim(dat)[1]
  # 
  # 
  # res.5y[i,1]<-sum(pseudo.obs.df$time.5[dat$trt==0]*dat$ps.IPTW.PO[dat$trt==0])/dim(dat)[1]
  # res.5y[i,2]<-sum(pseudo.obs.df$time.5[dat$trt==1]*dat$ps.IPTW.PO[dat$trt==1])/dim(dat)[1]
  # 
  # res.10y[i,1]<-sum(pseudo.obs.df$time.10[dat$trt==0]*dat$ps.IPTW.PO[dat$trt==0])/dim(dat)[1]
  # res.10y[i,2]<-sum(pseudo.obs.df$time.10[dat$trt==1]*dat$ps.IPTW.PO[dat$trt==1])/dim(dat)[1]

  # ################## 9. Pseudo-observations with centered weights 
  # #start.time<-Sys.time()
  # 
  i<-match("PO_adj_wt",methods) #Index number of this method


  #Centering the weights at 1 makes results unbiased
  weight.correction<-as.numeric(dat$trt==1)*dim(dat)[1]/sum(dat$psweight2[dat$trt==1]) +
    as.numeric(dat$trt==0)*dim(dat)[1]/sum(dat$psweight2[dat$trt==0])

    dat$ps.IPTW.PO<-dat$psweight2*weight.correction

  #dat$ps.IPTW.PO<-dat$ps.IPTW

  #Find median survival - need to use the whole curve

  #pseudo.obs<-pseudosurv(dat$time, dat$fail) #,tmax=cutoffs)
  cutoffs<-seq(from=1/12,to=last.obs.time,by=1/12)
  cutoffs.obs<-cutoffs[cutoffs>min(dat$time[dat$fail==1]) & cutoffs<max(dat$time[dat$fail==1])]
  
  #NOTE FOR DEMONSTRATION CODE ONLY: must select a subset of the observations unless run
  #on a high-memory system
  use<-rbinom(n=length(dat$time),size=1,p=0.1)
  #For original analysis
  #use<-rep(1, length(dat$time))
  trt.use<-dat$trt[use==1]


  
  
  cutoffs.obs<-cutoffs[cutoffs>min(dat$time[dat$fail==1 & use==1]) & cutoffs<max(dat$time[dat$fail==1& use==1])]
  pseudo.obs<-pseudosurv(dat$time[use==1], dat$fail[use==1],tmax=cutoffs.obs)


  #get the pseudo-observations
  pseudo.obs.mat<-as.matrix(pseudo.obs$pseudo)
  #Weight them by the propensity score weights
  pseudo.obs.mat.wt<-diag(as.numeric(dat$ps.IPTW.PO[use==1]))%*%pseudo.obs.mat
  #Survival at each time is the average of the weighted pseudo observations
  S_t_0<-colSums(pseudo.obs.mat.wt[(trt.use==0 ),])/dim(dat[use==1,])[1]
  S_t_1<-colSums(pseudo.obs.mat.wt[(trt.use==1),])/dim(dat[use==1,])[1]
  #Find the first time where the survival curve is <0.5
  res.median[i,1]<-min(pseudo.obs$time[S_t_0<0.5])
  res.median[i,2]<-min(pseudo.obs$time[S_t_1<0.5])

  #Pseudo mean - built in function
  pseudo.RMS = pseudomean(dat$time[use==1], dat$fail[use==1],tmax=last.obs.time)

  res.RMS[i,1]<-sum(pseudo.RMS[trt.use==0]*dat$ps.IPTW.PO[dat$trt==0 & use==1])/dim(dat[use==1,])[1]
  res.RMS[i,2]<-sum(pseudo.RMS[trt.use==1]*dat$ps.IPTW.PO[dat$trt==1 & use==1])/dim(dat[use==1,])[1]


  #Survival at specific timepoints
  #cutoffs <- c(1:14)
  #pseudo.obs<-pseudosurv(dat$time, dat$fail,tmax=cutoffs)

  pseudo.obs.df<-data.frame(pseudo.obs$pseudo)

  #Directly estimate the probability of surviving to time t
  #This is just the mean of the pseudo observations
  #Applying propensity score weights

  res.2y[i,1]<-sum(pseudo.obs.df$time.2[trt.use==0]*dat$ps.IPTW.PO[dat$trt==0 & use==1])/dim(dat[use==1,])[1]
  res.2y[i,2]<-sum(pseudo.obs.df$time.2[trt.use==1]*dat$ps.IPTW.PO[dat$trt==1& use==1])/dim(dat[use==1,])[1]


  res.5y[i,1]<-sum(pseudo.obs.df$time.5[trt.use==0]*dat$ps.IPTW.PO[dat$trt==0 & use==1])/dim(dat[use==1,])[1]
  res.5y[i,2]<-sum(pseudo.obs.df$time.5[trt.use==1]*dat$ps.IPTW.PO[dat$trt==1 & use==1])/dim(dat[use==1,])[1]

  res.10y[i,1]<-sum(pseudo.obs.df$time.10[trt.use==0]*dat$ps.IPTW.PO[dat$trt==0 & use==1])/dim(dat[use==1,])[1]
  res.10y[i,2]<-sum(pseudo.obs.df$time.10[trt.use==1]*dat$ps.IPTW.PO[dat$trt==1 & use==1])/dim(dat[use==1,])[1]
  # 
  # #end.time<-Sys.time()
  # #print(end.time-start.time)
  # 

  #Format the results to output
  
  diff.median<-res.median[,2]-res.median[,1]
  diff.RMS<-res.RMS[,2]-res.RMS[,1]
  diff.2y<-res.2y[,2]-res.2y[,1]
  diff.5y<-res.5y[,2]-res.5y[,1]
  diff.10y<-res.10y[,2]-res.10y[,1]
  
  res.median<-cbind(res.median,diff.median)
  res.RMS<-cbind(res.RMS,diff.RMS)
  res.2y<-cbind(res.2y,diff.2y)
  res.5y<-cbind(res.5y,diff.5y)
  res.10y<-cbind(res.10y,diff.10y)
  
  rownames(res.median)<-rownames(res.RMS)<-rownames(res.2y)<-rownames(res.5y)<-rownames(res.10y)<-methods  
  colnames(res.median)<-colnames(res.RMS)<-colnames(res.2y)<-colnames(res.5y)<-colnames(res.10y)<-c("Untreated","Treated","Difference")  
  
  return(list(res.median,res.RMS,res.2y,res.5y,res.10y))
  
}


#Helper functions
median.ggam <- function(mu, sigma, Q) {qgengamma(0.5, mu = mu, sigma = sigma, 
                                                 Q=Q)}
median.weibull <- function(shape, scale) {    qweibull(0.5, shape = shape, scale = scale)}


