#Plots to display results of simulation analyses
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(wesanderson)
library(khroma)

setwd("./")

#### Bias graphs

restmp<-read.csv("Bias_results_graph_input.csv")

res<-melt(restmp, id.vars = c("Scenario", "Method"))
names(res)<-c("Scenario","Method","Outcome","Bias")

#Remove unadjusted
res<-res[res$Method!="NUA",]

#Absolute value of bias
res$Bias<-abs(res$Bias)


#Change order of variables
res$Outcome<-factor(res$Outcome, levels=c("X2y","X5y","X10y","median","RMS"),
                    labels=c("2 years","5 years","10 years","Median","RMS"))

res$Method<-factor(res$Method, levels=c("NCM","CTV.LT","CTV.LPW","AFT.GG","AFT.WBL.LS","PO","WKM")
                   , labels=c("Cox","CTV LT","CTV PWC","AFT GG","AFT WBL LS","Pseudo","Wtd KM"))

res$Scenario<-factor(res$Scenario, levels=c("Base", "LPW" , "Modest NPH", "Modest TE" , "Small SS" ),
                     labels=c("Base", "PWC" , "Modest NPH", "Modest TE" , "Small SS" ))

library(khroma)
#All outcomes at once

p<-ggplot(data=res, aes(x=Scenario, y=Bias, fill=Scenario)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_grid(cols = vars(Method), rows=vars(Outcome),scales="free_y") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  scale_fill_manual(values=wes_palette(n=5, name="Cavalcanti1"))

print(p)

muted <- colour("muted")

p<-ggplot(data=res, aes(x=Method, y=Bias, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_grid(cols = vars(Scenario), rows=vars(Outcome),scales="free_y") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank())+
  scale_fill_manual(values=as.vector(muted(9)))

print(p)
##### Variance graphs


restmp<-read.csv("Var_results_graph_input.csv")

res<-melt(restmp, id.vars = c("Scenario", "Method"))
names(res)<-c("Scenario","Method","Outcome","Variance")

#Remove unadjusted
res<-res[res$Method!="NUA",]

#Standard error
res$SE<-sqrt(res$Variance)



#Change order of variables
res$Outcome<-factor(res$Outcome, levels=c("X2y","X5y","X10y","median","RMS"),
                    labels=c("2 years","5 years","10 years","Median","RMS"))

res$Method<-factor(res$Method, levels=c("NCM","CTV.LT","CTV.LPW","AFT.GG","AFT.WBL.LS","PO","WKM")
                   , labels=c("Cox","CTV LT","CTV PWC","AFT GG","AFT WBL LS","Pseudo","Wtd KM"))

res$Scenario<-factor(res$Scenario, levels=c("Base", "LPW" , "Modest NPH", "Modest TE" , "Small SS" ),
                     labels=c("Base", "PWC" , "Modest NPH", "Modest TE" , "Small SS" ))


#All outcomes at once

p<-ggplot(data=res, aes(x=Scenario, y=SE, fill=Scenario)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_grid(cols = vars(Method), rows=vars(Outcome),scales="free_y") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  scale_fill_manual(values=wes_palette(n=5, name="Cavalcanti1"))
print(p)


p<-ggplot(data=res, aes(x=Method, y=SE, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_grid(cols = vars(Scenario), rows=vars(Outcome),scales="free_y") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank())+
  scale_fill_manual(values=as.vector(muted(9)))

print(p)

