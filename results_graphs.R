#Plots to display results of simulation analyses
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
#library(wesanderson)
#library(khroma)
library(NatParksPalettes)

setwd("./")


#### Bias graphs

restmp<-read.csv("Tables_graphs/Fig1_input_bias.csv")

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

res$Scenario<-factor(res$Scenario, levels=c("Base", "LPW" , "Modest NPH", "Modest TE" , "Small SS","GG lambda -0.5"
                                            , "GG lambda 2.5" ),
                     labels=c("1 Base", "2 PWC" , "3 Mod NPH", "4 Mod TE" , "5 Sm SS", "6 GG -0.5"
                              , "7 GG 2.5" ))

#library(khroma)
#All outcomes at once

p<-ggplot(data=res, aes(x=Scenario, y=Bias, fill=Scenario)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_grid(cols = vars(Method), rows=vars(Outcome),scales="free_y") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  scale_fill_manual(values=natparks.pals("DeathValley", 7))

print(p)



#muted <- colour("muted")

p<-ggplot(data=res, aes(x=Method, y=Bias, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_grid(cols = vars(Scenario), rows=vars(Outcome),scales="free_y") +
  #facet()+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position="bottom")+
  scale_fill_manual(values=natparks.pals("DeathValley", 7))

print(p)

png("Tables_graphs/Figure1.png")
print(p)
dev.off()
##### Variance graphs


restmp<-read.csv("Tables_graphs/FigA1_input_var.csv")

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

res$Scenario<-factor(res$Scenario, levels=c("Base", "LPW" , "Modest NPH", "Modest TE" , "Small SS","GG lambda -0.5"
                                            , "GG lambda 2.5" ),
                     labels=c("1 Base", "2 PWC" , "3 Mod NPH", "4 Mod TE" , "5 Sm SS", "6 GG -0.5"
                              , "7 GG 2.5" ))

#All outcomes at once

p<-ggplot(data=res, aes(x=Scenario, y=SE, fill=Scenario)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_grid(cols = vars(Method), rows=vars(Outcome),scales="free_y") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  scale_fill_manual(values=natparks.pals("DeathValley", 7))
print(p)


p<-ggplot(data=res, aes(x=Method, y=SE, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge()) +
  facet_grid(cols = vars(Scenario), rows=vars(Outcome),scales="free_y") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position="bottom")+
  scale_fill_manual(values=natparks.pals("DeathValley", 7))

print(p)

png("Tables_graphs/FigureA1.png")
print(p)
dev.off()
