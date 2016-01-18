#    Mixed Effects Models and Extensions in Ecology with R (2009)
#    Zuur, Ieno, Walker, Saveliev and Smith.    Springer
#    This file was produced by Alain Zuur (highstat@highstat.com)
#    www.highstat.com

#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.




library(AED); data(DeerEcervi)

DeerEcervi$Ecervi.01 <- DeerEcervi$Ecervi
DeerEcervi$Ecervi.01[DeerEcervi$Ecervi>0]<-1
DeerEcervi$fSex <- factor(DeerEcervi$Sex)
DeerEcervi$CLength <- DeerEcervi$Length - mean(DeerEcervi$Length)
DeerEcervi$fFarm <- factor(DeerEcervi$Farm)

DE.glm<-glm(Ecervi.01 ~ CLength * fSex+fFarm, data = DeerEcervi,
         family = binomial)
drop1(DE.glm, test = "Chi")
summary(DE.glm)


plot(DeerEcervi$CLength,DeerEcervi$Ecervi.01,
     xlab="Length",ylab="Probability of \
      presence of E. cervi L1",main="Male data")
I<-order(DeerEcervi$CLength)
AllFarms<-unique(DeerEcervi$Farm)
for (j in AllFarms){
  mydata<-data.frame(CLength=DeerEcervi$CLength,fSex="1",
                     fFarm=AllFarms[j])
  n<-dim(mydata)[1]
  if (n>10){
     P.DE2<-predict(DE.glm,mydata,type="response")
     lines(mydata$CLength[I],P.DE2[I])
}}



library(MASS)
DE.PQL<-glmmPQL(Ecervi.01 ~ CLength * fSex,
      random = ~ 1 | fFarm, family = binomial, data = DeerEcervi)
summary(DE.PQL)

eta <- fitted(DE.PQL)
fv  <- exp(eta)/(1+exp(eta))
res.raw <- DeerEcervi$Ecervi.01 - fv
res.P <- (DeerEcervi$Ecervi.01 - fv) / sqrt(fv * (1-fv))
sd(res.P)


intervals(DE.PQL,which="var-cov")


DE.PQL1<-glmmPQL(Ecervi.01 ~ CLength * fSex,
      random = ~ 1 | fFarm,
      family = quasi(link=logit,variance="mu(1-mu)"),
      data = DeerEcervi)
summary(DE.PQL1)




g <- 0.8883697 + 0.0378608 * DeerEcervi$CLength
p.averageFarm1<-exp(g)/(1+exp(g))
I1<-order(DeerEcervi$CLength)              #Avoid spaghetti plot
plot(DeerEcervi$CLength,DeerEcervi$Ecervi.01,xlab="Length",
     ylab="Probability of presence of E. cervi L1")
lines(DeerEcervi$CLength[I1],p.averageFarm1[I],lwd=3)
p.Upp<-exp(g+1.96*1.462108)/(1+exp(g+1.96*1.462108))
p.Low<-exp(g-1.96*1.462108)/(1+exp(g-1.96*1.462108))
lines(DeerEcervi$CLength[I1],p.Upp[I1])
lines(DeerEcervi$CLength[I1],p.Low[I1])


library(lme4)
DE.lme4<-lmer(Ecervi.01 ~ CLength * fSex +(1|fFarm),
           family = binomial, data = DeerEcervi)
summary(DE.lme4)


library(glmmML)
DE.glmmML<-glmmML(Ecervi.01 ~ CLength * fSex,
             cluster = fFarm,family=binomial, data = DeerEcervi)
summary(DE.glmmML)






#Owls

library(AED) ; data(Owls)
library(nlme)
Owls$NCalls<-Owls$SiblingNegotiation
Owls$LBroodSize<-log(Owls$BroodSize)
Owls$fNest<-factor(Owls$Nest)


#GLMM
library(lme4)

O1.lmer<-lmer(NCalls~offset(LBroodSize)+
              SexParent*FoodTreatment+
              SexParent*ArrivalTime+(1|fNest),data=Owls,
              family=poisson)
summary(O1.lmer)

O2.lmer<-lmer(NCalls~offset(LBroodSize)+
              SexParent*FoodTreatment+
              SexParent+ArrivalTime+(1|fNest),data=Owls,
              family=poisson)
anova(O1.lmer,O2.lmer,test="F")



O3.lmer<-lmer(NCalls~offset(LBroodSize)+
              FoodTreatment+
              ArrivalTime+(1|fNest),data=Owls,
              family=poisson)
summary(O3.lmer)





library(mgcv)
O4.gamm<-gamm(NCalls~offset(LBroodSize)+
              FoodTreatment+s(ArrivalTime),
              random=list(fNest=~1),data=Owls,
              family=poisson)

summary(O4.gamm$gam)
anova(O4.gamm$gam)
plot(O4.gamm$gam)

summary(O4.gamm$lme)

E4<-resid(O4.gamm$lme,type="normalized")

O4.gamm<-gamm(NCalls~offset(LBroodSize)+
              FoodTreatment+s(ArrivalTime),
              random=list(fNest=~1),data=Owls,
              family=quasi(link=log,variance="mu"))
              
              
summary(O4.gamm$gam)

intervals(O4.gamm$lme,which="var-cov")



eta <- fitted(O4.gamm$lme)+Owls$LBroodSize
fv  <- exp(eta)
res.raw <- Owls$NCalls - fv
res.P <- (Owls$NCalls - fv) / sqrt(fv )
sd(res.P)



#Did not make it into the book

#####Californian birds

library(AED); data(RiceFieldBirds)
RFBirds<-RiceFieldBirds
RFBirds$Richness<-rowSums(RFBirds[,8:56] > 0)
RFBirds$fField <-factor(RFBirds$FIELD)
RFBirds$LA<-log(RFBirds$AREA)
RFBirds$fSptreat <- factor(RFBirds$SPTREAT)
RFBirds$DEPTH2 <- RFBirds$DEPTH^2
M0 <- glm(Richness~offset(LA)+fSptreat+DEPTH+
             DEPTH2,family=quasipoisson,
             data=RFBirds)
summary(M0)



library(lme4)
M.lme<-lmer(Richness~offset(LA)+fSptreat+DEPTH+DEPTH2+(1|fField),
            data = RFBirds,
            family=poisson)
summary(M.lme)



