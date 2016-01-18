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



library(AED); data(Boar)
#Linear regression:
B0=lm(Tb ~LengthCT, data = Boar)
MyData <- data.frame(LengthCT = seq(from = 46.5, to = 165, by = 1))
Pred <- predict(B0,newdata = MyData, type = "response")
#plot(x=Boar$LengthCT, y = Boar$Tb,xlab="Length", ylab="Probability of Tb")
plot(x=Boar$LengthCT, y = Boar$Tb,xlab="Length", ylab="Tb")
lines(MyData$LengthCT,Pred)



B1=glm(Tb ~LengthCT, family = binomial, data = Boar)
summary(B1)

MyData <- data.frame(LengthCT = seq(from = 46.5, to = 165, by = 1))
Pred <- predict(B1,newdata = MyData, type = "response")
plot(x=Boar$LengthCT, y = Boar$Tb,xlab="Length", ylab="Tb")
lines(MyData$LengthCT,Pred)



#Other link functions
B1.A=glm(Tb ~LengthCT, family = binomial(link="probit"), data = Boar,)
B1.B=glm(Tb ~LengthCT, family = binomial(link="cloglog"), data = Boar,)
Pred.A <- predict(B1.A,newdata = MyData, type = "response")
Pred.B <- predict(B1.B,newdata = MyData, type = "response")
lines(MyData$LengthCT,Pred.A,col=2,lty=2)
lines(MyData$LengthCT,Pred.B,col=4,lty=3)

legend("topleft",legend=c("logit","probit","clog-log"),
       lty=c(1,2,3), col=c(1,2,4),cex=0.6)
       


par(mfrow=c(2,2))
plot(B1)


summary(B1)
drop1(B1, test="Chi")



#Parasites in cod
library(AED); data(ParasiteCod)

ParasiteCod$fArea <- factor(ParasiteCod$Area)
ParasiteCod$fYear <- factor(ParasiteCod$Year)
P1 <- glm(Prevalence ~ fArea * fYear + Length,
               family = binomial, data = ParasiteCod)



library(AED); data(Tbdeer)
Z <- cbind(Tbdeer$OpenLand,Tbdeer$ScrubLand,Tbdeer$QuercusPlants,Tbdeer$QuercusTrees,
        Tbdeer$ReedDeerIndex,Tbdeer$EstateSize,Tbdeer$Fenced)
corvif(Z)
DeerNegCervi <- Tbdeer$DeerSampledCervi - Tbdeer$DeerPosCervi

Tbdeer$fFenced <- factor(Tbdeer$Fenced)
Deer1=glm(cbind(DeerPosCervi,DeerNegCervi)~
       OpenLand+ScrubLand+QuercusPlants+QuercusTrees+
       ReedDeerIndex+ EstateSize+fFenced,
       family=binomial, data = Tbdeer)
summary(Deer1)


Tbdeer$DeerPosProp <- Tbdeer$DeerPosCervi / Tbdeer$DeerSampledCervi
Deer2 <- glm(DeerPosProp~OpenLand+ScrubLand+QuercusPlants+
        QuercusTrees+ReedDeerIndex+EstateSize+
        fFenced,
        family=binomial,weights=DeerSampledCervi,data = Tbdeer)


Deer4 <- glm(cbind(DeerPosCervi,DeerNegCervi) ~ OpenLand,
       family=quasibinomial, data = Tbdeer)
drop1(Deer4,test="F")




MyData <- data.frame(OpenLand = seq(from = min(Tbdeer$OpenLand),
                                to = max(Tbdeer$OpenLand),by=0.01))
P1 <- predict(Deer4, newdata = MyData, type = "link", se = TRUE)
plot(MyData$OpenLand,exp(P1$fit)/(1+exp(P1$fit)),
     type="l",ylim=c(0,1),
     xlab="Percentage open land",
     ylab="Probability on E. cervi")
lines(MyData$OpenLand,exp(P1$fit+1.96*P1$se.fit)/
       (1+exp(P1$fit+1.96*P1$se.fit)),lty=2)
lines(MyData$OpenLand,exp(P1$fit-1.96*P1$se.fit)/
       (1+exp(P1$fit-1.96*P1$se.fit)),lty=2)
points(Tbdeer$OpenLand,Tbdeer$DeerPosProp)






library(AED); data(ParasiteCod)
library(mgcv)
ParasiteCod$fArea <- factor(ParasiteCod$Area)
ParasiteCod$fYear <- factor(ParasiteCod$Year)
P2 = gam(Prevalence~factor(Area)*factor(Year)+
           s(Length),family=binomial, data = ParasiteCod)
