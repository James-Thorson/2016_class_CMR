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



library(AED);data(RIKZ)

Beta<-vector(length=9)
for (i in 1:9){
 tmpout<-summary(lm(Richness~NAP,subset = (Beach==i),data=RIKZ))
 Beta[i]<-tmpout$coefficients[2,1]
}

ExposureBeach <- c("a","a","b","b","a",
                     "b","b","a","a")
tmp2 <- lm(Beta ~ factor(ExposureBeach),data=RIKZ)


library(nlme)
RIKZ$fBeach <- factor(RIKZ$Beach)
Mlme1 <- lme(Richness ~ NAP, random = ~1 | fBeach,data=RIKZ)
summary(Mlme1)





F0<-fitted(Mlme1,level=0)
F1<-fitted(Mlme1,level=1)
I<-order(RIKZ$NAP)
NAPs<-sort(RIKZ$NAP)
plot(NAPs,F0[I],lwd=4,type="l",ylim=c(0,22),
     ylab="Richness",xlab="NAP")
for (i in 1:9){
   x1<-RIKZ$NAP[RIKZ$Beach==i]
   y1<-F1[RIKZ$Beach==i]
   K<-order(x1)
   lines(sort(x1),y1[K])
}
text(RIKZ$NAP,RIKZ$Richness,RIKZ$Beach,cex=0.9)



Mlme2 <- lme(Richness ~ NAP,
             random = ~1 + NAP | fBeach, data = RIKZ)
summary(Mlme2)




Mlme3 <- lme(Richness ~ 1, random = ~1 | fBeach,
           data = RIKZ)
summary(Mlme3)




M.mixed <- lme(Richness ~ NAP, random = ~1 | fBeach,
                 method = "REML", data = RIKZ)
M.gls <- gls(Richness ~ NAP, method = "REML",
          correlation = corCompSymm(form =~1 | fBeach),
          data = RIKZ)





RIKZ$fExp<-RIKZ$Exposure
RIKZ$fExp[RIKZ$fExp==8]<-10
RIKZ$fExp<-factor(RIKZ$fExp,levels=c(10,11))
M0.ML <- lme(Richness ~ NAP, data = RIKZ,
              random = ~1| fBeach, method = "ML")
M0.REML <-lme(Richness ~ NAP, random = ~1|fBeach,
              method = "REML", data = RIKZ)
M1.ML <- lme(Richness ~ NAP+fExp, data = RIKZ,
              random = ~1| fBeach, method = "ML")
M1.REML <- lme(Richness ~NAP+fExp, data = RIKZ,
              random = ~1|fBeach, method = "REML")





Wrong1 <- gls(Richness ~ 1 + NAP, method = "REML",
               data = RIKZ)
Wrong2 <- lme(Richness ~ 1 + NAP, random = ~1|fBeach,
               method = "REML", data = RIKZ)
Wrong3 <- lme(Richness ~ 1 + NAP, method = "REML",
               random = ~1 + NAP | fBeach, data = RIKZ)



AIC(Wrong1,Wrong2,Wrong3)

RIKZ$fExp<-RIKZ$Exposure
RIKZ$fExp[RIKZ$fExp==8]<-10
RIKZ$fExp<-factor(RIKZ$fExp,levels=c(10,11))

Wrong4 <- lme(Richness ~1 + NAP * fExp,
             random = ~1 + NAP | fBeach,
             method = "REML", data = RIKZ)
anova(Wrong4)
summary(Wrong4)

#Drop the interaction

Wrong4 <- lme(Richness ~1 + NAP + fExp,
             random = ~1 + NAP | fBeach,
             method = "REML", data = RIKZ)

summary(Wrong4)









 lmc <- lmeControl(niterEM = 5200, msMaxIter = 5200)
 Wrong4A <- lme(Richness ~1 + NAP, method="ML",
             control = lmc, data = RIKZ,
               random = ~1+NAP|fBeach)
 Wrong4B <- lme(Richness ~ 1 + NAP + fExp,
               random = ~1 + NAP | fBeach, method="ML",
               data = RIKZ,control = lmc)
 Wrong4C <- lme(Richness ~1 + NAP * fExp,
               random = ~1 + NAP | fBeach, data = RIKZ,
               method = "ML", control = lmc)
 anova(Wrong4A, Wrong4B, Wrong4C)


 B1 <- gls(Richness ~ 1 + NAP * fExp,
            method = "REML", data = RIKZ)
 B2 <- lme(Richness ~1 + NAP * fExp, data = RIKZ,
        random = ~1 | fBeach, method = "REML")
 B3 <- lme(Richness ~ 1 + NAP * fExp,data = RIKZ,
        random = ~1 + NAP | fBeach, method = "REML")

summary(B2)

#Drop interaction
B2 <- lme(Richness ~1 + NAP + fExp, data = RIKZ,
        random = ~1 | fBeach, method = "REML")
summary(B2)



#Owls
library(AED) ; data(Owls)
names(Owls)

# "FoodTreatment"      "SexParent"
#[4] "ArrivalTime"        "SiblingNegotiation" "BroodSize"
#[7] "NegPerChick"



boxplot(NegPerChick~Nest,data=Owls)


boxplot(NegPerChick~FoodTreatment,data=Owls)
boxplot(NegPerChick~SexParent,data=Owls)

plot(x=Owls$ArrivalTime,y=Owls$NegPerChick)


M.lm=lm(NegPerChick~SexParent*FoodTreatment+SexParent*ArrivalTime,data=Owls)
plot(M.lm,select=c(1))

Owls$LogNeg<-log10(Owls$NegPerChick+1)
M2.lm=lm(LogNeg~SexParent*FoodTreatment+SexParent*ArrivalTime,data=Owls)
E=rstandard(M2.lm)

op<-par(mar=c(3,3,2,2))
boxplot(E~Nest,data=Owls,axes=FALSE,ylim=c(-3,3))
abline(0,0)
axis(2)
text(1:27,-2.5, levels(Owls$Nest),cex=0.75,srt=65)
par(op)



#Step 2 of protocol
library(nlme)
Form<-formula(LogNeg~SexParent*FoodTreatment+SexParent*ArrivalTime)
M.gls=gls(Form,data=Owls)

M1.lme=lme(Form,random=~1|Nest,method="REML",data=Owls)

anova(M.gls,M1.lme)

E2<-resid(M1.lme,type="normalized")
F2<-fitted(M1.lme)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=F2,y=E2,xlab="Fitted values",ylab=MyYlab)
boxplot(E2~SexParent,data=Owls,main="Sex of parent",ylab=MyYlab)
boxplot(E2~FoodTreatment,data=Owls,main="Food treatment",ylab=MyYlab)
plot(x=Owls$ArrivalTime,y=E,main="Arrival time",ylab=MyYlab,xlab="Time (hours)")
par(op)


#step 7 of the protocol
M1.lme=lme(Form,random=~1|Nest,method="REML",data=Owls)
summary(M1.lme)
anova(M1.lme)

M1.Full=lme(Form,random=~1|Nest,method="ML",data=Owls)
M1.A=update(M1.Full,.~.-SexParent:FoodTreatment)
M1.B=update(M1.Full,.~.-SexParent:ArrivalTime)
anova(M1.Full,M1.A)
anova(M1.Full,M1.B)


Form2<-formula(LogNeg~SexParent+FoodTreatment+SexParent*ArrivalTime)
M2.Full=lme(Form2, random= ~1| Nest, method = "ML", data = Owls)
M2.A=update(M2.Full, .~. -FoodTreatment)
M2.B=update(M2.Full, .~. -SexParent:ArrivalTime)
anova(M2.Full,M2.A)
anova(M2.Full,M2.B)



Form3 <- formula(LogNeg~SexParent+FoodTreatment+ArrivalTime)
M3.Full <- lme(Form3, random= ~1| Nest, method = "ML", data = Owls)
M3.A <- update(M3.Full, .~. -FoodTreatment)
M3.B <- update(M3.Full, .~. -SexParent)
M3.C <- update(M3.Full, .~. -ArrivalTime)
anova(M3.Full,M3.A)
anova(M3.Full,M3.B)
anova(M3.Full,M3.C)




Form4 <- formula(LogNeg ~ FoodTreatment + ArrivalTime)
M4.Full <- lme(Form4, random= ~1| Nest, method = "ML", data = Owls)
M4.A <- update(M4.Full, .~. -FoodTreatment)
M4.B <- update(M4.Full, .~. -ArrivalTime)
anova(M4.Full,M4.A)
anova(M4.Full,M4.B)


M5 <- lme(LogNeg ~ FoodTreatment + ArrivalTime, random= ~1| Nest, method = "REML", data = Owls)
summary(M5)


E2<-resid(M5,type="normalized")
F2<-fitted(M5)
op<-par(mfrow=c(2,2),mar=c(4,4,3,2))
MyYlab="Residuals"

plot(x=F2,y=E2,xlab="Fitted values",ylab=MyYlab)
boxplot(E2~SexParent,data=Owls,main="Sex of parent",ylab=MyYlab)
boxplot(E2~FoodTreatment,data=Owls,main="Food treatment",ylab=MyYlab)
plot(x=Owls$ArrivalTime,y=E,main="Arrival time",ylab=MyYlab,xlab="Time (hours)")
par(op)





library(lattice)

xyplot(E2~ArrivalTime|SexParent*FoodTreatment,data=Owls,
  ylab="Residuals",xlab="Arrival time (hours)",
      panel=function(x,y){
    panel.grid(h=-1, v= 2)
    panel.points(x,y,col=1)
    panel.loess(x,y,span=0.5,col=1,lwd=2)})
    
    


library(mgcv)
M6 <- gamm(LogNeg ~ FoodTreatment + +s(ArrivalTime),
        random=list(Nest=~1),data=Owls)
        
summary(M6)
plot(M6$gam)
anova(M6$gam)
summary(M6$gam)


M7=gamm(NegPerChick~FoodTreatment+
       s(ArrivalTime,by=as.numeric(FoodTreatment=="Deprived"))+
       s(ArrivalTime,by=as.numeric(FoodTreatment=="Satiated")),
        random=list(Nest=~1),data=Owls)

M8=gamm(NegPerChick~FoodTreatment+
       s(ArrivalTime,by=as.numeric(SexParent=="Female"))+
       s(ArrivalTime,by=as.numeric(SexParent=="Male")),
        random=list(Nest=~1),data=Owls)

AIC(M6$lme,M7$lme,M8$lme)

############################################

