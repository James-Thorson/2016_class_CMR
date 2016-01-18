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




library(AED); data(Boreality)
Boreality$Bor<-sqrt(1000*(Boreality$nBor+1)/(Boreality$nTot))
B.lm<-lm(Bor~Wet,data=Boreality)
summary(B.lm)

E<-rstandard(B.lm)
library(gstat)
mydata<-data.frame(E,Boreality$x,Boreality$y)
coordinates(mydata)<-c("Boreality.x","Boreality.y")
bubble(mydata,"E",col=c("black","grey"),
       main="Residuals",xlab="X-coordinates",
       ylab="Y-coordinates")

Vario1 = variogram(E ~ 1, mydata)
plot(Vario1)

Vario2 <- variogram(E ~ 1, mydata, alpha = c(0, 45, 90,135) )
plot(Vario2)


library(nlme)
D<-seq(0,1,by=0.1)
Mydata2<-data.frame(D=D)
cs1C <- corSpher(c(0.8,0.1), form = ~ D,nugget=T)
cs1C <- Initialize(cs1C,data=Mydata2)
v1C<-Variogram(cs1C)
plot(v1C,smooth=F,type="l",col=1)


library(nlme)
f1 <- formula(Bor ~ Wet)
B1.gls<-gls(f1, data = Boreality)
var1<-Variogram(B1.gls,form=~x+y,robust=TRUE,maxDist=2000,
             resType="pearson")
plot(var1,smooth=T)



B1A<-gls(f1,correlation=corSpher(form=~x+y,nugget=T),data=Boreality)
B1B<-gls(f1,correlation=corLin(form=~x+y,nugget=T),data=Boreality)
B1C<-gls(f1,correlation=corRatio(form=~x+y,nugget=T),data=Boreality)
B1D<-gls(f1,correlation=corGaus(form=~x+y,nugget=T),data=Boreality)
B1E<-gls(f1,correlation=corExp(form=~x+y,nugget=T),data=Boreality)
AIC(B1.gls,B1A,B1B,B1C,B1D,B1E)


Vario1E <- Variogram(B1E, form =~ x + y, robust = TRUE, maxDist = 2000,
                 resType = "pearson")
plot(Vario1E,smooth=FALSE)

Vario2E <- Variogram(B1E, form =~ x + y,
                     robust =  TRUE, maxDist = 2000,
                     resType = "normalized")
plot(Vario2E, smooth = FALSE)






#Section 7.3: Hawaii birds
library(AED); data(Hawaii)
Birds <- c(Hawaii$Stilt.Oahu, Hawaii$Stilt.Maui,
             Hawaii$Coot.Oahu, Hawaii$Coot.Maui)
Time <- rep(Hawaii$Year, 4)
Rain <- rep(Hawaii$Rainfall, 4)
ID <- factor(rep(c("Stilt.Oahu", "Stilt.Maui",
                     "Coot.Oahu", "Coot.Maui"),
                   each=length(Hawaii$Year)))

library(mgcv); library(nlme)

#In R versions 2.7.0 and earlier, use:
f1<-formula(Birds~Rain+ID+
        s(Time,by=as.numeric(ID=="Stilt.Oahu"))+
        s(Time,by=as.numeric(ID=="Stilt.Maui"))+
        s(Time,by=as.numeric(ID=="Coot.Oahu"))+
        s(Time,by=as.numeric(ID=="Coot.Maui")))

#In later R versions, use:
f1<-formula(Birds ~ Rain + ID + s(Time,by=ID))


        
HawA<-gamm(f1,method="REML",
        correlation=corAR1(form=~Time |ID ),
        weights=varIdent(form=~1|ID))
HawB<-gamm(f1,method="REML",
        correlation=corLin(form=~Time | ID,nugget=T),
        weights=varIdent(form=~1|ID))
HawC<-gamm(f1,method="REML",
        correlation=corGaus(form=~Time | ID,nugget=T),
        weights=varIdent(form=~1|ID))
HawD<-gamm(f1,method="REML",
        correlation=corExp(form=~Time | ID,nugget=T),
        weights=varIdent(form=~1|ID))
HawE<-gamm(f1,method="REML",
        correlation=corSpher(form=~Time | ID,nugget=T),
        weights=varIdent(form=~1|ID))
#Compare the models
AIC(HawA$lme,HawB$lme,HawC$lme,HawD$lme,
      HawE$lme)



#Section 7.4: Moby
library(AED); data(TeethNitrogen)
TN<-TeethNitrogen
N.Moby <-   TN$X15N[TN$Tooth == "Moby"]
Age.Moby <-    TN$Age[TN$Tooth == "Moby"]

library(mgcv); library(nlme)
f<-formula(N.Moby~s(Age.Moby))
Mob0<-gamm(f,method="REML")
Mob1<-gamm(f,cor=corSpher(form=~Age.Moby,nugget=T),
           method="REML")
Mob2<-gamm(f,cor=corLin(form=~Age.Moby,nugget=T),
           method="REML")
Mob3<-gamm(f,cor=corGaus(form=~Age.Moby,nugget=T),
           method="REML")
Mob4<-gamm(f,cor=corExp(form=~Age.Moby,nugget=T),
           method="REML")
Mob5<-gamm(f,cor=corRatio(form=~Age.Moby,nugget=T),
           method="REML")
Mob6<-gamm(f,cor=corAR1(form=~Age.Moby),method="REML")
AIC(Mob0$lme,Mob1$lme,Mob2$lme,Mob4$lme,Mob5$lme,Mob6$lme)

E0=resid(Mob0$lme,type="normalized")
E5=resid(Mob5$lme,type="normalized")





#All whales

lmc <- lmeControl(niterEM = 5200, msMaxIter = 5200)

#In R versions 2.7.0 and earlier, use:
 AllWhales.0<-gamm( X15N ~
s(Age,by=as.numeric(Tooth=="M2679/93"),bs="cr")+
s(Age,by=as.numeric(Tooth=="M2683/93"),bs="cr")+
s(Age,by=as.numeric(Tooth=="M2583/94(1)"),bs="cr")+
s(Age,by=as.numeric(Tooth=="M2583/94(7)"),bs="cr")+
s(Age,by=as.numeric(Tooth=="M2583/94(10)"),bs="cr")+
s(Age,by=as.numeric(Tooth=="M546/95"),bs="cr")+
s(Age,by=as.numeric(Tooth=="M143/96E"),bs="cr")+
s(Age,by=as.numeric(Tooth=="Moby"),bs="cr")+
s(Age,by=as.numeric(Tooth=="M447/98"),bs="cr")+
s(Age,by=as.numeric(Tooth=="I1/98"),bs="cr")+
factor(Tooth), method="REML",control=lmc,
data=TN)



AllWhales.corGaus <- gamm( X15N ~
  s(Age,by=as.numeric(Tooth=="M2679/93"),bs="cr")+
  s(Age,by=as.numeric(Tooth=="M2683/93"),bs="cr")+
  s(Age,by=as.numeric(Tooth=="M2583/94(1)"),bs="cr")+
  s(Age,by=as.numeric(Tooth=="M2583/94(7)"),bs="cr")+
  s(Age,by=as.numeric(Tooth=="M2583/94(10)"),bs="cr")+
  s(Age,by=as.numeric(Tooth=="M546/95"),bs="cr")+
  s(Age,by=as.numeric(Tooth=="M143/96E"),bs="cr")+
  s(Age,by=as.numeric(Tooth=="Moby"),bs="cr")+
  s(Age,by=as.numeric(Tooth=="M447/98"),bs="cr")+
  s(Age,by=as.numeric(Tooth=="I1/98"),bs="cr")+
  factor(Tooth), control=lmc,method="REML",
  correlation = corGaus(form=~Age|Tooth,nugget=T),
  data=TN)




#In later R versions, use:
AllWhales.0 <- gamm( X15N ~ s(Age, by=Tooth,bs="cr") + factor(Tooth),
                     method="REML",control=lmc,
                     data=TN)



AllWhales.corGaus <- gamm( X15N ~ s(Age, by=Tooth,bs="cr") + factor(Tooth),
                           control=lmc,method="REML",
                           correlation = corGaus(form=~Age|Tooth,nugget=T),
                           data=TN)




AIC(AllWhales.0$lme,AllWhales.corGaus$lme)
anova(AllWhales.0$gam)
anova(AllWhales.corGaus$gam)



#Section 7.5 Irish SDI data
library(AED); data(SDI2003)
library(lattice)
MyPch<-vector(length=dim(SDI2003)[1])
MyPch[SDI2003$Forested==1]<-16
MyPch[SDI2003$Forested==2]<-1

xyplot(Northing~Easting,aspect="iso",col=1,pch=MyPch,data=SDI2003)

library(geoR)
coords<-matrix(0,length(SDI2003$pH),2)
coords[,1]<-SDI2003$Easting;
coords[,2]<-SDI2003$Northing
gb<-list(data=SDI2003$pH,coords=coords)
plot(variog(gb,max.dist=200000))


library(nlme)
SDI2003$fForested <- factor(SDI2003$Forested)
SDI2003$LAltitude <- log(SDI2003$Altitude)
M1<-gls(pH~SDI*fForested*LAltitude, data =  SDI2003)
Vario1<-Variogram(M1,form =~ Easting + Northing, data=SDI2003,
                   nugget=T,maxDist=200000)
plot(Vario1)


M1C<-gls(pH ~ SDI * fForested * LAltitude,
           correlation=corRatio(form=~
           Easting+Northing,nugget=TRUE),data=SDI2003)
M1E<-gls(pH~SDI*Forested*LAltitude,
           correlation=corExp(form=~
           Easting+Northing,nugget=TRUE),data=SDI2003)

Vario1C<-Variogram(M1C,form =~ Easting + Northing, data=SDI2003,
                   nugget=T,maxDist=200000,
                   resType="normalized")
plot(Vario1C,smooth=FALSE)






library(gstat)
E <- resid(M1,type="normalized")
mydata<-data.frame(E,SDI2003$Easting,SDI2003$Northing)
coordinates(mydata)<-c("SDI2003.Easting","SDI2003.Northing")
bubble(mydata,"E",col=c("black","grey"),
       main="Normalised residuals",
       xlab="X-coordinates",ylab="Y-coordinates")


################################
#Limosa data

library(AED); data(Limosa)
names(Limosa)
Limosa$fID <- factor(Limosa$ID)
Limosa$fPeriod <-factor(Limosa$Period,levels=c(0,1,2),
              labels=c("Summer","LSummer.Fall","Winter"))
Limosa$fSex <-factor(Limosa$Sex, levels = c(0,1,2),
              labels=c("Unk","F","M"))

coplot(IntakeRate~Time|fPeriod*fSex,data=Limosa,xlab=c("Time (hours)"))
library(lattice)
xyplot(IntakeRate ~ Time | fID, data = Limosa,
  panel=function(x, y){
    panel.xyplot(x, y, col = 1, cex = 0.5, pch = 1)
    panel.grid(h = -1,v = 2)
    panel.abline(v = 0,lty = 2)
    if (length(x) > 5) panel.loess(x, y, span=0.9, col=1, lwd=2)
  })


Limosa$Time2 <-Limosa$Time^2-mean(Limosa$Time^2)
M.lm <- lm(IntakeRate~Time+Time2+fPeriod+fSex, data= Limosa)
drop1(M.lm,test="F")


library(nlme)

M1.gls <- gls(IntakeRate ~ Time + Time2 + fPeriod +
             fSex, data= Limosa)
E <- resid(M1.gls)

op<-par(mfrow=c(2,2))
boxplot(E~Limosa$fPeriod,main="Period")
abline(0,0)
boxplot(E~Limosa$fSex,main="Sex")
abline(0,0)
boxplot(E~Limosa$fSex * Limosa$fPeriod,main="Sex & Period")
abline(0,0)
boxplot(E~Limosa$ID,main="Day")
abline(0,0)
par(op)



M1.lme <- lme(IntakeRate ~ Time + Time2 + fPeriod +
              fSex, data= Limosa,
              weights = varIdent(form=~ 1 | fSex * fPeriod),
              random =~ 1 | fID, method = "REML")




M1.lme <- lme(IntakeRate ~ Time + Time2 + fPeriod +
        fSex, data= Limosa,
        weights = varIdent(form=~ 1 | fSex * fPeriod),
        random =~ 1 | fID, method = "ML")
M1.lmeA <- update(M1.lme, .~. -Time2)
M1.lmeB <- update(M1.lme, .~. -fPeriod)
M1.lmeC <- update(M1.lme, .~. -fSex)

anova(M1.lme,M1.lmeA)
anova(M1.lme,M1.lmeB)
anova(M1.lme,M1.lmeC)


M2.lme <- lme(IntakeRate ~ Time + Time2 + fSex, data= Limosa,
        weights = varIdent(form=~ 1 | fSex * fPeriod),
        random =~ 1 | fID, method = "ML")
M2.lmeA <- update(M2.lme, .~. -Time2)
M2.lmeB <- update(M2.lme, .~. -fSex)

anova(M2.lme,M2.lmeA)
anova(M2.lme,M2.lmeB)



M3.lme <- lme(IntakeRate ~ Time + fSex, data= Limosa,
        weights = varIdent(form=~ 1 | fSex * fPeriod),
        random =~ 1 | fID, method = "ML")
M3.lmeA <- update(M3.lme, .~. -Time)
M3.lmeB <- update(M3.lme, .~. -fSex)

anova(M3.lme,M3.lmeA)
anova(M3.lme,M3.lmeB)





M4.lme <- lme(IntakeRate ~ fSex, data= Limosa,
        weights = varIdent(form=~ 1 | fSex * fPeriod),
        random =~ 1 | fID, method = "ML")

M4.lmeA <- update(M4.lme, .~. -fSex)

anova(M4.lme,M4.lmeA)


M4.lme <- lme(IntakeRate ~ fSex, data= Limosa,
        weights = varIdent(form=~ 1 | fSex * fPeriod),
        random =~ 1 | fID, method = "REML")
summary(M4.lme)


0.064^2 / (0.064^2 + (0.136)^2)

0.064^2 / (0.064^2 + (0.436*0.136)^2)




M2.lm <- lm(IntakeRate ~ fSex, data= Limosa)
drop1(M2.lm, test = "F")


M5A.gls <- gls(IntakeRate ~ fSex, data= Limosa,
        weights = varIdent(form=~ 1 | fSex * fPeriod),
        method = "ML")

M5B.gls <- gls(IntakeRate ~ 1, data= Limosa,
        weights = varIdent(form=~ 1 | fSex * fPeriod),
        method = "ML")
anova(M5A.gls,M5B.gls)


        

summary(M4.lme)
