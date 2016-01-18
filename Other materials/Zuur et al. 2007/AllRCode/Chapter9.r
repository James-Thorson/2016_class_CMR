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




#Code for Chapter 9
library(AED); data(RoadKills)
RK <- RoadKills
names(RK)

#"Sector"       "X"            "Y"            "BufoCalamita"
#"TOT_N"
# [6] "S_RICH"       "OPEN_L"       "OLIVE"        "MONT_S"
#"MONT"
#[11] "POLIC"        "SHRUB"        "URBAN"        "WAT.RES"
#"L.WAT.C"
#[16] "L.D.ROAD"     "L.P.ROAD"     "D.WAT.RES"    "D.WAT.COUR"
#"D.PARK"
#[21] "N.PATCH"      "P.EDGE"       "L.SDI"


plot(RK$D.PARK,RK$TOT.N,xlab="Distance to park",
     ylab="Road kills")
M1<-glm(TOT.N~D.PARK,family=poisson,data=RK)
summary(M1)

MyData=data.frame(D.PARK=seq(from=0,to=25000,by=1000))
G<-predict(M1,newdata=MyData,type="link",se=T)
F<-exp(G$fit)
FSEUP<-exp(G$fit+1.96*G$se.fit)
FSELOW<-exp(G$fit-1.96*G$se.fit)
lines(MyData$D.PARK,F,lty=1)
lines(MyData$D.PARK,FSEUP,lty=2)
lines(MyData$D.PARK,FSELOW,lty=2)

RK$SQ.POLIC<-sqrt(RK$POLIC)
RK$SQ.WATRES<-sqrt(RK$WAT.RES)
RK$SQ.URBAN<-sqrt(RK$URBAN)
RK$SQ.OLIVE<-sqrt(RK$OLIVE)
RK$SQ.LPROAD<-sqrt(RK$L.P.ROAD)
RK$SQ.SHRUB<-sqrt(RK$SHRUB)
RK$SQ.DWATCOUR<-sqrt(RK$D.WAT.COUR)
M2<-glm(TOT.N~OPEN.L+MONT.S+SQ.POLIC+
         SQ.SHRUB+SQ.WATRES+L.WAT.C+SQ.LPROAD+
         SQ.DWATCOUR+D.PARK,family=poisson,data=RK)
summary(M2)

drop1(M2,test="Chi")

M3 <- glm(TOT.N ~ MONT.S + SQ.POLIC + D.PARK +
          SQ.SHRUB + SQ.WATRES + L.WAT.C + SQ.LPROAD +
          SQ.DWATCOUR, family = poisson, data = RK)
anova(M2, M3, test = "Chi")


M4<- glm(TOT.N ~ OPEN.L + MONT.S + SQ.POLIC+
         SQ.SHRUB + SQ.WATRES + L.WAT.C + SQ.LPROAD+
         SQ.DWATCOUR + D.PARK, family = quasipoisson, data = RK)


summary(M4)
drop1(M4,test="F")


M5<- glm(TOT.N ~ D.PARK, family = quasipoisson, data = RK)


summary(M5)
drop1(M5,test="F")



M5<-glm(TOT.N~D.PARK,family=quasipoisson,data=RK)
EP=resid(M5,type="pearson")
ED=resid(M5,type="deviance")
mu=predict(M5,type="response")
E=RK$TOT.N-mu
EP2=E/sqrt(7.630148*mu)
op <- par(mfrow=c(2,2))
plot(x=mu,y=E,main="Response residuals")
plot(x=mu,y=EP,main="Pearson residuals")
plot(x=mu,y=EP2,main="Pearson residuals scaled")
plot(x=mu,y=ED,main="Deviance residuals")
par(op)


library(MASS)
M6<-glm.nb(TOT.N~OPEN.L+MONT.S+SQ.POLIC+
         SQ.SHRUB+SQ.WATRES+L.WAT.C+SQ.LPROAD+
         SQ.DWATCOUR+D.PARK,link="log",data=RK)

M8<-glm.nb(TOT.N~OPEN.L+D.PARK,link="log",data=RK)
summary(M8)
drop1(M8,test="Chi")
par(mfrow=c(2,2))
plot(M8)


M9 <- glm(TOT.N ~ OPEN.L + D.PARK, family = poisson, data = RK)


llhNB = logLik(M8)
llhPoisson  =logLik(M9)
d <- 2 * (llhNB - llhPoisson)
pval <- pchisq(as.numeric(d), df=1, lower.tail=FALSE)/2



#Sea lice data

library(AED); data(Lice)
Lice$LVol <- log(Lice$Volume)
Lice$fStation <- factor(Lice$Station)
L0 <- glm(Copepod ~ offset(LVol) + fStation,
       family=poisson,data=Lice)



library(mgcv)
Lice$PW <- Lice$Production_week
Lice$fDepth <- factor(Lice$Depth)

#This code was used in R versions 2.7.0.
L1<-gam(Copepod~offset(LVol)+
      s(PW,by=as.numeric(Depth=="0m" & Station=="A"))+
      s(PW,by=as.numeric(Depth=="0m" & Station=="C"))+
      s(PW,by=as.numeric(Depth=="0m" & Station=="E"))+
      s(PW,by=as.numeric(Depth=="0m" & Station=="F"))+
      s(PW,by=as.numeric(Depth=="0m" & Station=="G"))+
      s(PW,by=as.numeric(Depth=="5m" & Station=="A"))+
      s(PW,by=as.numeric(Depth=="5m" & Station=="C"))+
      s(PW,by=as.numeric(Depth=="5m" & Station=="E"))+
      s(PW,by=as.numeric(Depth=="5m" & Station=="F"))+
      s(PW,by=as.numeric(Depth=="5m" & Station=="G"))+
     fDepth*fStation,
      family=negative.binomial(1),data=Lice)

#For higher R versions, use:

Lice$DS <- vector(length=nrow(Lice))
Lice$DS[Lice$Depth=="0m" & Lice$Station=="A"] <- "0mA"
Lice$DS[Lice$Depth=="0m" & Lice$Station=="C"] <- "0mC"
Lice$DS[Lice$Depth=="0m" & Lice$Station=="E"] <- "0mE"
Lice$DS[Lice$Depth=="0m" & Lice$Station=="F"] <- "0mF"
Lice$DS[Lice$Depth=="0m" & Lice$Station=="G"] <- "0mG"
Lice$DS[Lice$Depth=="5m" & Lice$Station=="A"] <- "5mA"
Lice$DS[Lice$Depth=="5m" & Lice$Station=="C"] <- "5mC"
Lice$DS[Lice$Depth=="5m" & Lice$Station=="E"] <- "5mE"
Lice$DS[Lice$Depth=="5m" & Lice$Station=="F"] <- "5mF"
Lice$DS[Lice$Depth=="5m" & Lice$Station=="G"] <- "5mG"
Lice$fDS <- factor(Lice$DS)


L1<-gam(Copepod~offset(LVol) + s(PW,by=fDS) + fDepth*fStation,
      family=negbin(c(1,10),link = log),
      data=Lice)

#Take a cup of coffee
#See also the helpfile of negbin
## unknown theta via outer iteration and AIC search...
## how to retrieve Theta...
L1$family$getTheta()




#In older R versions, use:
L3 <- gam(Copepod ~ offset(LVol) +
       s(PW, by = as.numeric(Depth=="0m")) +
       s(PW, by = as.numeric(Depth=="5m")) +
       fDepth + fStation, data = Lice,
       family = negative.binomial(1), gamma = 1.4)

L4 <- gam(Copepod ~ offset(LVol) +
       s(PW,by = as.numeric(Depth=="0m")) +
       s(PW,by = as.numeric(Depth=="5m")) +
       fDepth + fStation, data = Lice,
       family = poisson, gamma = 1.4)


#And in more recent R versions, use:
L3 <- gam(Copepod ~ offset(LVol) +
       s(PW, by = fDepth) + fDepth + fStation, data = Lice,
       family=negbin(c(1,10),link = log),
       gamma = 1.4)

L4 <- gam(Copepod ~ offset(LVol) +
       s(PW, by = fDepth) + fDepth + fStation, data = Lice,
       family = poisson, gamma = 1.4)




llhNB <- logLik(L3)
llhPoisson <- logLik(L4)
d <- 2 * (llhNB - llhPoisson)
pval <- 0.5 * pchisq(as.numeric(d), df = 1, lower.tail = FALSE)





#Figure 9.2
library(scatterplot3d)

x <- seq(0, 100)
y <- exp(0.01+0.03*x)
y
z <- 0*x

ymeas=rpois(length(y),lambda=y)
plot(x,ymeas,type="p",xlab="Covariate",ylab="Observed values")
lines(x,y)

rr=scatterplot3d(x, y, z, highlight.3d=TRUE, col.axis="black",
      col.grid="black", pch=20,zlim=c(0,0.4),type="l",lwd=3,
      xlab="Covariate",ylab="Possible values",zlab="Probability")


MyX=c(2,15,30,50,75)
for (i in 1:5){
  xi=MyX[i]
  yi=exp(0.01+0.03*xi)
  #if ( xi > 10) {yseq=round(seq(yi-0.8*yi,yi+0.8*yi,step=0.01))}
  if (xi <=100) {yseq=round(seq(0,20,by=0.1))}
  zi=dpois(yseq, lambda=yi)
  rb=cbind(xi,yseq,zi)
  rr$points3d(rb, col = 1,type="h",pch=30)
  }

