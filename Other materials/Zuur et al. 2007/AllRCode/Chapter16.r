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



library(AED); data(RoadKills)
RK <- RoadKills
names(RK)

# [1] "Sector"       "X"            "Y"            "BufoCalamita" "TOT.N"
# [6] "S.RICH"       "OPEN.L"       "OLIVE"        "MONT.S"       "MONT"
#[11] "POLIC"        "SHRUB"        "URBAN"        "WAT.RES"      "L.WAT.C"
#[16] "L.D.ROAD"     "L.P.ROAD"     "D.WAT.RES"    "D.WAT.COUR"   "D.PARK"
#[21] "N.PATCH"      "P.EDGE"       "L.SDI"

library(lattice)
xyplot(Y ~ X, aspect = "iso", col=1, pch=16, data = RK)

RK$SQ.POLIC <- sqrt(RK$POLIC)
RK$SQ.WATRES <- sqrt(RK$WAT.RES)
RK$SQ.URBAN <- sqrt(RK$URBAN)
RK$SQ.OLIVE <- sqrt(RK$OLIVE)
RK$SQ.LPROAD <- sqrt(RK$L.P.ROAD)
RK$SQ.SHRUB <- sqrt(RK$SHRUB)
RK$SQ.DWATCOUR <- sqrt(RK$D.WAT.COUR)
Z<-cbind(RK$OPEN.L, RK$SQ.OLIVE, RK$MONT.S,RK$MONT, RK$SQ.POLIC,
     RK$SQ.SHRUB, RK$SQ.URBAN, RK$SQ.WATRES, RK$L.WAT.C,
     RK$L.D.ROAD, RK$SQ.LPROAD,RK$D.WAT.RES, RK$SQ.DWATCOUR,
     RK$D.PARK, RK$N.PATCH, RK$P.EDGE, RK$L.SDI)
     
corvif(Z)



X12<-c(RK$OPEN.L, RK$SQ.OLIVE, RK$MONT.S, RK$SQ.POLIC, RK$SQ.SHRUB,
     RK$SQ.WATRES, RK$L.WAT.C, RK$L.D.ROAD, RK$SQ.LPROAD, RK$D.WAT.RES,
     RK$SQ.DWATCOUR, RK$D.PARK)
Killings12<-rep(RK$TOT.N, 12)
I12<-rep(c("OPEN.L","OLIVE","MONT.S","POLIC","SHRUB",
           "WATRES","L.WAT.C","L.D.ROAD","L.P.ROAD",
           "D.WAT.RES","D.WAT.COUR","D.PARK"),each=52)
ID12<-rep(I12,12)
library(lattice)
xyplot(Killings12~X12|ID12,col=1,
  strip = function(bg='white', ...)
    strip.default(bg='white', ...),
  scales = list(alternating = T,
     x = list(relation = "free"),
     y = list(relation = "same")),
  xlab="Explanatory variables",
  ylab="Number of roadkillings",
  panel=function(x,y){
    panel.grid(h=-1, v= 2)
    panel.points(x,y,col=1)
    panel.loess(x,y,col=1,lwd=2)})
    
    
X11<-c(RK$OPEN.L, RK$SQ.OLIVE, RK$MONT.S, RK$SQ.POLIC, RK$SQ.SHRUB,
       RK$SQ.WATRES, RK$L.WAT.C, RK$L.D.ROAD, RK$SQ.LPROAD, RK$D.WAT.RES,
       RK$SQ.DWATCOUR)

I11<-rep(c("OPEN.L","OLIVE","MONT.S","POLIC","SHRUB","WATRES","L.WAT.C",
       "L.D.ROAD","L.P.ROAD","D.WAT.RES","D.WAT.COUR"),each=52)
ID11<-rep(I11, 11)

RK$D.PARK.KM<-RK$D.PARK/1000
Space11<-rep(RK$D.PARK.KM, 11)

library(lattice)
xyplot(X11~Space11|ID11,col=1,
  strip = function(bg='white', ...) strip.default(bg='white', ...),
  scales = list(alternating = T, x = list(relation = "same"),y = list(relation = "free")),
  xlab="Distance from first point",ylab="Explanatory variables",
  panel=function(x,y){
    panel.grid(h=-1, v= 2)
    panel.points(x,y,col=1)
    panel.loess(x,y,col=1,lwd=2)
    })



library(mgcv)
library(MASS)


#For older R versions (2.7.0 and before), use the construction:
#family=negative.binomial(1)

#For newer R versions, use the construction:
#family=negbin(c(1,10),link = log)

#This is for older R versions
M1<-gam(TOT.N~s(OPEN.L)+s(MONT.S)+s(SQ.POLIC)+
         s(SQ.SHRUB)+s(SQ.WATRES)+s(L.WAT.C)+
         s(SQ.LPROAD)+s(SQ.DWATCOUR)+s(D.PARK),
         family=negative.binomial(1), data = RK)


M2 <- gam(TOT.N ~ s(OPEN.L) + s(D.PARK),
          family=negative.binomial(1), data = RK)



#And this is for newer versions
M1<-gam(TOT.N~s(OPEN.L)+s(MONT.S)+s(SQ.POLIC)+
         s(SQ.SHRUB)+s(SQ.WATRES)+s(L.WAT.C)+
         s(SQ.LPROAD)+s(SQ.DWATCOUR)+s(D.PARK),
         family=negbin(c(1,10),link = log), data = RK)
#Too many covariates..hence the error message

M2 <- gam(TOT.N ~ s(OPEN.L) + s(D.PARK),
          family=negbin(c(1,10),link = log), data = RK)




anova(M2)


E<-resid(M2,type="pearson")
I<-vector(length=length(E))
I[E<0]<-1
I[E>=0]<-16
library(lattice)
xyplot(Y~X,cex=2*abs(E)/max(abs(E)),pch=I,col=1,data=RK)




X12<-c(RK$OPEN.L, RK$SQ.OLIVE, RK$MONT.S, RK$SQ.POLIC, RK$SQ.SHRUB,
     RK$SQ.WATRES, RK$L.WAT.C, RK$L.D.ROAD, RK$SQ.LPROAD, RK$D.WAT.RES,
     RK$SQ.DWATCOUR, RK$D.PARK)
Killings12<-rep(E, 12)
I12<-rep(c("OPEN.L","OLIVE","MONT.S","POLIC","SHRUB",
           "WATRES","L.WAT_C","L.D.ROAD","L.P.ROAD",
           "D.WAT.RES","D.WAT.COUR","D.PARK"),each=52)
ID12<-rep(I12,12)
library(lattice)
xyplot(Killings12~X12|ID12,col=1,
  strip = function(bg='white', ...)
    strip.default(bg='white', ...),
  scales = list(alternating = T,
     x = list(relation = "free"),
     y = list(relation = "same")),
  xlab="Explanatory variables",
  ylab="Pearson residuals",
  panel=function(x,y){
    panel.grid(h=-1, v= 2)
    panel.points(x,y,col=1)
    panel.loess(x,y,col=1,lwd=2)})



#For older R versions:
M3<-gam(TOT.N~ s(D.PARK),
           family=negative.binomial(1), data = RK)

#For newer R versions:
M3<-gam(TOT.N~ s(D.PARK),
           family=negbin(c(1,15),link = log), data = RK)



           
           

M3Pred<-predict(M3, se = TRUE, type = "response")
plot(RK$D.PARK, RK$TOT.N, cex = 1.1, pch = 16,
      main = "Negative binomial GAM",
      xlab="Distance to park", ylab="Number of road killings")
I <- order(RK$D.PARK)
lines(RK$D.PARK[I], M3Pred$fit[I], lwd = 2)
lines(RK$D.PARK[I], M3Pred$fit[I] + 2 * M3Pred$se.fit[I],
      lty=2,lwd=2)
lines(RK$D.PARK[I], M3Pred$fit[I] - 2 * M3Pred$se.fit[I],
      lty = 2, lwd = 2)
for (i in 1:52){
   y <- rnbinom(100, size = 11.8, mu = M3Pred$fit[i])
   points(rep(RK$D.PARK[i], 100), y, cex = 0.5)
}


#Figure 16.9
MyI<-seq(1,52,by=5)
x<-seq(0,80)

YY<-1
for (i in MyI){
   y<-dnbinom(x,size=11.8,mu=M3Pred$fit[i])
   YY<-c(YY,y)
}

YY<-YY[-1]
XX<-rep(x,11)
ID<-rep(RK$D.PARK[MyI],each=80)
ID2<-factor(paste("D.PARK = ",ID,sep=""),
 levels=c(
 "D.PARK = 250.214",   "D.PARK = 2724.089",  "D.PARK = 5202.328",
 "D.PARK = 7668.833",  "D.PARK = 10047.63",  "D.PARK = 12470.968" ,
 "D.PARK = 14904.995", "D.PARK = 17235.045", "D.PARK = 19645.717" ,
 "D.PARK = 22119.102", "D.PARK = 24444.874")

)

xyplot(YY~XX|ID2,type="h",col=1,xlab="Number of road killings",
      ylab="Probabilities")



library(nlme)
RK$D.PARK.KM <- RK$D.PARK / 1000




#For older R versions:
M4 <- gamm(TOT.N ~ s(OPEN.L)+s(D.PARK),
            family = negative.binomial(theta=11.8), data = RK)
M4Var<-Variogram(M4$lme, form =~ D.PARK.KM, nugget=TRUE, data = RK)
plot(M4Var,col=1, smooth=FALSE)

M5 <- gamm(TOT.N ~ s(OPEN.L) + s(D.PARK), data = RK,
           family = negative.binomial(theta=11.8),
           correlation = corGaus(form =~ D.PARK.KM, nugget = TRUE))





#For newer R versions
M4 <- gamm(TOT.N ~ s(OPEN.L)+s(D.PARK),
            family=negbin(theta=14.0,link = log), data = RK)
            
            
M4Var<-Variogram(M4$lme, form =~ D.PARK.KM, nugget=TRUE, data = RK)
plot(M4Var,col=1, smooth=FALSE)

M5 <- gamm(TOT.N ~ s(OPEN.L) + s(D.PARK), data = RK,
           family=negbin(theta=14.0,link = log),
           correlation = corGaus(form =~ D.PARK.KM, nugget = TRUE))






