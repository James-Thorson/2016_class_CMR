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



library(AED); data(Antarcticbirds)
ABirds <- Antarcticbirds #saves some space
library(lattice)
Birds<-c(ABirds$ArrivalAP,ABirds$LayingAP,ABirds$ArrivalCP,ABirds$LayingCP,
         ABirds$ArrivalEP,ABirds$LayingEP)
AllYears<-rep(ABirds$Year,6)
MyNames<-c("Arrival Adelie Penguin" ,
    "Laying Adelie penguin","Arrival Cape Petrel",
    "Laying Cape Petrel","Arrival Emperor Penguin",
    "Laying Emperor Penguin")
ID1<-factor(rep(MyNames,each=length(ABirds$Year)),
levels=c(MyNames[1], MyNames[3], MyNames[5],
         MyNames[2], MyNames[4], MyNames[6]))
xyplot(Birds~AllYears|ID1, xlab="Years",ylab="Day",
      layout=c(3,2),
      strip = function(bg='white', ...)
      strip.default(bg='white', ...),
      scales = list(alternating = T,
                   x = list(relation = "same"),
                   y = list(relation = "free")),
      panel=function(x,y){
       panel.xyplot(x,y,col=1)
       panel.loess(x,y,col=1,span=0.5)
       panel.grid(h=-1,v=2)})

L.AP <- acf(ABirds$LayingAP, lag.max = 10,
      na.action =na.pass,
                   main="Laying dates Adelie Penguin")




ABirds$DifAP<-ABirds$LayingAP-ABirds$ArrivalAP
ABirds$DifCP<-ABirds$LayingCP-ABirds$ArrivalCP
ABirds$DifEP<-ABirds$LayingEP-ABirds$ArrivalEP
AllDif<-c(ABirds$DifAP,ABirds$DifCP,ABirds$DifEP,ABirds$MSA,ABirds$SOI)
AllYear<-rep(ABirds$Year,5)
IDDif<-rep(c("Difference AP","Difference CP",
              "Difference EP", "MSA","SOI"), each=55)
xyplot(AllDif ~ AllYear | IDDif, xlab = "Years", ylab = "Day",
      layout=c(3, 2),
      strip = function(bg='white', ...)
      strip.default(bg='white', ...),
      scales = list(alternating = T,
                    x = list(relation = "same"),
                    y = list(relation = "free")),
      panel=function(x,y){
        panel.xyplot(x,y,col=1)
        panel.loess(x,y,col=1,span=0.5)
        panel.grid(h=-1,v=2)})


library(mgcv)
library(nlme)
B1<-gamm(ArrivalAP~s(Year), correlation =
                     corARMA(form=~Year,p=1,q=0),
                     data = ABirds)
AIC(B1$lme)


M1<-lm(LayingAP~MSA, data = ABirds)
M2<-lm(LayingCP~MSA, data = ABirds)
M3<-lm(LayingEP~MSA, data = ABirds)
M4<-lm(ArrivalCP~MSA, data = ABirds)
summary(M1);summary(M2)
summary(M3);summary(M4)






Bird4<-c(ABirds$LayingAP,ABirds$LayingCP,ABirds$LayingEP,ABirds$ArrivalCP)
MSA4<-rep(ABirds$MSA,4)
ID4<-rep(c("Laying Adelie penguin",
      "Laying Cape Petrel","Laying Emperor Penguin",
      "Arrival Cape Petrel"),each=55)
xyplot(Bird4~MSA4|ID4, xlab="MSA",ylab="Day",
       layout=c(2,2),
     strip = function(bg='white', ...) strip.default(bg='white', ...),
     scales = list(alternating = T,
                 x = list(relation = "same"),
                y = list(relation = "free")),
     panel=function(x,y,subscripts,...){
       panel.xyplot(x,y,col=1)
       panel.grid(h=-1,v=2)
       I<-!is.na(y) & !is.na(x)
       tmp<-lm(y[I]~x[I]);  x1<-x[I];
       y1<-fitted(tmp);  I2<-order(x1)
       panel.lines(x1[I2],y1[I2],col=1,span=1) })



AP <- c(ABirds$ArrivalAP, ABirds$LayingAP)
SOI2 <- c(ABirds$SOI, ABirds$SOI)
Y2 <- c(ABirds$Year, ABirds$Year)
ID <- factor(rep(c("Arrival", "Laying"), each = 55))
library(nlme)
vf2 <- varIdent(form =~ 1 | ID)
M5 <- gls(AP ~ SOI2 + ID + SOI2:ID, weights = vf2,
        na.action = na.omit)
M6 <- gls(AP ~ SOI2 + ID + SOI2:ID, weights = vf2,
        na.action = na.omit,
        correlation = corAR1(form =~ Y2 | ID))
anova(M5, M6)




M7<-gls(AP ~ SOI2 + ID + SOI2:ID, weights=vf2,
                    na.action = na.omit,
        method = "ML", correlation = corAR1(form=~Y2|ID))
M8<-gls(AP~SOI2+ID,weights=vf2,na.action=na.omit,
         method="ML",correlation=corAR1(form=~Y2|ID))
anova(M7, M8)


M9 <- gls(AP ~ ID, weights = vf2, na.action = na.omit,
         method = "ML", correlation=corAR1(form=~Y2|ID))
summary(M9)











CP<-c(ABirds$ArrivalCP,ABirds$LayingCP)
SOI2<-c(ABirds$SOI,ABirds$SOI)
Y2<-c(ABirds$Year,ABirds$Year)
ID<-factor(rep(c("Arrival","Laying"),each=55))
vf2<-varIdent(form= ~ 1|ID)


M10<-gls(CP~SOI2+ID+SOI2:ID,weights=vf2,na.action=na.omit,method="ML")


plot(ABirds$SOI, ABirds$ArrivalCP,ylim=c(195,270),type="n",
   ylab="Arrival & laying dates")
points(ABirds$SOI,ABirds$ArrivalCP,pch=1)
points(ABirds$SOI,ABirds$LayingCP,pch=2)
MyX=data.frame(SOI2=seq(from=min(ABirds$SOI),to=max(ABirds$SOI),
                        length=20),
               ID="Arrival")
Pred1<-predict(M10,newdata=MyX)
lines(MyX$SOI2,Pred1)
MyX=data.frame(SOI2=seq(from=min(ABirds$SOI),to=max(ABirds$SOI),
                        length=20),
               ID="Laying")
Pred2<-predict(M10,newdata=MyX)
lines(MyX$SOI2,Pred2)


