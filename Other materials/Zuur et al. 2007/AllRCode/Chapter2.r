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


                               



library(AED)
data(Nereis)
dotchart(Nereis$concentration, xlab = "Concentration",
    ylab = "Order of observations",
    main = "Cleveland dotplot")


boxplot(concentration ~ factor(nutrient),
  varwidth = TRUE, xlab = "nutrient",
    main = "Boxplot of concentration conditional on\
    nutrient", ylab = "concentration", data = Nereis)



TeethNitrogen<-read.table("C:/Bookdata/TeethNitrogen.txt",
       header=T)
library(lattice)
xyplot(X15N ~ Age | factor(Tooth), type = "l", col = 1,
  xlab = "Estimated age",
  ylab = expression(paste(delta^{15}, "N")),
  strip = function(bg = 'white', ...)
  strip.default(bg = 'white', ...),
  data = TeethNitrogen)


library(AED); data(Clams)
 Clams$LNAFD <- log(Clams$AFD)
 Clams$LNLENGTH <- log(Clams$LENGTH)
 Clams$fMONTH <- factor(Clams$MONTH)
 library(lattice)
 coplot(LNAFD ~ LNLENGTH | fMONTH,data = Clams)
 M1 <- lm(LNAFD ~ LNLENGTH * fMONTH,
          data = Clams)
 drop1(M1, test = "F")







op <- par(mfrow=c(2,2),mar=c(5,4,1,2))
#par(mfrow=c(2,2),mar=c(5,4,1,2))
plot(M1,add.smooth=F,which=1)
E<-resid(M1)
hist(E,xlab="Residuals",main="")
plot(Clams$LNLENGTH,E,xlab="Log(Length)",ylab="Residuals")
plot(Clams$fMONTH,E,xlab="Month",ylab="Residuals")
par(op)

E1<-E[Clams$LNLENGTH<=2.75]
E2<-E[Clams$LNLENGTH>2.75]
var.test(E1,E2)



bartlett.test(E,Clams$fMONTH)




library(AED) ; data(TeethNitrogen)
TN<-TeethNitrogen
tmp<-lm(X15N~Age,subset=TN$Tooth=="Moby",data=TN)

op <- par(mfrow=c(2,2))
plot(tmp,add.smooth=F)
par(op)

N.Moby<-TN$X15N[TN$Tooth=="Moby"]
Age.Moby<-TN$Age[TN$Tooth=="Moby"]

plot(y=N.Moby,x=Age.Moby,xlab="Estimated age Moby",
     ylab=expression(paste(delta^{15},"N Moby")))
abline(tmp)						#Figure 2.9

summary(tmp)




library(AED); data(Nereis)
Nereis$fbiomass=factor(Nereis$biomass)
Nereis$fnutrient=factor(Nereis$nutrient)
M3<-lm(concentration~fbiomass*fnutrient,data=Nereis)
drop1(M3, test="F")

op <- par(mfrow=c(1,2))
plot(resid(M3)~Nereis$fbiomass,xlab="Biomass",
     ylab="Residuals")
plot(resid(M3)~Nereis$fnutrient,xlab="Nutrient",
     ylab="Residuals")
par(op)



library(AED); data(ISIT)
ISIT$fStation<-factor(ISIT$Station)
library(lattice)
xyplot(Sources~SampleDepth|fStation,
 data=ISIT,
 xlab="Sample Depth",ylab="Sources",
 strip = function(bg='white', ...)
 strip.default(bg='white', ...),
 panel = function(x, y) {
                panel.grid(h=-1, v= 2)
                I<-order(x)
                llines(x[I], y[I],col=1)})



