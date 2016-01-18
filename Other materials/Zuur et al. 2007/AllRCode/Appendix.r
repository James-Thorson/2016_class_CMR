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




library(AED); data(Loyn)
Loyn$fGRAZE <- factor(Loyn$GRAZE)

op<- par(mfrow=c(4,2),mar=c(3,3,3,1))
dotchart(Loyn$ABUND,main="ABUND",group=Loyn$fGRAZE)
plot(0,0,type="n",axes=F)
dotchart(Loyn$AREA,main="AREA",group=Loyn$fGRAZE)
dotchart(Loyn$DIST,main="DIST",group=Loyn$fGRAZE)
dotchart(Loyn$LDIST,main="LDIST",group=Loyn$fGRAZE)
dotchart(Loyn$YR.ISOL,main="YR.ISOL",group=Loyn$fGRAZE)
dotchart(Loyn$ALT,main="ALT",group=Loyn$fGRAZE)
dotchart(Loyn$GRAZE,main="GRAZE",group=Loyn$fGRAZE)
par(op)



Loyn$L.AREA<-log10(Loyn$AREA)
Loyn$L.DIST<-log10(Loyn$DIST)
Loyn$L.LDIST<-log10(Loyn$LDIST)


Z<-cbind(Loyn$ABUND,Loyn$L.AREA,Loyn$L.DIST,Loyn$L.LDIST,Loyn$YR.ISOL,Loyn$ALT,Loyn$GRAZE)
colnames(Z)<-c("ABUND","L.AREA","L.DIST","L.LDIST","YR.ISOL","ALT","GRAZE")
pairs(Z, lower.panel=panel.smooth2,
         upper.panel=panel.cor,diag.panel=panel.hist)


corvif(Z[,c(-1,-7)])



M1 <- lm(ABUND ~ L.AREA + L.DIST + L.LDIST + YR.ISOL + ALT +
     fGRAZE, data = Loyn)

summary(M1)
anova(M1)
drop1(M1,test="F")

#Verify drop1

M1A <- lm(ABUND ~ L.AREA + L.DIST + L.LDIST + YR.ISOL + ALT +
     fGRAZE, data = Loyn)

M1B <- lm(ABUND ~ L.DIST + L.LDIST + YR.ISOL + ALT +
     fGRAZE, data = Loyn)

anova(M1A,M1B)


#Model 1: ABUND ~ L.AREA + L.DIST + L.LDIST + YR.ISOL + ALT + fGRAZE
#Model 2: ABUND ~ L.DIST + L.LDIST + YR.ISOL + ALT + fGRAZE
#  Res.Df     RSS Df Sum of Sq      F   Pr(>F)
#1     46 1714.43
#2     47 2484.44 -1   -770.01 20.660 3.97e-05 ***

((2484.44-1714.43)/1) / (1714.43/(56-9))




#Verify anova table
M1A <- lm(ABUND ~ L.AREA + L.DIST + L.LDIST + YR.ISOL + ALT +
          fGRAZE, data = Loyn)

M1B <- lm(ABUND ~ L.AREA + L.DIST + L.LDIST + YR.ISOL ,
          data = Loyn)
anova(M1A,M1B)

M1C <- lm(ABUND ~ 1, data = Loyn)
anova(M1C)


M1C <- lm(ABUND ~ fGRAZE+L.AREA + L.DIST + L.LDIST + YR.ISOL + ALT, data = Loyn)
anova(M1A)
anova(M1C)


#Model 1: ABUND ~ L.AREA
#Model 2: ABUND ~ L.AREA + L.DIST
#  Res.Df     RSS Df Sum of Sq      F Pr(>F)
#1     54 2866.94
#2     53 2801.47  1     65.48 1.2388 0.2707

((2866.94-2801.47)/1) / (2801.47/(56-3))
#ok



#Analysis of Variance Table
#
#Response: ABUND
#          Df Sum Sq Mean Sq F value    Pr(>F)
#L.AREA     1 3471.0  3471.0 93.1303 1.247e-12 ***
#L.DIST     1   65.5    65.5  1.7568  0.191565





M2 <- lm(ABUND ~ L.AREA + L.DIST + L.LDIST + YR.ISOL + ALT, data = Loyn)
anova(M1, M2)
step(M1)


M3 <- lm(ABUND ~ L.AREA + fGRAZE, data = Loyn)
op <- par(mfrow = c(2, 2))
plot(M3) #standard graphical output
win.graph()
op <- par(mfrow = c(2, 2))
#Check for normality
E <- rstandard(M3)
hist(E)
#qqnorm(E)
#Check for independence: residuals versus individual #explanatory variables
plot(y = E, x = Loyn$L.AREA, xlab = "AREA", ylab = "Residuals")
abline(0,0)
plot(E ~ Loyn$fGRAZE, xlab = "GRAZE", ylab = "Residuals")
abline(0, 0)
par(op)




plot(Loyn$L.AREA,Loyn$ABUND)
#Loyn$SL.AREA<-sort(Loyn$L.AREA)
D1<-data.frame(L.AREA=Loyn$L.AREA[Loyn$GRAZE==1],fGRAZE="1")
D2<-data.frame(L.AREA=Loyn$L.AREA[Loyn$GRAZE==2],fGRAZE="2")
D3<-data.frame(L.AREA=Loyn$L.AREA[Loyn$GRAZE==3],fGRAZE="3")
D4<-data.frame(L.AREA=Loyn$L.AREA[Loyn$GRAZE==4],fGRAZE="4")
D5<-data.frame(L.AREA=Loyn$L.AREA[Loyn$GRAZE==5],fGRAZE="5")

P1<-predict(M3,newdata=D1)
P2<-predict(M3,newdata=D2)
P3<-predict(M3,newdata=D3)
P4<-predict(M3,newdata=D4)
P5<-predict(M3,newdata=D5)


D1<-data.frame(L.AREA = Loyn$L.AREA, fGRAZE = "1")
P1<-predict(M3, newdata = D1)



lines(D1$L.AREA,P1,lty=1)
lines(D2$L.AREA,P2,lty=2)
lines(D3$L.AREA,P3,lty=3)
lines(D4$L.AREA,P4,lty=4)
lines(D5$L.AREA,P5,lty=5)




library(mgcv)
AM1<-gam(ABUND~s(L.AREA)+s(L.DIST)+s(L.LDIST)+
    s(YR.ISOL)+s(ALT)+fGRAZE, data = Loyn)
anova(AM1)


AM2<-gam(ABUND ~ s(L.AREA, bs = "cs") + s(L.DIST, bs = "cs") +
                 s(L.LDIST,bs = "cs") + s(YR.ISOL, bs = "cs") +
                 s(ALT, bs = "cs") + fGRAZE, data = Loyn)
anova(AM2)



AM3 <- gam(ABUND ~ s(L.AREA, bs = "cs") + fGRAZE, data = Loyn)
plot(AM3)


E.AM3 <- resid(AM3)
Fit.AM3 <- fitted(AM3)
plot(x = Fit.AM3, y = E.AM3, xlab = "Fitted values",
       ylab = "Residuals")




M3<-lm(ABUND ~ L.AREA + fGRAZE, data = Loyn)
AM3<-gam(ABUND ~ s(L.AREA, bs = "cs") + fGRAZE, data = Loyn)
anova(M3, AM3)

