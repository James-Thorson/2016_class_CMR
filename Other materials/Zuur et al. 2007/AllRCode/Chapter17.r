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



library(AED); data(ISIT)
ISIT$fMonth <- factor(ISIT$Month)
ISIT$fStation <- factor(ISIT$Station)
ISIT$fYear <- factor(ISIT$Year)

ISIT2 <- ISIT[ISIT$fStation != "4" & ISIT$fStation != "5" & ISIT$fStation != "10",]
library(lattice)
MyLines <- function(xi, yi, ...){
   I <- order(xi)
   panel.lines(xi[I], yi[I], col = 1)}
xyplot(Sources ~ SampleDepth | fMonth,
 groups = fStation,
 data = ISIT2, xlab = "Depth", ylab = "Sources",
 panel = panel.superpose,
 panel.groups = MyLines)




 xyplot(Sources~SampleDepth / 1000 |
  fMonth * fYear,
  groups = fStation, data = ISIT2,
  strip = function(bg = 'white', ...)
  strip.default(bg = 'white', ...),
  scales = list(alternating = TRUE, x = list(relation =
  "same"),y = list(relation = "free")),
  xlab = "Depth (km)", ylab = expression(paste(Sources,
  " m"^{-3}, "")),
  panel = panel.superpose, panel.groups = MyLines)




#Approach 2
  #Extract coordinates
Xcoord<-vector(length=16)
Ycoord<-vector(length=16)
UStation<-unique(ISIT2$Station)
for (i in 1:16) {
 Xcoord[i]<-ISIT2$Xkm[UStation[i]==ISIT2$Station][1]
 Ycoord[i]<-ISIT2$Ykm[UStation[i]==ISIT2$Station][1]
}
#Calculate a distance matrix between the 16 stations
#using Pythagoras
D<-matrix(nrow=16,ncol=16)
for (i in 1:16){
  for (j in 1:16){
   D[i,j]<-sqrt((Ycoord[i]-Ycoord[j])^2+
        (Xcoord[i]-Xcoord[j])^2)}}

colnames(D)<-unique(ISIT2$Station)
rownames(D)<-unique(ISIT2$Station)
MyNames<-unique(ISIT2$Station)
#Apply clustering
Dist<-as.dist(D)
hc <- hclust(Dist, "ave")
plot(hc,labels=MyNames)



#Approach 3
library(mgcv)

MyDepth<-seq(from = min(ISIT2$SampleDepth),
             to = max(ISIT2$SampleDepth), by = 25)
NEWSOURCES <- matrix(nrow = 175,ncol = 16)
NEWSOURCES[] <- NA
j <- 1
for (i in MyNames){
  Mi <- gam(Sources ~ s(SampleDepth), subset = (ISIT2$Station==i), data = ISIT2)
  I <- MyDepth > min(ISIT2$SampleDepth[ISIT2$Station == i]) &
       MyDepth < max(ISIT2$SampleDepth[ISIT2$Station==i])
  mynewXdata <- data.frame(SampleDepth = MyDepth[I])
  tmp.pred <- predict(Mi, newdata = mynewXdata)
  NEWSOURCES[I,j] <- tmp.pred
  j <- j + 1
}

D<-cor(NEWSOURCES,use="pairwise.complete.obs")
colnames(D)<-unique(ISIT2$Station)
rownames(D)<-unique(ISIT2$Station)
Dist<-as.dist(1-D)
hc <-hclust(Dist, "ave")
plot(hc,labels=MyNames)



#Section 17.4
library(mgcv)
M1 <- gam(Sources ~ fStation + s(SampleDepth) +
     fMonth, data =ISIT2)
E <- resid(M1)
F <- fitted(M1)
op <- par(mfrow = c(2, 1), mar = c(5, 4, 1, 1))
plot(M1)
plot(F, E, xlab = "Fitted values", ylab = "Residuals")
par(op)


library(nlme)
ISIT2$Depth1000 <- ISIT2$SampleDepth / 1000
lmc <- lmeControl(niterEM = 5000, msMaxIter = 1000)
f1 <- formula(Sources ~ s(Depth1000) + fMonth)
M17.4A <- gamm(f1,random = list(fStation =~ 1),
      method = "REML", control = lmc, data = ISIT2)
M17.4B <- gamm(f1, random = list(fStation =~ 1),
      method = "REML", control = lmc, data = ISIT2,
      weights = varIdent(form =~ 1 | fStation))
M17.4C <- gamm(f1, random=list(fStation=~1), data = ISIT2,
     method = "REML", control = lmc,
      weights = varPower(form =~ Depth1000))
M17.4D<-gamm(f1, random = list(fStation =~ 1), data = ISIT2,
      method = "REML", control = lmc,
      weights=varComb(varIdent(form =~ 1 | fStation),
                    varPower(form =~ Depth1000)))
M17.4E<-gamm(f1, random=list(fStation =~ 1), data = ISIT2,
      method = "REML",control=lmc,
      weights = varComb(varIdent(form =~ 1 | fStation),
                    varPower(form =~ Depth1000 | fStation)))
                    
AIC(M17.4A$lme,M17.4B$lme,M17.4C$lme,M17.4D$lme,M17.4E$lme)



Vario17.4E<-Variogram(M17.4E$lme,form =~ Depth1000 | fStation,
   robust=TRUE, data = ISIT2)
plot(Vario17.4E)






#For older R versions (<2.7.0), use:
f1 <- formula(Sources ~
      s(Depth1000, by = as.numeric(Month == 3))+
      s(Depth1000, by = as.numeric(Month == 4))+
      s(Depth1000, by = as.numeric(Month == 8))+
      s(Depth1000, by = as.numeric(Month == 10))+
      fMonth)

#For more recent versions of R, use:
f1 <- formula(Sources ~ s(Depth1000, by = fMonth) + fMonth)







M17.5A <- gamm(f1,random = list(fStation =~ 1),
      method = "REML", control = lmc, data = ISIT2)

M17.5B <- gamm(f1, random = list(fStation =~ 1),
      method = "REML", control = lmc, data = ISIT2,
      weights = varIdent(form =~ 1 | fMonth))

M17.5C <- gamm(f1, random=list(fStation=~1), data = ISIT2,
     method = "REML", control = lmc,
      weights = varPower(form =~ Depth1000))
      
M17.5D<-gamm(f1, random = list(fStation =~ 1), data = ISIT2,
      method = "REML", control = lmc,
      weights=varComb(varIdent(form =~ 1 | fMonth),
                    varPower(form =~ Depth1000)))
                    
M17.5E<-gamm(f1, random=list(fStation =~ 1), data = ISIT2,
      method = "REML",control=lmc,
      weights = varComb(varIdent(form =~ 1 | fMonth),
                    varPower(form =~ Depth1000 | fMonth)))

AIC(M17.5A$lme,M17.5B$lme,M17.5C$lme,M17.5D$lme,M17.5E$lme)


anova(M17.5E$gam)
plot(M17.5E$gam)


E<-resid(M17.5E$lme,type="normalized")
F<-fitted(M17.5E$gam)
P<-predict(M17.5E$gam,se=T)
par(mfrow=c(2,1),mar=c(5,4,1,1))
plot(M17.5E$gam)
plot(F,E,xlab="Fitted values",ylab="Residuals")



Vario17.5E<-Variogram(M17.5E$lme,form=~Depth1000|Cruise,robust=T,resType="n")
plot(Vario17.5E)



#Figure 17.9
library(lattice)
MyLines2<-function(xi,yi,subscripts,...){
     Pi<-P$fit[subscripts]
     SEi<-P$se.fit[subscripts]
     I<-order(xi)
     panel.lines(xi[I], Pi[I],col=1)
     panel.lines(xi[I], Pi[I]+2*SEi[I],col=1,lty=2)
     panel.lines(xi[I], Pi[I]-2*SEi[I],col=1,lty=2)
     }


xyplot(Sources~Depth1000|fMonth, data =ISIT2,
  type="n",groups=Station,
  strip = function(bg='white', ...) strip.default(bg='white', ...),
  scales = list(alternating = T, x = list(relation = "same"),y = list(relation = "free")),
  xlab="Depth (km)",ylab=expression(paste(Sources," m"^{-3},"")),
  panel=panel.superpose,
  panel.groups = MyLines2)



#Figure 17.10:



xyplot(E~Depth1000|fStation, data = ISIT2,
type="n",col=1,xlab="Sample Depth",ylab="Normalized residuals",
strip = function(bg='white', ...) strip.default(bg='white', ...),
scales = list(alternating = T, x = list(relation = "free"),y = list(relation = "same")),
panel = function(x, y,subscripts) {
                panel.grid(h=-1, v= 2)
                I<-order(x)
                panel.points(x[I], y[I],col=1,cex=0.5)
                tmp<-gam(y~s(x))
                tmp2<-predict(tmp)
                print(unique(ISIT2$Station[subscripts]))
                print(summary(tmp))
                panel.lines(x[I],tmp2[I],col=1,lty=1)})

##################################################

#Based on geographical distances

G1<-ISIT2$Station==1 | ISIT2$Station==2  | ISIT2$Station==3
G2<-ISIT2$Station==6 | ISIT2$Station==9
G3<-ISIT2$Station==7 | ISIT2$Station==8  | ISIT2$Station==11
G4<-ISIT2$Station==12| ISIT2$Station==13 | ISIT2$Station==14 |
    ISIT2$Station==15| ISIT2$Station==19
G5<-ISIT2$Station==16| ISIT2$Station==17
G6<-ISIT2$Station==18


#For older R versions, use:
f1<-formula(Sources~
         s(Depth1000,by=as.numeric(G1))+
         s(Depth1000,by=as.numeric(G2))+
         s(Depth1000,by=as.numeric(G3))+
         s(Depth1000,by=as.numeric(G4))+
         s(Depth1000,by=as.numeric(G5))+
         s(Depth1000,by=as.numeric(G5))+
         fMonth)
         
         
         
         
         
#For more recent versions of R, use:
G.16 <- vector(length = nrow(ISIT2))

G.16[ISIT2$Station==1 | ISIT2$Station==2  | ISIT2$Station==3] <- 1
G.16[ISIT2$Station==6 | ISIT2$Station==9] <- 2
G.16[ISIT2$Station==7 | ISIT2$Station==8  | ISIT2$Station==11] <- 3
G.16[ISIT2$Station==12| ISIT2$Station==13 | ISIT2$Station==14 |
    ISIT2$Station==15| ISIT2$Station==19] <- 4
G.16[ISIT2$Station==16| ISIT2$Station==17] <- 5
G.16[ISIT2$Station==18] <- 6
fG.16<-factor(G.16)

f1<-formula(Sources~ s(Depth1000,by=fG.16) + fMonth)



M.GeoA <- gamm(f1,random = list(fStation =~ 1),
      method = "REML", control = lmc, data = ISIT2)

M.GeoB <- gamm(f1, random = list(fStation =~ 1),
      method = "REML", control = lmc, data = ISIT2,
      weights = varIdent(form =~ 1 | fMonth))

M.GeoC <- gamm(f1, random=list(fStation=~1), data = ISIT2,
     method = "REML", control = lmc,
      weights = varPower(form =~ Depth1000))

M.GeoD<-gamm(f1, random = list(fStation =~ 1), data = ISIT2,
      method = "REML", control = lmc,
      weights=varComb(varIdent(form =~ 1 | fMonth),
                    varPower(form =~ Depth1000)))

M.GeoE<-gamm(f1, random=list(fStation =~ 1), data = ISIT2,
      method = "REML",control=lmc,
      weights = varComb(varIdent(form =~ 1 | fMonth),
                    varPower(form =~ Depth1000 | fMonth)))

AIC(M.GeoA$lme,M.GeoB$lme,M.GeoC$lme,M.GeoD$lme,M.GeoE$lme)






#Make a vector G that defines the grouping structure.

#For older R versions, use:
G1<-ISIT2$Station==1 | ISIT2$Station==2 | ISIT2$Station==3 | ISIT2$Station==12
G2<-ISIT2$Station==6 | ISIT2$Station==7 | ISIT2$Station==8 | ISIT2$Station==9
G3<-ISIT2$Station==11 | ISIT2$Station==19
G4<-ISIT2$Station==13 | ISIT2$Station==14 | ISIT2$Station==15
G5<-ISIT2$Station==16 | ISIT2$Station==17 | ISIT2$Station==18
f1<-formula(Sources~
       s(Depth1000,by=as.numeric(G1))+
       s(Depth1000,by=as.numeric(G2))+
       s(Depth1000,by=as.numeric(G3))+
       s(Depth1000,by=as.numeric(G4))+
       s(Depth1000,by=as.numeric(G5))+
       factor(Month))

M.cor4A<-gamm(f1, random=list(Station=~1),
      method="REML",control=lmc,data=ISIT2)








#For newer R versions, use:
G.15 <- vector(length = nrow(ISIT2))
G.15[ISIT2$Station==1 | ISIT2$Station==2 | ISIT2$Station==3 | ISIT2$Station==12] <- 1
G.15[ISIT2$Station==6 | ISIT2$Station==7 | ISIT2$Station==8 | ISIT2$Station==9]  <- 2
G.15[ISIT2$Station==11 | ISIT2$Station==19] <- 3
G.15[ISIT2$Station==13 | ISIT2$Station==14 | ISIT2$Station==15] <- 4
G.15[ISIT2$Station==16 | ISIT2$Station==17 | ISIT2$Station==18] <- 5

f1<-formula(Sources ~ s(Depth1000,by=G.15) + factor(Month))

M.cor4A<-gamm(f1, random=list(Station=~1),
      method="REML",control=lmc,data=ISIT2)
