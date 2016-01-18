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

#Figure 3.1
op <- par(mfrow=c(2,2),mar=c(5,4,1,2))
Sources16<-ISIT$Sources[ISIT$Station==16]
Depth16<-ISIT$SampleDepth[ISIT$Station==16]
plot(Depth16,Sources16,type="p")

library(gam)
M2<-gam(Sources16~lo(Depth16,span=0.5))
plot(M2,se=T)					#Figure 3.1B


#It may be better to predict along an equidistant gradient
P2 <- predict(M2, se = TRUE)
plot(Depth16, Sources16, type = "p")
I1 <- order(Depth16)
lines(Depth16[I1], P2$fit[I1], lty = 1)
lines(Depth16[I1], P2$fit[I1] + 2 * P2$se[I1], lty = 2)
lines(Depth16[I1], P2$fit[I1] - 2 * P2$se[I1], lty = 2)

par(op)





#Figure 3.2
Data<-read.table("C://Bookdata/ISIT.txt",header=T)
library(lattice)

Sources<-Data$Sources[Data$Station==16]
Depth<-Data$SampleDepth[Data$Station==16]

plot(Depth,Sources,type="n")
TargetVal<-1500
Bin<-500
B1<-TargetVal-Bin
B2<-TargetVal+Bin

abline(v=B1,lty=2)
abline(v=B2,lty=2)
points(TargetVal,-1,pch=17,cex=2)
BinInterval<-vector(length=length(Depth))
BinInterval[1:length(Depth)]<-1
BinInterval[Depth>=B1 & Depth <=B2]<-16
points(Depth,Sources,pch=BinInterval)

S1<-Sources[BinInterval==16]
D1<-Depth[BinInterval==16]
tmp1<-lm(S1~D1)

DD1<-seq(B1,B2,length=100)
NewData<-data.frame(D1=DD1)
pred1<-predict(tmp1,newdata=NewData)

lines(DD1,pred1,lwd=2)

NewData<-data.frame(D1=1500)
pred1<-predict(tmp1,newdata=NewData)






#Figure 3.3
#Code to illustrate effect of span width
Data<-read.table("C://Bookdata/ISIT.txt",header=T)
library(lattice)

Sources<-Data$Sources[Data$Station==16]
Depth<-Data$SampleDepth[Data$Station==16]
I=order(Depth)
Sources1<-Sources[I]
Depth1=Depth[I]

library(gam)


M1=gam(Sources1~lo(Depth1,span=0.1))
M2=gam(Sources1~lo(Depth1,span=0.15))
M3=gam(Sources1~lo(Depth1,span=0.2))
M4=gam(Sources1~lo(Depth1,span=0.25))
M5=gam(Sources1~lo(Depth1,span=0.3))
M6=gam(Sources1~lo(Depth1,span=0.5))
M7=gam(Sources1~lo(Depth1,span=0.7))
M8=gam(Sources1~lo(Depth1,span=0.8))
M9=gam(Sources1~lo(Depth1,span=0.9))
M10=gam(Sources1~lo(Depth1,span=1.0))


Mp1=predict(M1,se=T)
Mp2=predict(M2,se=T)
Mp3=predict(M3,se=T)
Mp4=predict(M4,se=T)
Mp5=predict(M5,se=T)
Mp6=predict(M6,se=T)
Mp7=predict(M7,se=T)
Mp8=predict(M8,se=T)
Mp9=predict(M9,se=T)
Mp10=predict(M10,se=T)

xall=rep(Depth1,10)
yall=rep(Sources1,10)
id=rep(c("span = 0.1","span = 0.15","span = 0.2",
 "span = 0.25","span = 0.3","span = 0.5",
 "span = 0.7","span = 0.8","span = 0.9",
 "span = 1.0"),
  each=length(Depth1))

Pall=c(Mp1$fit,Mp2$fit,Mp3$fit,Mp4$fit,Mp5$fit,
       Mp6$fit,Mp7$fit,Mp8$fit,Mp9$fit,Mp10$fit)
SEall=c(Mp1$se.fit,Mp2$se.fit,Mp3$se.fit,Mp4$se.fit,Mp5$se.fit,
       Mp6$se.fit,Mp7$se.fit,Mp8$se.fit,Mp9$se.fit,Mp10$se.fit)
library(lattice)

xyplot(yall~xall|id,xlab="Depth",ylab="Sources",
panel=function(x,y,subscripts,...){
   panel.points(x,y,col=1,cex=0.5)
   panel.lines(xall[subscripts],Pall[subscripts],col=1,lwd=1)
    panel.lines(xall[subscripts],Pall[subscripts]+2*SEall[subscripts],col=1,lwd=1,lty=2)
    panel.lines(xall[subscripts],Pall[subscripts]-2*SEall[subscripts],col=1,lwd=1,lty=2)
   }
)

AIC(M1,M2,M3,M4,M5,M6,M7,M8,M9,M10)



#Figure 3.4
M2<-gam(Sources16~lo(Depth16,span=0.5))
E2 <- resid(M2)
F2 <- fitted(M2)
par(mfrow=c(2,2), mar=c(5,4,2,2))
plot(x = F2, y=E2, xlab="Fitted values", ylab="Residuals")
plot(x = Depth16, y=E2, xlab="Depth", ylab="Residuals")
hist(E2,main="", xlab="Residuals")



#Quit and restart R



#Figure 3.5
library(AED); data(ISIT)
library(mgcv)
op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
Sources16 <- ISIT$Sources[ISIT$Station == 16]
Depth16 <- ISIT$SampleDepth[ISIT$Station == 16]

plot(Depth16, Sources16, type = "p")
M3 <- gam(Sources16 ~ s(Depth16, fx = FALSE, k=-1,bs = "cr"))
plot(M3, se = TRUE)
M3pred <- predict(M3, se = TRUE, type = "response")
plot(Depth16, Sources16, type = "p")
I1 <- order(Depth16)
lines(Depth16[I1], M3pred$fit[I1], lty=1)
lines(Depth16[I1], M3pred$fit[I1]+2*M3pred$se[I1],lty=2)
lines(Depth16[I1], M3pred$fit[I1]-2*M3pred$se[I1],lty=2)




#Figure 3.6

E2 <- resid(M3)
F2 <- fitted(M3)
par(mfrow=c(2,2), mar=c(5,4,2,2))
plot(x = F2, y=E2, xlab="Fitted values", ylab="Residuals")
plot(x = Depth16, y=E2, xlab="Depth", ylab="Residuals")
hist(E2,main="", xlab="Residuals")






#Figue 3.7
#Code to illustrate additive modelling with mgcv
library(lattice)
x<-seq(0,1,length=25)

co<-matrix(nrow=6,ncol=4)
co[1:6,1:4]<-rnorm(mean=0,sd=5,24)
f1<-co[1,1]+co[1,2]*x+co[1,3]*x^2+co[1,4]*x^3
f2<-co[2,1]+co[2,2]*x+co[2,3]*x^2+co[2,4]*x^3
f3<-co[3,1]+co[3,2]*x+co[3,3]*x^2+co[3,4]*x^3
f4<-co[4,1]+co[4,2]*x+co[4,3]*x^2+co[4,4]*x^3
f5<-co[5,1]+co[5,2]*x+co[5,3]*x^2+co[5,4]*x^3
f6<-co[6,1]+co[6,2]*x+co[6,3]*x^2+co[6,4]*x^3
xall<-rep(x,6)
ID<-rep(c(1,2,3,4,5,6),each=25)
f<-c(f1,f2,f3,f4,f5,f6)
xyplot(f~xall|factor(ID),col=1,type="l",xlab="X values",ylab="Function f")







#Figure 3.8
#Quit and restart R
Data<-read.table("C://Bookdata/ISIT.txt",header=T)
library(lattice)
Sources<-Data$Sources[Data$Station==16]
Depth<-Data$SampleDepth[Data$Station==16]

Range<-max(Depth)-min(Depth)
Bins<-Range/4

B1<-min(Depth)+Bins
B2<-min(Depth)+2*Bins
B3<-min(Depth)+3*Bins

plot(Depth,Sources)
abline(v=B1,lty=2)
abline(v=B2,lty=2)
abline(v=B3,lty=2)

S1<-Sources[Depth<=B1]
D1<-Depth[Depth<B1]
tmp<-lm(S1~D1+I(D1^2) +I(D1^3))
F1<-fitted(tmp)
I1<-order(D1)
lines(D1[I1],F1[I1])

S1<-Sources[Depth > B1 & Depth<=B2]
D1<-Depth[Depth > B1 & Depth<=B2]
tmp<-lm(S1~D1+I(D1^2) +I(D1^3))
F1<-fitted(tmp)
I1<-order(D1)
lines(D1[I1],F1[I1])

S1<-Sources[Depth > B2 & Depth<=B3]
D1<-Depth[Depth > B2 & Depth<=B3]
tmp<-lm(S1~D1+I(D1^2) +I(D1^3))
F1<-fitted(tmp)
I1<-order(D1)
lines(D1[I1],F1[I1])

S1<-Sources[Depth > B3 ]
D1<-Depth[Depth > B3 ]
tmp<-lm(S1~D1+I(D1^2) +I(D1^3))
F1<-fitted(tmp)
I1<-order(D1)
lines(D1[I1],F1[I1])







#Figure 3.9
#Quit and restart R
Data<-read.table("C://Bookdata/ISIT.txt",header=T)
library(lattice)
Sources19<-Data$Sources[Data$Station==19]
Depth19<-Data$SampleDepth[Data$Station==19]

Depth01<-Depth19
Depth01<-Depth01-min(Depth01)
Depth01<-Depth01/max(Depth01)
I<-order(Depth01)
Depth01<-Depth01[I]
Sources19<-Sources19[I]




rk<-function(x,z){
((z-0.5)^2-1/12)*((x-0.5)^2-1/12)/4 -((abs(x-z)-0.5)^4 -0.5*(abs(x-z)-0.5)^2 +7/240)/24
}



spl.X<-function(x,xk){
 q<-length(xk)+2
 n<-length(x)
 X<-matrix(1,n,q)
 X[,2]<-x
 X[,3:q]<-outer(x,xk,FUN=rk)
 X
}


YALL<-vector(length=100*8)
XALL<-vector(length=100*8)
IDALL<-vector(length=100*8)
a1<-1
a2<-100

XkALL<-0
IDxk<-0
for (knots in 1: 8){
    xk<-1:knots/(knots+1)
    XkALL<-c(XkALL,xk)
    IDxk<-c(IDxk,rep(knots,length(xk)))

    X<-spl.X(Depth01,xk)
    tmp1<-lm(Sources19~X-1)
    xp<-1:100/100
    Xp<-spl.X(xp,xk)

    plot(Depth01,Sources19)
    lines(xp,Xp%*%coef(tmp1))
    YALL[a1:a2]<-Xp%*%coef(tmp1)
    XALL[a1:a2]<-xp
    IDALL[a1:a2]<-knots
    a1<-a1+100
    a2<-a2+100

    n<-length(xk)
    for (i in 1:n){
     abline(v=xk[i],lty=2)
     }
 }
XkALL<-XkALL[-1]
IDxk<-IDxk[-1]



library(lattice)

IDALL2<-IDALL+2
xyplot(YALL~XALL|factor(IDALL2),type="l",col=1,xlab="Depth (scaled between 0 and 1)",
        ylab="Sources",a1=1,
        panel=function(x,y,subscripts,...){
        panel.lines(x,y,col=1,lwd=2)
        panel.points(Depth01,Sources19,col=1,cex=0.5)#
        a<-IDALL2[subscripts]
        if (a[1] == 3){for (i in 1:1) {panel.abline(v=XkALL[i],lty=2);panel.abline(v=0,lty=2);panel.abline(v=1,lty=2)} }
        if (a[1] == 4){for (i in 1:2) {panel.abline(v=XkALL[i+1],lty=2);panel.abline(v=0,lty=2);panel.abline(v=1,lty=2)} }
        if (a[1] == 5){for (i in 1:3) {panel.abline(v=XkALL[i+1+2],lty=2);panel.abline(v=0,lty=2);panel.abline(v=1,lty=2)} }
        if (a[1] == 6){for (i in 1:4) {panel.abline(v=XkALL[i+1+2+3],lty=2);panel.abline(v=0,lty=2);panel.abline(v=1,lty=2)} }
        if (a[1] == 7){for (i in 1:5) {panel.abline(v=XkALL[i+1+2+3+4],lty=2);panel.abline(v=0,lty=2);panel.abline(v=1,lty=2)} }
        if (a[1] == 8){for (i in 1:6) {panel.abline(v=XkALL[i+1+2+3+4+5],lty=2);panel.abline(v=0,lty=2);panel.abline(v=1,lty=2)} }
        if (a[1] == 9){for (i in 1:7) {panel.abline(v=XkALL[i+1+2+3+4+5+6],lty=2);panel.abline(v=0,lty=2);panel.abline(v=1,lty=2)} }
        if (a[1] == 10){for (i in 1:8) {panel.abline(v=XkALL[i+1+2+3+4+5+6+7],lty=2);panel.abline(v=0,lty=2);panel.abline(v=1,lty=2)} }
        }
        )




#Figure 3.10
#Show penalized spline in action
spl.S<-function(xk){
  q<-length(xk)+2
  S<-matrix(0,q,q)
  S[3:q,3:q]<-outer(xk,xk,FUN=rk)
  S
}

mat.sqrt<-function(S){
  d<-eigen(S,symmetric=TRUE)
  d$values[d$values<0]<-0
  rS<-d$vectors%*%diag(d$values^0.5)%*%t(d$vectors)
  rS
}


prs.fit<-function(y,x,xk,lambda){
print(lambda)
  q<-length(xk)+2
  n<-length(x)
  Xa<-rbind(spl.X(x,xk),mat.sqrt(spl.S(xk))*sqrt(lambda))
  y[(n+1):(n+q)]<-0
  lm(y~Xa-1)
}




xk<-1:7/8
xp<-1:100/100


lambda<-1e-6


#Graph showing the effect of lambda
YALL<-vector(length=100*8)
XALL<-vector(length=100*8)
IDALL<-vector(length=100*8)
a1<-1
a2<-100
lambda<-10e-7
for (i in 1:8){
  a1
  mod.2<-prs.fit(Sources19,Depth01,xk,lambda)
  Xp<-spl.X(xp,xk)
  YALL[a1:a2]<-Xp%*%coef(mod.2)
  XALL[a1:a2]<-xp
  IDALL[a1:a2]<-rep(lambda,100)
  a1<-a1+100
  a2<-a2+100
  lambda<-lambda*10
 }



xyplot(YALL~XALL|factor(IDALL),type="l",col=1,xlab="Depth (scaled between 0 and 1)",
        ylab="Sources",a1=1,
        panel=function(x,y,subscripts,...){
        panel.lines(x,y,col=1,lwd=2)
        panel.points(Depth01,Sources19,col=1,cex=0.5)}
        )





#Figure 3.11
#Graph showing the effect of lambda
lambda<-10e-6
n<-length(Depth01)
V<-rep(0,60)
Vl<-rep(0,60)
for (i in 1:60){
  mod<-prs.fit(Sources19,Depth01,xk,lambda)
  trA<-sum(influence(mod)$hat[1:n])
  rss<-sum((Sources19-fitted(mod)[1:n])^2)
  V[i]<-n*rss/(n-trA)^2
  Vl[i]<-lambda
  lambda<-lambda*1.09
 }
 plot(x=Vl,y=V,type="l",xlab="lambda",ylab="GCV")





#Quit and restart R
#Figure 3.12
library(AED); data(ISIT)
S8 <- ISIT$Sources[ISIT$Station == 8]
D8 <- ISIT$SampleDepth[ISIT$Station == 8]
S13 <- ISIT$Sources[ISIT$Station == 13]
D13 <- ISIT$SampleDepth[ISIT$Station == 13]
So <- c(S8, S13); De <- c(D8, D13)
ID <- rep(c(8, 13), c(length(S8), length(S13)))
mi <- max(min(D8), min(D13))
ma <- min(max(D8), max(D13))
I1 <- De > mi & De < ma
op <- par(mfrow = c(1, 2))
plot(D8[I1], S8[I1], pch = 16, xlab = "Depth",
    ylab = "Sources", col = 1, main = "Station 8",
    xlim = c(500, 3000), ylim = c(0, 40))
plot(D13[I1], S13[I1], pch = 16, xlab = "Depth",
    ylab = "Sources", col = 1, main = "Station 13",
    xlim = c(500, 3000), ylim = c(0, 40))
par(op)



#Page 57
library(mgcv)
fID <-factor(ID)  #This construction is needed for the vis.gam
M4 <- gam(So ~ s(De) + fID, subset = I1)
summary(M4)
anova(M4)

#Figure 3.13
plot(M4)



#Figure 3.14
par(mar=c(2,2,2,2))
vis.gam(M4,theta=120,color="heat")



#Figure 3.15
gam.check(M4)


#page 60
M5<-gam(So ~ s(De)+
             s(De, by = as.numeric(ID == 13)),
             subset = I1)
anova(M5)



#Figure 3.16
plot(M5)

#Figure 3.17
gam.check(M5)


anova(M4,M5, test="F")
AIC(tmp)
AIC(tmp1)



M6 <- gam(So~s(De,by = as.numeric(ID== 8))+
             s(De,by = as.numeric(ID== 13))-1, subset=I1)


#Or:
M6 <- gam(So~s(De,by = ID) + factor(ID),subset=I1)




#Section 3.5
Data<-read.table("C:\\BookData\\VegetationIndex.txt",header=T)

attach(Data)
names(Data)

#"SAMPLEYR" "Time"     "Transect" "Richness" "ROCK"     "LITTER"
#"BARESOIL"
#[8] "FallPrec" "SprTmax"

Z<-Data[,c(-1,-2,-3)]

pairs(Z)

source("D:\\applicat\\HighlandStatistics\\consultancy\\2007\\RIKZ2007\\supportroutines.R")

#Figure 3.18
pairs(Z, upper.panel=panel.smooth2,
         lower.panel=panel.cor,cex.labels=1.5)


cor(Z,use="pairwise.complete.obs")




library(mgcv)
M7<-gam(Richness~s(ROCK,bs="cs")+s(LITTER,bs="cs")+
                 s(BARESOIL,bs="cs")+
                 s(FallPrec,bs="cs")+
                 s(SprTmax,bs="cs"))
anova(M7)
summary(M7)



#Figure 3.19
par(mfrow=c(2,3),mar=c(4,4,2,2))
plot(M7)






tmp1<-gam(Richness~s(ROCK,bs="cs")+
                 s(BARESOIL,bs="cs")+
                 s(FallPrec,bs="cs")+
                 s(SprTmax,bs="cs"))
anova(tmp1)



tmp2<-gam(Richness~s(ROCK,bs="cs")+
                 s(BARESOIL,bs="cs")+
                 s(SprTmax,bs="cs"))
anova(tmp2)


anova(tmp)
anova(tmp1)
anova(tmp2)


tmp3<-gam(Richness~s(ROCK,BARESOIL)+
                 s(SprTmax,bs="cs"))






#Figure 3.20
tmp1<-gam(Richness~s(ROCK,bs="cs")+
                 s(BARESOIL,bs="cs")+
                 s(SprTmax,bs="cs"))
anova(tmp1)


par(mfrow=c(2,2),mar=c(4,4,2,2))
plot(tmp1)








library(AED); data(ISIT)
library(mgcv)
op <- par(mfrow=c(2,2),mar=c(5,4,1,2))
Sources16<-ISIT$Sources[ISIT$Station==16]
Depth16<-ISIT$SampleDepth[ISIT$Station==16]
plot(Depth16,Sources16,type="p")
tmp<-gam(Sources16~s(Depth16,fx=F,k=-1,bs="cr"))
plot(tmp,se=T)
tmppred<-predict(tmp,se=T,type="response")
plot(Depth16,Sources16,type="p")
I<-order(Depth16)
lines(Depth16[I],tmppred$fit[I],lty=1)
lines(Depth16[I],tmppred$fit[I]+2*tmppred$se[I],lty=2)
lines(Depth16[I],tmppred$fit[I]-2*tmppred$se[I],lty=2)







library(AED); data(ISIT)
S8<-ISIT$Sources[ISIT$Station==8]
D8<-ISIT$SampleDepth[ISIT$Station==8]
S13<-ISIT$Sources[ISIT$Station==13]
D13<-ISIT$SampleDepth[ISIT$Station==13]
So<-c(S8,S13)
De<-c(D8,D13)
ID<-rep(c(8,13),c(length(S8),length(S13)))
mi<-max(min(D8),min(D13))
ma<-min(max(D8),max(D13))
I<-De>mi & De < ma
par(mfrow=c(1,2))
plot(D8[I],S8[I],pch=16,xlab="Depth",
     ylab="Sources",col=1,main="Station 8",
     xlim=c(500,3000),ylim=c(0,40))
plot(D13[I],S13[I],pch=16,xlab="Depth",
     ylab="Sources",col=1,main="Station 13",
     xlim=c(500,3000),ylim=c(0,40))

library(mgcv)
M4 <-gam(So~s(De)+factor(ID),subset=I)
summary(M4)
anova(M4)



library(AED); data(Vegetation)
library(mgcv)
tmp<-gam(Richness~s(ROCK,bs="cs")+s(LITTER,bs="cs")+
           s(BARESOIL,bs="cs")+s(FallPrec,bs="cs")+
         s(SprTmax,bs="cs"),data=Vegetation)
