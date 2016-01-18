#    Mixed Effects Models and Extensions with R (2009)
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



#The R code for Chapter 18 is rather complicated. We copied and pasted
#code from various tinn-R files into this file. It needs some polishing.
#Give a shout if you have error messages (highstat@highstat.com)
#Latest update: 25 August 2009. Alain Zuur



library(AED); data(RIKZDATAEnv)

RIKZ2 <- RIKZDATAEnv
RIKZ1 <- RIKZ2[RIKZ2$Year>1990,]
I <- !is.na(RIKZ1$DIN)
RIKZ <- RIKZ1[I,]

RIKZ$LDIN <- log(RIKZ$DIN)
RIKZ$fStation <- factor(RIKZ$Station)
library(lattice)
RIKZ$MyTime <- RIKZ$Year + RIKZ$dDay3 / 365



#Figure 18.3
xyplot(LDIN ~ MyTime | Station, data = RIKZ,
       xlab = "Time", col = 1, type = "h",
       strip = function(bg = 'white', ...)
       strip.default(bg = 'white', ...))




#Figure 18.4
boxplot(LDIN~fStation, data = RIKZ, xaxt = "n")
text(1:31, par("usr")[3] - 0.25, srt = 45, adj = 1,
   labels = levels(RIKZ$fStation), xpd = TRUE, cex = 0.75)




#Figure 18.5
RIKZ$fMonth <- factor(RIKZ$Month)
bwplot(LDIN ~ fMonth | Area, data = RIKZ,
   xlab="Month",
   strip = function(bg = 'white', ...)
   strip.default(bg = 'white', ...), col = 1,
 scales = list(rot = 45, cex = .6))








RIKZ$X <- RIKZ$X31UE_ED50  #spatial coordinates
RIKZ$Y <- RIKZ$X31UN_ED50  #spatial coordinates
library(mgcv)
M1 <- gamm(LDIN ~ s(Year) + s(dDay3) + s(X, Y),
           random = list(fStation =~ 1), data = RIKZ )


#Figure 18.6
op <- par(mfrow = c(2, 2))
plot(M1$gam, select = c(1))
plot(M1$gam, select = c(2))
plot(M1$gam, select = c(3))
E <- resid(M1$lme, type = "normalized")
F <- fitted(M1$lme)
plot(x = F, y = E, xlab = "Fitted values", ylab = "Residuals",
     cex = 0.3)
par(op)




#For older R versions (<2.7.0), use:
ID<-c("WZ","GM","VD","ED","OS","WS","KZ","NZ","NC","VM")
M2<-gamm(LDIN~
      s(Year,by=as.numeric(Area=="WZ"),bs="cr")+
      s(Year,by=as.numeric(Area=="GM"),bs="cr")+
      s(Year,by=as.numeric(Area=="VD"),bs="cr")+
      s(Year,by=as.numeric(Area=="ED"),bs="cr")+
      s(Year,by=as.numeric(Area=="OS"),bs="cr")+
      s(Year,by=as.numeric(Area=="WS"),bs="cr")+
      s(Year,by=as.numeric(Area=="KZ"),bs="cr")+
      s(Year,by=as.numeric(Area=="NZ"),bs="cr")+
      s(Year,by=as.numeric(Area=="NC"),bs="cr")+
      s(Year,by=as.numeric(Area=="VM"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="WZ"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="GM"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="VD"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="ED"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="OS"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="WS"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="KZ"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="NZ"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="NC"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="VM"),bs="cr")+
      s(X, Y),random=list(fStation=~1), data = RIKZ)


#For more recent R versions, use:
#M2<-gamm(LDIN~ s(Year,by=Area,bs="cr")+s(X, Y),random=list(fStation=~1), data = RIKZ)
ID<-levels(RIKZ$Area)
M2<-gamm(LDIN~ Area + s(Year,  by = Area, bs="cr") +
               s(dDay3, by = Area, bs="cr") +
               s(X, Y),
               random=list(fStation=~1), data = RIKZ)
#Perhaps remove Area as fStation is used as a random effect?





#Figure 18.7
#Change the path name of this file to where ever you stored it.
source("D:\\applicat\\HighlandStatistics\\Book2\\AllRCode\\supportroutines4.R")


#Run this first..before the rest of the code is pasted
out8D<-mygamplot2(M2$gam)


OUTTrend<-out8D[out8D[,5]<=10,]

#This out8D contains the actual values that is used by plot.gam to
#make the graphs that you normally get with the plot.gam function
#The first 10 blocks contain the fitted long-term trends


#time,fit,ul,ll,id

#But for new R version, use:

IDFull<-rep(ID,each=100)

library(lattice)
time2<-OUTTrend[,1]
fit2<-OUTTrend[,2]
ul2<-OUTTrend[,3]
ll2<-OUTTrend[,4]
id2<-OUTTrend[,5]

timespecial<-1991+time2/365

xyplot(fit2~time2|IDFull,type="l",col=1,xlab="Time (years)",
       ylab="Trends",
  strip = function(bg='white', ...) strip.default(bg='white', ...),
  panel = function(x, y,subscripts) {
                panel.grid(h=-1, v= 2)
                I<-order(x)
                llines(x[I], y[I],col=1)
                zup<-ul2[subscripts]
                zlow<-ll2[subscripts]
                llines(x[I], zup[I],col=1,lty=2)
                llines(x[I], zlow[I],col=1,lty=2)},
  scales = list(alternating = T,
                x = list(relation = "same"),
                y = list(relation = "same")))



#Figure 18.8
OUTSS<-out8D[out8D[,5]>=11,]
library(lattice)
time2<-OUTSS[,1]
fit2<-OUTSS[,2]
ul2<-OUTSS[,3]
ll2<-OUTSS[,4]
id2<-OUTSS[,5]

xyplot(fit2~time2|IDFull,type="l",col=1,xlab="Time (days)",
       ylab="Seasonal patterns",
  strip = function(bg='white', ...) strip.default(bg='white', ...),
  panel = function(x, y,subscripts) {
                panel.grid(h=-1, v= 2)
                I<-order(x)
                llines(x[I], y[I],col=1)
                zup<-ul2[subscripts]
                zlow<-ll2[subscripts]
                llines(x[I], zup[I],col=1,lty=2)
                llines(x[I], zlow[I],col=1,lty=2)},
  scales = list(alternating = T,
                x = list(relation = "same"),
                y = list(relation = "same")))



#Figure 18.9

E2<-resid(M2$lme,type="n")
plot(E2 ~ RIKZ$fMonth, xlab = "Month", ylab = "Normalised residuals")





n<-length(RIKZ$Month)
RIKZ$M14<-vector(length=n)
RIKZ$M14[1:n]<-0
RIKZ$M14[RIKZ$Month >= 1 & RIKZ$Month  <= 3]<-1
RIKZ$M14[RIKZ$Month >= 4 & RIKZ$Month  <= 6]<-2
RIKZ$M14[RIKZ$Month >= 7 & RIKZ$Month  <= 9]<-3
RIKZ$M14[RIKZ$Month >= 10 & RIKZ$Month <= 12]<-4
RIKZ$fM14 <- factor(RIKZ$M14)

ffM14<-RIKZ$fM14


#Final model....for older R versions
M3<-gamm(LDIN~
      s(Year,by=as.numeric(Area=="WZ"),bs="cr")+
      s(Year,by=as.numeric(Area=="GM"),bs="cr")+
      s(Year,by=as.numeric(Area=="VD"),bs="cr")+
      s(Year,by=as.numeric(Area=="ED"),bs="cr")+
      s(Year,by=as.numeric(Area=="OS"),bs="cr")+
      s(Year,by=as.numeric(Area=="WS"),bs="cr")+
      s(Year,by=as.numeric(Area=="KZ"),bs="cr")+
      s(Year,by=as.numeric(Area=="NZ"),bs="cr")+
      s(Year,by=as.numeric(Area=="NC"),bs="cr")+
      s(Year,by=as.numeric(Area=="VM"),bs="cr")+
      fMonth * Area + s(X,Y),
      random=list(fStation=~1), data = RIKZ,
      weights=varIdent(form=~1|ffM14))


#For newer R versions
M3<-gamm(LDIN~ s(Year,by=Area,bs="cr") + fMonth * Area + s(X,Y),
      random=list(fStation=~1), data = RIKZ,
      weights=varIdent(form=~1|ffM14))
#Make some coffee, or go for lunch!

E3<-resid(M3$lme,type="n")


#Figure 18.10
plot(E3~RIKZ$fMonth,xlab="Month",ylab="Normalized residuals")




#Figure 18.11
MyTime <- RIKZ$Year + RIKZ$dDay3/365
library(lattice)
xyplot(E3~MyTime|Station,data=RIKZ,
        type="h",col=1,ylab="Residuals",xlab="Time (year)")



#Figure 18.12
ESpace<-tapply(E3,RIKZ$Station,FUN=mean)
X1<-tapply(RIKZ$X,RIKZ$Station,FUN=mean)
Y1<-tapply(RIKZ$Y,RIKZ$Station,FUN=mean)
LDIN1<-tapply(RIKZ$LDIN,RIKZ$Station,FUN=mean)

library(gstat)
mydata<-data.frame(ESpace,LDIN1,X1,Y1)

coordinates(mydata)<-c("X1","Y1")
bubble(mydata,"ESpace",col=c("black","grey"),main="Normalized residuals",
  xlab="X-coordinates",ylab="Y-coordinates")


#Save output to a textfile
rownames(out8D) <- seq(1:nrow(out8D))
MyFile<-"D:\\applicat\\HighlandStatistics\\Book2\\Data\\RIKZ\\DINRESULTS.txt"
write.table(out8D,file=MyFile)









####################################################
#Section 18.4 Temperature
####################################################
library(AED); data(RIKZDATAEnv)
Data2 <- RIKZDATAEnv
Data1<-Data2[Data2$Year>1990,]
I<-!is.na(Data1$T)
Data<-Data1[I,]
names(Data)

# [1] "Sample"     "Date"       "DateNr"     "dDay1"      "dDay2"
# [6] "dDay3"      "Station"    "Area"       "X31UE_ED50" "X31UN_ED50"
#[11] "Year"       "Month"      "Season"     "DIN"        "DIP"
#[16] "SIL"        "TN"         "TP"         "SAL"        "T"
#[21] "ZICHT"      "ZS"         "CHLFa"

attach(Data)
TEMP<-T


#Fi
library(lattice)
MyTime<-Year+dDay3/365
#xyplot(TEMP~MyTime|Station,xlab="Time",
#   strip = function(bg='white', ...) strip.default(bg='white', ...),
#   col=1,type="h")
#


#Figure 18.13
bwplot(TEMP~factor(Month)|Area,xlab="Month",
   strip = function(bg='white', ...) strip.default(bg='white', ...),
   col=1)





X<-X31UE_ED50
Y<-X31UN_ED50

library(mgcv)


#For older (<2.7.0) use:
ID<-c("WZ","GM","VD","ED","OS","WS","KZ","NZ","NC","VM")
tmp2<-gamm(TEMP~
      s(Year,by=as.numeric(Area=="WZ"),bs="cr")+
      s(Year,by=as.numeric(Area=="GM"),bs="cr")+
      s(Year,by=as.numeric(Area=="VD"),bs="cr")+
      s(Year,by=as.numeric(Area=="ED"),bs="cr")+
      s(Year,by=as.numeric(Area=="OS"),bs="cr")+
      s(Year,by=as.numeric(Area=="WS"),bs="cr")+
      s(Year,by=as.numeric(Area=="KZ"),bs="cr")+
      s(Year,by=as.numeric(Area=="NZ"),bs="cr")+
      s(Year,by=as.numeric(Area=="NC"),bs="cr")+
      s(Year,by=as.numeric(Area=="VM"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="WZ"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="GM"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="VD"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="ED"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="OS"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="WS"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="KZ"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="NZ"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="NC"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="VM"),bs="cr")+s(X,Y),
      random=list(Station=~1))



#For newer R versions, use:
ID<-levels(Area)
tmp2<-gamm(TEMP~
      s(Year, by = Area, bs = "cr") +
      s(dDay3, by = Area, bs = "cr") + s(X,Y) + Area,
      random=list(Station=~1))





#Figures 18.14 and 18.15
library(AED)
#source("D:\\applicat\\HighlandStatistics\\consultancy\\2007\\RIKZ2007\\supportroutines4.R")
#The code in this file has been added to the AED package.


out8D<-mygamplot2(tmp2$gam)



IDFull<-rep(ID,each=100)


OUTTrend<-out8D[out8D[,5]<=10,]
library(lattice)
time2<-OUTTrend[,1]
fit2<-OUTTrend[,2]
ul2<-OUTTrend[,3]
ll2<-OUTTrend[,4]
id2<-OUTTrend[,5]

timespecial<-1991+time2/365

xyplot(fit2~time2|IDFull,type="l",col=1,xlab="Time (years)",
       ylab="Trends",
  strip = function(bg='white', ...) strip.default(bg='white', ...),
  panel = function(x, y,subscripts) {
                panel.grid(h=-1, v= 2)
                I<-order(x)
                llines(x[I], y[I],col=1)
                zup<-ul2[subscripts]
                zlow<-ll2[subscripts]
                llines(x[I], zup[I],col=1,lty=2)
                llines(x[I], zlow[I],col=1,lty=2)},
  scales = list(alternating = T,
                x = list(relation = "same"),
                y = list(relation = "same")))



#Seasonal part
OUTTrend<-out8D[out8D[,5]>=11,]
library(lattice)
time2<-OUTTrend[,1]
fit2<-OUTTrend[,2]
ul2<-OUTTrend[,3]
ll2<-OUTTrend[,4]
id2<-OUTTrend[,5]

xyplot(fit2~time2|IDFull,type="l",col=1,xlab="Time (days)",
       ylab="Seasonal patterns",
  strip = function(bg='white', ...) strip.default(bg='white', ...),
  panel = function(x, y,subscripts) {
                panel.grid(h=-1, v= 2)
                I<-order(x)
                llines(x[I], y[I],col=1)
                zup<-ul2[subscripts]
                zlow<-ll2[subscripts]
                llines(x[I], zup[I],col=1,lty=2)
                llines(x[I], zlow[I],col=1,lty=2)},
  scales = list(alternating = T,
                x = list(relation = "same"),
                y = list(relation = "same")))




#Figure 18.16
plot(tmp2$gam,select=21)




#Save the results for later
rownames(out8D) <- seq(1:nrow(out8D))
MyFile<-"D:\\applicat\\HighlandStatistics\\Book2\\Data\\RIKZ\\TEMPRESULTS.txt"
write.table(out8D,file=MyFile)





##############################################################
# DIAT1
##############################################################


Data1<-read.table(file="D:\\applicat\\HighlandStatistics\\consultancy\\2007\\RIKZ2007\\RIKZdata3aAggregatedJune17_2007.txt",header=T)

I<-!is.na(Data1$DIAT1)
Data<-Data1[I,]
names(Data)




# [1] "Sample"     "Date"       "DateNr"     "dDay1"      "dDay2"
# [6] "dDay3"      "Station"    "Area"       "X31UE_ED50" "X31UN_ED50"
#[11] "Lab"        "Year"       "Month"      "Season"     "AMPHIDOM"
#[16] "CERATIAC"   "DIAT1"      "DIAT2"      "DIAT3"      "DINOPHYC"
#[21] "FLAG1"      "FLAG2"      "GLENODIN"   "GONYAULA"   "GYMNODIN"
#[26] "HETEROCA"   "PERIDINC"   "PHAEOCYST"  "PROROCEN"   "PYROCYST"
#[31] "PYROPHAC"


attach(Data)


LDIAT1<-log(1+DIAT1)


#Start with DIAT 1

library(lattice)

#Figure 18.17
bwplot(LDIAT1~factor(Lab)|factor(Area),xlab="Month",
   strip = function(bg='white', ...) strip.default(bg='white', ...),
   col=1)





X<-X31UE_ED50
Y<-X31UN_ED50
library(mgcv)



#For older R versions, use:
ID<-c("WZ","GM","VD","ED","OS","WS","KZ","NZ","NC","VM")
tmp2<-gamm(LDIAT1~
      s(Year,by=as.numeric(Area=="WZ"),bs="cr")+
      s(Year,by=as.numeric(Area=="GM"),bs="cr")+
      s(Year,by=as.numeric(Area=="VD"),bs="cr")+
      s(Year,by=as.numeric(Area=="ED"),bs="cr")+
      s(Year,by=as.numeric(Area=="OS"),bs="cr")+
      s(Year,by=as.numeric(Area=="WS"),bs="cr")+
      s(Year,by=as.numeric(Area=="KZ"),bs="cr")+
      s(Year,by=as.numeric(Area=="NZ"),bs="cr")+
      s(Year,by=as.numeric(Area=="NC"),bs="cr")+
      s(Year,by=as.numeric(Area=="VM"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="WZ"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="GM"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="VD"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="ED"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="OS"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="WS"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="KZ"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="NZ"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="NC"),bs="cr")+
      s(dDay3,by=as.numeric(Area=="VM"),bs="cr")+s(X,Y),
      random=list(Station=~1))




#For newer R versions, use:
ID<-levels(Area)
tmp2<-gamm(LDIAT1~
      s(Year,  by=Area, bs="cr")+
      s(dDay3, by=Area, bs="cr")+s(X,Y) + Area,
      random=list(Station=~1))





###############
#plot smoothers
#Run supportroutines3.r


library(mgcv)
library(AED)


out8D<-mygamplot2(tmp2$gam)
IDFull<-rep(ID,each=100)


OUTTrend<-out8D[out8D[,5]<=10,]
library(lattice)
time2<-OUTTrend[,1]
fit2<-OUTTrend[,2]
ul2<-OUTTrend[,3]
ll2<-OUTTrend[,4]
id2<-OUTTrend[,5]

timespecial<-1991+time2/365


#Figure 18.18
xyplot(fit2~time2|IDFull,type="l",col=1,xlab="Time (years)",
       ylab="Trends",
  strip = function(bg='white', ...) strip.default(bg='white', ...),
  panel = function(x, y,subscripts) {
                panel.grid(h=-1, v= 2)
                I<-order(x)
                llines(x[I], y[I],col=1)
                zup<-ul2[subscripts]
                zlow<-ll2[subscripts]
                llines(x[I], zup[I],col=1,lty=2)
                llines(x[I], zlow[I],col=1,lty=2)},
  scales = list(alternating = T,
                x = list(relation = "same"),
                y = list(relation = "same")))





#Seasonal, if fitted
OUTTrend<-out8D[out8D[,5]>=11,]
library(lattice)
time2<-OUTTrend[,1]
fit2<-OUTTrend[,2]
ul2<-OUTTrend[,3]
ll2<-OUTTrend[,4]
id2<-OUTTrend[,5]

#Figure 18.19
xyplot(fit2~time2|IDFull,type="l",col=1,xlab="Time (days in the year)",
       ylab="Seasonal patterns",
  strip = function(bg='white', ...) strip.default(bg='white', ...),
  panel = function(x, y,subscripts) {
                panel.grid(h=-1, v= 2)
                I<-order(x)
                llines(x[I], y[I],col=1)
                zup<-ul2[subscripts]
                zlow<-ll2[subscripts]
                llines(x[I], zup[I],col=1,lty=2)
                llines(x[I], zlow[I],col=1,lty=2)},
  scales = list(alternating = T,
                x = list(relation = "same"),
                y = list(relation = "same")))




#Figure 18.20
plot(tmp2$gam)


rownames(out8D) <- seq(1:nrow(out8D))
MyFile<-"D:\\applicat\\HighlandStatistics\\Book2\\Data\\RIKZ\\DIATRESULTS.txt"
write.table(out8D,file=MyFile)


#####################################################################
# Compare all results
#####################################################################




#File to read all trends

MyFile2<-"D:\\applicat\\HighlandStatistics\\Book2\\Data\\RIKZ\\DINRESULTS.txt"
DIN<-read.table(file=MyFile2,header=T)
colnames(DIN)<-c("time","fit","ul","ll","id")



MyFile8<-"D:\\applicat\\HighlandStatistics\\Book2\\Data\\RIKZ\\TEMPRESULTS.txt"
TEMP<-read.table(file=MyFile8,header=T)
colnames(TEMP)<-c("time","fit","ul","ll","id")

################################

MyFile9<-"D:\\applicat\\HighlandStatistics\\Book2\\Data\\RIKZ\\DIATRESULTS.txt"
DIAT1<-read.table(file=MyFile9,header=T)
colnames(DIAT1)<-c("time","fit","ul","ll","id")




#Old R versions:
ID<-c("WZ","GM","VD","ED","OS","WS","KZ","NZ","NC","VM")


#New R versions:
library(AED); data(RIKZDATAEnv)
RIKZ2 <- RIKZDATAEnv
RIKZ1 <- RIKZ2[RIKZ2$Year>1990,]
I <- !is.na(RIKZ1$DIN)
RIKZ <- RIKZ1[I,]
ID<-levels(RIKZ$Area)



ToMatrix<-function(x,MyName,ID){
  Mat<-matrix(nrow=100,ncol=10)
  for (i in 1:10){
     Mat[,i]<-x[x$id==i,2]
  }
  colnames(Mat)<-paste(MyName,ID,sep="")
  Mat
}


DIN10<-ToMatrix(DIN,"DIN",ID)
TEMP10<-ToMatrix(TEMP,"TEM",ID)

DIAT110<-ToMatrix(DIAT1,"DIAT1",ID)




MyTime<-DIN[DIN$id==1,1]

Z<-cbind(DIN10,TEMP10)
Y<-cbind(DIAT110)

#rownames(Z)<-MyTime


dim(Y)
dim(Z)

Res<-matrix(nrow=10,ncol=2)
for (i in 1:10){
  for (j in 1:2){
     j1<-i+10*(j-1)
     Res[i,j]<-cor(Y[,i],Z[,j1])

}}

rownames(Res)<-colnames(Y)
colnames(Res)<-c("DIN","TEMP")

#Note the difference in the order of the stations.
#Other differences are due to small changes in mgcv and coding
#of the smoothers.
options(digits=2)
Res



#These are not significant!
ToPlot<-matrix(1,nrow=2,ncol=10*2)
#for (i in 1:6){
#  ToPlot[2,8+10*(i-1)]<-0
#  ToPlot[2,4+10*(i-1)]<-0
#  ToPlot[2,9+10*(i-1)]<-0
#}




plot(0,0,xlim=c(1,10),ylim=c(-1,1),type="n",axes=F,xlab="",ylab="Pearson correlation")
axis(2)
axis(1,at=seq(1:10),labels=colnames(Y),cex.axis=1,srt = 45)
for (i in 1:10){
  a1<-c(i,i)
  a2<-c(i+0.05,i+0.05)
  b1<-c(Res[i,1],0)
  b2<-c(Res[i,2],0)
  if (ToPlot[1,i]==1) lines(a1,b1,col=1,lty=1,lwd=1)
  if (ToPlot[2,i]==1) lines(a2,b2,col=1,lty=2,lwd=2)
  }
  abline(0,0)


leg.txt <- c("DIN","TEMP")
y.col <- c(1,1)


legend(x="topleft", leg.txt,col=y.col,lty=c(1,2),lwd=c(1,2),cex=0.7)













ID<-c("WZ","GM","VD","ED","OS","WS","KZ","NZ","NC","VM")


s1<-c(644618,	5899537)
s2<-c(568549,	5732245)
s3<-c(544866,	5734293)
s4<-c(762713,	5918656)
s5<-c(566123,	5713095)
s6<-c(546986,	5694684)
s7<-c(611852,	5873694)
s8<-c(566854,	5891403)
s9<-c(564500,	6058290)
s10<-c(553452,	5711596)


Coord<-rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)
Coord10<-rbind(Coord)/10000


Dat<-cbind(Res,Coord/10000)

library(lattice)



#
k<-1

a1<-1+10*(k-1)
a10<-10+10*(k-1)
X<-rep(Dat[a1:a10,3],2)
Y<-rep(Dat[a1:a10,4],2)
Block<-factor(rep(c("DIN","TEMP"),each=10))
Cor<-c(Res[a1:a10,1],Res[a1:a10,2])
Loc<-rep(ID,2)

MyVar<-c("DIAT1")

a1<-1 +10*(k-1)
a2<-10+10*(k-1)
ToPlotVec<-as.vector(t(ToPlot[1:2,a1:a2]))

library(lattice)
xyplot(Y~X|Block,type="n",col=1,
   xlab="X coordinate",ylab="Y coordinate", main=MyVar[k],
   strip = function(bg='white', ...) strip.default(bg='white', ...),
   panel = function(x, y,subscripts,...) {
                panel.grid(h=-1, v= 2)
                x1<-X[subscripts]
                y1<-Y[subscripts]
                tp1<-ToPlotVec[subscripts]
                z1<-abs(Cor[subscripts])
                Loc1<-Loc[subscripts]
                panel.text(x1,y1,Loc1,col=1*tp1,cex=1*z1)})




#gplot stuff


library(gplots)
#DAT=read.table(file("clipboard"), header=TRUE, dec=".")
DAT<-Res
Y=as.matrix(round(DAT,2))

heatmap.2(Y, main="Correlation matrix",
  cellnote=Y, notecol="black",                      #plot correlations
  col=colorpanel(4, "purple", "grey", "red"),
  trace="none",                                     #no trace lines
  density.info="histogram", keysize=1, denscol="white",     #legend
  cexCol=1, cexRow=1, margins=c(10,10),      #label size & space
          )



heatmap.2(Y, main="Correlation matrix",
          cellnote=Y, notecol="black",
          col=colorpanel(10, "blue", "green", "red"),trace="none",
          density.info="histogram", keysize=1, denscol="white",
          cexCol=1, cexRow=1, margins=c(10,10),
          Rowv=FALSE, Colv=FALSE
          )

