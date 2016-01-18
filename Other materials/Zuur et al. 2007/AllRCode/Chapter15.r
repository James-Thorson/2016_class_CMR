#   Mixed Effects Models and Extensions in Ecology with R (2009)
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




library(AED); data(Ythan)
library(lattice)
Birds <- as.vector(as.matrix(Ythan[,2:8]))
X <- as.vector(as.matrix(Ythan[,9:15]))
YX14 <- c(Birds,X)
Year14 <- rep(Ythan$Year,14)
N <- length(Ythan$Year)
ID14<-factor(rep(names(Ythan[,2:15]),each = N),
     levels=c("wheat",
    "barley","oats","cattle","sheep","pigs",
    "nitrate","Oystercatcher","Turnstone","Curlew",
    "BartailedGodwit","Redshank","Knot","Dunlin"))


xyplot(YX14~Year14|ID14,xlab="Time",ylab="Variable",
   layout=c(4,4),
   scales=list(alternating=T,x=list(relation = "same"),
               y = list(relation = "free")),
   panel=function(x,y){
     panel.xyplot(x,y,col=1)
     panel.grid(h=-1,v=2)
     panel.loess(x,y,col=1,span=0.5)
     I2<-is.na(y)
     panel.text(x[I2],min(y,na.rm=T),'|',cex=0.5)  })


##########  Section 15.3 Estimation of trends


Birds7 <- as.vector(as.matrix(Ythan[,2:8]))
BirdNames<-c("Oystercatcher","Turnstone",
   "Curlew","BartailedGodwit","Redshank","Knot",
   "Dunlin")
ID7<-factor(rep(BirdNames, each = N), levels = BirdNames)
Year7<-rep(Ythan$Year,7)




Oyst.01 <- as.numeric(ID7=="Oystercatcher")
Turn.01 <- as.numeric(ID7=="Turnstone")
Curl.01 <- as.numeric(ID7=="Curlew")
Bart.01 <- as.numeric(ID7=="BartailedGodwit")
Reds.01 <- as.numeric(ID7=="Redshank")
Knot.01 <- as.numeric(ID7=="Knot")
Dunl.01 <- as.numeric(ID7=="Dunlin")

#For older R versions (2.7.0 and before), use:
f7 <- formula(Birds7 ~ ID7 +
    s(Year7, by = Oyst.01, bs = "cr") +
    s(Year7, by = Turn.01, bs = "cr") +
    s(Year7, by = Curl.01, bs = "cr") +
    s(Year7, by = Bart.01, bs = "cr") +
    s(Year7, by = Reds.01, bs = "cr") +
    s(Year7, by = Knot.01, bs = "cr") +
    s(Year7, by = Dunl.01, bs = "cr"))
    
    
#For later R versions, use:
f7 <- formula(Birds7 ~ ID7 + s(Year7, by = ID7, bs = "cr"))



library(mgcv);library(nlme)
lmc<-lmeControl(niterEM=5000,msMaxIter=1000)

M0<-gamm(f7,
      control=lmc,method="REML",
      weights=varIdent(form=~1|ID7))

AIC(M0$lme)


p1<-plot(M0$lme,resid(.,type="n")~fitted(.),abline=0,col=1)
p2<-plot(M0$lme,resid(.,type="n")~Year7,abline=0,col=1,xlab="Year")
p3<-plot(M0$lme,resid(.,type="n")~fitted(.)|ID7,
         abline=0,col=1,par.strip.text=list(cex = 0.75))


print(p1, position = c(0,0,1,1), split=c(1,1,2,2), more=TRUE)
print(p2, position = c(0,0,1,1), split=c(2,1,2,2), more=TRUE)
print(p3, position = c(0,0,2,1), split=c(1,2,2,2), more=FALSE)


plot(M0$lme,resid(.,type="n")~Year7|ID7,abline=0,col=1,xlab="Year")
########################################################


###########   Failed approach 2  ##########
Birds7 <- as.vector(as.matrix(scale(Ythan[,2:8])))
BirdNames<-c("Oystercatcher","Turnstone",
   "Curlew","BartailedGodwit","Redshank","Knot",
   "Dunlin")
ID7<-factor(rep(BirdNames, each = N), levels = BirdNames)
Year7<-rep(Ythan$Year,7)


Oyst.01 <- as.numeric(ID7=="Oystercatcher")
Turn.01 <- as.numeric(ID7=="Turnstone")
Curl.01 <- as.numeric(ID7=="Curlew")
Bart.01 <- as.numeric(ID7=="BartailedGodwit")
Reds.01 <- as.numeric(ID7=="Redshank")
Knot.01 <- as.numeric(ID7=="Knot")
Dunl.01 <- as.numeric(ID7=="Dunlin")

#For R version 2.7.0 and earlier:
f7 <- formula(Birds7 ~ ID7 +
    s(Year7, by = Oyst.01, bs = "cr") +
    s(Year7, by = Turn.01, bs = "cr") +
    s(Year7, by = Curl.01, bs = "cr") +
    s(Year7, by = Bart.01, bs = "cr") +
    s(Year7, by = Reds.01, bs = "cr") +
    s(Year7, by = Knot.01, bs = "cr") +
    s(Year7, by = Dunl.01, bs = "cr"))

#For later versions, use:
f7 <- formula(Birds7 ~ ID7 + s(Year7, by = ID7, bs = "cr"))




M2<-gamm(f7,
      control=lmc,method="REML",
      weights = varPower(form =~ fitted(.) | ID7))




#######################################################
# Adding a correlation structure
#15.4 Dealing with independence
#####Section 15.4

Birds7 <- as.vector(as.matrix(Ythan[,2:8]))
BirdNames<-c("Oystercatcher","Turnstone",
   "Curlew","BartailedGodwit","Redshank","Knot",
   "Dunlin")
ID7<-factor(rep(BirdNames, each = N), levels = BirdNames)
Year7<-rep(Ythan$Year,7)

Oyst.01 <- as.numeric(ID7=="Oystercatcher")
Turn.01 <- as.numeric(ID7=="Turnstone")
Curl.01 <- as.numeric(ID7=="Curlew")
Bart.01 <- as.numeric(ID7=="BartailedGodwit")
Reds.01 <- as.numeric(ID7=="Redshank")
Knot.01 <- as.numeric(ID7=="Knot")
Dunl.01 <- as.numeric(ID7=="Dunlin")

#For R version 2.7.0 and before:
f7 <- formula(Birds7 ~ ID7 +
    s(Year7, by = Oyst.01, bs = "cr") +
    s(Year7, by = Turn.01, bs = "cr") +
    s(Year7, by = Curl.01, bs = "cr") +
    s(Year7, by = Bart.01, bs = "cr") +
    s(Year7, by = Reds.01, bs = "cr") +
    s(Year7, by = Knot.01, bs = "cr") +
    s(Year7, by = Dunl.01, bs = "cr"))



#For later R versions, use:
f7 <- formula(Birds7 ~ ID7 + s(Year7, by = ID7, bs = "cr"))





library(mgcv);library(nlme)
lmc<-lmeControl(niterEM=5000,msMaxIter=1000)

M0<-gamm(f7,
      control=lmc,method="REML",
      weights=varIdent(form=~1|ID7))




Vario2 <- Variogram(M0$lme, form =~ Year7 | ID7,
                  maxDist = 10, robust = TRUE,
                  resType="normalized")

plot(Vario2, lwd=2,pch=16,cex=1.5,smooth=FALSE)


#Add spatial correlations


M0A <- gamm(f7, method = "REML",
     control = lmc,
     weights=varIdent(form=~1|ID7),
     correlation = corSpher(form =~ Year7 | ID7,
                   nugget = TRUE, fixed = FALSE))

M0B <- gamm(f7, method = "REML",
     control = lmc,
     weights=varIdent(form=~1|ID7),
     correlation = corLin(form =~ Year7 | ID7,
                   nugget = TRUE, fixed = FALSE))

M0C <- gamm(f7, method = "REML",
     control = lmc,
     weights=varIdent(form=~1|ID7),
     correlation = corGaus(form =~ Year7 | ID7,
                   nugget = TRUE, fixed = FALSE))

M0D <- gamm(f7, method = "REML",
     control = lmc,
     weights=varIdent(form=~1|ID7),
     correlation = corRatio(form =~ Year7 | ID7,
                   nugget = TRUE, fixed = FALSE))

M0E <- gamm(f7, method = "REML",
     control = lmc,
     weights=varIdent(form=~1|ID7),
     correlation = corExp(form =~ Year7 | ID7,
                   nugget = TRUE, fixed = FALSE))



M0F <- gamm(f7, method = "REML",
     control = lmc,
     weights=varIdent(form=~1|ID7),
     correlation = corAR1(form =~ Year7 | ID7))


AIC(M0$lme,M0A$lme,M0B$lme,M0C$lme,M0D$lme,M0E$lme,M0F$lme)



anova(M0$gam)


#######################################################

#Three trends versus 7 trends

ORB.01<- as.numeric(ID7 =="Oystercatcher" |
                    ID7 =="Redshank" |
                    ID7 =="BartailedGodwit")

TCK.01 <- as.numeric(ID7=="Turnstone" |
                     ID7=="Curlew" |
                     ID7=="Knot")
D.01<- as.numeric(ID7=="Dunlin")


M0.7<-gamm(f7,
      control=lmc,method="ML",
      weights=varIdent(form=~1|ID7))


#For older R versions, use:
M0.3<-gamm(Birds7 ~ ID7 +
    s(Year7, by = ORB.01, bs = "cr") +
    s(Year7, by = TCK.01, bs = "cr") +
    s(Year7, by = D.01, bs = "cr"),
    method = "ML", control = lmc,
    weights=varIdent(form=~1|ID7))


#For R versons 2.7.1 and above, use:
ID3 <- vector(length=length(ID7))
ID3[ID7 == "Oystercatcher" |
    ID7 == "Redshank" |
    ID7 == "BartailedGodwit"] <- "ORB"
ID3[ID7 == "Turnstone" |
    ID7 == "Curlew" |
    ID7 == "Knot"] <- "TCK"
ID3[ID7 == "Dunlin"] <- "D"
fID3 <- factor(ID3)

M0.3<-gamm(Birds7 ~ ID7 +
    s(Year7, by = fID3, bs = "cr"),
    method = "ML", control = lmc,
    weights=varIdent(form=~1|ID7))




AIC(M0.7$lme,M0.3$lme)



P0<-predict(M0.7$gam, se = TRUE)
Isna<-is.na(Birds7)

F <- P0$fit
Fup <-P0$fit + 1.96* P0$se.fit
Flow <-P0$fit - 1.96* P0$se.fit

xyplot(F+Fup+Flow~ Year7[!Isna] | ID7[!Isna], xlab = "Time",
   ylab = "Fitted values",    lty=c(1,2,2),col=1,
   type=c("l","l","l"),
   scales = list(alternating = TRUE,
               x = list(relation = "same"),
               y = list(relation = "free")))






############################################################
######### To transform or not to transform
#Show that square root transformed data removes trends


Birds7 <- as.vector(as.matrix(sqrt(Ythan[,2:8])))

#For older R versions (<2.7.0), use:
f7 <- formula(Birds7 ~ ID7 +
    s(Year7, by = Oyst.01, bs = "cr") +
    s(Year7, by = Turn.01, bs = "cr") +
    s(Year7, by = Curl.01, bs = "cr") +
    s(Year7, by = Bart.01, bs = "cr") +
    s(Year7, by = Reds.01, bs = "cr") +
    s(Year7, by = Knot.01, bs = "cr") +
    s(Year7, by = Dunl.01, bs = "cr"))

#For more recent R versions (>2.7.0), use:
f7 <- formula(Birds7 ~ ID7 +
    s(Year7, by = ID7, bs = "cr"))



M0<-gamm(f7,method="REML",
       control = lmc,
       weights=varIdent(form=~1|ID7))



p1<-plot(M0$lme,resid(.,type="n")~fitted(.),abline=0,col=1)
p2<-plot(M0$lme,resid(.,type="n")~Year7,abline=0,col=1,xlab="Year")
p3<-plot(M0$lme,resid(.,type="n")~fitted(.)|ID7,
         abline=0,col=1,par.strip.text=list(cex = 0.75))


print(p1, position = c(0,0,1,1), split=c(1,1,2,2), more=TRUE)
print(p2, position = c(0,0,1,1), split=c(2,1,2,2), more=TRUE)
print(p3, position = c(0,0,2,1), split=c(1,2,2,2), more=FALSE)



plot(M0$gam)

#Does adding an auto-correlation help here?
M0.plain<-gamm(f7,method="REML",
       control = lmc)


M0<-gamm(f7,method="REML",
       control = lmc,
       weights=varIdent(form=~1|ID7))
       
M0E <- gamm(f7, method = "REML",
       control = lmc,
       weights=varIdent(form=~1|ID7),
       correlation = corExp(form =~ Year7 | ID7,
                   nugget = TRUE, fixed = FALSE))

M0F <- gamm(f7, method = "REML",
       control = lmc,
       weights=varIdent(form=~1|ID7),
       correlation = corAR1(form =~ Year7 | ID7))

AIC(M0$lme,M0E$lme,M0F$lme)
AIC(M0$lme,M0.plain$lme)

plot(M0$lme)


P0<-predict(M0$gam, se = TRUE)
Isna<-is.na(Birds7)

F <- P0$fit
Fup <-P0$fit + 1.96* P0$se.fit
Flow <-P0$fit - 1.96* P0$se.fit

xyplot(F+Fup+Flow~ Year7[!Isna] | ID7[!Isna], xlab = "Time",
   ylab = "Fitted values",    lty=c(1,2,2),col=1,
   type=c("l","l","l"),
   scales = list(alternating = TRUE,
               x = list(relation = "same"),
               y = list(relation = "free")))
               

###########################################################


###########################################################





#Figure 15.6




X<-as.vector(as.matrix(Ythan[,9:15]))
IDX<-factor(rep(names(Ythan[,9:15]), each = 27))
Xnames<-names(Ythan[,9:15])


FF<-0
YY<-0
SSlow<-0
SSup<-0
II<-0
for (i in 1:7){
 print(i)
 x1<-seq(from=1969,to=1995)
 y1<-X[IDX==Xnames[i]]
 i1<-!is.na(y1)
 y2<-scale(y1[i1])
 x2<-x1[i1]
 if(length(x2)>5) {
     tmp<-loess(y2~x2,span=0.7)
     z2<-seq(min(x2),max(x2))
     mydata=data.frame(x2=z2)
     tmp2<-predict(tmp,mydata,se=T)
     n<-length(z2)
     FF<-c(FF,tmp2$fit)
     YY<-c(YY,z2)
     SSlow<-c(SSlow,tmp2$fit-2*tmp2$se.fit)
     SSup<-c(SSup,tmp2$fit+2*tmp2$se.fit)
     II<-c(II,rep(Xnames[i],each=length(z2)))
     }}


n<-length(FF)
FF1<-as.vector(FF[2:n])
YY1<-as.vector(YY[2:n])
SSlow1<-as.vector(SSlow[2:n])
SSup1<-as.vector(SSup[2:n])
II1<-II[2:n]

Isna<-is.na(Birds7)

FF2<-c(FF1,F)
SSup2<-c(SSup1,Fup)



II2<-c(II1,rep(BirdNames, each = N)[!Isna])
II2<-factor(II2,
   levels=c("wheat","nitrate","sheep","cattle","oats","pigs",
           "barley","Oystercatcher","Turnstone",
   "Curlew","BartailedGodwit","Redshank","Knot",
   "Dunlin"))

YY2<-c(YY1, Year7[!Isna])
SSlow2<-c(SSlow1,Flow)

xyplot(FF2~YY2|II2,xlab="Time",ylab="Smoothers",xlim=c(1968,1996),
     scales = list(alternating = T,
          x = list(relation = "same",tick.number=5),
                   y = list(relation = "same")),
     panel=function(x,y,subscripts){
       p1<-SSlow2[subscripts]
       p2<-SSup2[subscripts]
       panel.lines(x,y,col=1)
       panel.lines(x,p1,col=1,lty=2)
       panel.lines(x,p2,col=1,lty=2)
       panel.grid(h=-1,v=2)
       })


