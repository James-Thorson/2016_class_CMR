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



library(AED) ; data(Hawaii)

Hawaii$Birds<-sqrt(Hawaii$Moorhen.Kauai)
plot(Hawaii$Year,Hawaii$Birds,xlab="Year",
     ylab="Moorhen abundance on Kauai")


library(nlme)
M0 <- gls(Birds ~ Rainfall + Year, na.action = na.omit, data = Hawaii)
summary(M0)


E<-residuals(M0,type="normalized")
I<-!is.na(Hawaii$Birds)
Efull<-vector(length=length(Hawaii$Birds))
Efull<-NA
Efull[I]<-E
acf(Efull,na.action = na.pass,
    main="Auto-correlation plot for residuals")

M1<-gls(Birds ~ Rainfall + Year, na.action = na.omit,
        correlation = corCompSymm(form =~ Year),
        data=Hawaii)

M2<-gls(Birds ~ Rainfall + Year, na.action = na.omit,
      correlation = corAR1(form =~ Year), data = Hawaii)
summary(M2)




cs1 <- corARMA(c(0.2), p = 1, q = 0)
cs2 <- corARMA(c(0.3, -0.3), p = 2, q = 0)
M3arma1<-gls(Birds~Rainfall+Year,na.action=na.omit,
            correlation=cs1,data = Hawaii)
M3arma2<-gls(Birds~Rainfall+Year,na.action=na.omit,
            correlation=cs2,data = Hawaii)
AIC(M3arma1,M3arma2)




### Multivariate analysis

library(AED) ; data(Hawaii)

Birds<-c(Hawaii$Stilt.Oahu,Hawaii$Stilt.Maui,Hawaii$Coot.Oahu,Hawaii$Coot.Maui)
Time<-rep(Hawaii$Year,4)
Rain<-rep(Hawaii$Rainfall,4)
ID<-factor(rep(c("Stilt.Oahu","Stilt.Maui",
         "Coot.Oahu","Coot.Maui"),each=length(Hawaii$Year)))
library(lattice)
xyplot(Birds~Time|ID,col=1)
library(mgcv)
#This code was fitted using R version 2.6.
BM1<-gamm(Birds~Rain+ID+
        s(Time,by=as.numeric(ID=="Stilt.Oahu"))+
        s(Time,by=as.numeric(ID=="Stilt.Maui"))+
        s(Time,by=as.numeric(ID=="Coot.Oahu"))+
        s(Time,by=as.numeric(ID=="Coot.Maui")),
      weights=varIdent(form=~1|ID))

#For more recent R versions, you need to change to code:
#See also the help file of  ?gam.models

vf1 <- varIdent(form=~1|ID)
BM1<-gamm(Birds ~ Rain + ID + s(Time,by=ID),
      weights=varIdent(form=~1|ID))


summary(BM1$gam)

BM2<-gamm(Birds~Rain+ID+
       s(Time,by=as.numeric(ID=="Stilt.Oahu"))+
       s(Time,by=as.numeric(ID=="Stilt.Maui"))+
       s(Time,by=as.numeric(ID=="Coot.Oahu"))+
       s(Time,by=as.numeric(ID=="Coot.Maui")),
     correlation=corAR1(form=~Time |ID ),
     weights=varIdent(form=~1|ID))
AIC(BM1$lme, BM2$lme)

#Or recent R code:
BM2<-gamm(Birds~Rain+ID+
       s(Time,by=ID),
     correlation=corAR1(form=~Time |ID ),
     weights=varIdent(form=~1|ID))



anova(BM2$gam)


E2 <- resid(BM2$lme, type="normalized")
EAll<-vector(length = length(Birds))
EAll[] <- NA
I <- !is.na(Birds)
EAll[I] <- E2
library(lattice)
xyplot(EAll ~ Time | ID, col = 1, ylab = "Residuals")



E1<-EAll[ID=="Stilt.Oahu"]
E2<-EAll[ID=="Stilt.Maui"]
E3<-EAll[ID=="Coot.Oahu"]
E4<-EAll[ID=="Coot.Maui"]
par(mfrow=c(2,2))
acf(E1,na.action = na.pass)
acf(E2,na.action = na.pass)
acf(E3,na.action = na.pass)
acf(E4,na.action = na.pass)
D<-cbind(E1,E2,E3,E4)
cor(D,use="pairwise.complete.obs")





################


library(AED) ; data(Owls)
names(Owls)

# "FoodTreatment"      "SexParent"
#[4] "ArrivalTime"        "SiblingNegotiation" "BroodSize"
#[7] "NegPerChick"

library(AED) ; data(Owls)
library(nlme)
Owls$LogNeg<-log10(Owls$NegPerChick+1)
Form<-formula(LogNeg~SexParent*FoodTreatment+SexParent*ArrivalTime)
M2.gls<-gls(Form,correlation=corCompSymm(form=~1|Nest),
            method="REML",
            data=Owls)


library(lattice)
xyplot(LogNeg~ArrivalTime|Nest,data=Owls,type="h",col=1,
    subset=(FoodTreatment=="Deprived"),main="Deprived")




M1.lme=lme(Form,random=~1|Nest,method="REML",data=Owls)
M2.gls<-gls(Form,correlation=corCompSymm(form=~1|Nest),
            method="REML",
            data=Owls)
            
M3.gls<-gls(Form,correlation=corAR1(form=~1|Nest/FoodTreatment),
            method="REML",
            data=Owls)


M4 <- lme(Form,random=~1|Nest,method="REML",data=Owls,
     correlation=corAR1(form=~1|Nest/FoodTreatment),)

summary(M2.gls)
summary(M3.gls)



xyplot(LogNeg~ArrivalTime|Nest,data=Owls,type="h",col=1,
    subset=(FoodTreatment=="Satiated"))


