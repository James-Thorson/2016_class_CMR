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




library(AED); data(Squid)
Squid$fMONTH=factor(Squid$MONTH)
M1 <- lm(Testisweight ~ DML * fMONTH,data=Squid)
op <- par(mfrow = c(2,2), mar=c(4,4,2,2))
plot(M1, which=c(1), col=1, add.smooth=F, caption="")
plot(Squid$fMONTH, resid(M1), xlab="Month",
     ylab="Residuals")
plot(Squid$DML, resid(M1),xlab="DML",ylab="Residuals")
par(op)


library(nlme)
M.lm<-gls(Testisweight~DML*fMONTH,data=Squid)
vf1Fixed<-varFixed(~DML)
M.gls1<-gls(Testisweight~DML*fMONTH,
              weights=vf1Fixed,data=Squid)
anova(M.lm,M.gls1)


vf2 <- varIdent(form= ~ 1|fMONTH)
M.gls2 <- gls(Testisweight ~ DML*fMONTH,
              weights=vf2, data =Squid)
anova(M.lm,M.gls1,M.gls2)


anova(M.lm,M.gls2)
summary(M.gls2)

library(lattice)
E <- resid(M.lm)
coplot(E~DML|fMONTH,data=Squid)


vf3 <- varPower(form =~ DML)
M.gls3 <- gls(Testisweight ~ DML * fMONTH,
              weights = vf3,data=Squid)

vf4 <- varPower(form=~ DML | fMONTH)
M.gls4<-gls(Testisweight ~ DML * fMONTH, data = Squid,
              weights = vf4)


vf5 <- varExp(form =~ DML)
M.gls5 <- gls(Testisweight ~ DML * fMONTH,
              weights = vf5, data = Squid)


vf6<-varConstPower(form =~ DML)
M.gls6<-gls(Testisweight ~ DML * fMONTH,
            weights = vf6, data = Squid)


vf7 <- varConstPower(form =~ DML | fMONTH)
M.gls7<-gls(Testisweight ~ DML * fMONTH,
              weights = vf7, data = Squid)

vf8 <- varComb(varIdent(form= ~ 1 | fMONTH) ,
                 varExp(form =~ DML) )
M.gls8<-gls(Testisweight ~ DML * fMONTH,
              weights = vf8, data = Squid)

anova(M.lm,M.gls1,M.gls2,M.gls3,M.gls4,
        M.gls5,M.gls6,M.gls7,M.gls8)



anova(M.lm,M.gls4)


E1 <- resid(M.gls4)
coplot(E1 ~ DML | fMONTH,ylab="Ordinary residuals",
      data = Squid)

E2 <- resid(M.gls4, type = "normalized")
coplot(E2 ~ DML | fMONTH, data = Squid,
       ylab = "Normalised residuals")


anova(M.gls4)


#####################################################
#Second example

library(AED); data(Biodiversity)
Biodiv <- Biodiversity #saves some space
Biodiv$fTreatment <- factor(Biodiv$Treatment)
Biodiv$fNutrient  <- factor(Biodiv$Nutrient)

boxplot(Concentration~factor(Treatment)*factor(Nutrient),data=Biodiv)
M0<-lm(Concentration~Biomass*factor(Treatment)*factor(Nutrient),data=Biodiv)
plot(M0,which=c(1),add.smooth=FALSE)





M0<-gls(Concentration~Biomass*fTreatment*fNutrient,
      data = Biodiv)
M1A<-gls(Concentration~Biomass*fTreatment*
    fNutrient,weights=varIdent(
    form=~1|fTreatment*fNutrient),
    data = Biodiv)
M1B<-gls(Concentration~Biomass*fTreatment*
    fNutrient,
    weights=varIdent(form=~1|fNutrient),
    data=Biodiv)
M1C<-gls(Concentration~Biomass*fTreatment*
    fNutrient,
    weights=varIdent(form=~1|fTreatment),
    data = Biodiv)


anova(M0,M1A,M1B,M1C)
anova(M1A)




library(nlme)
M2A1<-gls(Concentration ~ Biomass + fTreatment +
           fNutrient +
           Biomass:fTreatment +
           Biomass:fNutrient +
           fTreatment:fNutrient +
           Biomass:fTreatment:fNutrient,
           weights = varIdent(form =~ 1 | fTreatment *
                                         fNutrient),
           method = "ML", data = Biodiv)

M2A2<-gls(Concentration ~ Biomass + fTreatment +
           Nutrient +
           Biomass:fTreatment +
           Biomass:fNutrient +
           fTreatment:fNutrient,
           weights=varIdent(form =~ 1 |
                                      fTreatment*
                                      fNutrient),
          method="ML", data = Biodiv)

anova(M2A1,M2A2)



vfOptim<-varIdent(form=~1|fTreatment*fNutrient)
#Assess significance of all three two-way interactions
#Full model
M3.Full<-gls(Concentration~Biomass+fTreatment+fNutrient+
          Biomass:fTreatment+
          Biomass:fNutrient+
          fTreatment:fNutrient,
          weights=vfOptim,
          method="ML",data=Biodiv)
          

#Drop Biomass:fTreatment
M3.Drop1<-gls(Concentration~Biomass+fTreatment+fNutrient+
          Biomass:fNutrient+
          fTreatment:fNutrient,
          weights=vfOptim,
          method="ML",data=Biodiv)
anova(M3.Full,M3.Drop1)

#Drop Biomass:fNutrient
M3.Drop2<-gls(Concentration~Biomass+fTreatment+fNutrient+
          Biomass:fTreatment+
          fTreatment:fNutrient,
          weights=vfOptim,
          method="ML",data=Biodiv)
anova(M3.Full,M3.Drop2)


#factor(Treatment):fNutrient
M3.Drop3<-gls(Concentration~Biomass+fTreatment+fNutrient+
          Biomass:fTreatment+
          Biomass:fNutrient,
          weights=vfOptim,
          method="ML",data=Biodiv)
anova(M3.Full,M3.Drop3)




#Conclusion: drop Biomass:fTreatment



############################################
#Alternative coding
fFull<-formula(Concentration~Biomass+fTreatment+fNutrient+
          Biomass:fTreatment+
          Biomass:fNutrient+
          fTreatment:fNutrient)
          
M3.Full<-gls(fFull,
          weights=vfOptim,
          method="ML",data=Biodiv)


#Drop Biomass:fTreatment
M3.Drop1<-update(M3.Full,.~.-Biomass:fTreatment)
anova(M3.Full,M3.Drop1)

#Drop Biomass:fNutrient
M3.Drop2<-update(M3.Full,.~.-Biomass:fNutrient)
anova(M3.Full,M3.Drop2)


#fTreatment:fNutrient
M3.Drop3<-update(M3.Full,.~.-fTreatment:fNutrient)
anova(M3.Full,M3.Drop3)




#########  New round
#Full model
M4.Full<-gls(Concentration~Biomass+fTreatment+fNutrient+
          Biomass:fNutrient+
          fTreatment:fNutrient,
          weights=vfOptim,
          method="ML",data=Biodiv)

#Drop Biomass:fNutrient
M4.Drop1 <- update(M4.Full, .~. -Biomass:fNutrient )
anova(M4.Full,M4.Drop1)

#Drop fTreatment:fNutrient
M4.Drop2 <- update(M4.Full, .~. -fTreatment:fNutrient )
anova(M4.Full,M4.Drop2)





##############################
#New full model
M5.Full<-gls(Concentration~Biomass+fTreatment+fNutrient+
          fTreatment:fNutrient,
          weights=vfOptim,
          method="ML",data=Biodiv)

#Drop fTreatment:fNutrient
M5.Drop1 <- update(M5.Full, .~. -fTreatment:fNutrient)
anova(M5.Full,M5.Drop1)

#Drop Biomass
M5.Drop2 <- update(M5.Full, .~. -Biomass)
anova(M5.Full,M5.Drop2)




##############################
#New full model
M6.Full<-gls(Concentration~fTreatment+fNutrient+
          fTreatment:fNutrient,
          weights=vfOptim,
          method="ML",data=Biodiv)

M6.Drop1 <- update(M6.Full, .~. -fTreatment:fNutrient)
anova(M6.Full,M6.Drop1)



#The aftermath
MFinal<-gls(Concentration~fTreatment+fNutrient+
          fTreatment:fNutrient,
          weights=vfOptim,
          method="REML",data=Biodiv)


E<-resid(MFinal,type="normalized")
Fit=fitted(MFinal)

op<-par(mfrow=c(1,2))
plot(x=Fit,y=E,xlab="Fitted values",ylab="Residuals",
main="Residuals versus fitted values")
identify(Fit,E)
hist(E,nclass=15)
par(op)

summary(MFinal)


boxplot(E~fTreatment*fNutrient,data=Biodiv)


PV=predict(MFinal)

boxplot(E~fTreatment*fNutrient,data=Biodiv)


boxplot(predict(MFinal)~fTreatment*fNutrient,data=Biodiv)










