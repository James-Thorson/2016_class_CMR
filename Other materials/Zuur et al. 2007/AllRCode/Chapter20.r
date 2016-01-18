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



library(AED); data(Cetaceans)

Cetaceans$fSpecies <- factor(Cetaceans$Species)
Cetaceans$fDolphinID <- factor(Cetaceans$DolphinID)


boxplot(Age~fSpecies, data = Cetaceans)
boxplot(Age~fDolphinID, data = Cetaceans)
I<-Cetaceans$Sex==0
Cetaceans2<-Cetaceans[!I,]


library(nlme)
Cetaceans2$fSex <- factor(Cetaceans2$Sex)
Cetaceans2$fLocation <- factor(Cetaceans2$Location)
Cetaceans2$fStain <- factor(Cetaceans2$Stain)
f1<-formula(Age ~ fSex * fStain * fLocation)
M1<-gls(f1,method="REML", data = Cetaceans2)




M2 <- lme(f1, random =~1 | fSpecies/fDolphinID, data = Cetaceans2, method="REML")
anova(M1,M2)


M3<-lme(f1,random=~1|fSpecies/fDolphinID,
           weights=varIdent(form=~1|fLocation),data = Cetaceans2)

f1<-formula(Age ~ fStain + fLocation + fStain : fLocation)

M3<-lme(f1,random=~1|fSpecies/fDolphinID, method = "REML",
           weights=varIdent(form=~1|fLocation),data = Cetaceans2)
options(digits=4)
summary(M3)



