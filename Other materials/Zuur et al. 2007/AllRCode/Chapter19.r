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



library(AED); data(Bees)

Bees$fHive <- factor(Bees$Hive)
Bees$LSpobee <- log10(Bees$Spobee + 1)

op<- par(mfrow = c(1, 2), mar = c(3, 4, 1, 1))
dotchart(Bees$Spobee, groups = Bees$fHive)
dotchart(Bees$LSpobee, groups = Bees$fHive)
par(op)


Bees$Infection01 <- Bees$Infection
Bees$Infection01[Bees$Infection01 > 0] <- 1
Bees$fInfection01 <- factor(Bees$Infection01)

boxplot(LSpobee ~ fInfection01, data = Bees, varwidth = TRUE)


M1 <- lm(LSpobee ~ fInfection01 * BeesN, data = Bees)
E1 <- rstandard(M1)
plot(E1 ~ Bees$fHive, ylab = "Standardised residuals", xlab = "Hives")
abline(0, 0)



library(nlme)
M2<-gls(LSpobee ~ fInfection01 * BeesN, data = Bees)
M3<-lme(LSpobee ~ fInfection01 * BeesN,
   random =~ 1 | fHive, data = Bees)
anova(M2,M3)



M4<-lme(LSpobee ~ fInfection01 * BeesN,
          random =~ 1 + BeesN | fHive, data = Bees)
M5<-lme(LSpobee ~ fInfection01 * BeesN,
          random =~ 1 + fInfection01 | fHive, data = Bees)
anova(M2,M3,M4,M5)



plot(M3, col = 1)

M6<-lme(LSpobee ~ fInfection01 * BeesN,
          random =~ 1 | fHive,
          weights = varIdent(form =~ 1 | fInfection01), data =  Bees)
anova(M3,M6)




M7full<-lme(LSpobee ~ fInfection01 * BeesN,
         random=~1|fHive,
         weights=varIdent(form =~ 1 | fInfection01),
         method="ML", data = Bees)
M7sub<-update(M7full, .~. -fInfection01 : BeesN )
anova(M7full,M7sub)




M8full<-lme(LSpobee ~ fInfection01 + BeesN,
          random =~ 1|fHive, method="ML", data = Bees,
          weights = varIdent(form =~ 1 | fInfection01))

M8sub1<-update(M8full, .~. -fInfection01 )
M8sub2<-update(M8full, .~. -BeesN )
anova(M8full,M8sub1)
anova(M8full,M8sub2)


M9full<-lme(LSpobee ~ fInfection01,
          random =~ 1|fHive, method="ML", data = Bees,
          weights = varIdent(form =~ 1 | fInfection01))

M9sub1<-update(M9full, .~. -fInfection01 )
anova(M9full,M9sub1)


Mfinal<-lme(LSpobee ~ fInfection01,
          random =~ 1|fHive, data = Bees,
          weights = varIdent(form =~ 1 | fInfection01),
          method="REML")
summary(Mfinal)


plot(Mfinal)
qqnorm(Mfinal)
qqnorm(Mfinal,~ranef(.),col=1)


intervals(Mfinal)