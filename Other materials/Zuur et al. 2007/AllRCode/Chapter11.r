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




library(AED); data(Snakes)


plot(table(Snakes$N_days),ylab="Frequencies")

library(MASS)
M1 <- glm.nb(N_days~Size_cm+PDayRain+Tot_Rain+Road+
            Road_Loc+Size_cm:PDayRain+Size_cm:Tot_Rain+
            Size_cm:Road+Size_cm:Road_Loc+
            PDayRain:Tot_Rain+ PDayRain:Road+
            PDayRain:Road_Loc+Tot_Rain:Road,
            data=Snakes)
            
M2A <- glm.nb(N_days ~ PDayRain + Tot_Rain + Road_Loc +
         PDayRain:Tot_Rain, data = Snakes)

library(VGAM)
M2B <- vglm(N_days ~ PDayRain + Tot_Rain + Road_Loc +
        PDayRain:Tot_Rain, family = negbinomial, data = Snakes)
summary(M2B)

M3A <- vglm(N_days~PDayRain+Tot_Rain+Road_Loc+
           PDayRain:Tot_Rain,
           family=posnegbinomial,
           control = vglm.control(maxit = 100),
           data=Snakes)
           
           
Z=cbind(coef(M2A),coef(M3A)[-2])
ZSE=cbind(sqrt(diag(vcov(M2A))),
          sqrt(diag(vcov(M3A))[-1]))
Comp= cbind(Z[,1],Z[,2],ZSE[,1],ZSE[,2])
Comb=round(Comp,digits=3)
colnames(Comb)=
          c("NB","Trunc.NB","SE NB","SE Trunc.NB")
Comb





#Parasite data
library(AED); data(ParasiteCod)

ParasiteCod$fArea <- factor(ParasiteCod$Area)
ParasiteCod$fYear <- factor(ParasiteCod$Year)

I=is.na(ParasiteCod$Intensity) | is.na(ParasiteCod$fArea) |
  is.na(ParasiteCod$fYear) | is.na(ParasiteCod$Length)

ParasiteCod2 <- ParasiteCod[!I,]


plot(table(ParasiteCod2$Intensity),ylab="Frequencies",
xlab="Observed intensity values")

library(pscl)
f1 <- formula(Intensity ~ fArea*fYear +
               Length | fArea*fYear + Length)
Zip1 <- zeroinfl(f1, dist = "poisson", link = "logit",
              data = ParasiteCod2)

f1 <- formula(Intensity ~ fArea * fYear +
               Length | fArea*fYear + Length)
               

Nb1 <- zeroinfl(f1, dist = "negbin", link = "logit",
              data = ParasiteCod2)
library(lmtest)
lrtest(Zip1,Nb1)


summary(Nb1)



#Full model
f1A <-formula(Intensity ~ fArea * fYear |
               fArea * fYear + Length)
#Drop interaction from count model
f1B <-formula(Intensity ~ fArea + fYear+
        Length | fArea * fYear + Length)
#Drop length from binomial model
f1C<-formula(Intensity ~ fArea * fYear+
        Length | fArea * fYear)
#Drop interaction from binomial model
f1D<-formula(Intensity ~ fArea * fYear+
        Length | fArea + fYear + Length)
        
Nb1A = zeroinfl(f1A, dist = "negbin", link = "logit", data = ParasiteCod2)
Nb1B = zeroinfl(f1B, dist = "negbin", link = "logit", data = ParasiteCod2)
Nb1C = zeroinfl(f1C, dist = "negbin", link = "logit", data = ParasiteCod2)
Nb1D = zeroinfl(f1D, dist = "negbin", link = "logit", data = ParasiteCod2)

lrtest(Nb1,Nb1A) ;lrtest(Nb1,Nb1B)
lrtest(Nb1,Nb1C) ;lrtest(Nb1,Nb1D)
summary(Nb1B)

EstPar <- coef(Nb1B,model="zero")
Z <- model.matrix(Nb1B,model="zero")
g <- Z %*% EstPar
p <- exp(g) / (1 + exp(g))


EstPar2 <- coef(Nb1B, model = "count")
X <- model.matrix(Nb1B, model = "count")
g <- X %*% EstPar2
mu1 <- exp(g)


mu <- (1 - p) * mu1
k <- Nb1B$theta
VarY <- ((mu^2) / k + mu) * (1 - p) +
           (mu^2) * (p^2 + p)
 EP <- (ParasiteCod2$Intensity - mu) / sqrt(VarY)


H1A <- hurdle(f1, dist = "poisson", link = "logit", data =ParasiteCod2)
H1B <- hurdle(f1, dist = "negbin", link = "logit", data =ParasiteCod2)


summary(H1B)

fFinal <- formula(Intensity ~ fArea + fYear +
               Length | fArea*fYear )

HFinal <- hurdle(f1, dist = "negbin", link = "logit", data =ParasiteCod2)
summary(HFinal)









##############################
library(MASS)
#library(VGAM)
library(pscl)

#######################################

#Code for ZIP and hurdle, and comaprisons

##############################################################################
# function to calculate the performance statistics

MyStatistics = function(Atmp1,Atmp2,Atmp3,Atmp4,Atmp5,AIntensity,
                       WhatToDo=1,P1=Null,P2=Null,P3=Null,P4=Null,P5=Null){
        if (WhatToDo==1){
            #Steps 2 to 5
            Z1=stats::predict(Atmp1,type="response")
            Z2=stats::predict(Atmp2,type="response")
            Z3=stats::predict(Atmp3,type="response")
            Z4=stats::predict(Atmp4,type="response")
            Z5=stats::predict(Atmp5,type="response")
        }
        if (WhatToDo==2){
            #Step 6
            Z1=P1
            Z2=P2
            Z3=P3
            Z4=P4
            Z5=P5
        }

# Pearson prediction correlation with the true data

        Pois.cor  = cor(AIntensity,Z1)
        QPois.cor = cor(AIntensity,Z2)
        NB.cor    = cor(AIntensity,Z3)
        ZIP.cor   = cor(AIntensity,Z4)
        Hur.cor   = cor(AIntensity,Z5)
        Cor5=c(Pois.cor,QPois.cor,NB.cor,ZIP.cor,Hur.cor)

# Spearman prediction correlation with the true data

        Pois.Rcor  = cor(AIntensity,Z1,method = "spearman")
        QPois.Rcor = cor(AIntensity,Z2,method = "spearman")
        NB.Rcor    = cor(AIntensity,Z3,method = "spearman")
        ZIP.Rcor   = cor(AIntensity,Z4,method = "spearman")
        Hur.Rcor   = cor(AIntensity,Z5,method = "spearman")
        RCor5=c(Pois.Rcor,QPois.Rcor,NB.Rcor,ZIP.Rcor,Hur.Rcor)

# Regression coefficients (bias and shrinkage)

        regPois  = coef(lm(AIntensity~Z1))
        regQPois = coef(lm(AIntensity~Z2))
        regNB    = coef(lm(AIntensity~Z3))
        regZIP   = coef(lm(AIntensity~Z4))
        regHur   = coef(lm(AIntensity~Z5))

# Root mean square errors

        n=length(AIntensity)
        RMSE1 = sqrt(sum((Z1-AIntensity)^2)/n)
        RMSE2 = sqrt(sum((Z2-AIntensity)^2)/n)
        RMSE3 = sqrt(sum((Z3-AIntensity)^2)/n)
        RMSE4 = sqrt(sum((Z4-AIntensity)^2)/n)
        RMSE5 = sqrt(sum((Z5-AIntensity)^2)/n)
        RMSE=c(RMSE1,RMSE2,RMSE3,RMSE4,RMSE5)

# Mean absolute error (MAE was used to avoid the zero divide in weigths calculation)

        AVE1=sum(abs(Z1-AIntensity))/n
        AVE2=sum(abs(Z2-AIntensity))/n
        AVE3=sum(abs(Z3-AIntensity))/n
        AVE4=sum(abs(Z4-AIntensity))/n
        AVE5=sum(abs(Z5-AIntensity))/n
        AVE=c(AVE1,AVE2,AVE3,AVE4,AVE5)

        list(Cor.app=Cor5,RCor.app=RCor5,
             reg.app=cbind(regPois,regQPois,regNB,regZIP,regHur),
             RMSE.app=RMSE,AVE.app=AVE)
}
##############################################################################

#Parasite data
library(AED); data(ParasiteCod)

I=is.na(Data$Intensity) | is.na(Data$Area) |
  is.na(Data$Year) | is.na(Data$Length)

Data2=Data[!I,]

names(Data2)
attach(Data2)


library(MASS)
#library(VGAM)
library(pscl)

# Fit model on the full data set to calculate the apparent value of the statistics.

tmp1 = glm(Intensity~factor(Area)*factor(Year)+Length,family=poisson)
tmp2 = glm(Intensity~factor(Area)*factor(Year)+Length,family=quasipoisson)
tmp3 = glm.nb(Intensity~factor(Area)*factor(Year)+Length)
tmp4 = zeroinfl(Intensity~factor(Area)*factor(Year)+Length|
               factor(Area)*factor(Year)+Length,
               zero.dist="binomial",dist="negbin",link="logit")

tmp5 = hurdle(Intensity~factor(Area)*factor(Year)+Length|
              factor(Area)*factor(Year)+Length,
              zero.dist="binomial",dist="negbin",link="logit")


# Prepare bootstap selections and weigths (Harrell FE. Design library. http://lib.stat.cmu.edu/S/Harrell/)

n.boot <- 200
n.data <- nrow(Data2)

boot.Sel.ready<-FALSE
for(i in 1:10) { # 10 attempts to create proper bootstrapp selections ...
  boot.Sel <- matrix(integer(1),nrow=n.data,ncol=n.boot)
  boot.Wts <- matrix(TRUE,      nrow=n.data,ncol=n.boot)
  for(i in 1:n.boot) {
    boot.Sel[,       i] <- (cur.sel <- sample(n.data,replace=TRUE))
    boot.Wts[cur.sel,i] <- FALSE  # selected samples
  }
  n.omit<-apply(boot.Wts,1,sum)
  if(all(n.omit > 0)) { boot.Sel.ready<-TRUE; break; }
}
if(!boot.Sel.ready) {
  stop("not every observation omitted at least once in bootstrap samples.\nRe--run with larger n.boot")
}
boot.Wts <- apply(boot.Wts/n.omit, 2, sum)/n.data
hist(boot.Wts)

##############################################################################
# Bootstrap approach Potts and Elith

#-----------------------------------------------------------------------------
# statistics for the weights calculation from Steyerberg et al. (2001)

Stat.no_inf<-list(Cor.app =rep(0,5),
                  RCor.app=rep(0,5),
                  reg.app =matrix(rep(0,2*5),2),
                  RMSE.app=rep(0,5),
                  AVE.app=rep(0,5) )
Stat.test  <-list(Cor.app =rep(0,5),
                  RCor.app=rep(0,5),
                  reg.app =matrix(rep(0,2*5),2),
                  RMSE.app=rep(0,5),
                  AVE.app =rep(0,5) )
test.wts   <-rep(0,n.data)
test.fault <-rep(FALSE,n.boot)

#-----------------------------------------------------------------------------
#Step 2:

Stat.App=MyStatistics(tmp1,tmp2,tmp3,tmp4,tmp5,Intensity)

#-----------------------------------------------------------------------------
#calculate weights from Steyerberg et al. (2001)

for(i.boot in 1:n.boot) {
  cat("*")
  print(i.boot)
  Iv = seq(1:n.data)
  Is = boot.Sel[,i.boot]

# Remember which ones have been selected

  IsSelected     = rep(FALSE,n.data)
  IsSelected[Is] = TRUE #IsSelected == TRUE if sample was selected selected

  Data3=Data2[Is,]

# Fit models on the selected data

  Intensity.b = Data3$Intensity
  Area.b      = Data3$Area
  Year.b      = Data3$Year
  Length.b    = Data3$Length

  Btmp1 = glm(Intensity.b~factor(Area.b)*factor(Year.b)+Length.b,family=poisson)
  Btmp2 = glm(Intensity.b~factor(Area.b)*factor(Year.b)+Length.b,family=quasipoisson)
  Btmp3 = glm.nb(Intensity.b~factor(Area.b)*factor(Year.b)+Length.b)
  Btmp4 = try(zeroinfl(Intensity.b~factor(Area.b)*factor(Year.b)+Length.b|
                factor(Area.b)*factor(Year.b)+Length.b,
                zero.dist="binomial",dist="poisson",link="logit"),silent=TRUE)
  if(class(Btmp4)=="try-error") { test.fault[i.boot]<-T; next; }

  Btmp5 = try(hurdle(Intensity.b~factor(Area.b)*factor(Year.b)+Length.b|
                factor(Area.b)*factor(Year.b)+Length.b,
                zero.dist="binomial",dist="negbin",link="logit"),silent=TRUE)
  if(class(Btmp5)=="try-error") { test.fault[i.boot]<-T; next; }

# Not selected data

  DataNS=Data2[!IsSelected,]
  MyNewData=data.frame(Area.b=DataNS$Area,Year.b=DataNS$Year,Length.b=DataNS$Length)
  p1=stats::predict(Btmp1,newdata=MyNewData,type="response")
  p2=stats::predict(Btmp2,newdata=MyNewData,type="response")
  p3=stats::predict(Btmp3,newdata=MyNewData,type="response")
  p4=stats::predict(Btmp4,newdata=MyNewData,type="response")
  p5=stats::predict(Btmp5,newdata=MyNewData,type="response")

  Stat.excl=MyStatistics(Btmp1,Btmp2,Btmp3,Btmp4,Btmp5,DataNS$Intensity,
                          WhatToDo=2,P1=p1,P2=p2,P3=p3,P4=p4,P5=p5)

# Test performance evaluation

  Stat.test$Cor.app <-Stat.test$Cor.app  + Stat.excl$Cor.app *boot.Wts[i.boot]
  Stat.test$RCor.app<-Stat.test$RCor.app + Stat.excl$RCor.app*boot.Wts[i.boot]
  Stat.test$reg.app <-Stat.test$reg.app  + Stat.excl$reg.app *boot.Wts[i.boot]
  Stat.test$RMSE.app<-Stat.test$RMSE.app + Stat.excl$RMSE.app*boot.Wts[i.boot]
  Stat.test$AVE.app <-Stat.test$AVE.app  + Stat.excl$AVE.app *boot.Wts[i.boot]

# No information performance evaluation ...

  Intensity.p=sample(Data2$Intensity)
  Area.p     =Data2$Area
  Year.p     =Data2$Year
  Length.p   =Data2$Length

  Ptmp1 = glm(Intensity.p~factor(Area.p)*factor(Year.p)+Length.p,family=poisson)
  Ptmp2 = glm(Intensity.p~factor(Area.p)*factor(Year.p)+Length.p,family=quasipoisson)
  Ptmp3 = glm.nb(Intensity.p~factor(Area.p)*factor(Year.p)+Length.p)
  Ptmp4 = try(zeroinfl(Intensity.p~factor(Area.p)*factor(Year.p)+Length.p|
                 factor(Area.p)*factor(Year.p)+Length.p,
                 zero.dist="binomial",dist="negbin",link="logit"),silent=TRUE)
  if(class(Ptmp4)=="try-error") { test.fault[i.boot]<-T; next; }

  Ptmp5 = try(hurdle(Intensity.p~factor(Area.p)*factor(Year.p)+Length.p|
                factor(Area.p)*factor(Year.p)+Length.p,
                zero.dist="binomial",dist="negbin",link="logit"),silent=TRUE)
  if(class(Ptmp5)=="try-error") { test.fault[i.boot]<-T; next; }

  Pp1=stats::predict(Ptmp1,type="response")
  Pp2=stats::predict(Ptmp2,type="response")
  Pp3=stats::predict(Ptmp3,type="response")
  Pp4=stats::predict(Ptmp4,type="response")
  Pp5=stats::predict(Ptmp5,type="response")
  Stat.perm=MyStatistics(Ptmp1,Ptmp2,Ptmp3,Ptmp4,Ptmp5,Data2$Intensity,
                          WhatToDo=2,P1=Pp1,P2=Pp2,P3=Pp3,P4=Pp4,P5=Pp5)
  Stat.no_inf$Cor.app <-Stat.no_inf$Cor.app  + Stat.perm$Cor.app
  Stat.no_inf$RCor.app<-Stat.no_inf$RCor.app + Stat.perm$RCor.app
  Stat.no_inf$reg.app <-Stat.no_inf$reg.app  + Stat.perm$reg.app
  Stat.no_inf$RMSE.app<-Stat.no_inf$RMSE.app + Stat.perm$RMSE.app
  Stat.no_inf$AVE.app <-Stat.no_inf$AVE.app  + Stat.perm$AVE.app

}
cat("\n")
c.num<-sum(!test.fault)
if(c.num / n.boot < 0.5) {
  stop("Too few successes in weights calculation ...")
}
Stat.no_inf$Cor.app <-Stat.no_inf$Cor.app  / c.num
Stat.no_inf$RCor.app<-Stat.no_inf$RCor.app / c.num
Stat.no_inf$reg.app <-Stat.no_inf$reg.app  / c.num
Stat.no_inf$RMSE.app<-Stat.no_inf$RMSE.app / c.num
Stat.no_inf$AVE.app <-Stat.no_inf$AVE.app  / c.num

# Store ave(Stat.test) for the use at step 7

Stat.test.ave<-Stat.test

# Room for weigths at step 7,8

Stat.R<-list(Cor.app =rep(0,5),
             RCor.app=rep(0,5),
             reg.app =matrix(rep(0,2*5),2),
             RMSE.app=rep(0,5),
             AVE.app=rep(0,5))
Stat.w<-list(Cor.app =rep(0,5),
             RCor.app=rep(0,5),
             reg.app =matrix(rep(0,2*5),2),
             RMSE.app=rep(0,5),
             AVE.app=rep(0,5))
Stat.perf<-list(Cor.app =rep(0,5),
             RCor.app=rep(0,5),
             reg.app =matrix(rep(0,2*5),2),
             RMSE.app=rep(0,5),
             AVE.app=rep(0,5))
Stat.O<-list(Cor.app =rep(0,5),
             RCor.app=rep(0,5),
             reg.app =matrix(rep(0,2*5),2),
             RMSE.app=rep(0,5),
             AVE.app=rep(0,5))
boot.fault <-rep(FALSE,n.boot)

#-----------------------------------------------------------------------------
# All we need to calculate weights at step 7 is ready, run bootstrap ...

for(i.boot in 1:n.boot) {
  cat("*")
  print(i.boot)

  #-----------------------------------------------------------------------------
  # Step 3:

  Iv = seq(1:n.data)
  Is = boot.Sel[,i.boot]

  #Remember which ones have been selected

  IsSelected     = rep(FALSE,n.data)
  IsSelected[Is] = TRUE #IsSelected == TRUE if sample was selected selected

  #-----------------------------------------------------------------------------
  # Step 4

  Data3=Data2[Is,]

  Intensity.b = Data3$Intensity
  Area.b      = Data3$Area
  Year.b      = Data3$Year
  Length.b    = Data3$Length

  Btmp1 = glm(Intensity.b~factor(Area.b)*factor(Year.b)+Length.b,family=poisson)
  Btmp2 = glm(Intensity.b~factor(Area.b)*factor(Year.b)+Length.b,family=quasipoisson)
  Btmp3 = glm.nb(Intensity.b~factor(Area.b)*factor(Year.b)+Length.b)
  Btmp4 = try(zeroinfl(Intensity.b~factor(Area.b)*factor(Year.b)+Length.b|
                factor(Area.b)*factor(Year.b)+Length.b,
                zero.dist="binomial",dist="negbin",link="logit"),silent=TRUE)
  if(class(Btmp4)=="try-error") { boot.fault[i.boot]<-T; next; }
  Btmp5 = try(hurdle(Intensity.b~factor(Area.b)*factor(Year.b)+Length.b|
                factor(Area.b)*factor(Year.b)+Length.b,
                zero.dist="binomial",dist="negbin",link="logit"),silent=TRUE)
  if(class(Btmp5)=="try-error") { boot.fault[i.boot]<-T; next; }

  Stat.boot=MyStatistics(Btmp1,Btmp2,Btmp3,Btmp4,Btmp5,Intensity.b)

  #-----------------------------------------------------------------------------
  # Step 5: Randomise the response variable

  Is=sample(Iv,n.data,replace=FALSE)
  Stat.perm=MyStatistics(Btmp1,Btmp2,Btmp3,Btmp4,Btmp5,Intensity[Is])

  #-----------------------------------------------------------------------------
  # Step 6
  #Not selected data

  DataNS=Data2[!IsSelected,]
  MyNewData=data.frame(Area.b=DataNS$Area,Year.b=DataNS$Year,Length.b=DataNS$Length)
  p1=stats::predict(Btmp1,newdata=MyNewData,type="response")
  p2=stats::predict(Btmp2,newdata=MyNewData,type="response")
  p3=stats::predict(Btmp3,newdata=MyNewData,type="response")
  p4=stats::predict(Btmp4,newdata=MyNewData,type="response")
  p5=stats::predict(Btmp5,newdata=MyNewData,type="response")

  Stat.test=MyStatistics(Btmp1,Btmp2,Btmp3,Btmp4,Btmp5,DataNS$Intensity,
                         WhatToDo=2,P1=p1,P2=p2,P3=p3,P4=p4,P5=p5)

  #-----------------------------------------------------------------------------
  # Step 7

# calculate relative overfitting R

  Stat.R$Cor.app <-(Stat.test$Cor.app -Stat.App$Cor.app) /(Stat.no_inf$Cor.app -Stat.App$Cor.app)
  Stat.R$RCor.app<-(Stat.test$RCor.app-Stat.App$RCor.app)/(Stat.no_inf$RCor.app-Stat.App$RCor.app)
  Stat.R$reg.app <-(Stat.test$reg.app -Stat.App$reg.app) /(Stat.no_inf$reg.app -Stat.App$reg.app)
  Stat.R$RMSE.app<-(Stat.test$RMSE.app-Stat.App$RMSE.app)/(Stat.no_inf$RMSE.app-Stat.App$RMSE.app)
  Stat.R$AVE.app <-(Stat.test$AVE.app -Stat.App$AVE.app) /(Stat.no_inf$AVE.app -Stat.App$AVE.app)

# calculate weight w

  Stat.w$Cor.app <-0.632/(1.0-0.368*Stat.R$Cor.app)
  Stat.w$RCor.app<-0.632/(1.0-0.368*Stat.R$RCor.app)
  Stat.w$reg.app <-0.632/(1.0-0.368*Stat.R$reg.app)
  Stat.w$RMSE.app<-0.632/(1.0-0.368*Stat.R$RMSE.app)
  Stat.w$AVE.app <-0.632/(1.0-0.368*Stat.R$AVE.app)

# calculate a best estimate of predictive performance

  Stat.perf$Cor.app <-(1.0-Stat.w$Cor.app) *Stat.App$Cor.app  + Stat.w$Cor.app *Stat.test.ave$Cor.app
  Stat.perf$RCor.app<-(1.0-Stat.w$RCor.app)*Stat.App$RCor.app + Stat.w$RCor.app*Stat.test.ave$RCor.app
  Stat.perf$reg.app <-(1.0-Stat.w$reg.app) *Stat.App$reg.app  + Stat.w$reg.app *Stat.test.ave$reg.app
  Stat.perf$RMSE.app<-(1.0-Stat.w$RMSE.app)*Stat.App$RMSE.app + Stat.w$RMSE.app*Stat.test.ave$RMSE.app
  Stat.perf$AVE.app <-(1.0-Stat.w$AVE.app) *Stat.App$AVE.app  + Stat.w$AVE.app *Stat.test.ave$AVE.app

  #-----------------------------------------------------------------------------
  # Step 8 - Measure how optimistic the fit on the bootstrap sample was

  Stat.O$Cor.app  <-Stat.O$Cor.app +(Stat.boot$Cor.app -Stat.perf$Cor.app)
  Stat.O$RCor.app <-Stat.O$RCor.app+(Stat.boot$RCor.app-Stat.perf$RCor.app)
  Stat.O$reg.app  <-Stat.O$reg.app +(Stat.boot$reg.app -Stat.perf$reg.app)
  Stat.O$RMSE.app <-Stat.O$RMSE.app+(Stat.boot$RMSE.app-Stat.perf$RMSE.app)
  Stat.O$AVE.app  <-Stat.O$AVE.app +(Stat.boot$AVE.app -Stat.perf$AVE.app)
}
#-----------------------------------------------------------------------------
# Step 9 - Calculate an average optimism

c.num<-sum(!boot.fault)
if(c.num / n.boot < 0.5) {
  stop("Too few successes in bootdtrap ...")
}
Stat.O$Cor.app <-Stat.O$Cor.app  / c.num
Stat.O$RCor.app<-Stat.O$RCor.app / c.num
Stat.O$reg.app <-Stat.O$reg.app  / c.num
Stat.O$RMSE.app<-Stat.O$RMSE.app / c.num
Stat.O$AVE.app <-Stat.no_inf$AVE.app / c.num

# Correct for optimism

Stat.App$Cor.app  <-Stat.App$Cor.app  - Stat.O$Cor.app
Stat.App$RCor.app <-Stat.App$RCor.app - Stat.O$RCor.app
Stat.App$reg.app  <-Stat.App$reg.app  - Stat.O$reg.app
Stat.App$RMSE.app <-Stat.App$RMSE.app - Stat.O$RMSE.app
Stat.App$AVE.app  <-Stat.App$AVE.app  - Stat.O$AVE.app


Stat.App$Cor.app
Stat.App$RCor.app
Stat.App$reg.app
Stat.App$RMSE.app
Stat.App$AVE.app