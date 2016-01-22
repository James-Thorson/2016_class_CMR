############################
#  Inference in R and JAGS
#  Concepcion, Chile
#  22 enero 2016
#  Prof: Noble Hendrix & Jim Thorson
#--------------------------
#  Day 5 - Tuna State Space Model
#--------------------------

#tuna data as presented in Meyer and Millar CJFAS 65:1078-1086
C<-c(15.9,25.7,28.5,23.7,25.0,33.3,28.2,19.7,17.5,19.3,21.6,23.1,22.5,22.5,23.6,29.1,14.4,13.2,28.4,34.6,37.5,25.9,25.3)

I<-c(61.89,78.98,55.59,44.61,56.89,38.27,33.84,36.13,41.95,36.63,36.33,38.82,34.32,37.64,34.01,32.16,26.88,36.61,30.07,30.75,23.36,22.36,21.91)
the.years<-1967:1989
par(mfrow = c(1,2))
plot(the.years, C, ylab = "Catch (1000 t)", xlab = "Year", type = 'l', lwd = 3)
plot(the.years, I, ylab = "CPUE (kg/100 hooks)", xlab = "Year", type = 'l', lwd = 3)

#Run model in JAGS
#example of using rjags rather than R2jags 
library(rjags)

N<-23
tuna.data<-list("N"=N,"C"=C,"I"=I)
tuna.inits<-function()
list(K = runif(1,10,1000) ,r = runif(1,0.01, 1.2), P = runif(N,0.5,1), iq = runif(1,10,100), isigma2 = runif(1,80,150 ), itau2 = runif(1,90, 120))
tuna.params<-c("r", "K", "MSP", "EMSP", "B1990", "sigma2", "tau2", "I.new") 
tuna.jags<-jags.model(file="tuna.bug", data=tuna.data, inits=tuna.inits, n.chains=3, n.adapt=1000)

#samplers being used for each Gibbs step
list.samplers(tuna.jags)  

#burn-in
update(tuna.jags, n.iter=10000) 
#samples
tuna.sim<-coda.samples(tuna.jags, variable.names=tuna.params, n.iter=50000, thin = 5)

summary(tuna.sim)
gelman.diag(tuna.sim, multivariate = F)

#update another 50K
tuna.sim<-coda.samples(tuna.jags, variable.names=tuna.params, n.iter=50000, thin = 5)
gelman.diag(tuna.sim, multivariate = F)


I.ind<- which(varnames(tuna.sim)=="I.new[1]")

#Observed versus modeled
par(mfrow = c(1,1))
plot(the.years, I, ylab = "CPUE (kg/100 hooks)", xlab = "Year", pch = 15)
lines(the.years, summary(tuna.sim)[[2]][I.ind:(I.ind + 22), 3], lwd = 2)
lines(the.years, summary(tuna.sim)[[2]][I.ind:(I.ind + 22), 1], lwd = 2, lty = 2)
lines(the.years, summary(tuna.sim)[[2]][I.ind:(I.ind + 22), 5], lwd = 2, lty = 2)

