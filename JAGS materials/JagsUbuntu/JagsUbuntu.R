############################
#Example of running JAGS using rjags - can be used in Ubuntu
#--------------------------

#tuna data as presented in Meyer and Millar CJFAS 65:1078-1086
C<-c(15.9,25.7,28.5,23.7,25.0,33.3,28.2,19.7,17.5,19.3,21.6,23.1,22.5,22.5,23.6,29.1,14.4,13.2,28.4,34.6,37.5,25.9,25.3)

I<-c(61.89,78.98,55.59,44.61,56.89,38.27,33.84,36.13,41.95,36.63,36.33,38.82,34.32,37.64,34.01,32.16,26.88,36.61,30.07,30.75,23.36,22.36,21.91)
the.years<-1967:1989
par(mfrow = c(1,2))
plot(the.years, C, ylab = "Catch (1000 t)", xlab = "Year", type = 'l', lwd = 3)
plot(the.years, I, ylab = "CPUE (kg/100 hooks)", xlab = "Year", type = 'l', lwd = 3)

#Construct model in JAGS
#example of using rjags rather than R2jags 
library(rjags)

#write the model within R
sink("tuna.mod.jags")
cat("
model 
{

#time step [1] conditions
Pmed[1] <-0
P[1]~dlnorm(Pmed[1], isigma2)T(0.05,1.6)
 
#time steps of model  
for( t in 2 : N )
 	{
      Pmed[t] <- log(max(P[t - 1] + (r * P[t - 1]) * (1 - P[t - 1]) - C[t - 1] / K, 0.001) )	
	 P[t] ~ dlnorm(Pmed[t],isigma2)T(0.05,1.5)
 	}
for( t in 1 : N )
	{
	Imed[t] <- log((q * K) * P[t])
 	I[t] ~ dlnorm(Imed[t],itau2)
	
 	#posterior predictions
 	index[t]<- log(q*K*P[t])
	I.new[t]~dlnorm(index[t], itau2)
	}

 #priors 
 r ~ dlnorm( -1.38, 3.845)I(0.01,1.2)
 isigma2 ~ dgamma(3.785,0.0102)
 sigma2 <- 1/isigma2
 itau2 ~ dgamma(1.709,0.00861)
 tau2 <- 1/itau2
 iq ~ dgamma(0.001,0.001)I( 0.5,100)
  q <- 1/iq
  K ~ dlnorm(5.0429,3.7603)I(10,1000)

#additional parameters and preditions
MSP <-  r*K/4
EMSP <-  r/(2*q)
P1990 <-  P[N] + r*P[N]*(1-P[N]) - C[N]/K
B1990 <-  P1990*K
 }
",fill = TRUE)
sink()


N<-23
tuna.data<-list("N"=N,"C"=C,"I"=I)
tuna.inits<-function()
list(K = runif(1,10,1000) ,r = runif(1,0.01, 1.2), P = runif(N,0.5,1), iq = runif(1,10,100), isigma2 = runif(1,80,150 ), itau2 = runif(1,90, 120))
tuna.params<-c("r", "K", "MSP", "EMSP", "B1990", "sigma2", "tau2", "I.new") 
tuna.jags<-jags.model(file="tuna.mod.jags", data=tuna.data, inits=tuna.inits, n.chains=3, n.adapt=1000)

#samplers being used for each Gibbs step
list.samplers(tuna.jags)  

#burn-in
update(tuna.jags, n.iter=10000) 
#samples
tuna.sim<-coda.samples(tuna.jags, variable.names=tuna.params, n.iter=50000, thin = 5)

