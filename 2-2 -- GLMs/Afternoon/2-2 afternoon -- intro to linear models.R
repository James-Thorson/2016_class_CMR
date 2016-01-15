
setwd( "C:/Users/James.Thorson/Desktop/Project_git/2016_classes_private/CMR models/2-2 -- GLMs/Afternoon/" )

###### Example GLm code

# simulate example data
TrueMean = 3
TrueSize = 1
Counts = rnbinom(20, mu=TrueMean, size=TrueSize) # Var = mu + mu^2/size
Covariate = rnorm( length(Counts), mean=0, sd=1)
#Counts = rpois(10000, lambda=TrueMean) # Var = mu + mu^2/size

# Plot example data
png(file="Lab_1_Counts.png", width=4, height=4, res=200, units="in")
  par(mar=c(3,3,2,0), mgp=c(1.5,0.25,0))
  hist(Counts)
dev.off()

#### Method 1 -- Nonlinear optimization
NegLogLike_Fn = function(Par, Data){
  # Parameters
  Mean_hat = Par[1]
  Size_hat = Par[2]
  Slope_hat = Par[3]
  # Log-likelihood
  LogLike_i = dnbinom( Data$Counts, mu=exp(Mean_hat + Data$Covariate*Slope_hat), size=Size_hat, log=TRUE )
  NegLogLike = -1 * sum(LogLike_i)
  return( NegLogLike )
}
# Run model
Data = list( 'Counts'=Counts, 'Covariate'=Covariate )
Start = c(1,1, 1)
NegLogLike_Fn( Par=Start, Data=Data)
Opt = optim( par=Start, fn=NegLogLike_Fn, Data=Data, lower=c(0.01,0.01,-Inf), upper=Inf, method="L-BFGS-B", hessian=TRUE )
# Estimated parameters
print( Opt$par ) # Estimated parameters
# Estimated standard errors
print( sqrt(diag( solve(Opt$hessian) )) ) # square root of diagonal elements of the inverse-hessian matrix

#### Method 2 -- GLM function in R
library(MASS)
Glm = glm.nb( Counts ~ 1 )
# Estimated parameters and standard errors
summary(Glm)
  
#### Method 3 -- GLM in JAGS
library(R2jags)

# Define JAGS model
NegBin = function(){
  # Prior distributions
  Size ~ dunif(0.001,10)
  Ln_Mean ~ dunif(0.001,10)
  # Derived quantities
  VarInf <- exp(Ln_Mean) / Size
  Mean <- exp(Ln_Mean)
  # Change to JAGS parameterization
  p <- 1 / (VarInf + 1)     # WIKIPEDIA: p <- 1 / (1-VarInf)
  r <- Mean * p / (1-p)
  # Sampling declarations
  for(i in 1:Nobs){
    Counts[i] ~ dnegbin(p,r)
  }
}

# Generate inputs for JAGS
Nsim = Nburnin = 5e2
Data = list(Counts=Counts, Nobs=length(Counts))
# Run jags
Jags <- jags(model.file=NegBin, working.directory=NULL, data=Data, parameters.to.save=c("Ln_Mean","Size"), n.chains=3, n.thin=1, n.iter=Nsim+Nburnin, n.burnin=Nburnin)
# Look at estimates
Jags$BUGSoutput$summary

####### Convergence diagnostics
# Look at effective sample sizes
Jags$BUGSoutput$summary
# Plot MCMC trace
png(file="Lab_1_Trace.png", width=6, height=4, res=200, units="in")
  par(mar=c(3,3,2,0), mgp=c(1.5,0.25,0), xaxs="i", yaxs="i")
  traceplot(Jags, mfrow=c(1,3), ask=FALSE)
dev.off()

####### Check goodness-of-fit
png(file="Lab_1_Posterior_predictive_check.png", width=6, height=4, res=200, units="in")
  par(mfrow=c(1,2), mar=c(3,3,2,0), mgp=c(2,0.5,0))
  hist( Jags$BUGSoutput$sims.list$Pred, freq=FALSE, xlim=c(0,max(Data$Counts)), breaks=seq(0,100,by=1), main="Predictive distribution", xlab="Counts" ) 
  hist( Data$Counts, freq=FALSE, xlim=c(0,max(Data$Counts)), breaks=seq(0,100,by=1), main="Available data", xlab="Counts" ) 
dev.off()

####### Interpret results
# Plot posteriors
png(file="Lab_1_Posteriors.png", width=6, height=4, res=200, units="in")
  par(mfrow=c(1,2), mar=c(3,3,2,0), mgp=c(1.5,0.25,0), xaxs="i", yaxs="i")
  hist(exp(Jags$BUGSoutput$sims.list[["Ln_Mean"]]), breaks=25, main="Mean", xlab="Value", ylab="Density")
    abline( v=TrueMean, col="red"  )
  hist(Jags$BUGSoutput$sims.list[["Size"]], breaks=25, main="Size", xlab="Value", ylab="Density")
    abline( v=TrueSize, col="red" )
dev.off()




  
