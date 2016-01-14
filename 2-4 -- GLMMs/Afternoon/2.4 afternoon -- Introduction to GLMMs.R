
############
# Example 1 - Linear mixed model
############
library(R2jags)
library(lme4)

###### Simulate data
# Parameters
Nsite = 10
Nobs_per_site = 10
Site_Mean = 10
Site_Sd = 2
Obs_Sd = 2

# Bookkeeping
Site_i = rep( 1:Nsite, each=Nobs_per_site)

# Simulation
Mean_s = rnorm( Nsite, mean=Site_Mean, sd=Site_Sd )
Obs_i = rnorm( Nsite*Nobs_per_site, mean=Mean_s[Site_i], sd=Obs_Sd )

# Plot data
library(lattice)
#histogram( ~ Obs_i | factor(Site_i), breaks=seq( min(Obs_i), max(Obs_i), length=10) )      #
histogram( ~ Obs_i | factor(Site_i), breaks=seq( min(Obs_i), max(Obs_i), length=10), type="density", panel=function(x,...){ panel.histogram(x, ...); panel.mathdensity(dmath=dnorm, col="black", args = list(mean=mean(x),sd=sd(x))) } )      #

# Plot theoretical distributions
X = seq(0,20, length=1000)
hist(Obs_i, freq=FALSE, breaks=40)
lines( x=X, y=dnorm(X, mean=Site_Mean, sd=Site_Sd), col="red")
for(i in 1:Nsite){
  Y = dnorm(X, mean=Mean_s[i], sd=Obs_Sd)
  lines( x=X, y=Y, col="blue")
}

###### Fit using R
# No site level (Not recommended)
LM = lm( Obs_i ~ 1 )
print( summary(LM) )

# Using fixed effects (Not recommended)
LM = lm( Obs_i ~ 0 + factor(Site_i) )
print( summary(LM) )

# Using mixed effects (Recommended)
LMM = lmer( Obs_i ~ 1 + (1 | factor(Site_i)) )
print( summary(LMM) )

# jim's attempt at replicating this from first principles
Mu = mean(Obs_i)
  Mu_s = tapply( Obs_i, INDEX=Site_i, FUN=mean)
Sigma = sd( Mu_s )
  Sigma_s = sd( Obs_i - Mu_s[Site_i] )
Weights_hat = c( 1/Sigma^2, Nobs_per_site/Sigma_s^2 )
  Weights_hat = Weights_hat / sum(Weights_hat)
# Predictions
Mu_s_hat = ( Mu*Weights_hat[1] + Mu_s*Weights_hat[2] )
cbind( Mu_s_hat, ranef(LMM)[['factor(Site_i)']]+fixef(LMM)['(Intercept)'] )

###### Fit using JAGS
# Define JAGS model
LMM_Jags = function(){
  # Prior distributions
  Mean ~ dunif(0,20)
  Sd_Site ~ dunif(0,10)
  Sd_Obs ~ dunif(0,10)
  # Derived parameters
  Tau_Site <- pow(Sd_Site,-2)
  Tau_Obs <- pow(Sd_Obs,-2)
  # Probability of random effects
  for(i in 1:Nsite){
    Mean_s[i] ~ dnorm( Mean, Tau_Site)
  }
  # Probability of data
  for(i in 1:Nobs){
    Obs_i[i] ~ dnorm(Mean_s[Site_i[i]], Tau_Obs)
  }
}

# Generate inputs for JAGS
Nsim = Nburnin = 1e3
Data = list(Obs_i=Obs_i, Nobs=length(Obs_i), Site_i=Site_i, Nsite=Nsite)
# Run jags
Jags <- jags(model.file=LMM_Jags, working.directory=NULL, data=Data, parameters.to.save=c("Mean","Sd_Site","Sd_Obs"), n.chains=3, n.thin=1, n.iter=Nsim+Nburnin, n.burnin=Nburnin)
# Look at estimates
print( Jags$BUGSoutput$summary )
# Plot MCMC trace
traceplot(Jags, mfrow=c(2,2), ask=FALSE)
# Plot posteriors
par(mfrow=c(2,2), mar=c(3,3,2,0), mgp=c(1.5,0.25,0), xaxs="i", yaxs="i")
  for(i in 1:4) hist(Jags$BUGSoutput$sims.matrix[,i], breaks=25, main=colnames(Jags$BUGSoutput$sims.matrix)[i], xlab="Value", ylab="Density")


############
# Example 2 -- Generalized linear mixed model
############

###### Simulate data
# Parameters
Nsite = 10
Nobs_per_site = 10
Site_logMean = log(10)
Site_logSd = 1

# Bookkeeping
Site_i = rep( 1:Nsite, each=Nobs_per_site)

# Simulation
Mean_s = rlnorm( Nsite, meanlog=Site_logMean, sdlog=Site_logSd )
Obs_i = rpois( Nsite*Nobs_per_site, lambda=Mean_s[Site_i] )

# Plot data
library(lattice)
histogram( ~ Obs_i | factor(Site_i), breaks=seq( min(Obs_i), max(Obs_i), length=10), type="density", panel=function(x,...){ panel.histogram(x, ...); panel.mathdensity(dmath=dnorm, col="black", args = list(mean=mean(x),sd=sd(x))) } )      #

###### Fit using R
# No site level (Not recommended)
GLM = glm( Obs_i ~ 0, family="poisson" )
print( summary(GLM) )

# Using fixed effects (Not recommended)
GLM = glm( Obs_i ~ 0 + factor(Site_i), family="poisson" )
print( summary(GLM) )

# Using mixed effects (Recommended)
library(lme4)
GLMM = glmer( Obs_i ~ 1 + (1 | factor(Site_i)), family="poisson" )
print( summary(GLMM) )

###### Fit using JAGS
library(R2jags)
# Define JAGS model
LMM_Jags = function(){
  # Prior distributions
  logMean ~ dunif(0,20)
  logSd_Site ~ dunif(0,10)
  # Derived parameters
  logTau_Site <- pow(logSd_Site,-2)
  # Probability of random effects
  for(i in 1:Nsite){
    logMean_s[i] ~ dnorm( logMean, logTau_Site)
  }
  # Derived statistics
  for(i in 1:Nsite){
    Mean_s[i] <- exp(logMean_s[i])
  }
  # Probability of data
  for(i in 1:Nobs){
    Obs_i[i] ~ dpois(Mean_s[Site_i[i]])
  }
}

# Generate inputs for JAGS
Nsim = Nburnin = 1e3
Data = list(Obs_i=Obs_i, Nobs=length(Obs_i), Site_i=Site_i, Nsite=Nsite)
# Run jags
Jags <- jags(model.file=LMM_Jags, working.directory=NULL, data=Data, parameters.to.save=c("logMean","logSd_Site"), n.chains=3, n.thin=1, n.iter=Nsim+Nburnin, n.burnin=Nburnin)
# Look at estimates
print( Jags$BUGSoutput$summary )
# Plot MCMC trace
traceplot(Jags, mfrow=c(2,2), ask=FALSE)
# Plot posteriors
par(mfrow=c(2,2), mar=c(3,3,2,0), mgp=c(1.5,0.25,0), xaxs="i", yaxs="i")
  for(i in 1:3) hist(Jags$BUGSoutput$sims.matrix[,i], breaks=25, main=colnames(Jags$BUGSoutput$sims.matrix)[i], xlab="Value", ylab="Density")

