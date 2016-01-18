
library(R2jags)

############
# Example -- State-space model
############

###### Simulate data
# Parameters
Nobs = 100
Obs_Sd = 1
Process_Sd = 1
Mu = 0.2
N1 = 10

# Simulation
N_t = rep(NA, Nobs)
N_t[1] = N1
for(t in 2:Nobs) N_t[t] = rnorm(1, mean=N_t[t-1]+Mu, sd=Process_Sd)
B_t = rnorm(Nobs, mean=N_t, sd=Obs_Sd)

# Plot simulated data
plot(B_t)
points( N_t, col="blue")
lines( N1 + Mu*(1:Nobs-1), col="red" )

###### Fit using JAGS
# Define JAGS model
StateSpace_Jags = function(){
  # Prior distributions
  Obs_Sd ~ dunif(0,20)
  Process_Sd ~ dunif(0,10)
  Mu ~ dunif(-10,10)
  N_t[1] ~ dunif(0,100)
  # Derived parameters
  Obs_Tau <- pow(Obs_Sd,-2)
  Process_Tau <- pow(Process_Sd,-2)
  # Probability of random effects
  for(t in 2:Nobs){
    N_t_hat[t] <- N_t[t-1] + Mu
    N_t[t] ~ dnorm( N_t_hat[t], Process_Tau )
  }
  # Probability of data
  for(t in 1:Nobs){
    B_t[t] ~ dnorm( N_t[t], Obs_Tau )
  }
}

# Generate inputs for JAGS
Nsim = Nburnin = 5e3
Data = list(Nobs=length(B_t), B_t=B_t)
# Run jags
Jags <- jags(model.file=StateSpace_Jags, working.directory=NULL, data=Data, parameters.to.save=c("Obs_Sd","Process_Sd","Mu","N_t"), n.chains=3, n.thin=1, n.iter=Nsim+Nburnin, n.burnin=Nburnin)
# Look at estimates
print( Jags$BUGSoutput$summary )
# Plot MCMC trace
par(mar=c(2,2,0,0))
traceplot(Jags, mfrow=c(5,5), ask=TRUE)
