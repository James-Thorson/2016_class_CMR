
library(R2jags)

############
# Version 1 -- Original model
############

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

Data = read.csv( paste0("C:/Users/James.Thorson/Desktop/Project_git/2016_class_CMR/3-1 -- Review and state-space models/Exercise/state_space_count_data.csv"))

# Generate inputs for JAGS
Nsim = Nburnin = 1e3
Data = list(Nyears=max(Data$Year)-min(Data$Year)+1, Nobs=nrow(Data), B_i=Data$Simulated_counts, Year_i=Data$Year-min(Data$Year)+1)
# Run jags
Jags <- jags(model.file=StateSpace_Jags, working.directory=NULL, data=Data, parameters.to.save=c("Obs_Sd","Process_Sd","Mu","N_t"), n.chains=3, n.thin=1, n.iter=Nsim+Nburnin, n.burnin=Nburnin)

############
# Version 2 -- change indexing so that we can use multiple observations in each year
############

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
  for(t in 2:Nyears){
    N_t_hat[t] <- N_t[t-1] + Mu
    N_t[t] ~ dnorm( N_t_hat[t], Process_Tau )
  }
  # Probability of data
  for(i in 1:Nobs){
    B_i[i] ~ dnorm( N_t[ Year_i[i] ], Obs_Tau )
  }
}

Data_Orig = read.csv( paste0("C:/Users/James.Thorson/Desktop/Project_git/2016_class_CMR/3-1 -- Review and state-space models/Exercise/state_space_count_data.csv"))

# Generate inputs for JAGS
Nsim = Nburnin = 1e3
Data = list("Nyears"=max(Data_Orig$Year)-min(Data_Orig$Year)+1,
     "Nobs"=nrow(Data_Orig),
     "B_i"=Data_Orig$Simulated_counts,
     "Year_i"=Data_Orig$Year-min(Data_Orig$Year)+1 )

# Run jags
Jags <- jags(model.file=StateSpace_Jags, working.directory=NULL, data=Data, parameters.to.save=c("Obs_Sd","Process_Sd","Mu","N_t"), n.chains=3, n.thin=1, n.iter=Nsim+Nburnin, n.burnin=Nburnin)

############
# Version 3 -- Change to Gompertz production funcction and Poisson distribution for data
############

# Define JAGS model
StateSpace_Jags = function(){
  # Prior distributions
  Process_Sd ~ dunif(0,10)
  Mu ~ dunif(-10,10)
  N_t[1] ~ dunif(0,100)
  rho ~ dunif(-1,1)
  # Derived parameters
  Process_Tau <- pow(Process_Sd,-2)
  Carrying_capacity <- Mu / (1-rho)
  # Probability of random effects
  for(t in 2:Nyears){
    log_N_t_hat[t] <- rho*log(N_t[t-1]) + Mu
    N_t[t] ~ dlnorm( log_N_t_hat[t], Process_Tau )
  }
  # Probability of data
  for(i in 1:Nobs){
    B_i[i] ~ dpois( N_t[Year_i[i]] )
  }
}

Data = read.csv( paste0("C:/Users/James.Thorson/Desktop/Project_git/2016_class_CMR/3-1 -- Review and state-space models/Exercise/state_space_count_data.csv"))

# Generate inputs for JAGS
Nsim = Nburnin = 5e3
Data = list(Nyears=max(Data$Year)-min(Data$Year)+1, Nobs=nrow(Data), B_i=Data$Simulated_counts, Year_i=Data$Year-min(Data$Year)+1)
# Run jags
Jags <- jags(model.file=StateSpace_Jags,
working.directory=NULL, data=Data,
parameters.to.save=c("Carrying_capacity","rho","Process_Sd","Mu","N_t"),
n.chains=3, n.thin=1, n.iter=Nsim+Nburnin, n.burnin=Nburnin)


