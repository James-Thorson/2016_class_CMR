
library(R2jags)

# Function to generate data
Sim_Fn = function( n_per_group, n_groups=10, logSD=1, logmean=2, beta_u=1){
  # Simulate measured and unmeasured random components
  u_s = rnorm(n_groups)
  r_s = rnorm(n_groups, mean=logmean, sd=logSD)

  # Simulate samples
  s_i = rep( 1:n_groups, each=n_per_group)
  mean_s = exp( beta_u*u_s + r_s )
  c_i = rpois( n=n_per_group*n_groups, lambda=mean_s[s_i] )

  # Bundle and return
  DF = data.frame( "site"=s_i, "count"=c_i, "covariate"=u_s[s_i])
  return( DF )
}

# Define JAGS model
LMM_Jags_1 = function(){
  # Prior distributions
  log_Mean ~ dunif(0,20)
  Sd_Site ~ dunif(0,10)
  Beta ~ dunif(-10,10)
  # Derived parameters
  Tau_Site <- pow(Sd_Site,-2)
  # Probability of random effects
  for(i in 1:Nsite){
    log_Mean_s[i] ~ dnorm( log_Mean, Tau_Site)
  }
  # Probability of data
  for(i in 1:Nobs){
    Obspred_i[i] <- exp( log_Mean_s[Site_i[i]] + Beta*Covariate_i[i] )
    Obs_i[i] ~ dpois( Obspred_i[i] )
  }
}
LMM_Jags_2 = function(){
  # Prior distributions
  log_Mean ~ dunif(0,20)
  #Sd_Site ~ dunif(0,10)
  Beta ~ dunif(-10,10)
  # Derived parameters
  #Tau_Site <- pow(Sd_Site,-2)
  # Probability of random effects
  #for(i in 1:Nsite){
  #  log_Mean_s[i] ~ dnorm( log_Mean, Tau_Site)
  #}
  # Probability of data
  for(i in 1:Nobs){
    #Obspred_i[i] <- exp( log_Mean_s[Site_i[i]] + Beta*Covariate_i[i] )
    Obspred_i[i] <- exp(log_Mean + Beta*Covariate_i[i])
    Obs_i[i] ~ dpois( Obspred_i[i] )
  }
}

Results = array(NA, dim=c(5,2,3))
for( repI in 1:5){
  # Illustrate function
  SimData = Sim_Fn( n_per_group=10, n_groups=10, beta_u=1)

  # Generate inputs for JAGS
  Nsim = Nburnin = 1e4
  # Run Model 1
  Data = list(Obs_i=SimData$count, Nobs=nrow(SimData), Site_i=SimData$site, Nsite=length(unique(SimData$site)), Covariate_i=SimData$covariate)
  Jags_1 <- jags(model.file=LMM_Jags_1, working.directory=NULL, data=Data, parameters.to.save=c("log_Mean","Sd_Site","Beta"), n.chains=3, n.thin=1, n.iter=Nsim+Nburnin, n.burnin=Nburnin)
  Results[repI,1,] = Jags_1$BUGSoutput$summary['Beta',c('mean','25%','75%')]
  # Run Model 1
  Data = list(Obs_i=SimData$count, Nobs=nrow(SimData), Covariate_i=SimData$covariate)
  Jags_2 <- jags(model.file=LMM_Jags_2, working.directory=NULL, data=Data, parameters.to.save=c("log_Mean","Beta"), n.chains=3, n.thin=1, n.iter=Nsim+Nburnin, n.burnin=Nburnin)
  Results[repI,2,] = Jags_2$BUGSoutput$summary['Beta',c('mean','25%','75%')]
}
