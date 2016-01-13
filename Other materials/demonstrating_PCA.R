
# Simulate data set
library(mvtnorm)
V = matrix(c(2,1,1,1), 2,2)
Y = rmvnorm( n=1000, sigma=V)

# Built-in PCA
PCA = princomp(Y)
eigenscores = PCA$scores

# By-hand PCA
eigenvec = eigen( cov(Y) )$vectors
eigenscores = t(solve(eigenvec) %*% t(Y))

##################
# Explore fitting in JAGS
##################
library(R2jags)

# Define JAGS model
PCA = function(){
  # Prior distributions
  Sigma ~ dwish( R, k )
  for(p in 1:Nvar){
    Mean[p] ~ dunif(-10,10)
  }
  # Sampling declarations
  for(i in 1:Nobs){
    Y[i,] ~ dmnorm( Mean, Sigma )
  }
}

# Generate inputs for JAGS
Nsim = Nburnin = 1e3
Data = list(Y=Y, R=diag(2), k=2, Nobs=nrow(Y), Nvar=ncol(Y))
# Run jags
Jags <- jags(model.file=PCA, working.directory=NULL, data=Data, parameters.to.save=c("Mean","Sigma"), n.chains=3, n.thin=1, n.iter=Nsim+Nburnin, n.burnin=Nburnin)

