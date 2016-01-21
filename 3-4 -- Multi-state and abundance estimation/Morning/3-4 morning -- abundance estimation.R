#############################
# Statistical Inference in R and JAGS
#  Prof: Noble Hendrix & Jim Thorson
#  noblehendrix@gmail.com
#--------------------------
#  Day 3 Lab code - Cormack-Jolly-Seber
#--------------------------


library(R2jags)
#######  OPEN POPULATION MODEL ###############

################
#Cormack-Jolly-Seber CJS
################
# Models with constant parameters
# Define parameter values
n_occasions <- 6                   # Number of capture occasions
n_marked <- rep(50, n_occasions-1)   # Annual number of newly marked individuals
survival_prob = 0.6
detection_prob = 0.4

# Define function to simulate a capture-history (CH) matrix
simul.cjs <- function(survival_prob, detection_prob, n_marked){
  n_occasions <- length(n_marked) + 1
  z_it = y_it = array(0, dim=c(sum(n_marked),n_occasions) )
  # Define a vector with the occasion of marking
  mark_occasion <- rep(1:length(n_marked), times=n_marked)
  # Fill the CH matrix
  for (i in 1:nrow(y_it)){
    z_it[i, mark_occasion[i]] = 1      # Write an 1 at the release occasion
    y_it[i, mark_occasion[i]] = 1
    if( mark_occasion[i]!=n_occasions ){
      for(t in (mark_occasion[i]+1):n_occasions){
        # Bernoulli trial: does individual survive occasion?
        z_it[i,t] = rbinom(1, size=1, prob=survival_prob*z_it[i,t-1])
        # Bernoulli trial: is individual recaptured?
        y_it[i,t] = rbinom(1, size=1, prob=detection_prob*z_it[i,t])
      } #t
    }
  } #i
  Return = list( "z_it"=z_it, "y_it"=y_it)
  return(Return)
}

# Execute function
SimList <- simul.cjs(survival_prob, detection_prob, n_marked)
y_it <- SimList$y_it

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
first_mark <- apply(y_it, 1, get.first)

# Specify model in BUGS language
CJS = function(){

  survival_prob ~ dunif(0, 1)         # Prior for mean survival
  detection_prob ~ dunif(0, 1)           # Prior for mean recapture

  # Likelihood
  for (i in 1:n_indiv){
    # Define latent state at first capture
    z_it[i,first_mark[i]] <- 1
    for (t in (first_mark[i]+1):n_occasions){
      # State process
      zstar_it[i,t] <- survival_prob * z_it[i,t-1]
      z_it[i,t] ~ dbern(zstar_it[i,t])
      # Observation process
      ystar_it[i,t] <- detection_prob * z_it[i,t]
      y_it[i,t] ~ dbern(ystar_it[i,t])
    } #t
  } #i
}

# Bundle data
jags_data <- list(y_it=SimList$y_it, first_mark=first_mark, n_indiv=nrow(SimList$y_it), n_occasions=ncol(SimList$y_it))

# Initial values
# In JAGS we have to give good initial values for the latent state z. At all occasions when an individual was observed, its state is z = 1 for sure. In addition, if an individual was not observed at an occasion, but was alive for sure, because it was observed before and thereafter (i.e. has a capture history of e.g. {101} or {10001}), then we know that the individual was alive at all of these occasions, and thus z = 1. Therefore, we should provide initial values of z = 1 at these positions as well. The following function provides such initial values from the observed capture histories:
known.state.cjs <- function(ch){
  state <- array(NA, dim=dim(ch))
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]==1))+1
    n2 <- max(which(ch[i,]==1))
    if(n2>=n1) state[i,n1:n2] = 1
  }
  return(state)
}

# (Note that the function known.state.cjs is used in section 7.3.1 as well for another purpose) 
inits <- function(){list(z_it=known.state.cjs(SimList$y_it))}

# Parameters monitored
parameters <- c("survival_prob", "detection_prob")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 1 min)
Jags <- jags(data=jags_data, inits=inits, parameters=parameters, model.file=CJS, n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory=getwd())

# Summarize posteriors
print(cjs.c.c, digits = 3)


