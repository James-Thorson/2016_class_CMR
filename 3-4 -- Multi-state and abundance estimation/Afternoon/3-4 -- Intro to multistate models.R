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
juvenile_survival = 0.7
adult_survival = 0.9
juvenile_growth = 0.5
juvenile_detection = 0.7
adult_detection = 0.7

Transition_Prob = cbind( c(1,0,0), c(juvenile_survival*(1-juvenile_growth),juvenile_survival*juvenile_growth,1-juvenile_survival), c(1-adult_survival,0,adult_survival))
Detection_Prob = cbind( c(1,0,0), c(1-juvenile_detection,juvenile_detection,0), c(1-adult_detection,0,adult_detection))
rownames(Transition_Prob) = colnames(Transition_Prob) = c("Dead", "Juvenile", "Adult")
rownames(Detection_Prob) = colnames(Detection_Prob) = c("Unobserved", "Juvenile", "Adult")

# Define function to simulate a capture-history (CH) matrix
simulate_multistate_cjs <- function(Transition_Prob, Detection_Prob, n_marked){
  n_occasions <- length(n_marked) + 1
  d_its = z_its = array(0, dim=c(sum(n_marked),n_occasions,nrow(Transition_Prob)) )
  d_it = d_its[,,1]
  z_it = array(NA, dim=dim(d_it))
  # Define a vector with the occasion of marking
  mark_occasion <- rep(1:length(n_marked), times=n_marked)
  # Fill the CH matrix
  for (i in 1:dim(d_its)[1]){
    z_its[i,mark_occasion[i],] = rmultinom(n=1, size=1, prob=c(0,rep(1,dim(z_its)[3]-1)))[,1]      # Write an 1 at the release occasion
    z_it[i,mark_occasion[i]] = which(z_its[i,mark_occasion[i],]==1)
    d_its[i,mark_occasion[i],] = z_its[i,mark_occasion[i],]
    d_it[i,mark_occasion[i]] = which(d_its[i,mark_occasion[i],]==1)
    if( mark_occasion[i]!=n_occasions ){
      for(t in (mark_occasion[i]+1):n_occasions){
        # Bernoulli trial: does individual survive occasion?
        z_its[i,t,] = rmultinom(1, size=1, prob=Transition_Prob %*% z_its[i,t-1,])
        z_it[i,t] = which(z_its[i,t,]==1)
        # Bernoulli trial: is individual recaptured?
        d_its[i,t,] = rmultinom(1, size=1, prob=Detection_Prob %*% z_its[i,t,])
        d_it[i,t] = which(d_its[i,t,]==1)
      } #t
    }
  } #i
  Return = list("d_its"=d_its, "z_its"=z_its, "d_it"=d_it, "z_it"=z_it)
  return(Return)
}

# Execute function
SimList = simulate_multistate_cjs(Transition_Prob=Transition_Prob, Detection_Prob=Detection_Prob, n_marked=n_marked)
CH = SimList$d_it
CH2 = SimList$d_its

# Create vector with occasion of marking
first_fn = function(x) min(which(x!=0))
first_mark = apply(CH, MARGIN=1, FUN=first_fn)

# Bundle data
jags_data <- list(d_it=CH, first_mark=first_mark, n_indiv=nrow(CH), n_occasions=ncol(CH), dummy=1)

# (Note that the function known.state.cjs is used in section 7.3.1 as well for another purpose)

inits <- function(){
  Return = list(z_it=SimList$z_it)
  return(Return)
}

# Parameters monitored
parameters <- c("juvenile_survival", "adult_survival", "juvenile_growth", "juvenile_detection", "adult_detection")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

# Specify model in BUGS language
Multistate_CJS = function(){

  # Priors
  juvenile_survival ~ dunif(0,1)
  adult_survival ~ dunif(0,1)
  juvenile_growth ~ dunif(0,1)
  juvenile_detection ~ dunif(0,1)
  adult_detection ~ dunif(0,1)

  # Transition matrix
  Trans_Prob[1,1] <- 1
  Trans_Prob[2,1] <- 0
  Trans_Prob[3,1] <- 0
  Trans_Prob[1,2] <- (1-juvenile_survival)
  Trans_Prob[2,2] <- juvenile_survival*(1-juvenile_growth)
  Trans_Prob[3,2] <- juvenile_survival*juvenile_growth
  Trans_Prob[1,3] <- (1-adult_survival)
  Trans_Prob[2,3] <- 0
  Trans_Prob[3,3] <- adult_survival

  # Detection matrix
  Detect_Prob[1,1] <- 1
  Detect_Prob[2,1] <- 0
  Detect_Prob[3,1] <- 0
  Detect_Prob[1,2] <- (1-juvenile_detection)
  Detect_Prob[2,2] <- juvenile_detection
  Detect_Prob[3,2] <- 0
  Detect_Prob[1,3] <- (1-adult_detection)
  Detect_Prob[2,3] <- 0
  Detect_Prob[3,3] <- adult_detection

  # Likelihood
  Beta[1] <- 0
  Beta[2] ~ dunif(0,1)
  Beta[3] <- 1-Beta[2]
  for (i in 1:n_indiv){
    # Define latent state at first capture
    z_it[i,first_mark[i]] ~ dcat( Beta[1:3] )
    for(t in (first_mark[i]+1):n_occasions){
      # State process
      zstar_its[i,t,1:3] = Trans_Prob[1:3,z_it[i,t-1]]
      z_it[i,t] ~ dcat( zstar_its[i,t,1:3] )
      # Observation process
      dstar_its[i,t,1:3] = Detect_Prob[1:3,z_it[i,t]]
      d_it[i,t] ~ dcat( dstar_its[i,t,1:3] )
    } #t
  } #i

  dummy ~ dnorm( Trans_Prob[1,3], 1 )
}

# Call JAGS from R (BRT 1 min)  # inits=inits,
Jags = jags(data=jags_data, parameters=parameters, inits=inits, model.file=Multistate_CJS, n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory=getwd())

# Summarize posteriors
print(Jags, digits = 3)



