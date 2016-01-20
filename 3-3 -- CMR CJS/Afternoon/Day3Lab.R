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
n.occasions <- 6                   # Number of capture occasions
marked <- rep(50, n.occasions-1)   # Annual number of newly marked individuals
phi <- rep(0.65, n.occasions-1)
p <- rep(0.4, n.occasions-1)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked))
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Define function to simulate a capture-history (CH) matrix
simul.cjs <- function(PHI, P, marked){
   n.occasions <- dim(PHI)[2] + 1
   CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))
   # Define a vector with the occasion of marking
   mark.occ <- rep(1:length(marked), marked[1:length(marked)])
   # Fill the CH matrix
   for (i in 1:sum(marked)){
      CH[i, mark.occ[i]] <- 1       # Write an 1 at the release occasion
      if (mark.occ[i]==n.occasions) next
      for (t in (mark.occ[i]+1):n.occasions){
         # Bernoulli trial: does individual survive occasion?
         sur <- rbinom(1, 1, PHI[i,t-1])
         if (sur==0) break		# If dead, move to next individual 
         # Bernoulli trial: is individual recaptured? 
         rp <- rbinom(1, 1, P[i,t-1])
         if (rp==1) CH[i,t] <- 1
         } #t
      } #i
   return(CH)
   }

# Execute function
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in BUGS language
sink("cjs-c-c.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- mean.phi
      p[i,t] <- mean.p
      } #t
   } #i

mean.phi ~ dunif(0, 1)         # Prior for mean survival
mean.p ~ dunif(0, 1)           # Prior for mean recapture

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2])


# Initial values
# In JAGS we have to give good initial values for the latent state z. At all occasions when an individual was observed, its state is z = 1 for sure. In addition, if an individual was not observed at an occasion, but was alive for sure, because it was observed before and thereafter (i.e. has a capture history of e.g. {101} or {10001}), then we know that the individual was alive at all of these occasions, and thus z = 1. Therefore, we should provide initial values of z = 1 at these positions as well. The following function provides such initial values from the observed capture histories:
known.state.cjs <- function(ch){
   state <- ch
   for (i in 1:dim(ch)[1]){
      n1 <- min(which(ch[i,]==1))
      n2 <- max(which(ch[i,]==1))
      state[i,n1:n2] <- 1
      state[i,n1] <- NA
      }
   state[state==0] <- NA
   return(state)
   }

# (Note that the function known.state.cjs is used in section 7.3.1 as well for another purpose) 
library(R2jags)
inits <- function(){list(mean.phi = runif(1, 0, 1), mean.p = runif(1, 0, 1), z = known.state.cjs(CH))}

# Parameters monitored
parameters <- c("mean.phi", "mean.p")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 1 min)
cjs.c.c <- jags(jags.data, inits, parameters, "cjs-c-c.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(cjs.c.c, digits = 3)

#######################
# Random time effects in CJS model
######################

# Define parameter values
n.occasions <-   6                # Number of capture occasions
marked <- rep(30, n.occasions-1)   # Annual number of newly marked individuals
mean.phi <- 0.65
var.phi <- 1                       # Temporal variance of survival
p <- rep(0.4, n.occasions-1)

# Determine annual survival probabilities
mu.phi<- qlogis(mean.phi)
logit.phi <- rnorm(n.occasions-1, mu.phi , var.phi^0.5)
phi <- plogis(logit.phi)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Simulate capture-histories
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Specify model in BUGS language
sink("cjs-temp-raneff.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      logit(phi[i,t]) <- mu + epsilon[t]
      p[i,t] <- mean.p
      } #t
   } #i
for (t in 1:(n.occasions-1)){
   epsilon[t] ~ dnorm(0, tau)
   }

#mu ~ dnorm(0, 0.1)                      # Prior for logit of mean survival
#mean.phi <- 1 / (1+exp(-mu))            # Logit transformation
mean.phi ~ dunif(0, 1)                   # Prior for mean survival
mu <- log(mean.phi / (1-mean.phi))       # Logit transformation
sigma ~ dunif(0, 10)                     # Prior for standard deviation
tau <- pow(sigma, -2)
sigma2 <- pow(sigma, 2)                  # Temporal variance
mean.p ~ dunif(0, 1)                     # Prior for mean recapture

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH))

# Function to create a matrix of initial values for latent state z
cjs.init.z <- function(ch,f){
   for (i in 1:dim(ch)[1]){
      if (sum(ch[i,])==1) next
      n2 <- max(which(ch[i,]==1))
      ch[i,f[i]:n2] <- NA
      }
   for (i in 1:dim(ch)[1]){
   ch[i,1:f[i]] <- NA
   }
   return(ch)
   }


# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), mean.phi = runif(1, 0, 1), sigma = runif(1, 0, 10), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.p", "sigma2")

# MCMC settings
ni <- 10000
nt <- 6
nb <- 5000
nc <- 3

# Call JAGS from R (BRT 4 min)
cjs.ran <- jags(jags.data, inits, parameters, "cjs-temp-raneff.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(cjs.ran, digits = 3)

# Produce histogram
hist(cjs.ran$BUGSoutput$sims.list$sigma2, col = "gray", las = 1, xlab = expression(sigma^2), main = "", breaks = 100)
abline(v = var.phi, col = "red", lwd = 2)

###################
#Model with temporal covariates in CJS and random effects
####################

n.occasions <- 8                  # Number of capture occasions
marked <- rep(15, n.occasions-1)   # Annual number of newly marked individuals
mean.phi <- 0.65
p <- rep(0.4, n.occasions-1)
beta <- -0.3                       # Slope of survival-winter relationship	
r.var <- 0.2                       # Residual temporal variance

# Draw annual survival probabilities
winter <- rnorm(n.occasions-1, 0, 1^0.5)
logit.phi <- qlogis(mean.phi) + beta*winter + rnorm(n.occasions-1, 0, r.var^0.5)
phi <- plogis(logit.phi)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Simulate capture-histories
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)


# Specify model in BUGS language
sink("cjs-cov-raneff.jags")
cat("
model {

# Priors and constraints
for (i in 1:nind){
   for (t in f[i]:(n.occasions-1)){
      logit(phi[i,t]) <- mu + beta*x[t] + epsilon[t]
      p[i,t] <- mean.p
      } #t
   } #i
for (t in 1:(n.occasions-1)){
   epsilon[t] ~ dnorm(0, tau)
   phi.est[t] <- 1 / (1+exp(-mu-beta*x[t]-epsilon[t])) # Yearly survival
   }
mu ~ dnorm(0, 0.001)                     # Prior for logit of mean survival
mean.phi <- 1 / (1+exp(-mu))             # Logit transformation
beta ~ dnorm(0, 0.001)I(-10, 10)         # Prior for slope parameter
sigma ~ dunif(0, 10)                     # Prior on standard deviation
tau <- pow(sigma, -2)
sigma2 <- pow(sigma, 2)                  # Residual temporal variance
mean.p ~ dunif(0, 1)                     # Prior for mean recapture

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- 1
   for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
      } #t
   } #i
}
",fill = TRUE)
sink()

# Bundle data
jags.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH), x = winter)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), mu = rnorm(1), sigma = runif(1, 0, 5), beta = runif(1, -5, 5), mean.p = runif(1, 0, 1))}  

# Parameters monitored
parameters <- c("mean.phi", "mean.p", "phi.est", "sigma2", "beta")

# MCMC settings
ni <- 20000
nt <- 6
nb <- 10000
nc <- 3

# Call JAGS from R (BRT 12 min)
cjs.cov <- jags(jags.data, inits, parameters, "cjs-cov-raneff.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(cjs.cov, digits = 3)

# Produce graph
par(mfrow = c(1, 2), las = 1)
hist(cjs.cov$BUGSoutput$sims.list$beta, nclass = 25, col = "gray", main = "", xlab = expression(beta), ylab = "Frequency")
abline(v = -0.3, col = "red", lwd = 2)
hist(cjs.cov$BUGSoutput$sims.list$sigma2, nclass = 50, col = "gray", main = "", xlab = expression(sigma^2), ylab = "Frequency", xlim=c(0, 3))
abline(v = 0.2, col = "red", lwd = 2)



