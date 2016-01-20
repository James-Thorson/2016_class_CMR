#############################
# Statistical Inference in R and JAGS
#  Prof: Noble Hendrix & Jim Thorson
#  19 January 2016
#--------------------------
#  Desafio 2 - solution
#--------------------------


#1. SIMULATE DATA
#
data.fn <- function(R = 250, J = 3, K = 10, psi1 = 0.4, phi.alpha = -1, phi.beta = 1.5, range.p = c(0.2, 0.4), range.gamma = c(0, 0.1)) {
# Function to simulate detection/nondetection data for dynamic site-occupancy model
# Annual variation in probabilities of patch survival, colonization and detection 

# Function arguments:
# R - Number of sites
# J - Number of replicate surveys
# K - Number of years
# psi1 - occupancy probability in first year
# range.p - bounds of uniform distribution from which annual p drawn 
# range.psi and range.gamma - same for survival and colonization probability

   # Set up some required arrays
   site <- 1:R					# Sites
   year <- 1:K					# Years
   psi <- rep(NA, K)				# Occupancy probability
   muZ <- z <- array(dim = c(R, K))	# Expected and realized occurrence
   y <- array(NA, dim = c(R, J, K))	# Detection histories
   
###Code for constructing covariable here - variability in each site
	X<-  rnorm(K-1, mean = 0, sd = 1)
	
   # Determine initial occupancy and demographic parameters
   psi[1] <- psi1				# Initial occupancy probability
   p <- runif(n = K, min = range.p[1], max = range.p[2])
   
#### Write the probability of persistence as a function of the covariate and a transformation
   phi <- plogis(phi.alpha + phi.beta*X)
   gamma <- runif(n = K-1, min = range.gamma[1], max = range.gamma[2])

   # Generate latent states of occurrence
   # First year
   z[,1] <- rbinom(R, 1, psi[1])		# Initial occupancy state
   # Later years
   for(i in 1:R){				# Loop over sites
      for(k in 2:K){				# Loop over years
      # Prob for occupied = phi[k-1] if z[i, k-1] = 1 (persistence) or gamma[k-1] if z[i,k-1] = 0 (colonization)
      
      #Have to change the probability of persistence for each site and year
         muZ[k] <- z[i, k-1]*phi[k-1] + (1-z[i, k-1])*gamma[k-1] 
         #Actual state state given the probability muZ
         z[i,k] <- rbinom(1, 1, muZ[k])
         }
      }

# Plot realised occupancy
   par(mfrow = c(1,1))
   plot(year, apply(z, 2, mean), type = "l", xlab = "Year", ylab = "Occupancy or Detection prob.", col = "red", xlim = c(0,K+1), ylim = c(0,1.1), lwd = 2, lty = 1, frame.plot = FALSE, las = 1)
   lines(year, p , type = "l", col = "blue", lwd = 2, lty = 1)


   # Generate detection/nondetection (observation) data
   for(i in 1:R){
      for(k in 1:K){
         prob <- z[i,k] * p[k]
         for(j in 1:J){
            y[i,j,k] <- rbinom(1, 1, prob)
            }
         }
      }

   # Compute annual population occupancy
   for (k in 2:K){
      psi[k] <- psi[k-1]*phi[k-1] + (1-psi[k-1])*gamma[k-1]
      }

   # Plot apparent occupancy
   psi.app <- apply(apply(y, c(1,3), max), 2, mean)
   lines(year, psi.app, type = "l", col = "black", lwd = 2)
   legend("bottomleft", c("true occupancy", "detection", "observed occupancy"), lty = c(1,1,1), col = c(2,4,1), lwd = c(2,2,2), cex = 0.8 )
   
   
   # Return data
   return(list(R = R, J = J, K = K, psi = psi, psi.app = psi.app, z = z, phi = phi, gamma = gamma, p = p, y = y, X = X, phi.alpha = phi.alpha, phi.beta = phi.beta))
}

data <- data.fn(R = 250, J = 3, K = 10, psi1 = 0.6, range.p = c(0.1, 0.9),phi.alpha = -1, phi.beta = 1.5 , range.gamma = c(0.1, 0.5))

attach(data)
str(data)

# Specify model in BUGS language
sink("Dynocc2.jags")
cat("
model {

# Specify priors
psi1 ~ dunif(0, 1)
alpha.phi~dnorm(0, 0.1)
beta.phi~dnorm(0, 0.1)
for (k in 1:(nyear-1)){
   logit(phi[k]) <- alpha.phi + beta.phi*X[k]  # phi is the probability of persistence if z = 1
   gamma[k] ~ dunif(0, 1)  #gamma is the probability of colonization if z = 0
   p[k] ~ dunif(0, 1)   #p is the probability of detection
   }
p[nyear] ~ dunif(0, 1)

# Ecological submodel: Define state conditional on parameters
for (i in 1:nsite){  #the first year
   z[i,1] ~ dbern(psi1)  
   for (k in 2:nyear){    #additional years after the first
      muZ[i,k]<- z[i,k-1]*phi[k-1] + (1-z[i,k-1])*gamma[k-1]
      z[i,k] ~ dbern(muZ[i,k])
      } #k
   } #i

# Observation model
for (i in 1:nsite){
   for (j in 1:nrep){
      for (k in 1:nyear){
         muy[i,j,k] <- z[i,k]*p[k]
         y[i,j,k] ~ dbern(muy[i,j,k])
         } #k
      } #j
   } #i

# Derived parameters: Sample and population occupancy, growth rate and turnover
psi[1] <- psi1
n.occ[1]<-sum(z[1:nsite,1])
for (k in 2:nyear){
   psi[k] <- psi[k-1]*phi[k-1] + (1-psi[k-1])*gamma[k-1]
   n.occ[k] <- sum(z[1:nsite,k])
   growthr[k-1] <- psi[k]/psi[k-1]       
   #turnover = Pr(z=0 at time t)*Pr(colonization from t to t+1)/Pr(z = 1 at time t+1)
   turnover[k-1] <- (1 - psi[k-1]) * gamma[k-1]/psi[k]
   }
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = data$y, nsite = dim(data$y)[1], nrep = dim(data$y)[2], nyear = dim(data$y)[3], X=data$X)

# Initial values
zst <- apply(data$y, c(1, 3), max)	# Observed occurrence as inits for z
inits <- function(){ list(z = zst)}

# Parameters monitored
params <- c("psi", "alpha.phi", "beta.phi", "gamma", "p", "n.occ", "growthr", "turnover") 


# MCMC settings
ni <- 2500
nt <- 4
nb <- 500
nc <- 3

# Call JAGS from R (BRT 3 min)
out <- jags(win.data, inits, params, "Dynocc2.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, dig = 2)
psiall <- paste("psi[", 1:K, "]", sep="")
print(cbind(data$psi, out$BUGSoutput$summary[psiall, c(1, 2, 3, 7)]), dig = 3)
print( c( data$phi.alpha, out$BUGSoutput$summary["alpha.phi", c(1, 2, 3, 7)]), dig = 3 )
print(c(data$phi.beta, out$BUGSoutput$summary["beta.phi", c(1, 2, 3, 7)]), dig = 3)
gammaall <- paste("gamma[", 1:(K-1), "]", sep="")
print(cbind(data$gamma, out$BUGSoutput$summary[gammaall, c(1, 2, 3, 7)]), dig = 3)
pall <- paste("p[", 1:K, "]", sep="")
print(cbind(data$p, out$BUGSoutput$summary[pall, c(1, 2, 3, 7)]), dig = 3)

#Plot the estimates and the true values

plot(1:K, data$psi, type = "l", xlab = "Year", ylab = "Occupancy probability", col = "red", xlim = c(0,K+1), ylim = c(0,1), lwd = 2, lty = 1, frame.plot = FALSE, las = 1)
lines(1:K, data$psi.app, type = "l", col = "black", lwd = 2)
points(1:K, out$BUGSoutput$mean$psi, type = "l", col = "blue", lwd = 2)
#Plot of 95% posterior interval 
segments(1:K, out$BUGSoutput$summary[psiall,3], 1:K, out$BUGSoutput$summary[psiall,7], col = "blue", lwd = 1)




   


