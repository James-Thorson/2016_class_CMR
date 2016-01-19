#############################
# Statistical Inference in R and JAGS
#  Prof: Noble Hendrix & Jim Thorson
#  noblehendrix@gmail.com
#--------------------------
#  Day 2 Lab code - Patch occupancy models
#--------------------------

#############################
#	PART I: CLOSED PATCH OCCUPANCY
#############################
#R code based on workshop by Marc Kery, Richard Chandler and Andy Royle 2012



#1. SIMULATE DATA 

###### PROCESS MODEL

# Create a covariate called vegHtnSites <- 100set.seed(443)                   # so that we all get the same values of vegHtvegHt <- sort(runif(nSites, 1, 3)) # sort for convenience# Suppose that occupancy probability increases with vegHt# The relationship is described by an intercept of -3 and#    a slope parameter of 2 on the logit scalepsi <- plogis(-3 + 2*vegHt)# Now go to 100 sites and observe their occurrence state (perfectly)z <- rbinom(nSites, 1, psi)

# We can fit a model that relates abundance to vegHt using the glm() function#  with "family=binomial"
fm.glm1 <- glm(z ~ vegHt, family=binomial)
summary(fm.glm1)

# Plot the resultsplot(vegHt, z, xlab="Vegetation height", ylab="Occurrence (z)")glm1.est <- coef(fm.glm1)plot(function(x) plogis(-3 + 2*x), 1, 3, add=TRUE, lwd=3)plot(function(x) plogis(glm1.est[1] + glm1.est[2]*x), 1, 3, add=TRUE,lwd=3, col="blue")legend(2.5, 0.2, c("Truth", "Estimate"), col=c("black", "blue"), lty=1,lwd=3)

#Comment: discrepancy between truth and realized occurrence is due to the discrete nature of binary data.

#### OBSERVATION MODEL #########
# Introduce measurement error - so that the occupancy is a latent variable that is observed imperfectly as a function of a covariate, wind

nVisits <- 3wind <- array(rnorm(nSites * nVisits), dim = c(nSites, nVisits))p <- plogis(1 - 2*wind)# plot(p ~ wind)y <- matrix(NA, nSites, nVisits)for(i in 1:nSites) {    y[i,] <- rbinom(n= nVisits, size = z[i], prob = p[i,])}# Look at the datacbind(z=z, y1=y[,1], y2=y[,2], y3=y[,3])

#Fit the model in unmarked()
# Load library, format data and summarizelibrary(unmarked)
#unmarked has its own unique framework for including observed data and covariates - not that we have covariate for state process and covariate for observation processumf <- unmarkedFrameOccu(y = y, siteCovs = data.frame(vegHt = vegHt), obsCovs = list(wind = wind))  summary(umf)

# Fit a model and extract estimates# Detection covariates follow first tilde, then come occupancy covariatesfm.occ <- occu(~wind ~vegHt, data=umf)  
summary(fm.occ)
# All estimates are on logit-scalepar(mfrow = c(1,2))beta1 <- coef(fm.occ)plot(function(x) plogis(beta1[1] + beta1[2]*x), 1, 3, xlab="Vegetationheight", ylab="Occupancy probability", ylim = c(0, 1))plot(function(x) plogis(beta1[3] + beta1[4]*x), -3, 3, xlab="Wind",ylab="Detection probability", ylim = c(0, 1))


#What are the random effects?
ranef(fm.occ)  #each site posterior (empirical - bayes) estimate of occupancy

#What is the finite-sample occupancy - that is what is the number of sites occupied in the sample of sites actually studied?
sum(ranef(fm.occ)@post[,2,])

# Predictions of probability occupied for new values of vegHt, say 1.2 and 3.1newdat <- data.frame(vegHt=c(1.2, 3.1))predict(fm.occ, type="state", newdata=newdat)

#####
#Fit model in JAGS
library(R2jags)
# Specify model in JAGS  languagesink("model.txt")cat("model {# Priors - alpha.occ ~ dnorm(0, 0.4)beta.occ ~ dnorm(0, 0.4)alpha.p ~ dnorm(0, 0.4)beta.p ~ dnorm(0, 0.4)# Likelihood

for (i in 1:R) {   #process portion of the model
   logit(psi[i]) <- alpha.occ + beta.occ * vegHt[i]
   #psi <- plogis(-3 + 2*vegHt)  - equation used to simulate data
   # True state model for the partially observed true state   z[i] ~ dbern(psi[i])             # True occupancy z at site i
   #z <- rbinom(nSites, 1, psi)   - equation used to simulate data      #Observation portion of the model
   for (j in 1:T) {
   	  #Equation for covariates affecting observation process      logit(p[i,j]) <- alpha.p + beta.p * wind[i,j]
      #p <- plogis(1 - 2*wind) -  equation used to simulate data
      
      #can only see animals if they were truly there, (z[i]=1)
      p.eff[i,j] <- z[i] * p[i,j]
      
      # Actual observations      y[i,j] ~ dbern(p.eff[i,j])    # Detection-nondetection at i and j      #y[i,] <- rbinom(n= nVisits, size = z[i], prob = p[i,]) - equation used to simulate data
            } #j} #i# Derived quantitiesocc.fs <- sum(z[])}",fill = TRUE)sink()

# Bundle datawin.data <- list(y = y, vegHt = vegHt, wind = wind, R = nrow(y), T = ncol(y))

# Initial values# Number of occupied sites among those studiedzst <- apply(y, 1, max)inits <- function(){list(z = zst, alpha.occ = runif(1, -3, 3), beta.occ =runif(1, -3, 3), alpha.p = runif(1, -3, 3), beta.p = runif(1, -3, 3))}# Parameters monitoredparams <- c("alpha.occ", "beta.occ", "alpha.p", "beta.p", "occ.fs")

# add "z" for estimates of the random effects# MCMC settingsni <- 5000nt <- 2nb <- 2000nc <- 3

out2 <- jags(win.data, inits, params, "model.txt", n.chains = nc,   n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd() ) traceplot(out2)# Summarize posteriorsprint(out2, dig = 3)

#Put JAGS into format for reading into CODA for MCMC diagnostics
out2.mcmc<- as.mcmc(out2)
gelman.diag(out2.mcmc, multivariate = F)


#SOME PLOTS TO COMPARE ESTIMATES TO THE TRUE RELATIONSHIPS
# Graphical comparison between truth and estimatespredHt <- seq(1,3,,100)predW <- seq(-3,3,,100)
mle <- coef(fm.occ)bayesJ <- out2$BUGSoutput$summary[1:4,1]

#Occupancy probability
par(mfrow = c(1,2), mar = c(5,4,4,4), cex.lab = 1.2, cex.main = 1.2)
plot(vegHt, z, xlab="Vegetation height", ylab="", main="Occupancy probability", las = 1, lty = 1)lines(vegHt, psi, lwd=3, col="black")lines(predHt, plogis(mle[1] + mle[2]*predHt), lwd = 3, col="blue", lty = 1) lines(predHt, plogis(bayesJ[1] + bayesJ[3]*predHt), lwd = 3, col="red", lty = 3)
 legend(1.5, 0.25, c("Truth", "unmarked",  "JAGS"), col=c("black", "blue", "red"), lty=c(1, 1,  3), lwd=3)

#Detection probability
plot(wind, y, xlab="Wind", ylab="", main="Detection probability", las = 1, lty = 1)lines(sort(wind), p[order(wind)], lwd=3, col="black")lines(predW, plogis(mle[3] + mle[4]*predW), lwd = 3, col="blue", lty = 1) lines(predW, plogis(bayesJ[2] + bayesJ[4]*predW), lwd = 3, col="red", lty = 3)

#############################
#	PART 2: DYNAMIC PATCH OCCUPANCY
#############################
#R code from Bayesian Population Modeling Ch 13 by Kery and Schaub 2012

#1. SIMULATE DATA

data.fn <- function(R = 250, J = 3, K = 10, psi1 = 0.4, range.p = c(0.2, 0.4), range.phi = c(0.6, 0.8), range.gamma = c(0, 0.1)) {
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

   # Determine initial occupancy and demographic parameters
   psi[1] <- psi1				# Initial occupancy probability
   p <- runif(n = K, min = range.p[1], max = range.p[2])
   phi <- runif(n = K-1, min = range.phi[1], max = range.phi[2])
   gamma <- runif(n = K-1, min = range.gamma[1], max = range.gamma[2])

   # Generate latent states of occurrence
   # First year
   z[,1] <- rbinom(R, 1, psi[1])		# Initial occupancy state
   # Later years
   for(i in 1:R){				# Loop over sites
      for(k in 2:K){				# Loop over years
      # Prob for occupied = phi[k-1] if z[i, k-1] = 1 (persistence) or gamma[k-1] if z[i,k-1] = 0 (colonization)
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
   return(list(R = R, J = J, K = K, psi = psi, psi.app = psi.app, z = z, phi = phi, gamma = gamma, p = p, y = y))
}

data <- data.fn(R = 250, J = 3, K = 10, psi1 = 0.6, range.p = c(0.1, 0.9), range.phi = c(0.7, 0.9), range.gamma = c(0.1, 0.5))

attach(data)
str(data)

# Specify model in BUGS language
sink("Dynocc.jags")
cat("
model {

# Specify priors
psi1 ~ dunif(0, 1)
for (k in 1:(nyear-1)){
   phi[k] ~ dunif(0, 1)  # phi is the probability of persistence if z = 1
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
win.data <- list(y = data$y, nsite = dim(data$y)[1], nrep = dim(data$y)[2], nyear = dim(data$y)[3])

# Initial values
zst <- apply(data$y, c(1, 3), max)	# Observed occurrence as inits for z
inits <- function(){ list(z = zst)}

# Parameters monitored
params <- c("psi", "phi", "gamma", "p", "n.occ", "growthr", "turnover") 


# MCMC settings
ni <- 2500
nt <- 4
nb <- 500
nc <- 3

# Call JAGS from R (BRT 3 min)
out <- jags(win.data, inits, params, "Dynocc.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, dig = 2)
psiall <- paste("psi[", 1:K, "]", sep="")
print(cbind(data$psi, out$BUGSoutput$summary[psiall, c(1, 2, 3, 7)]), dig = 3)
phiall <- paste("phi[", 1:(K-1), "]", sep="")
print(cbind(data$phi, out$BUGSoutput$summary[phiall, c(1, 2, 3, 7)]), dig = 3)
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




   


