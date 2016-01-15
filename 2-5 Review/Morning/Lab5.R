#  Inference with R and JAGS
#  Bayesian Models
#  15 January 2015
#  Noble Hendrix & Jim Thorson
#--------------------------
#  Day 5 Lab code
#--------------------------

####################
# Poisson GLM
# Example taken from Kery and Schaub Bayesian Population Analysis 2012

data.fn <- function(n = 40, alpha = 3.5576, beta1 = -0.0912, beta2 = 0.0091, beta3 = -0.00014){
# n: Number of years
# alpha, beta1, beta2, beta3: coefficients of a 
#    cubic polynomial of count on year

# Generate values of time covariate
year <- 1:n

# Signal: Build up systematic part of the GLM
log.expected.count <- alpha + beta1 * year + beta2 * year^2 + beta3 * year^3
expected.count <- exp(log.expected.count)

# Noise: generate random part of the GLM: Poisson noise around expected counts
C <- rpois(n = n, lambda = expected.count)

# Plot simulated data
plot(year, C, type = "b", lwd = 2, col = "black", main = "", las = 1, ylab = "Population size", xlab = "Year", cex.lab = 1.2, cex.axis = 1.2)
lines(year, expected.count, type = "l", lwd = 3, col = "red")

return(list(n = n, alpha = alpha, beta1 = beta1, beta2 = beta2, beta3 = beta3, year = year, expected.count = expected.count, C = C))
}

data <- data.fn()

fm <- glm(C ~ year + I(year^2) + I(year^3), family = poisson, data = data)
summary(fm)


##  JAGS Model ##
library(R2jags)

# Specify model in BUGS language
sink("GLM_Poisson.jags")
cat("
model {

# Priors
alpha ~ dunif(-20, 20)
beta1 ~ dunif(-10, 10)
beta2 ~ dunif(-10, 10)
beta3 ~ dunif(-10, 10)

# Likelihood: Note key components of a GLM on one line each
for (i in 1:n){
	# 1. Linear predictor
   log.lambda[i] <- alpha + beta1 * year[i] + beta2 * pow(year[i],2) + beta3 * pow(year[i],3)  					
   log(lambda[i]) <- log.lambda[i]  # 2. Link function
   C[i] ~ dpois(lambda[i])          # 3. Distribution for measurement error                    
   } #i
   
   #make prediction for next year
   log(lambda.new) <- alpha + beta1 * year.new + beta2 * pow(year.new,2) + beta3 * pow(year.new,3)  #predicted signal or mean level
   C.new ~ dpois(lambda.new) #distribution to add random part
}
",fill = TRUE)
sink()

#Need to center covariate data to avoid numerical overflow issues
# Bundle data
mean.year <- mean(data$year)             # Mean of year covariate
sd.year <- sd(data$year)                 # SD of year covariate
s.year<- (data$year - mean.year) / sd.year

win.data <- list(C = data$C, n = length(data$C), year = s.year,
year.new = max(s.year) + 1/sd.year)


# Initial values
inits <- function() list(alpha = runif(1, -2, 2), beta1 = runif(1, -3, 3))

# Parameters monitored
params <- c("alpha", "beta1", "beta2", "lambda", "beta3", "lambda.new", "C.new")

# MCMC settings
ni <- 2000
nt <- 2
nb <- 1000
nc <- 3


# Call JAGS from R 
pois.jags <- jags(data = win.data, inits = inits, parameters.to.save = params, model.file = "GLM_Poisson.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)# working.directory = getwd())

#look at estimates and Rhat values
pois.jags


#############################
# Poisson GLMM - Mixed Model
############################
#Falcon example and add random effects
#-------------------------------------
data2.fn <- function(n = 40, alpha = 3.5576, beta1 = -0.0912, beta2 = 0.0091, beta3 = -0.00014, sd = 0.1){
   # n: Number of years
   # alpha, beta1, beta2, beta3: coefficients of a 
   #    cubic polynomial of count on year
   # sd: standard deviation of normal distribution assumed for year effects

   # Generate values of time covariate
   year <- 1:n

   # First level of noise: generate random year effects
   eps <- rnorm(n = n, mean = 0, sd = sd)

   # Signal (plus first level of noise): build up systematic part of the GLM and add the random year effects
   log.expected.count <- alpha + beta1 * year + beta2 * year^2 + beta3 * year^3 + eps
   expected.count <- exp(log.expected.count)

   # Second level of noise: generate random part of the GLM: Poisson noise around expected counts
   C <- rpois(n = n, lambda = expected.count)

   # Plot simulated data
   plot(year, C, type = "b", lwd = 2, main = "", las = 1, ylab = "Population size", xlab = "Year", ylim = c(0, 1.1*max(C)))
   lines(year, expected.count, type = "l", lwd = 3, col = "red")

   return(list(n = n, alpha = alpha, beta1 = beta1, beta2 = beta2, beta3 = beta3, year = year, sd = sd, expected.count = expected.count, C = C))
   }

data2 <- data2.fn()

#add lines from data1 to see how the process parts compare
lines(data$expected.count, col = 4, lwd = 3)
legend("topleft", c("No random effects", "Random effects"), col = c(4,2), lwd = c(3,3))


#################
#Fit GLMM in R - 
################
library(lme4)
# Create a factor 'yr' for using in random effects call 
yr <- factor(data2$year)         

#standardize covariates
mny <- mean(data2$year)
sdy <- sd(data2$year)
cov1 <- (data2$year - mny) / sdy
cov2 <- cov1 * cov1
cov3 <- cov1 * cov1 * cov1
glmm.fit <- glmer(C ~ (1 | yr) + cov1 + cov2 + cov3, family = poisson, data = data2)
glmm.fit

# Plot simulated data
  plot(data$year, data$C, type = "b", lwd = 2, main = "", las = 1, ylab = "Population size", xlab = "Year", ylim = c(0, 1.1*max(data$C)))
   lines(data$year, data$expected.count, type = "l", lwd = 3, col = "red")
#Make predictions from R glmm.fit object
R.predictions <- exp(fixef(glmm.fit)[1] + fixef(glmm.fit)[2]*cov1 + fixef(glmm.fit)[3]*cov2 + fixef(glmm.fit)[4]*cov3 + unlist(ranef(glmm.fit)))
#Plot predictions from R glmm model
lines(data$year, R.predictions, col = "green", lwd = 2, type = "l")


##########
##GLMM in JAGS 
##########
library(R2jags)

sink("GLMM_Poisson.jags")
cat("
model {

# Priors
alpha ~ dunif(-20, 20)
beta1 ~ dunif(-10, 10)
beta2 ~ dunif(-10, 10)
beta3 ~ dunif(-10, 10)
sd ~ dunif(0, 5)
tau <- 1 / (sd*sd)

# Likelihood: note key components of a GLMM 
for (i in 1:n){
    # 1. Linear predictor including random year effect
   log.lambda[i] <- alpha + beta1 * year[i] + beta2 * pow(year[i],2) + beta3 * pow(year[i],3) + eps[i]     
   # 2. Distribution for random effect
   eps[i] ~ dnorm(0, tau) 			
   # 3. Link function
   log(lambda[i]) <- log.lambda[i]  
   # 4. Distribution for random part
   C[i] ~ dpois(lambda[i])             
   }
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = data2$C, n = length(data2$C), year = cov1)

# Initial values
inits <- function() list(alpha = runif(1, -2, 2), beta1 = runif(1, -3, 3), sd = runif(1, 0,1))

# Parameters monitored
params <- c("alpha", "beta1", "beta2", "beta3", "lambda", "sd", "eps")

# MCMC settings
ni <- 30000
nt <- 10
nb <- 20000
nc <- 3

# Call JAGS from R (BRT <1 min)
glmm1<- jags(win.data, inits, params, "GLMM_Poisson.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)#, working.directory = getwd())
	
# Summarize posteriors
print(glmm1, dig = 2)

glmm1.mcmc<- as.mcmc(glmm1)

JAGS.predictions <- glmm1$BUGSoutput$mean$lambda
lines(data2$year, JAGS.predictions, type = "l", lwd = 3, col = "blue", lty = 2)

#95% credible intervals:
attach.jags(glmm1)
lambda.q<- apply(lambda, 2, function(x) quantile(x,c(0.025, 0.975)) )
detach.jags()

lines(data2$year, lambda.q[1,], col = "blue", lty = 2)
lines(data2$year, lambda.q[2,], col = "blue", lty = 2)


##########
#Adding spatial random effects with multiple populations
#########

#5 populations of falcons with random effects for each population and year
# Example 4.3.1 from Kery & Schaub 2012: Mixed models with random effects for variability among groups (site and year effects)
# 
data3.fn <- function(nsite = 5, nyear = 40, alpha = 4.18456, beta1 = 1.90672, beta2 = 0.10852, beta3 = -1.17121, sd.site = 0.5, sd.year = 0.2){
   # nsite: Number of populations
   # nyear: Number of years
   # alpha, beta1, beta2, beta3: cubic polynomial coefficients of year
   # sd.site: standard deviation of the normal distribution assumed for the population intercepts alpha
   # sd.year: standard deviation of the normal distribution assumed for the year effects
   # We standardize the year covariate so that it runs from about -1 to 1

   # Generate data structure to hold counts and log(lambda)
   C <- log.expected.count <- array(NA, dim = c(nyear, nsite))

   # Generate covariate values
   year <- 1:nyear
   yr <- (year-20)/20	# Standardize
   site <- 1:nsite

   # Draw two sets of random effects from their respective distribution
   alpha.site <- rnorm(n = nsite, mean = alpha, sd = sd.site)
   eps.year <- rnorm(n = nyear, mean = 0, sd = sd.year)

   # Loop over populations
   for (j in 1:nsite){
      # Signal (plus first level of noise): build up systematic part of the GLM including random site and year effects
      log.expected.count[,j] <- alpha.site[j] + beta1 * yr + beta2 * yr^2 + beta3 * yr^3 + eps.year
      expected.count <- exp(log.expected.count[,j])

      # Second level of noise: generate random part of the GLM: Poisson noise around expected counts
      C[,j] <- rpois(n = nyear, lambda = expected.count)
      }

   # Plot simulated data
   matplot(year, C, type = "l", lty = 1, lwd = 2, main = "", las = 1, ylab = "Population size", xlab = "Year")

   return(list(nsite = nsite, nyear = nyear, alpha.site = alpha.site, beta1 = beta1, beta2 = beta2, beta3 = beta3, year = year, sd.site = sd.site, sd.year = sd.year, expected.count = expected.count, C = C))
   }

data3 <- data3.fn(nsite = 6, nyear = 40, sd.site = 0.3, sd.year = 0.2)


# Specify model in BUGS language
sink("GLMM_Poisson2.jags")
cat("
model {

# Priors
for (j in 1:nsite){
   alpha[j] ~ dnorm(mu.alpha, tau.alpha)		# 4. Random site effects
   }
mu.alpha ~ dnorm(0, 0.01)				# Hyperparameter 1
tau.alpha <- 1 / (sd.alpha*sd.alpha)	        #Hyperparameter 2
sd.alpha ~ dunif(0, 2)
for (p in 1:3){
   beta[p] ~ dnorm(0, 0.01)
   }

tau.year <- 1 / (sd.year*sd.year)
sd.year ~ dunif(0, 1)				# Hyperparameter 3

# Likelihood
for (i in 1:nyear){
   eps[i] ~ dnorm(0, tau.year)                # 1. Random year effects
   for (j in 1:nsite){
   	log.lambda[i,j] <- alpha[j] + beta[1] * year[i] + beta[2] * pow(year[i],2) + beta[3] * pow(year[i],3) + eps[i]    # 2. Linear predictor including random site and random year effects
      lambda[i,j] <- exp(log.lambda[i,j])     # 3. Link function
      C[i,j] ~ dpois(lambda[i,j])             # 4. Distribution for random part
      
      
      }  #j
   }  #i
}
",fill = TRUE)
sink()


# Bundle data
win.data <- list(C = data3$C, nsite = ncol(data3$C), nyear = nrow(data3$C), year = (data3$year-20) / 20) # Note year standardized

# Initial values
inits <- function() list(mu.alpha = runif(1, 0, 2), alpha = runif(data3$nsite, -1, 1), beta = runif(3, -1, 1), sd.alpha = runif(1, 0, 0.1), sd.year = runif(1, 0, 0.1))

# Parameters monitored (may want to add "lambda")
params <- c("mu.alpha", "alpha", "beta", "sd.alpha", "sd.year")

# MCMC settings (may have to adapt)
ni <- 10000
nt <- 50
nb <- 2000
nc <- 3

#NOTE LONG TIME TO RUN THIS MODEL IF HAVE LOTS OF SITES
# Call JAGS from R 

out <- jags(win.data, inits, params, "GLMM_Poisson2.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)#, working.directory = getwd())

# Summarize posteriors
print(out, dig = 3)

### What if we want to run the model for longer?
#replaces out with new draws (per chain), treating the values in the original fit as burn-in
out<- update(out, 20000)
print(out, dig = 3)




######  ADDITIONAL TO DO ON YOUR OWN:  
##### FIT GLM AND GLMM MODELS TO PEREGRIN FALCON DATA SETS
#####  KERY CH 4
######

#Fit to the real data set of peregrine abundances:

# Read data
peregrine <- read.table("falcons.txt", header = TRUE)
attach(peregrine)
plot(Year, Pairs, type = "b", lwd = 2, main = "", las = 1, ylab = "Pair count", xlab = "Year", ylim = c(0, 200), pch = 16, xlim = c(min(Year), max(Year)+1))

# Bundle data
mean.year <- mean(Year)        # Mean of year covariate
sd.year <- sd(Year)            # SD of year covariate
s.year<- (Year - mean.year) / sd.year

win.data <- list(C = Pairs, n = length(Pairs), year = s.year, year.new = max(s.year) + 1/sd.year)

# Initial values
inits <- function() list(alpha = runif(1, -2, 2), beta1 = runif(1, -3, 3))

# Parameters monitored
params <- c("alpha", "beta1", "beta2", "lambda", "beta3", "lambda.new", "C.new")

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call JAGS from R (BRT < 1 min)
out1 <- jags(data = win.data, inits = inits, parameters.to.save = params, model.file = "GLM_Poisson.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)#, working.directory = getwd())

# Summarize posteriors
print(out1, dig = 3) 


JAGS.predictions <- out1$BUGSoutput$mean$lambda
lines(Year, JAGS.predictions, type = "l", lwd = 3, col = "blue", lty = 2)
points(max(Year) + 1, out1$BUGSoutput$mean$lambda.new, col = "red")

#95% credible intervals:
attach.jags(out1)
lambda.q<- apply(lambda, 2, function(x) quantile(x,c(0.025, 0.975)) )
lambda.q.new<- quantile(lambda.new, c(0.025, 0.975))
detach.jags()

lines(Year, lambda.q[1,], col = "blue", lty = 2)
lines(Year, lambda.q[2,], col = "blue", lty = 2)
points(max(Year) +1, lambda.q.new[1], pch = "-", col = "red")
points(max(Year) +1, lambda.q.new[2], pch = "-", col = "red")

detach(peregrine)


# Analysis of peregrine data with random effects model:

peregrine <- read.table("falcons.txt", header = TRUE)

yr <- factor(peregrine$Year)
mny <- mean(peregrine$Year)
sdy <- sd(peregrine$Year)
cov1 <- (peregrine$Year - mny) / sdy
cov2 <- cov1 * cov1
cov3 <- cov1 * cov1 * cov1
glmm <- glmer(peregrine$Pairs ~ (1 | yr) + cov1 + cov2 + cov3, family = poisson, data = peregrine)
glmm

# Bundle data
win.data <- list(C = peregrine$Pairs, n = length(peregrine$Pairs), year = cov1)

# Initial values
inits <- function() list(alpha = runif(1, -2, 2), beta1 = runif(1, -3, 3), sd = runif(1, 0,1))

# Parameters monitored
params <- c("alpha", "beta1", "beta2", "beta3", "lambda", "sd", "eps")

# MCMC settings (may have to adapt)
ni <- 30000
nt <- 10
nb <- 20000
nc <- 3

# Call JAGS from R (BRT < 1 min)
falc.jags<- jags(win.data, inits, params, "GLMM_Poisson.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)#, working.directory = getwd())

falc.mcmc<- as.mcmc(falc.jags)
#Gelman-Rubin statistic
gelman.diag(falc.mcmc, multivariate = F)
geweke.diag(falc.mcmc)

# Summarize posteriors
print(falc.jags, dig = 3)