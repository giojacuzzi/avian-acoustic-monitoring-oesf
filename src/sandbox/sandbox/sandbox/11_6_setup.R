# 11.3 Metacommunity data from the Swiss breeding bird survey MHB
# ------------------------------------------------------------------------

## Code modified to use the built-in data set MHB2014 instead of the file "MHB_2014.csv"
library(AHMbook)
data(MHB2014)
?MHB2014
str(MHB2014)
# NB some of the data preprocessing on p.644 has already been done.

# Check the detection data in 3D array MHB2014$count: site x rep x species
( nsite <- nrow(MHB2014$sites) )    # number of sites in Swiss MHB
nrep <- 3                           # maximum number of replicate surveys per season
( nspec <- nrow(MHB2014$species) )  # 158 species occur in the 2014 data
dim(MHB2014$count) == c(nsite, nrep, nspec) # check

# Create the detection/nondetection (1/0) array
y <- MHB2014$count ; y[y > 1] <- 1  ## 'Y' replaced with 'y'
str(y)

# Check data for one species, here chaffinch, and pull them out from 3D array
(tmp <- y[, , "Common Chaffinch"])

# Frequency distribution of number of surveys actually carried out per site in 2014
# NB MHB2014$sites$nsurvey gives the number of surveys *planned*.
table(nsurveys <- apply(!is.na(y[,,1]), 1, sum))

# Which site has all NA data in 2014 ?
(NAsites <- which(nsurveys == 0) )

# Observed number of occupied sites
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
# For the 'all NA' site, max returns -Inf with a warning
tmp[tmp == -Inf] <- NA         # Change -Inf to NA
sort(obs.occ <- apply(tmp, 2, sum, na.rm = TRUE))

# Plot species 'occurrence frequency' distribution (not shown)
plot(sort(obs.occ), xlab = "Species number", ylab = "Number of quads with detections")

# Drop data from species that were not observed in 2014
toss.out <- which(obs.occ == 0)
y <- y[,,-toss.out]
obs.occ <- obs.occ[-toss.out]

# Redefine nspec as the number of species observed in 2014: 145
( nspec <- dim(y)[3] )

str(y)

# Get observed number of species per site
tmp <- apply(y, c(1,3), max, na.rm = TRUE)
tmp[tmp == "-Inf"] <- NA
sort(C <- apply(tmp, 1, sum))     # Compute and print sorted species counts

plot(table(C), xlim = c(0, 60), xlab = "Observed number of species", ylab = "Number of quadrats", frame = FALSE)
abline(v = mean(C, na.rm = TRUE), col = "blue", lwd = 3)

# 11.5.1 Simple Poisson regression for the observed community size
# ------------------------------------------------------------------------
# Get covariates and standardise them
# Quadrat elevation and forest cover
data <- MHB2014
orig.ele <- data$sites$elev[1:nsite]
(mean.ele <- mean(orig.ele, na.rm = TRUE))
(sd.ele <- sd(orig.ele, na.rm = TRUE))
ele <- (orig.ele - mean.ele) / sd.ele
orig.forest <- data$sites$forest[1:nsite]
(mean.forest <- mean(orig.forest, na.rm = TRUE))
(sd.forest <- sd(orig.forest, na.rm = TRUE))
forest <- (orig.forest - mean.forest) / sd.forest

# Average date and duration of survey
tmp <- cbind(data$date)[1:nsite,] 
orig.mdate <- apply(tmp, 1, mean, na.rm = TRUE)
(mean.mdate <- mean(orig.mdate[-NAsites]))   # drop unsurved site
(sd.mdate <- sd(orig.mdate[-NAsites]))
mdate <- (orig.mdate - mean.mdate) / sd.mdate
mdate[NAsites] <- 0                 # impute mean for missing

tmp <- cbind(data$dur)[1:nsite,] 
orig.mdur <- apply(tmp, 1, mean, na.rm = TRUE)
(mean.mdur <- mean(orig.mdur[-NAsites]))
(sd.mdur <- sd(orig.mdur[-NAsites]))
mdur <- (orig.mdur - mean.mdur) / sd.mdur
mdur[NAsites] <- 0                  # impute mean for missing

# Get observed species richness per site and rep and plot
CC <- apply(y, c(1,2), sum, na.rm = TRUE)
CC[CC == 0] <- NA            # 0 means not surveyed
matplot(t(CC), type = 'l', lty = 1, lwd = 2, xlab = "First to third survey", ylab = "Number of species detected", frame = F)  # Fig. 11�6 right

# Get survey date and survey duration and standardise both
# Survey date (this is Julian date, with day 1 being April 1)
orig.DAT <- cbind(data$date)[1:nsite,]
(mean.date <- mean(orig.DAT, na.rm = TRUE))
(sd.date <- sd(c(orig.DAT), na.rm = TRUE))
DAT <- (orig.DAT - mean.date) / sd.date      # scale
DAT[is.na(DAT)] <- 0                         # impute missings
# Survey duration (in minutes)
orig.DUR <- cbind(data$dur)[1:nsite,]
(mean.dur <- mean(orig.DUR, na.rm = TRUE))
(sd.dur <- sd(c(orig.DUR), na.rm = TRUE))
DUR <- (orig.DUR - mean.dur) / sd.dur        # scale
DUR[is.na(DUR)] <- 0                         # mean impute missings

# 11.6.1 Simplest community occupancy model: n-fold single species 
#        occupancy model with species treated as fixed effects
# ------------------------------------------------------------------------
# In community occupancy models, in the absence of covariates that vary by replicate survey j,itis convenient to aggregate binary detection/nondetection data yijk into site- and occasion-specific counts by summing over replicates, and then model detection frequency ysumikdi.e., the number of detections of species k at site idas a binomial random variable with a binomial index given by the number of surveys (which may or may not be the same for all sites and/or species). Binomial versions of the models are computationally much more efficient to analyze in BUGS.

# Collapse 3D detection/nondetection data to 2D detection frequencies
ysum <- apply(y, c(1,3), sum, na.rm = T) # Collapse to detection frequency
ysum[NAsites,] <- NA                     # Have to NA out sites with NA data

# Bundle and summarize data set
str( win.data <- list(ysum = ysum, M = nrow(ysum), J = MHB2014$sites$nsurvey[1:nsite], nspec = dim(ysum)[2]) )


# NIMBLE ################### ################### ################### ################### ###################
library(nimble)
# Specify model in BUGS language
# "model5.txt" in book code
Section11p6p1_code <- nimbleCode( {
  
  # Priors
  for(k in 1:nspec){          # Loop over species
    psi[k] ~ dunif(0, 1)
    p[k] ~ dunif(0, 1)
  }
  
  # Ecological model for latent occurrence z (process model)
  for(k in 1:nspec){          # Loop over species
    for (i in 1:M) {         # Loop over sites
      z[i,k] ~ dbern(psi[k])
    }
  }
  
  # Observation model for observed data Y
  for(k in 1:nspec){          # Loop over species
    for (i in 1:M) {
      mup[i,k] <- z[i,k] * p[k]
      ysum[i,k] ~ dbin(mup[i,k], J[i])
    }
  }
  
  # Derived quantities
  for(k in 1:nspec){          # Loop over species
    Nocc.fs[k] <- sum(z[1:M,k]) # Add up number of occupied sites among the 267
  }
  for (i in 1:M) {            # Loop over sites
    Nsite[i] <- sum(z[i,1:nspec])   # Add up number of occurring species at each site
  }
}
)

# Initial values
zst <- apply(y, c(1,3), max) # Observed occurrence as inits for z
zst[is.na(zst)] <- 1
inits <- function() list(z = zst, 
                         psi = rep(0.4, nspec), 
                         p = rep(0.4, nspec))
# Parameters monitored
params <- c("psi", "p", "Nsite", "Nocc.fs")

# MCMC settings
ni <- 2500   ;   nt <- 2   ;   nb <- 500   ;   nc <- 3

m5 <- nimbleModel(
  code = Section11p6p1_code,
  constants = win.data,
  inits = inits(),
  calculate = FALSE) # Faster to skip calculation here

cm5 <- compileNimble(m5)

mcmc5 <- buildMCMC(m5, monitors = params, thin = nt, useConjugacy = FALSE) # Faster to skip checks for conjugacy

cmcmc5 <- compileNimble(mcmc5, project = m5)

out5 <- runMCMC(cmcmc5, niter = ni, samplesAsCodaMCMC = TRUE)

# Call nimble from R
out5 <- nimbleMCMC(
  code = Section11p6p1_code,
  constants = win.data,
  inits = inits,
  monitors = params,
  nburnin = nb,
  niter = 2*ni,
  samplesAsCodaMCMC = TRUE,
  check = FALSE
)

library(mcmcplots)
colnames(out5)
params_of_interest <- c('Nocc.fs[50]', 'Nsite[1]', 'p[1]', 'psi[1]') ## small subsample of parameters
mcmcplot(out5[, params_of_interest])

# BUGS ################### ################### ################### ################### ###################

# Specify model in BUGS language
sink("model5.txt")
cat("
model {

# Priors
for (k in 1:nspec) { # Loop over species
  psi[k] ~ dunif(0,1)
  p[k] ~ dunif(0,1)
}

# Ecological model for latent occurrence z (process model)
for (k in 1:nspec) { # loop over species
  for (i in 1:M) { # loop over sites
    z[i,k] ~ dbern(psi[k])
  }
}

# Observation model for observed data ysum
for (k in 1:nspec) {
  for (i in 1:M) {
    mup[i,k] <- z[i,k] * p[k]
    ysum[i,k] ~ dbin(mup[i,k], J[i])
  }
}

# Derived quantities
for (k in 1:nspec) {
  Nocc.fs[k] <- sum(z[,k]) # add up number of occupied sites
}
for (i in 1:M) {
  Nsite[i] <- sum(z[i,]) # add up number of occurring species at each site
}

}    
", fill = TRUE)
sink()

# Initial values
zst <- apply(y, c(1,3), max) # observed occurrence as inits for z
zst[is.na(zst)] <- 1
inits <- function() list(z = zst, psi = rep(0.4, nspec), p = rep(0.4, nspec))

# Parameters monitored
params <- c("psi", "p", "Nsite", "Nocc.fs")

# MCMC settings
ni <- 2500; nt <-2; nb <- 500; nc <- 3

# Call JAGS from R
out5_jags <- jags(win.data, inits, params, "model5.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
par(mfrow=c(4,4))
# traceplot(out5_jags)
print(out5_jags, dig = 3)

# Compare observed and estimated site species richness
par(cex = 1.3)
twice <- which(data$nsurvey[1:267] == 2)
plot(C[twice], out5_jags$summary[291:557,1][twice], xlab = "Observed number of species", ylab = "Estimated number of species", frame = F, xlim = c(0, 60), ylim = c(0, 70), col = "red", pch = 16)
segments(C[twice], out5_jags$summary[291:557,3][twice], C[twice], out5_jags$summary[291:557,7][twice], col = "red")
points(C[-twice], out5_jags$summary[291:557,1][-twice], col = "blue", pch = 16)
segments(C[-twice], out5_jags$summary[291:557,3][-twice], C[-twice], out5_jags$summary[291:557,7][-twice], col = "blue")

# Observed and estimated number of occupied sites for each species in a table and a plot
cbind(obs.occu = obs.occ, out5_jags$summary[558:702, c(1,3,7)])
plot(obs.occ, out5_jags$summary[558:702, 1], xlab = "Observed number of occupied sites", ylab = "Estimated version of quantity", ylim = c(0, 267), frame = F, pch = 16)
abline(0,1)
segments(obs.occ, out5_jags$summary[558:702,3], obs.occ, out5_jags$summary[558:702,7], col = "grey", lwd = 2)

# Estimated occupancy and detection probability for each species (model 5)
plot(out5_jags$summary[1:145,1], out5_jags$summary[146:290,1], xlab = "Occupancy estimate", ylab = "Detection estimate", xlim = c(0,1), ylim = c(0,1), frame = F, pch = 16)
segments(out5_jags$summary[1:145,3], out5_jags$summary[146:290,1], out5_jags$summary[1:145,7], out5_jags$summary[146:290,1], col = "grey", lwd = 2)
segments(out5_jags$summary[1:145,1], out5_jags$summary[146:290,3], out5_jags$summary[1:145,1], out5_jags$summary[146:290,7], col = "grey", lwd = 2)

# 11.6.2 Community occupancy model with bivariate species-specific random effects
#
# ------------------------------------------------------------------------
# Bundle and summarize data set
str( win.data <- list(ysum = ysum, M = nrow(ysum), J = MHB2014$sites$nsurvey[1:nsite], nspec = dim(ysum)[2], R = matrix(c(5,0,0,1), ncol = 2), df = 3) )

# Specify model in BUGS language
sink("model6.txt")
cat("
model {

# Priors
for (k in 1:nspec) { # group lpsi and lp together in array eta
  lpsi[k] <- eta[k,1]
  lp[k] <- eta[k,2]
  eta[k, 1:2] ~ dmnorm(mu.eta[], Omega[,])
}

# Hyperpriors
# Priors for mu.lpsi=mu.eta[1] and mu.lp=mu.eta[2]
# probs = community means of occupancy and detection probability
for (v in 1:2) {
  mu.eta[v] <- log(probs[v] / (1-probs[v]))
  probs[v] ~ dunif(0,1)
}
# Prior for variance-covariance matrix
Omega[1:2, 1:2] ~ dwish(R[,], df)
Sigma[1:2, 1:2] <- inverse(Omega[,])

# Ecological model for latent occurrence z (process model)
for (k in 1:nspec) {
  logit(psi[k]) <- lpsi[k] # must take outside of i loop (b/c only indexed k)
  for (i in 1:M) {
    z[i,k] ~ dbern(psi[k])
  }
}

# Observation model for observed data ysum
for (k in 1:nspec) {
  logit(p[k]) <- lp[k]
  for (i in 1:M) {
    mu.p[i,k] <- z[i,k] * p[k]
    ysum[i,k] ~ dbin(mu.p[i,k], J[i])
  }
}

# Derived quantities
rho <- Sigma[1,2] / sqrt(Sigma[1,1] * Sigma[2,2]) # Correlation coefficient
for (k in 1:nspec) {
  Nocc.fs[k] <- sum(z[,k]) # Number of occupied sites among the 267
}
for (i in 1:M) {
  Nsite[i] <- sum(z[i,]) # Number of occurring species
}

}", fill=TRUE)
sink()

# Initial values
zst <- apply(y, c(1,3), max) # observed occurrence as starting values for z
zst[is.na(zst)] <- 1
inits <- function() list(z = zst, Omega = matrix(c(1,0,0,1), ncol = 2), eta = matrix(0, nrow = nspec, ncol = 2))

# Parameters monitored
params <- c("mu.eta", "probs", "psi", "p", "Nsite", "Nocc.fs", "Sigma", "rho")

#MCMC settings
ni <- 20000 ; nt <- 15 ; nb <- 5000 ; nc <- 3

# Call JAGS
out6 <- jags(win.data, inits, params, "model6.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

# Summarize posteriors
print(out6, 3)

# Graphically compare some estimates between fixed- and random-effects model
par(mfrow = c(2,2))

# Species-specific occupancy (probability scale)
plot(out5_jags$summary[1:145,1], out6$summary[5:149,1], main = "Species-specific occupancy probability") ; abline(0,1)

# Species-specific detection (probability scale)
plot(out5_jags$summary[146:290,1], out6$summary[150:294,1], main = "Species-specific detection probability") ; abline(0,1)

# Site-specific species richness
plot(out5_jags$summary[291:557,1], out6$summary[295:561,1], main = "Site-specific species richness (conditional on list of 145 detected)") ; abline(0,1)

# Species-specific number of presences
plot(out5_jags$summary[558:702,1], out6$summary[562:706,1], main = "Species-specific number of presences (in 267 sites)") ; abline(0,1)

# Estimated occupancy and detection probability for each species
parmfrow=c(1,0)
plot(out6$summary[5:149,1], out6$summary[150:294,1], xlab = "Occupancy estimate", ylab = "Detection estimate", xlim = c(0,1), ylim = c(0,1), frame = F, pch = 16)
segments(out6$summary[5:149,3], out6$summary[150:294,1], out6$summary[5:149,7], out6$summary[150:294,1], col = "grey", lwd = 2)
segments(out6$summary[5:149,1], out6$summary[150:294,3], out6$summary[5:149,1], out6$summary[150:294,7], col = "grey", lwd = 2)

# Plot posterior distribution of site-specific species richness (Nsite)
par(mfrow = c(3,3), mar = c(5,4,3,2))
for(i in c(9, 32, 162, 12, 27, 30, 118, 159, 250)){
  plot(table(out6$sims.list$Nsite[,i]), main = paste("Quadrat", i), xlab = "Local species richness", ylab = "", frame = F, xlim = c((min(C[i], out6$sims.list$Nsite[,i], na.rm = T)-2), max(out6$sims.list$Nsite[,i]) ))
  abline(v = C[i], col = "grey", lwd = 4)
  browser()
}


# 11.6.2 Modeling species-specific effects in community occupancy models
#
# ------------------------------------------------------------------------
# Frequently, we want to compare groups of speciesdfor instance, guilds (e.g., herbivores vs carnivores), life-history categories (slow vs fast species), migration modes (resident vs migratory species), or other distinctions between species. We then treat species again as random effects, but specify a linear model for the hyperparameters that govern the species-specific occupancy and detection parameters. We illustrate this here. In our data file, there are data on three traits that characterize species: body length (cm), body mass (g), and wingspan (cm). As a first example for the modeling of species-specific effects, we will fit a community occupancy model with separate parameters for three species size groups that we define by body mass. We try to obtain three groups with roughly similar numbers of species, and take as cut points 1.3 and 2.6 for the log10 of body mass; this corresponds to about 20 and 398 g. One species with 11-kg body mass (the mute swan) is assigned to group 3, too.

# Look at distribution of body mass among 145 observed species
mass <- tapply(MHB2014$species$body.mass, MHB2014$species$specid, mean) # Get mean species mass (all 158 spec)
mass <- mass[-toss.out]                           # Only retain 145 observed species
hist(log10(mass), breaks = 40, col = "grey")      # Look at log10
gmass <- as.numeric(log10(mass) %/% 1.3 + 1)      # size groups 1, 2 and 3
gmass[gmass == 4] <- 3                            # Mute swan is group 3, too

# Bundle and summarize data set
str( win.data <- list(ysum = ysum, g = gmass, M = nrow(ysum), J = MHB2014$sites$nsurvey[1:nsite], nspec = dim(ysum)[2]) )

# Specify modle in BUGS
sink("model7.txt")
cat("
model {

# Priors: note group effects specified in this section
for (k in 1:nspec) {
  lpsi[k] ~ dnorm(mu.lpsi[g[k]], tau.lpsi[g[k]]) # note g-dependence now
  lp[k] ~ dnorm(mu.lp[g[k]], tau.lp[g[k]])
}

# Hyperpriors
for (g in 1:3) { # loop over 3 groups (g)
  mu.lpsi[g] <- logit(mu.psi[g]) # everythign is indexed g now
  mu.lp[g] <- logit(mu.p[g])
  mu.psi[g] ~ dunif(0,1)
  mu.p[g] ~ dunif(0,1)
  tau.lpsi[g] <- pow(sd.lpsi[g], -2)
  sd.lpsi[g] ~ dunif(0,5)
  tau.lp[g] <- pow(sd.lp[g], -2)
  sd.lp[g] ~ dunif(0,5)
}

# Ecological model for latent occurrence z (process model)
for (k in 1:nspec) { # no change at all down here in model
  logit(psi[k]) <- lpsi[k]
  for (i in 1:M) {
    z[i,k] ~ dbern(psi[k])
  }
}

# Observation model for observed data ysum
for (k in 1:nspec) {
  logit(p[k]) <- lp[k]
  for (i in 1:M) {
    mu.px[i,k] <- z[i,k] * p[k] # call mu.px to avoid conflict with above
    ysum[i,k] ~ dbin(mu.px[i,k], J[i])
  }
}

# Derived quantities
for (k in 1:nspec) {
  Nocc.fs[k] <- sum(z[,k])
}
for (i in 1:M) {
  Nsite[i] <- sum(z[i,])
}

}
", fill = TRUE)
sink()

# Initial values
zst <- apply(y, c(1,3), max)
zst[is.na(zst)] <- 1
inits <- function() list(z = zst)

# Parameters monitored
params <- c("mu.psi", "mu.lpsi", "sd.lpsi", "mu.p", "mu.lp", "sd.lp")

# MCMC settings
ni <- 6000 ; nt <- 2 ; nb <- 2000 ; nc <- 3

# Call JAGS from R (ART 6 min), look at convergence and summarize posteriors
out7 <- jags(win.data, inits, params, "model7.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
print(out7, dig = 3)

MCMCplot(out7, params = c("mu.psi"))

# Instead of binning body mass, treat it as a continuous covariate
# We next fit a linear regression of the hyperparameters governing species-specific values of occupancy and detection probability on the continuous mass covariate (which we log-transform). This is model 8, and it looks like this (the only change to the previous model is again in the part specifying the nature of species heterogeneity).

# Bundle and summarize data set
logmass <- as.numeric(log10(mass)) # Take log10 of body mass
str( win.data <- list(ysum = ysum, logmass = logmass, M = nrow(ysum), J = MHB2014$sites$nsurvey[1:nsite], nspec = dim(ysum)[2]) )


# Specify model in BUGS language
sink("model8.txt")
cat(" model {

# Priors
for(k in 1:nspec){ # loop over species
  lpsi[k] ~ dnorm(mu.lpsi[k], tau.lpsi[k]) # now all indexed by k, not g
  tau.lpsi[k] <- 1/var.lpsi[k]
  lp[k] ~ dnorm(mu.lp[k], tau.lp[k])
  tau.lp[k] <- 1/var.lp[k]
  mu.lpsi[k] <- delta0.lpsi + delta1.lpsi * logmass[k]
  mu.lp[k] <- delta0.lp + delta1.lp * logmass[k]
  log(var.lpsi[k]) <- phi0.lpsi + phi1.lpsi * logmass[k]
  log(var.lp[k]) <- phi0.lp + phi1.lp * logmass[k]
}

# Priors for regression params for means
delta0.lpsi ~ dnorm(0, 0.01)
delta1.lpsi ~ dnorm(0, 0.01)
delta0.lp ~ dnorm(0, 0.01)
delta1.lp ~ dnorm(0, 0.01)

# Priors for regression params for variances
phi0.lpsi ~ dnorm(0, 0.01)
phi1.lpsi ~ dnorm(0, 0.01)
phi0.lp ~ dnorm(0, 0.01)
phi1.lp ~ dnorm(0, 0.01)

# Ecological model for latent occurrence z (process model)
for(k in 1:nspec){
  logit(psi[k]) <- lpsi[k]
  for (i in 1:M) {
    z[i,k] ~ dbern(psi[k])
  }
}

# Observation model for observed data ysum
for(k in 1:nspec){ # Loop over species
  logit(p[k]) <- lp[k]
  for (i in 1:M) {
    mu.p[i,k] <- z[i,k] * p[k]
    ysum[i,k] ~ dbin(mu.p[i,k], J[i])
  }
}

# Derived quantities
for(k in 1:nspec){ # Loop over species
  Nocc.fs[k] <- sum(z[,k]) # Number of occupied sites among the 267
}
for (i in 1:M) { # Loop over sites
  Nsite[i] <- sum(z[i,]) # Number of occurring species at each site
}

}", fill = TRUE)
sink()

# Initial values
zst <- apply(y, c(1,3), max)
zst[is.na(zst)] <- 1
inits <- function() list(z = zst, delta0.lpsi = rnorm(1), delta1.lpsi = rnorm(1), delta0.lp = rnorm(1), delta1.lp = rnorm(1), phi0.lpsi = rnorm(1), phi1.lpsi = rnorm(1), phi0.lp = rnorm(1), phi1.lp = rnorm(1))

# Parameters monitored
params <- c("delta0.lpsi", "delta1.lpsi", "delta0.lp", "delta1.lp", "phi0.lpsi", "phi1.lpsi", "phi0.lp", "phi1.lp", "psi", "p", "Nocc.fs", "Nsite")

# MCMC settings
ni <- 12000 ; nt <- 2 ; nb <- 2000 ; nc <- 3

# Call JAGS from R, look at convergence and summarize posteriors
out8 <- jags(win.data, inits, params, "model8.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(out8, dig = 3)

# Get covariate values for prediction
predm <- seq(10, 10000,,500) # Predict for mass of 10g to 10 kg
pred.logm <- log10(predm)

# Compute predictions (all in one array)
tmp <- out8$sims.list # Grab simulation list
nsamp <- out8$mcmc.info$n.samples # Number of MCMC samples
pred <- array(NA, dim = c(500, nsamp, 4)) # Array for predictions

for(i in 1:nsamp){ # Fill array
  pred[,i,1] <- plogis(tmp$delta0.lpsi[i] + tmp$delta1.lpsi[i] * pred.logm)
  pred[,i,2] <- plogis(tmp$delta0.lp[i] + tmp$delta1.lp[i] * pred.logm)
  pred[,i,3] <- exp(tmp$phi0.lpsi[i] + tmp$phi1.lpsi[i] * pred.logm)
  pred[,i,4] <- exp(tmp$phi0.lp[i] + tmp$phi1.lp[i] * pred.logm)
}

# Plot posterior mean and a random sample of 100 from posterior of regression
selection <- sample(1:nsamp, 100) # Choose random sample of MCMC output
par(mfrow = c(2,2), mar = c(5,5,2,2))
matplot(predm, pred[,selection,1], ylab = "Occupancy mean", xlab = "Body mass (g)", type = "l", lty = 1, lwd = 1, col = "grey", ylim = c(0, 0.4), frame = F)
lines(predm, apply(pred[,,1], 1, mean), lwd = 3, col = "blue")
matplot(predm, pred[,selection,2], ylab = "Detection mean", xlab = "Body mass (g)", type = "l", lty = 1, lwd = 1, col = "grey", ylim = c(0, 0.8), frame = F)
lines(predm, apply(pred[,,2], 1, mean), lwd = 3, col = "blue")
matplot(predm, pred[,selection,3], ylab = "Occupancy variance", xlab = "Body mass (g)", type = "l", lty = 1, lwd = 1, col = "grey", ylim = c(0, 8), frame = F)
lines(predm, apply(pred[,,3], 1, mean), lwd = 3, col = "blue")
matplot(predm, pred[,selection,4], ylab = "Detection variance", xlab = "Body mass (g)", type = "l", lty = 1, lwd = 1, col = "grey", ylim = c(0, 8), frame = F)
lines(predm, apply(pred[,,4], 1, mean), lwd = 3, col = "blue")

