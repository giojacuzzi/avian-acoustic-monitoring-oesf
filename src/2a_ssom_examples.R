####################################################################################
# A basic single-season occupancy model with only occupancy and detection intercepts
# Model is fit with unmarked, jags, and nimble for comparison
#
# Input:
# - A 4-dimensional array with dimensions [unit × survey × season × species], with each element representing an observation (e.g. max confidence)
path_community_array = "data/cache/1_derive_community_array/community_array.rds"
#
# Output:
# - TODO

####################################################################################

library(nimble)
library(MCMCvis)
library(ggplot2)
library(unmarked)
library(jagsUI)
library(dplyr)

s = "Brown Creeper" # species s
t = "2020"          # season t

community_array = readRDS(path_community_array)
y = community_array[, , t, s] # Model observations for a single species from a single season
y = ifelse(is.na(y), NA, ifelse(y >= 0.95, 1, 0)) # Threshold confidence scores to obtain binary presence-absence data

# Clean array
unit_survey_counts = rowSums(!is.na(y))
if (sum(rowSums(!is.na(y)) == 0)) {
  warning("Discarding ", sum(rowSums(!is.na(y)) == 0), " units with no survey observations...")
  y = y[rowSums(!is.na(y)) > 0, ]
}

# DEBUG
# y = y[,1:10]
# DEBUG

nunit = nrow(y)   # unit i
nsurvey = ncol(y) # survey j

niter = 15000
nburnin = 7500
nchains = 2
nthin = 1

# unmarked ---------------------------------------------
summary((y_unmarked = unmarkedFrameOccu(y = as.matrix(y))))

ssom_unmarked = occu(formula = 
                     ~1 # detection formula
                     ~1, # occupancy formula 
                     data = y_unmarked)
summary_unmarked = round(predict(ssom_unmarked, newdata = data.frame(site = 1), type = "state"), 3) # Occupancy estimate (95% CI)
summary_unmarked = rbind(summary_unmarked, round(predict(ssom_unmarked, newdata = data.frame(site = 1), type = "det"), 3))   # Detection estimate (95% CI)
rownames(summary_unmarked) = c('psi', 'rho')

# jags -------------------------------------------------
jags_file = tempfile(fileext = ".txt")
sink(jags_file)
cat("
model {
  # Uninformative priors
  psi ~ dunif(0, 1)
  rho ~ dunif(0, 1)
  # Likelihood
  for (i in 1:I){ # For each site
    
    z[i] ~ dbern(psi) # Occupancy state model
    
    for (j in 1:J) { # For each survey
      
      y[i,j] ~ dbern(z[i] * rho) # Observation model
    }
  }
}
", fill = T)
sink()

init_func = function() {
  list( # Initial values to avoid data/model conflicts
    psi = runif(1, 0, 1),
    rho = runif(1, 0, 1),
    z = apply(y, 1, max, na.rm = T) # Initialize z[i] as 1 if a detection occurred at unit i, and 0 otherwise
  )
}

ssom_jags <- jags(
  data = list(
    y = y,
    I = nunit,
    J = nsurvey
  ),
  inits = init_func,
  parameters.to.save = c("psi", "rho"),
  model.file = jags_file,
  n.chains = nchains,
  n.thin = nthin,
  n.iter = niter,
  n.burnin = nburnin
)
summary(ssom_jags)

# nimble -----------------------------------------------
ssom_nimble_code = nimbleCode({
  # Uninformative priors
  psi ~ dunif(0, 1)
  rho ~ dunif(0, 1)
  # Likelihood
  for (i in 1:I){ # For each site
    
    z[i] ~ dbern(psi) # Occupancy state model
    
    for (j in 1:J) { # For each survey

        y[i,j] ~ dbern(z[i] * rho) # Observation model
    }
  }
})

init_func = function() {
  list( # Initial values to avoid data/model conflicts
    psi = runif(1, 0, 1),
    rho = runif(1, 0, 1),
    z = apply(y, 1, max, na.rm = T) # Initialize z[i] as 1 if a detection occurred at unit i, and 0 otherwise
  )
}

# Fit the model
ssom_nimble = nimbleMCMC(
  ssom_nimble_code,
  constants = list(
    I = nunit,
    J = nsurvey
  ),
  data = list(
    y = y
  ),
  inits = init_func,
  monitors = c("psi", "rho"),
  thin = nthin,
  niter = niter,
  nburnin = nburnin,
  nchains = nchains
)

# Summarize posteriors
MCMCsummary(ssom_nimble, digits = 2)
# MCMCtrace(ssom_nimble, exact = T, pdf = F)
# MCMCplot(ssom_nimble, params = c('psi','rho'), ci=c(50,95))

# Comparisons --------------------------------------------------------------------------------------------

summary_unmarked
MCMCsummary(ssom_jags)
MCMCsummary(ssom_nimble)
MCMCplot(ssom_jags, params = c('psi','rho'),ci=c(50,95))
MCMCplot(ssom_nimble, params = c('psi','rho'),ci=c(50,95))

# nimble, ignoring missing surveys ------------------------------------------------------------------------

y_aligned <- t(apply(y, 1, function(row) {
  non_na <- row[!is.na(row)]       # keep non-NA values
  c(non_na, rep(NA, length(row) - length(non_na)))  # pad with NAs to the right
}))
unit_effort = as.vector(rowSums(!is.na(y_shifted)))

y_aligned_zero <- y_aligned
y_aligned_zero[is.na(y_aligned_zero)] <- 0

NA_ssom_nimble_code = nimbleCode({
  psi ~ dunif(0, 1)
  rho ~ dunif(0, 1)
  for (i in 1:I){ # For each site
    z[i] ~ dbern(psi) # Occupancy state model
    for (j in 1:unit_effort[i]) { # For each survey
      y[i,j] ~ dbern(z[i] * rho) # Observation model
    }
  }
})

NA_init_func = function() {
  list( # Initial values to avoid data/model conflicts
    psi = runif(1, 0, 1),
    rho = runif(1, 0, 1),
    z = apply(y_aligned_zero, 1, max, na.rm = T) # Initialize z[i] as 1 if a detection occurred at unit i, and 0 otherwise
  )
}

NA_ssom_nimble = nimbleMCMC(
  NA_ssom_nimble_code,
  constants = list(
    I = nunit,
    unit_effort = unit_effort
  ),
  data = list(
    y = y_aligned_zero
  ),
  inits = NA_init_func,
  monitors = c("psi", "rho"),
  thin = nthin,
  niter = niter,
  nburnin = nburnin,
  nchains = nchains
)
