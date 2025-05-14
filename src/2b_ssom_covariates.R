####################################################################################
# A single-season occupancy model with covariates
#
# INPUT:
path_community_array = "data/cache/1_derive_community_array/community_array.rds"
path_community_survey_data = "data/cache/1_derive_community_array/community_survey_data.rds"
path_environmental_data = "data/environment/PAM_PreHarvest_Habitat_results_DD_WD_TM.xlsx"
path_unit_key = "data/unit_key.csv"
####################################################################################

library(MCMCvis)
library(ggplot2)
library(jagsUI)
library(tidyverse)
library(readxl)
library(unmarked)
library(coda)

s = "Brown Creeper" # species s
t = "2020"          # season t
threshold = 0.95

# community_array = readRDS(path_community_array)
# y = community_array[, , t, s] # Model observations for a single species from a single season
# y = ifelse(is.na(y), NA, ifelse(y >= threshold, 1, 0)) # Threshold confidence scores to obtain binary presence-absence data

community_survey_data = readRDS(path_community_survey_data)
y = community_survey_data[, , t, s] # Model observations for a single species from a single season
y = matrix( # Threshold confidence scores to obtain binary presence-absence data
  unlist(lapply(y, function(x) if (!is.null(x)) as.integer(any(x$confidence >= threshold, na.rm = TRUE)) else NA)),
  nrow = dim(y)[1],
  ncol = dim(y)[2],
  dimnames = dimnames(y)[1:2]
)

# Clean array
unit_survey_counts = rowSums(!is.na(y))
if (sum(rowSums(!is.na(y)) == 0)) {
  warning("Discarding ", sum(rowSums(!is.na(y)) == 0), " units with no survey observations...")
  y = y[rowSums(!is.na(y)) > 0, ]
}
survey_unit_counts = colSums(!is.na(y))
if (sum(survey_unit_counts == 0) > 0) {
  warning("Discarding ", sum(survey_unit_counts == 0), " surveys with no unit observations...")
  y = y[ , survey_unit_counts > 0]
}

nunit = nrow(y)   # unit i
nsurvey = ncol(y) # survey j

# Get environmental data
env_data = read_csv(path_unit_key)
env_data = env_data[env_data$season == t, ]
env_data = env_data[env_data$unit %in% rownames(y), ]
stratum = env_data %>% mutate(unit = factor(unit, levels = rownames(y))) %>% arrange(unit) %>% select(unit, stratum) # arrange to match y

# env_data_snags_deadwood = read_excel(path_environmental_data, sheet = 6, skip = 1) %>% rename(unit = 1, stratum = 2)
# env_data_snags_deadwood = env_data_snags_deadwood %>% mutate(unit = sapply(strsplit(as.character(unit), "_"), `[`, 1))
# env_data_snags_deadwood[,c('unit', 'vol_alldown_m3')]

# unmarked ---------------------------------------------
hab_levels = c('STAND INIT', 'COMP EXCL', 'THINNED', 'MATURE')
hab = factor(stratum$stratum, levels = hab_levels)

summary((
  y_unmarked = unmarkedFrameOccu(
    y = as.matrix(y),
    siteCovs = data.frame(hab = hab)
  )
))

ssom_unmarked = occu(formula = 
                     ~1      # constant detection probability across units and surveys
                     ~hab-1, # occupancy varies by habitat (no intercept, each habitat has a coefficient)
                     data = y_unmarked)
summary(ssom_unmarked)

# Predict occupancy probabilities for each habitat type
# Create a new data frame with each habitat type
new_data = data.frame(hab = factor(c("STAND INIT", "COMP EXCL", "THINNED", "MATURE"),
                                    levels = c("STAND INIT", "COMP EXCL", "THINNED", "MATURE")))
pred_psi = predict(ssom_unmarked, type = "state", newdata = new_data)
pred_psi$habitat = new_data$hab
ggplot(pred_psi, aes(x = habitat, y = Predicted)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  ylim(0.0, 1.0) +
  labs(title = paste(s, t, "(unmarked)"), x = "Habitat", y = "Occupancy probability") +
  theme_minimal()

# jags -------------------------------------------------
niter = 50000
nburnin = 7500
nchains = 2
nthin = 1

jags_file = tempfile(fileext = ".txt")
sink(jags_file)
cat("
model {
  # Uninformative prior for detection probability
  rho ~ dunif(0, 1)
  # Priors for each habitat's effect on occupancy
  for (h in 1:H) {
    beta[h] ~ dnorm(0, 1)
  }
  
  # Likelihood
  for (i in 1:I){ # For each site
    
    logit(psi[i]) <- beta[habitat[i]]  # no intercept
    z[i] ~ dbern(psi[i]) # Occupancy state model
    
    for (j in 1:J) { # For each survey
      
      y[i,j] ~ dbern(z[i] * rho) # Observation model
    }
  }
}
", fill = T)
sink()

init_func = function() {
  list( # Initial values to avoid data/model conflicts
    rho = runif(1, 0, 1),
    z = apply(y, 1, max, na.rm = T) # Initialize z[i] as 1 if a detection occurred at unit i, and 0 otherwise
  )
}

ssom_jags <- jags(
  data = list(
    y = y,
    I = nunit,
    J = nsurvey,
    H = length(unique(hab)),
    habitat = as.numeric(hab)
  ),
  inits = init_func,
  parameters.to.save = c("beta", "rho"),
  model.file = jags_file,
  n.chains = nchains,
  n.thin = nthin,
  n.iter = niter,
  n.burnin = nburnin
)

MCMCsummary(ssom_jags)

## Diagnostics
MCMCplot(ssom_jags, params = c('beta','rho'),ci=c(50,95))

# Check for convergence with trace and density plots
MCMCtrace(ssom_jags$samples, ISB = FALSE, pdf = F, exact = TRUE, post_zm = TRUE, type = 'trace', Rhat = TRUE, n.eff = TRUE)
MCMCtrace(ssom_jags$samples, ISB = FALSE, pdf = F, exact = TRUE, post_zm = TRUE, type = 'density', Rhat = TRUE, n.eff = TRUE, ind = TRUE)
gelman.diag(ssom_jags$samples) # Gelman-Rubin diagnostics

## Investigate posterior distributions for covariates

# Visualize posterior distributions for occupancy probability in relation to habitat covariates
df <- as.data.frame(plogis(ssom_jags$sims.list$beta))
colnames(df) <- hab_levels
df_long <- pivot_longer(df, cols = everything(), names_to = "Group", values_to = "Value")
ggplot(df_long, aes(x = Value, fill = Group)) +
  geom_histogram(binwidth = 0.01, alpha = 0.6, position = "identity") +
  labs(title = "Posterior distributions by habitat covariate", x = "Value", y = "Frequency") +
  theme_minimal()

# Predict the mean occupancy probability for each habitat (beta coefficient)
pred_summary = data.frame(habitat = factor(hab_levels, levels = hab_levels))
# Calculate the posterior mean probability for each beta coefficient by transforming
# each MCMC sample from log-odds to probability and averaging over the samples
pred_summary$mean  = apply(ssom_jags$sims.list$beta, 2, function(x) mean(plogis(x)))
# Compute 95% credible intervals for the posterior probability of each beta coefficient
pred_summary$lower = apply(ssom_jags$sims.list$beta, 2, function(x) quantile(plogis(x), 0.025))
pred_summary$upper = apply(ssom_jags$sims.list$beta, 2, function(x) quantile(plogis(x), 0.975))
ggplot(pred_summary, aes(x = habitat, y = mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  ylim(0.0, 1.0) +
  labs(title = paste(s, t, "(jags)"), x = "Habitat", y = "Occupancy probability") +
  theme_minimal()
