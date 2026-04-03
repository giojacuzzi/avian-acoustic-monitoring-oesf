####################################################################################
# Marginal SR and FDis prediction
#
# CONFIG:
path_msom = "data/cache/models/V4_msom_V4_nofp_nofp_all.rds"
#
# INPUT:
path_trait_data = "data/cache/2_traits/1_agg_traits/trait_data.csv"
####################################################################################

source("src/global.R")

# Load data for multi-species occupancy model --------------------------------------------------

message("Loading data for multi-species occupancy model ", path_msom)
model_data = readRDS(path_msom)

(msom_summary = model_data$msom_summary)
(msom = model_data$msom)
(groups = model_data$groups %>% arrange(common_name))
(sites = model_data$sites)
(species = model_data$species)
(stages = model_data$stages)
msom_results = model_data

message("Loading species trait data from ", path_trait_data)
species_traits = read_csv(path_trait_data, show_col_types = FALSE) %>% filter(common_name %in% species)

# ── Extract species-specific homerange scaling parameters ─────────────────────
# param_alpha_homerange_data$data[[a]] is a named list (length I),
# each element is a scaled vector with attributes "scaled:center" and "scaled:scale"

hr_data <- msom_results$param_alpha_data$param_alpha_homerange_data
# hr_data$name should be c("pcnt_standinit", "pcnt_thin", "pcnt_mature")
# confirming the mapping:
#   alpha_homerange1 = pcnt_standinit
#   alpha_homerange2 = pcnt_thin
#   alpha_homerange3 = pcnt_mature
#   compex is reference (pcnt_compex = 1 - sum of others, not modelled directly)

sims   <- msom$sims.list
n_iter <- dim(sims$u)[1]   # 6000
I      <- dim(sims$u)[2]   # 67 species

# For each homerange variable and species, extract center (mean) and scale (SD)
# used during standardisation — shape: [3 variables × I species]
hr_center <- matrix(NA_real_, nrow = 3, ncol = I,
                    dimnames = list(hr_data$name, msom_results$species))
hr_scale  <- matrix(NA_real_, nrow = 3, ncol = I,
                    dimnames = list(hr_data$name, msom_results$species))

for (a in 1:3) {
  for (i in seq_along(msom_results$species)) {
    sp <- msom_results$species[i]
    scaled_vec <- hr_data$data[[a]][[sp]]
    hr_center[a, i] <- attr(scaled_vec, "scaled:center")
    hr_scale[a, i]  <- attr(scaled_vec, "scaled:scale")
  }
}

# Helper: standardise a raw proportion value for each species using their
# individual scaling params. Returns a vector of length I.
std_for_species <- function(raw_value, var_idx) {
  (raw_value - hr_center[var_idx, ]) / hr_scale[var_idx, ]
}

# ── Define scenarios as raw (unstandardised) proportions ──────────────────────
# Each scenario: 100% of one stand type → other homerange proportions = 0%
# compex: pcnt_standinit=0, pcnt_thin=0, pcnt_mature=0  (all reference)
# Note: x_alpha_plot covariates = 0 (global standardised mean) in all scenarios,
# so stage index has no effect and is omitted.

scenarios <- list(
  compex    = list(pcnt_standinit = 0, pcnt_thin = 0, pcnt_mature = 0),
  standinit = list(pcnt_standinit = 1, pcnt_thin = 0, pcnt_mature = 0),
  thin      = list(pcnt_standinit = 0, pcnt_thin = 1, pcnt_mature = 0),
  mature    = list(pcnt_standinit = 0, pcnt_thin = 0, pcnt_mature = 1)
)

richness_pred <- matrix(NA_real_,
                        nrow     = n_iter,
                        ncol     = length(scenarios),
                        dimnames = list(NULL, names(scenarios)))

for (nm in names(scenarios)) {
  sc <- scenarios[[nm]]
  
  # Standardised homerange values — one value per species [length I]
  # Recycled across iterations via matrix multiplication below
  z_hr1 <- std_for_species(sc$pcnt_standinit, 1)  # [I]
  z_hr2 <- std_for_species(sc$pcnt_thin,      2)  # [I]
  z_hr3 <- std_for_species(sc$pcnt_mature,    3)  # [I]
  
  # logit(psi): [n_iter × I]
  # sims$alpha_homerangeN is [n_iter × I]; z_hrN is [I] → broadcast via sweep
  logit_psi <-
    sims$u +
    sweep(sims$alpha_homerange1, 2, z_hr1, `*`) +
    sweep(sims$alpha_homerange2, 2, z_hr2, `*`) +
    sweep(sims$alpha_homerange3, 2, z_hr3, `*`)
  # alpha_plot terms drop out (x_alpha_plot1/2/3 = 0, global mean)
  # alpha_point, alpha_season terms drop out (= 0, global mean)
  
  richness_pred[, nm] <- rowSums(plogis(logit_psi))
}

# ── Posterior summary ──────────────────────────────────────────────────────────
richness_summary <- data.frame(
  scenario = names(scenarios),
  mean  = colMeans(richness_pred),
  sd    = apply(richness_pred, 2L, sd),
  q2.5  = apply(richness_pred, 2L, quantile, 0.025),
  q50   = apply(richness_pred, 2L, quantile, 0.500),
  q97.5 = apply(richness_pred, 2L, quantile, 0.975),
  row.names = NULL
)

print(richness_summary, digits = 3)

# ── Plot ───────────────────────────────────────────────────────────────────────
library(ggplot2)

# Build long data frame directly — avoids melt() label scrambling
richness_long <- data.frame(
  scenario  = rep(colnames(richness_pred), each = n_iter),
  richness  = as.vector(richness_pred)
)

richness_long$scenario <- factor(richness_long$scenario,
                                 levels = c("standinit", "thin", "mature", "compex"))

# Quick sanity check — medians should match the table
tapply(richness_long$richness, richness_long$scenario, median)

ggplot(richness_long, aes(x = scenario, y = richness, fill = scenario)) +
  stat_summary(
    fun.data = \(x) data.frame(y    = median(x),
                               ymin = quantile(x, 0.025),
                               ymax = quantile(x, 0.975)),
    geom = "crossbar", width = 0.4
  ) +
  labs(
    title    = "Predicted marginal species richness by stand type",
    subtitle = "Homerange = 100% focal type; all other covariates at mean (0)",
    x = "Stand type", y = "Expected species richness"
  ) +
  theme_bw() + theme(legend.position = "none")

###################

# Create species-by-traits matrix
trait_matrix = species_traits %>% as.data.frame() %>%
  select(group_nest_ps,
         group_forage_substrate, # TODO: Consider more traits?
         group_diet,
         group_migrant,
         mass)
rownames(trait_matrix) = species_traits$common_name
all(rownames(trait_matrix) == species)  # must be TRUE

library(FD)

# ── 1. Fixed Gower distance matrix ────────────────────────────────────────────
trait_dist <- gowdis(trait_matrix)

# ── 2. Build all psi matrices first, rbind into one block ─────────────────────
# Single fdisp() call is faster than one per scenario
psi_list <- lapply(scenarios, function(sc) {
  z_hr1 <- std_for_species(sc$pcnt_standinit, 1)
  z_hr2 <- std_for_species(sc$pcnt_thin,      2)
  z_hr3 <- std_for_species(sc$pcnt_mature,    3)
  
  psi_mat <- plogis(
    sims$u +
      sweep(sims$alpha_homerange1, 2, z_hr1, `*`) +
      sweep(sims$alpha_homerange2, 2, z_hr2, `*`) +
      sweep(sims$alpha_homerange3, 2, z_hr3, `*`)
  )
  colnames(psi_mat) <- species
  rownames(psi_mat) <- paste0(nm, "_iter_", seq_len(n_iter))
  psi_mat
})

psi_all <- do.call(rbind, psi_list)  # [n_iter * 4 × I]

# ── 3. Single fdisp() call ────────────────────────────────────────────────────
fd_res <- fdisp(trait_dist, psi_all)

# Split result back into [n_iter × scenarios] matrix
fdis_pred <- matrix(fd_res$FDis,
                    nrow     = n_iter,
                    ncol     = length(scenarios),
                    dimnames = list(NULL, names(scenarios)))

# ── 4. Posterior summary ───────────────────────────────────────────────────────
fdis_summary <- data.frame(
  scenario = names(scenarios),
  mean  = colMeans(fdis_pred),
  sd    = apply(fdis_pred, 2L, sd),
  q2.5  = apply(fdis_pred, 2L, quantile, 0.025),
  q50   = apply(fdis_pred, 2L, quantile, 0.500),
  q97.5 = apply(fdis_pred, 2L, quantile, 0.975),
  row.names = NULL
)

print(fdis_summary, digits = 3)

# ── 5. Plot ────────────────────────────────────────────────────────────────────
fdis_long <- data.frame(
  scenario = rep(colnames(fdis_pred), each = n_iter),
  fdis     = as.vector(fdis_pred)
)
fdis_long$scenario <- factor(fdis_long$scenario,
                             levels = c("standinit", "thin", "mature", "compex"))

ggplot(fdis_long, aes(x = scenario, y = fdis, fill = scenario)) +
  stat_summary(
    fun.data = \(x) data.frame(y    = median(x),
                               ymin = quantile(x, 0.025),
                               ymax = quantile(x, 0.975)),
    geom = "crossbar", width = 0.4
  ) +
  labs(
    title    = "Posterior functional dispersion (FDis) by stand type",
    subtitle = "Species weighted by posterior occupancy probability; all other covariates at mean (0)",
    x = "Stand type", y = "FDis"
  ) +
  theme_bw() + theme(legend.position = "none")
