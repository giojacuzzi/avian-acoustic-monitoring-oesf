# INPUT:
path_msom = "data/cache/models/V4_msom_V4_nofp_nofp_all.rds" # "data/cache/models/msom_nofp_all_2026-02-12_19:37:00.rds"
path_trait_data = "data/cache/2_traits/1_agg_traits/trait_data.csv"
##########################################################################################################

source("src/global.R")

strata_cols = c(
  "standinit" = "orange",
  "compex"  = "forestgreen",
  "thin"    = "purple",
  "mature"     = "tan4"
)

# Load data for multi-species occupancy model --------------------------------------------------

message("Loading data for multi-species occupancy model ", path_msom)
model_data = readRDS(path_msom)

msom_summary = model_data$msom_summary
msom = model_data$msom
groups = model_data$groups %>% arrange(common_name)
sites = model_data$sites
species = model_data$species
seasons = model_data$seasons
stages = model_data$stages
strata = factor(stages$stratum_4, levels = c("standinit", "compex", "thin", "mature"))

message("Loading species trait data from ", path_trait_data)
species_traits = read_csv(path_trait_data, show_col_types = FALSE) %>% filter(common_name %in% species)

z = msom$sims.list$z
str(z)
n_iter    = dim(z)[1]
n_sites   = dim(z)[2]
n_seasons = dim(z)[3]
n_species = dim(z)[4]


compute_scenario_occupancy = function(species_to_plot, model_data) {
  
  msom    = model_data$msom
  species = model_data$species
  seasons = model_data$seasons
  stages  = model_data$stages
  
  param_alpha_plot_data      = model_data$param_alpha_data$param_alpha_plot_data
  param_alpha_homerange_data = model_data$param_alpha_data$param_alpha_homerange_data
  
  J       = length(model_data$sites)
  Tseason = length(seasons)
  sims    = msom$sims.list
  n_iter  = dim(sims$u)[1]
  
  stage_idx = stages$stage_idx
  
  scenarios = list(
    standinit = list(stage = 2L, hr = c(1, 0, 0)),
    compex    = list(stage = 1L, hr = c(0, 0, 0)),
    thin      = list(stage = 4L, hr = c(0, 1, 0)),
    mature    = list(stage = 3L, hr = c(0, 0, 1))
  )
  
  n_plot = nrow(param_alpha_plot_data)
  
  stage_plot_means = lapply(seq_len(n_plot), function(p) {
    mat = matrix(param_alpha_plot_data$data[[p]], nrow = J, ncol = Tseason)
    vapply(1:4, function(s) mean(mat[stage_idx == s, ], na.rm = TRUE), numeric(1))
  })
  
  # Accumulate per-species results and per-scenario [n_iter] richness vectors
  species_rows = list()
  # richness_iters[[scen]] will be an [n_iter] vector accumulated by summing cum_psi
  richness_iters = setNames(
    lapply(names(scenarios), function(.) rep(0, n_iter)),
    names(scenarios)
  )
  
  for (sp in species_to_plot) {
    
    sp_i = which(species == sp)
    if (length(sp_i) == 0) { warning("Species not found: ", sp); next }
    
    for (scen_name in names(scenarios)) {
      
      scen = scenarios[[scen_name]]
      s    = scen$stage
      
      x_plot_scen = vapply(stage_plot_means, `[`, numeric(1), s)
      
      x_hr_scen = vapply(seq_len(3L), function(p) {
        vec = param_alpha_homerange_data$data[[p]][[sp]]
        (scen$hr[p] - attr(vec, "scaled:center")) / attr(vec, "scaled:scale")
      }, numeric(1))
      
      # logit(psi): [n_iter] vector — season fixed at mean (= 0), so omitted
      logit_psi = sims$u[, sp_i] +
        sims$alpha_plot1[, s, sp_i] * x_plot_scen[1] +
        sims$alpha_plot2[, s, sp_i] * x_plot_scen[2] +
        sims$alpha_plot3[, s, sp_i] * x_plot_scen[3] +
        sims$alpha_homerange1[, sp_i] * x_hr_scen[1] +
        sims$alpha_homerange2[, sp_i] * x_hr_scen[2] +
        sims$alpha_homerange3[, sp_i] * x_hr_scen[3]
      # alpha_point* = 0 (standardized mean) → omitted
      # alpha_season * 0                     → omitted
      
      cum_psi = plogis(logit_psi)   # [n_iter]; single-season psi = cumulative psi
      
      species_rows = c(species_rows, list(tibble(
        species  = sp,
        scenario = scen_name,
        mean     = mean(cum_psi),
        lower95  = quantile(cum_psi, 0.025),
        upper95  = quantile(cum_psi, 0.975)
      )))
      
      # Accumulate expected richness: E[richness] = sum_i E[z_i]
      richness_iters[[scen_name]] = richness_iters[[scen_name]] + cum_psi
    }
  }
  
  species_df = bind_rows(species_rows) %>%
    mutate(scenario = factor(scenario, levels = names(scenarios)))
  
  richness_df = imap_dfr(richness_iters, function(rich_vec, scen_name) {
    tibble(
      scenario = scen_name,
      mean     = mean(rich_vec),
      lower95  = quantile(rich_vec, 0.025),
      upper95  = quantile(rich_vec, 0.975)
    )
  }) %>% mutate(scenario = factor(scenario, levels = names(scenarios)))
  
  list(species = species_df, richness = richness_df)
}

focal_species = sort(intersect(species, c(
  "vaux's swift", "marbled murrelet", "rufous hummingbird", "northern goshawk", "common nighthawk", "evening grosbeak", "golden-crowned kinglet", "olive-sided flycatcher", "pine siskin", "pileated woodpecker", "willow flycatcher"
)))
old_associates = sort(intersect(species, c(
  "vaux's swift", "marbled murrelet", "pacific-slope flycatcher", "hammond's flycatcher", "hermit warbler", "varied thrush", "hermit thrush", "townsend's warbler", "golden-crowned kinglet", "brown creeper", "chestnut-backed chickadee", "pacific wren", "red-breasted sapsucker", "northern goshawk", "pileated woodpecker"
)))
early_associates = sort(intersect(species, c(
  "rufous hummingbird", "willow flycatcher", "orange-crowned warbler", "lazuli bunting", "violet-green swallow", "black-throated gray warbler", "spotted towhee", "mountain quail", "olive-sided flycatcher", "macgillivray's warbler", "wilson's warbler", "wrentit"
)))

# ===============================================================================

results = compute_scenario_occupancy(focal_species, model_data)

# Species-specific marginal expected habitat use probability
ggplot(results$species, aes(x = scenario, y = mean, ymin = lower95, ymax = upper95, color = scenario)) +
  geom_linerange(linewidth = 0.8) + geom_point(size = 3) +
  scale_color_manual(values = strata_cols, guide = "none") +
  facet_wrap(~ species) +
  labs(x = "Stage scenario", y = "Marginal expected use probability")
# FINDINGS:
# - Standinit and mature promote the maximum probability of use across focal species, while compex does not enhance probability for any species.
# - Mature may function as long-term marginal habitat for RuHm, PnSk, OlSiFly
# - Thinned can substantially promote use for the OlSiFly
# - GCK and EG show no preference for any stage

# Marginal expected focal species richness
ggplot(results$richness, aes(x = scenario, y = mean, ymin = lower95, ymax = upper95, color = scenario)) +
  geom_linerange(linewidth = 0.8) + geom_point(size = 3) +
  scale_color_manual(values = strata_cols, guide = "none") +
  labs(x = "Stage scenario", y = "Marginal expected focal species richness")
# FINDINGS:
# - Highest richness of priority species stand init > mature > thin > compex
# - Competitive exclusion stage supports the fewest number of priority species
