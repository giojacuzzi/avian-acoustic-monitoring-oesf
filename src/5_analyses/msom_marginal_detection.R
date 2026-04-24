# INPUT:
path_msom = "data/cache/models/prefinal_msom_jags_nofp_all.rds" # "data/cache/models/msom_nofp_all_2026-02-12_19:37:00.rds"
path_trait_data = "data/cache/2_traits/1_agg_traits/trait_data.csv"
##########################################################################################################

source("src/global.R")

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

match_s = stages %>% distinct() %>% arrange(stage_idx)
match_g = groups %>% select(group, group_idx) %>% distinct() %>% arrange(group_idx)
match_i = tibble(species = species, species_idx = 1:length(species))

param_beta_point_data = model_data$param_beta_data$param_beta_point_data

mu_v = msom_summary %>%
  filter(str_starts(param, "mu.v")) %>%
  mutate(
    s = as.integer(str_match(param, "mu\\.v\\[(\\d+),(\\d+)\\]")[,2]),
    g = as.integer(str_match(param, "mu\\.v\\[(\\d+),(\\d+)\\]")[,3])
  ) %>%
  left_join(match_s, by = c("s" = "stage_idx"))
ggplot(mu_v, aes(x = mean, y = param, color = stratum_4)) +
  geom_vline(xintercept = 0) + geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0.1, position = position_dodge(width = 0.5))

v = msom_summary %>%
  filter(str_starts(param, "v")) %>%
  mutate(
    s = as.integer(str_match(param, "v\\[(\\d+),(\\d+)\\]")[,2]),
    i = as.integer(str_match(param, "v\\[(\\d+),(\\d+)\\]")[,3])
  ) %>%
  left_join(match_s, by = c("s" = "stage_idx")) %>% left_join(match_i, by = c("i" = "species_idx"))

mu_beta_point = msom_summary %>%
  filter(str_starts(param, "mu.beta_point")) %>%
  mutate(param = str_extract(param, "beta_point\\d+")) %>%
  mutate(
    p = as.integer(str_match(param, "mu\\.beta_point(\\d+)\\[(\\d+)\\]")[,2]),
    g = as.integer(str_match(param, "mu\\.beta_point(\\d+)\\[(\\d+)\\]")[,3])
  ) %>% left_join(param_beta_point_data %>% select(param, name), by = c("param")) %>% left_join(match_g, by = c("g" = "group_idx"))

beta_point = msom_summary %>% mutate(full_param = param) %>%
  filter(str_starts(full_param, "beta_point")) %>%
  mutate(param = str_extract(full_param, "^[^\\[]+")) %>%
  mutate(
    p = as.integer(str_match(full_param, "beta_point(\\d+)\\[(\\d+)\\]")[,2]),
    i = as.integer(str_match(full_param, "beta_point(\\d+)\\[(\\d+)\\]")[,3])
  ) %>% left_join(param_beta_point_data %>% select(param, name), by = c("param")) %>% left_join(match_i, by = c("i" = "species_idx")) %>% left_join(groups, by = c("species" = "common_name"))

beta_point %>% filter(param %in% c("beta_point3", "beta_point3_sq"))
param_beta_point_data

#############
# Species-specific

# ============================================================
# Marginal detection probability vs. yday (quadratic)
# ============================================================

library(tidyverse)

# ---- 1. Recover yday scaling parameters ------------------
yday_mat      <- param_beta_point_data$data[[3]]
yday_mean_raw <- attr(yday_mat, "scaled:center")  # 159.2048
yday_sd_raw   <- attr(yday_mat, "scaled:scale")   # 30.22745
yday_std_vals <- yday_mat[is.finite(yday_mat)]

# ---- 2. Pull posterior samples ---------------------------
beta3_sims    <- msom$sims.list$beta_point3     # [n_iter, n_species]
beta3_sq_sims <- msom$sims.list$beta_point3_sq  # [n_iter, n_species]
v_sims        <- msom$sims.list$v               # [n_iter, n_stages, n_species]

# Marginal intercept: average v over stages (equal weight)
v_marg_sims <- apply(v_sims, c(1, 3), mean)    # [n_iter, n_species]

# ---- 3. Prediction grid ----------------------------------
n_grid       <- 200
yday_std_seq <- seq(min(yday_std_vals), max(yday_std_vals), length.out = n_grid)
yday_raw_seq <- yday_std_seq * yday_sd_raw + yday_mean_raw

# ---- 4. Posterior predictive summaries per species -------
pred_df <- map_dfr(seq_len(n_species), function(i) {
  
  # logit(p) posterior matrix: [n_iter x n_grid]
  logit_p <-
    v_marg_sims[, i]    %o% rep(1, n_grid) +   # intercept
    beta3_sims[, i]     %o% yday_std_seq    +   # linear yday
    beta3_sq_sims[, i]  %o% yday_std_seq^2      # quadratic yday
  
  p_mat <- plogis(logit_p)  # [n_iter x n_grid]
  
  tibble(
    species_idx = i,
    yday_std    = yday_std_seq,
    yday        = yday_raw_seq,
    p_mean      = colMeans(p_mat),
    p_lo        = apply(p_mat, 2, quantile, 0.025),
    p_hi        = apply(p_mat, 2, quantile, 0.975)
  )
}) %>%
  left_join(match_i, by = "species_idx") %>%
  left_join(groups %>% select(common_name, group, grouping),
            by = c("species" = "common_name"))

# ---- 5. Plot: facet by group, one line per species -------
ggplot(pred_df, aes(x = yday, group = species)) +
  geom_ribbon(aes(ymin = p_lo, ymax = p_hi, fill = species),
              alpha = 0.15, show.legend = FALSE) +
  geom_line(aes(y = p_mean, color = species),
            linewidth = 0.7, show.legend = FALSE) +
  facet_wrap(~ group, scales = "free_y") +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(
    breaks = c(91, 121, 152, 182, 213),
    labels = c("Apr", "May", "Jun", "Jul", "Aug")
  ) +
  labs(
    x        = "Day of year",
    y        = "Detection probability p(detection | present)",
    title    = "Marginal species-specific detection probability vs. day of year",
    subtitle = "Posterior mean ± 95% CrI; marginalised over stage (equal weights)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    panel.grid.minor = element_blank()
  )

# ---- Optional: one panel per species ---------------------
ggplot(pred_df, aes(x = yday)) +
  geom_ribbon(aes(ymin = p_lo, ymax = p_hi),
              fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = p_mean),
            color = "steelblue", linewidth = 0.8) +
  facet_wrap(~ species, ncol = 8) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(
    breaks = c(121, 182),
    labels = c("May", "Jul")
  ) +
  labs(
    x     = NULL,
    y     = "p(detection | present)",
    title = "Species-specific marginal detection probability vs. yday"
  ) +
  theme_bw(base_size = 9) +
  theme(
    strip.text = element_text(size = 7),
    axis.text  = element_text(size = 6)
  )
