# 4_community_comp.R #############################################################################
# Community composition analyses, propagating posterior uncertainty
#
# INPUT:
path_msom = "data/cache/models/V4_msom_V4_nofp_nofp_all.rds" # "data/cache/models/msom_nofp_all_2026-02-12_19:37:00.rds"
path_trait_data = "data/cache/2_traits/1_agg_traits/trait_data.csv"
path_occurrence_predictor_plot_data = "data/cache/4_msom/1_assemble_msom_data/V3_occurrence_predictor_plot_data.rds"
##################################################################################################

source("src/global.R")

# Load data --------------------------------------------------

message("Loading data for multi-species occupancy model ", path_msom)
model_data = readRDS(path_msom)

msom_summary = model_data$msom_summary
groups  = model_data$groups
sites   = model_data$sites
species = model_data$species
seasons = model_data$seasons
stages  = factor(model_data$stages$stratum_4,  levels = c("standinit", "compex", "thin", "mature"))

z_raw = model_data$msom$sims.list$z
str(z_raw)
n_iter    = dim(z_raw)[1]
n_sites   = dim(z_raw)[2]
n_seasons = dim(z_raw)[3]
n_species = dim(z_raw)[4]

rm(model_data)
gc()

message("Loading species trait data from ", path_trait_data)
species_traits = read_csv(path_trait_data, show_col_types = FALSE) %>% filter(common_name %in% species)

message("Loading occurrence predictor plot scale data from ", path_occurrence_predictor_plot_data)
occurrence_predictor_plot_data = readRDS(path_occurrence_predictor_plot_data)[[1]] %>% filter(site %in% sites)
stopifnot(all(occurrence_predictor_plot_data$site == sites))

# "We retained posterior distributions of latent occurrence states averaged across years to propagate uncertainty through analyses of community composition."
# Average over years (dim 3) while retaining iterations, sites, species, using accumulation to avoid memory limits
z_tmean = array(0, dim = c(n_iter, n_sites, n_species))
for (i in seq_len(n_species)) {
  acc = matrix(0, nrow = n_iter, ncol = n_sites)
  for (t in seq_len(n_seasons)) {
    acc = acc + z_raw[, , t, i]
  }
  z_tmean[, , i] = acc / n_seasons
}
str(z_tmean)

# Posterior estimated richness (averaged across years) ----------------------------------
{
  # Posterior estimated richness
  post_tmean_rich = matrix(0, nrow = n_iter, ncol = n_sites)
  for (i in seq_len(n_species)) {
    post_tmean_rich = post_tmean_rich + z_tmean[, , i]
  }

  # Posterior mean estimated richness, for summary and visualization
  post_tmean_rich_stats = data.frame(
    site      = sites,
    stage     = stages,
    mean      = apply(post_tmean_rich, 2, mean),
    bci_2.5   = apply(post_tmean_rich, 2, quantile, probs = 0.025),
    bci_97.5  = apply(post_tmean_rich, 2, quantile, probs = 0.975)
  )
  ggplot(post_tmean_rich_stats, aes(x = stage, y = mean, color = stage)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.25, alpha = 0.35, aes(color = stage)) +
    scale_color_manual(values = unname(stage_colors))

  # ANOVA + Tukey on every posterior iteration
  anova_pvals = numeric(n_iter)
  tukey_post  = vector("list", n_iter)
  for (i in seq_len(n_iter)) {
    df_i = data.frame(
      richness = post_tmean_rich[i, ],
      stage    = stages
    )
    fit             = aov(richness ~ stage, data = df_i)
    anova_pvals[i]  = summary(fit)[[1]]["stage", "Pr(>F)"]
    tukey_post[[i]] = TukeyHSD(fit)$stage
    if (i %% 1000 == 0) message(i, "/", n_iter)
  }
  
  # Summarize Tukey posterior
  tukey_arr  = aperm(simplify2array(tukey_post), c(3, 1, 2))
  pair_names = rownames(tukey_post[[1]])
  
  contrasts = lapply(seq_along(pair_names), function(p) {
    diffs = tukey_arr[, p, 1]
    pvals = tukey_arr[, p, 4]
    data.frame(
      pair      = pair_names[p],
      mean      = mean(diffs),
      bci_2.5   = quantile(diffs, 0.025),
      bci_97.5  = quantile(diffs, 0.975),
      prob_gt0  = mean(diffs > 0),   # Probability that stage A is > stage B
      prop_sig  = mean(pvals < 0.05) # Proportion of draws with significant Tukey tests
    )
  }) %>% bind_rows()
  
  print(contrasts, digits = 3)
  message("Proportion of posterior iterations with significant ANOVA: ", mean(anova_pvals < 0.05))
  
  # TAKEAWAYS:
  # - Stand initiation has the highest richness (8-10 more species than closed-canopy stages, BCI 6.5-11.4)
  # - Compex has the lowest richness (1-2 fewer species on average than thinned and mature)
  # - Differences in richness are clear between thin and compex, but weaker between thin and mature
}

# TODO: Posterior estimated richness (each year)

# Principal coordinates analysis (PCoA) -----------------------------------------------------------
{
  # Average posterior draws of Bray-curtis dissimilarity matrix
  d_sum = matrix(0, n_sites, n_sites)
  for (i in seq_len(n_iter)) {
    if (i %% 1000 == 0) message(i, "/", n_iter)
    d_sum = d_sum + as.matrix(vegdist(z_tmean[i, , ], method = "bray"))
  }
  d_mean_post = as.dist(d_sum / n_iter)

  # PCoA
  pcoa_res = cmdscale(d_mean_post, eig = TRUE, k = 2)
  pct  = round(pcoa_res$eig / sum(pcoa_res$eig[pcoa_res$eig > 0]) * 100, 1)
  pcoa_df = as.data.frame(pcoa_res$points) %>% setNames(c("PCoA1", "PCoA2")) %>% mutate(site = sites, stage = stages)
  fig_PCoA_sites = ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = stage)) +
    geom_point(size = 2, alpha = 0.8) +
    stat_ellipse(aes(fill = stage), geom = "polygon", alpha = 0.1, level = 0.95) +
    scale_color_manual(values = unname(stage_colors)) +
    scale_fill_manual(values = unname(stage_colors)) +
    labs(x = paste0("PCoA1 (", pct[1], "%)"), y = paste0("PCoA2 (", pct[2], "%)")) +
    theme_bw(); print(fig_PCoA_sites)
  
  { # Species ordination correlations
    species_fit = envfit(pcoa_res$points, z_tmean, permutations = 999)
    species_scores = as.data.frame(species_fit$vectors$arrows) |>
      setNames(c("PCoA1", "PCoA2")) |>
      mutate(species = species, r2 = species_fit$vectors$r, p = species_fit$vectors$pvals) |>
      arrange(desc(r2))
    top_species = species_scores %>% filter(p < 0.05) %>% filter(r2 > 0.0)
    arrow_scale = 0.6 * max(abs(pcoa_df[, c("PCoA1", "PCoA2")])) / max(abs(top_species[, c("PCoA1", "PCoA2")]))
    fig_PCoA_species = ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = stage)) +
      stat_ellipse(aes(fill = stage), geom = "polygon", alpha = 0.05, level = 0.95) +
      geom_segment(data = top_species,
                   aes(x = 0, y = 0, xend = PCoA1 * r2 * arrow_scale, yend = PCoA2 * r2 * arrow_scale),
                   inherit.aes = FALSE, arrow = arrow(length = unit(0.2, "cm")), color = "black", alpha = 0.5) +
      geom_text(data = top_species,
                aes(x = PCoA1 * r2 * arrow_scale, y = PCoA2 * r2 * arrow_scale, label = species),
                inherit.aes = FALSE, color = "black", alpha = 0.5, size = 3, hjust = 0.5, vjust = -0.5) +
      scale_color_manual(values = alpha(unname(stage_colors), 0.5)) +
      scale_fill_manual(values = unname(stage_colors)) +
      labs(x = paste0("PCoA1 (", pct[1], "%)"), y = paste0("PCoA2 (", pct[2], "%)"), color = "Stratum", fill = "Stratum") +
      theme_bw(); print(fig_PCoA_species)
    print(species_scores)
    }

  { # Environmental ordination correlations
    env_data = occurrence_predictor_plot_data %>%
      select(-c(site, ces, treatment, wadnr_patch_stratum,
                thinning_status, thinning_treatment,
                stage_3, stage_4, stratum_4, stratum_5, scale),
             -ends_with("_median"), -ends_with("_sd"),
             -starts_with("sdi"), -starts_with("tree_acre_4_"), -starts_with("tree_acre_30_"), -starts_with("qmd_4_"),
             -starts_with("qmd_t100_"), -starts_with("ht"), -starts_with("pcnt_"), -starts_with("ba_"))
    env_fit = envfit(pcoa_res$points, env_data, permutations = 999, na.rm = TRUE)
    env_scores = as.data.frame(env_fit$vectors$arrows) %>%
      setNames(c("PCoA1", "PCoA2")) %>%
      mutate(variable = rownames(.), r2 = env_fit$vectors$r, p = env_fit$vectors$pvals) %>%
      arrange(desc(r2))
    top_env = env_scores %>% filter(p < 0.05) %>% filter(r2 > 0.0)
    arrow_scale_env = 0.6 * max(abs(pcoa_df[, c("PCoA1", "PCoA2")])) / max(abs(top_env[, c("PCoA1", "PCoA2")]))
    fig_PCoA_env = ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = stage)) +
      stat_ellipse(aes(fill = stage), geom = "polygon", alpha = 0.05, level = 0.95) +
      geom_segment(data = top_env, aes(x = 0, y = 0, xend = PCoA1 * r2 * arrow_scale_env, yend = PCoA2 * r2 * arrow_scale_env),
                   inherit.aes = FALSE, arrow = arrow(length = unit(0.2, "cm")), color = "black", alpha = 0.5) +
      geom_text(data = top_env, aes(x = PCoA1 * r2 * arrow_scale_env, y = PCoA2 * r2 * arrow_scale_env, label = variable),
                inherit.aes = FALSE, color = "black", alpha = 0.5, size = 3, hjust = 0.5, vjust = -0.5) +
      scale_color_manual(values = alpha(unname(stage_colors), 0.5)) +
      scale_fill_manual(values = unname(stage_colors)) +
      labs(x = paste0("PCoA1 (", pct[1], "%)"), y = paste0("PCoA2 (", pct[2], "%)"), color = "Stratum", fill = "Stratum") +
      theme_bw(); print(fig_PCoA_env)
    print(env_scores)
  }
  
  # Findings:
  # - Strong primary environmental filter is a canopy closure / disturbance gradient (PCoA axis 1)
  # - Weaker secondary elevational gradient (PCoA axis 2)
  # - Example PCoA1 species ordination correlations: brown creeper + (mature trees, closed), song sparrow - (edge, open habitat)
  # - Example PCoA2 species ordination correlations: gray jay + (montane), kingfisher/bald eagle - (valley bottom / riparian / aquatic)
}

# Community composition -----------------------------------------------------------
{
  ## Overall differences
  
  R2_post = numeric(n_iter) # Proportion of variance explained by stage
  F_post  = numeric(n_iter) # Larger pseudo-F indicates greater difference between groups relative to variation within groups
  p_post  = numeric(n_iter) # Significant frequentist difference exists
  
  for (i in seq_len(n_iter)) {
    if (i %% 1000 == 0) message(i, "/", n_iter)
    z_i        = z_tmean[i, , ]
    d          = vegdist(z_i, method = "bray")
    fit        = adonis2(d ~ stages, permutations = 50) # TODO: 999 permutations?
    R2_post[i] = fit$R2[1]
    F_post[i]  = fit$F[1]
    p_post[i]  = fit$`Pr(>F)`[1]
  }
  
  res_adonis2 = data.frame(
    stat      = c("R2", "F", "p"),
    mean      = c(mean(R2_post),             mean(F_post),             mean(p_post)),
    median    = c(median(R2_post),           median(F_post),           median(p_post)),
    bci_2.5   = c(quantile(R2_post, 0.025),  quantile(F_post, 0.025),  quantile(p_post, 0.025)),
    bci_97.5  = c(quantile(R2_post, 0.975),  quantile(F_post, 0.975),  quantile(p_post, 0.975))
  )
  print(res_adonis2, digits = 1)
  
  message("Proportion of posterior iterations with significant PERMANOVA: ", mean(p_post < 0.05))
}
