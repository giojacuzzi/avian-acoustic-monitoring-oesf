# 4_community_comp.R #############################################################################
# Community composition analyses, propagating posterior uncertainty
#
# CONFIG:
thin = 10 # Number of iterations to thin from the z_tmean posterior to speed up computation (0 or more)
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

if (thin > 0) {
  message("Thinning z_raw posterior to every ", thin, " iterations")
  thin_idx = seq(1, n_iter, by = thin)
  z_raw  = z_raw[thin_idx, , , ]
  n_iter   = length(thin_idx)
  str(z_raw)
}

rm(model_data)
gc()

message("Loading species trait data from ", path_trait_data)
species_traits = read_csv(path_trait_data, show_col_types = FALSE) %>% filter(common_name %in% species)

message("Loading occurrence predictor plot scale data from ", path_occurrence_predictor_plot_data)
occurrence_predictor_plot_data = readRDS(path_occurrence_predictor_plot_data)[[1]] %>% filter(site %in% sites)
stopifnot(all(occurrence_predictor_plot_data$site == sites))

# "We retained posterior distributions of latent occurrence states averaged across years to propagate uncertainty through analyses of community composition."
# Average over years (dim 3) while retaining iterations, sites, species, using accumulation to avoid memory limits
message("Averaging posterior occurrence states across years")
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
  pb = progress_bar$new(format = progress_bar_format, total = n_iter, clear = FALSE)
  for (i in seq_len(n_iter)) {
    df_i = data.frame(
      richness = post_tmean_rich[i, ],
      stage    = stages
    )
    fit             = aov(richness ~ stage, data = df_i)
    anova_pvals[i]  = summary(fit)[[1]]["stage", "Pr(>F)"]
    tukey_post[[i]] = TukeyHSD(fit)$stage
    pb$tick()
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
  }) %>% bind_rows() %>% remove_rownames()
  
  print(contrasts, digits = 2)
  message("Proportion of posterior iterations with significant ANOVA: ", mean(anova_pvals < 0.05))
  
  # Findings:
  # - Stand initiation has the highest richness (8-10 more species than closed-canopy stages, BCI 6.5-11.4)
  # - Compex has the lowest richness (1-2 fewer species on average than thinned and mature)
  # - Differences in richness are clear between thin and compex, but weaker between thin and mature
}

# TODO: Posterior estimated richness (each year)

# Principal coordinates analysis (PCoA) -----------------------------------------------------------
{
  # Average posterior draws of Bray-curtis dissimilarity matrix
  d_sum = matrix(0, n_sites, n_sites)
  pb = progress_bar$new(format = progress_bar_format, total = n_iter, clear = FALSE)
  for (i in seq_len(n_iter)) {
    d_sum = d_sum + as.matrix(vegdist(z_tmean[i, , ], method = "bray"))
    pb$tick()
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
    species_fit = envfit(pcoa_res$points, apply(z_tmean, c(2, 3), mean), permutations = 999)
    species_scores = as.data.frame(species_fit$vectors$arrows) |>
      setNames(c("PCoA1", "PCoA2")) |>
      mutate(species = species, r2 = species_fit$vectors$r, p = species_fit$vectors$pvals) |>
      arrange(desc(r2))
    top_species = species_scores %>% filter(p < 0.05) %>% filter(r2 > 0.25)
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
    top_env = env_scores %>% filter(p < 0.05) %>% filter(r2 > 0.25)
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
  # - Some "outlier" sites more similar to other stages -> are these located near stands of the other stage?
  # - Example PCoA1 species ordination correlations: brown creeper + (mature trees, closed), song sparrow - (edge, open habitat)
  # - Example PCoA2 species ordination correlations: gray jay + (montane), kingfisher/bald eagle - (valley bottom / riparian / aquatic)
}

# Community composition -----------------------------------------------------------
{
  ## Among-stage variation
  # Q: What proportion of variance is explained by stage, and are there overall differences among stages?
  R2_post = numeric(n_iter) # Proportion of variance explained by stage
  F_post  = numeric(n_iter) # Larger pseudo-F indicates greater difference between groups relative to variation within groups
  p_post  = numeric(n_iter) # Significant frequentist difference exists
  
  pb = progress_bar$new(format = progress_bar_format, total = n_iter, clear = FALSE)
  for (i in seq_len(n_iter)) {
    d          = vegdist(z_tmean[i, , ], method = "bray")
    fit        = adonis2(d ~ stages, permutations = 50) # TODO: 999 (for p value) or 0 permutations
    R2_post[i] = fit$R2[1]
    F_post[i]  = fit$F[1]
    p_post[i]  = fit$`Pr(>F)`[1]
    pb$tick()
  }
  
  res_adonis2 = data.frame(
    stat      = c("R2", "F", "p"),
    mean      = c(mean(R2_post),             mean(F_post),             mean(p_post)),
    median    = c(median(R2_post),           median(F_post),           median(p_post)),
    bci_2.5   = c(quantile(R2_post, 0.025),  quantile(F_post, 0.025),  quantile(p_post, 0.025)),
    bci_97.5  = c(quantile(R2_post, 0.975),  quantile(F_post, 0.975),  quantile(p_post, 0.975))
  )
  print(res_adonis2, digits = 1)
  # TODO: Findings:
  # - Substantial variation in taxonomic composition by stage
  # - Stage explains roughly 27% of variation in taxonomic composition
  
  message("Proportion of posterior iterations with significant PERMANOVA: ", mean(p_post < 0.05))
  
  ## Within-stage dispersion (mean distance to centroid)
  # Q: Are some stages more compositionally variable than others?
  dispersion_post = matrix(NA, nrow = n_iter, ncol = nlevels(stages))
  colnames(dispersion_post) = levels(stages)
  
  pb = progress_bar$new(format = progress_bar_format, total = n_iter, clear = FALSE)
  for (i in seq_len(n_iter)) {
    d  = vegdist(z_tmean[i, , ], method = "bray")
    bd = betadisper(d, stages)
    dispersion_post[i, ] = tapply(bd$distances, stages, mean)
    pb$tick()
  }
  
  dispersion_df = do.call(rbind, lapply(levels(stages), function(s) {
    data.frame(
      stage    = s,
      mean     = mean(dispersion_post[, s]),
      median   = median(dispersion_post[, s]),
      bci_2.5  = quantile(dispersion_post[, s], 0.025),
      bci_97.5 = quantile(dispersion_post[, s], 0.975)
    )
  })) %>% remove_rownames()
  print(dispersion_df, digits = 3)
  # TODO: Findings:
  # - Dispersion differs significantly among strata, so some of the PERMANOVA signal could reflect differences in within-stratum spread rather than centroid location alone
}

# Pairwise differences in community composition -----------------------------------------------------------
{
  ## Between-stage variation
  # Q: Which stages are different in their composition?
  stage_pairs = combn(levels(stages), 2, simplify = FALSE)
  pairwise_results = lapply(stage_pairs, function(pair) {
    message(pair[1], " vs ", pair[2])
    keep    = which(stages %in% pair)
    stg_sub = droplevels(stages[keep])
    
    R2_post = numeric(n_iter)
    F_post  = numeric(n_iter)
    p_post  = numeric(n_iter)
    
    pb = progress_bar$new(format = progress_bar_format, total = n_iter, clear = FALSE)
    for (i in seq_len(n_iter)) {
      d     = vegdist(z_tmean[i, keep, ], method = "bray")
      fit   = adonis2(d ~ stg_sub, permutations = how(nperm = 50)) # TODO: 999 (for p value) or 0
      R2_post[i] = fit$R2[1]
      F_post[i]  = fit$F[1]
      p_post[i]  = fit$`Pr(>F)`[1]
      pb$tick()
    }
    data.frame(
      pair     = paste(pair, collapse = " vs "),
      stat     = c("R2", "F", "p"),
      mean     = c(mean(R2_post),            mean(F_post),            mean(p_post)),
      median   = c(median(R2_post),          median(F_post),          median(p_post)),
      bci_2.5  = c(quantile(R2_post, 0.025), quantile(F_post, 0.025), quantile(p_post, 0.025)),
      bci_97.5 = c(quantile(R2_post, 0.975), quantile(F_post, 0.975), quantile(p_post, 0.975))
    )
  })
  
  pairwise_df = do.call(rbind, pairwise_results)
  print(pairwise_df %>% filter(stat == "R2"), digits = 2)
  print(pairwise_df %>% filter(stat == "F"),  digits = 2)
  print(pairwise_df %>% filter(stat == "p"),  digits = 2)
  # Interpretation:
  # - Small R2 -> similar assemblages; large R2 -> distinct assemblages
  # - A significant p with low R2 means detectable but minor difference
  # Findings:
  # - Standinit vs. everything else: by far the largest compositional differences, confirming standinit supports a fundamentally distinct bird community, not just a richer one.
  # - Among the remaining three strata: all significant but weak, suggesting compex, thin, and mature share broadly similar assemblages with modest compositional differences. Thin–mature is  most distinct of closed canopy, while compex–thin are most similar.
}

# TODO: Between- and within-stage mean dissimilarity
{
  ## Between- and within-stage mean dissimilarity
  # Q: How compositionally different are pairs of stages from each other?
  # Q: How compositionally homogeneous are sites within the same stage?
  stage_pairs     = combn(levels(stages), 2, simplify = FALSE)
  pair_names      = sapply(stage_pairs, paste, collapse = " vs ")
  between_post    = matrix(NA, nrow = n_iter, ncol = length(stage_pairs))
  colnames(between_post) = pair_names
  within_post     = matrix(NA, nrow = n_iter, ncol = nlevels(stages))
  colnames(within_post) = levels(stages)
  
  pb = progress_bar$new(format = progress_bar_format, total = n_iter, clear = FALSE)
  for (i in seq_len(n_iter)) {
    d_mat <- as.matrix(vegdist(z_tmean[i, , ], method = "bray"))
    for (j in seq_along(stage_pairs)) {
      pair <- stage_pairs[[j]]
      m    <- d_mat[stages == pair[1], stages == pair[2]]
      between_post[i, j] <- mean(m)
    }
    for (s in levels(stages)) {
      m <- d_mat[stages == s, stages == s]
      within_post[i, s] <- mean(m[upper.tri(m)])
    }
    pb$tick()
  }
  summarise_post <- function(mat) {
    do.call(rbind, lapply(colnames(mat), function(nm) {
      data.frame(
        group    = nm,
        mean     = mean(mat[, nm]),
        median   = median(mat[, nm]),
        bci_2.5  = quantile(mat[, nm], 0.025),
        bci_97.5 = quantile(mat[, nm], 0.975)
      )
    }))
  }
  
  between_df = summarise_post(between_post) %>% arrange(mean)
  within_df  = summarise_post(within_post)
  print(between_df, digits = 3) # Higher values -> greater compositional turnover between that pair of stages
  print(within_df,  digits = 3) # Higher values -> more compositionally heterogeneous sites within that stage
}

# TODO: Turnover vs nestedness decomposition

# TODO: Within-stage partition




