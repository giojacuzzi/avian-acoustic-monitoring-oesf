# 4_community_comp.R #############################################################################
# Community composition analyses, propagating posterior uncertainty
#
# CONFIG:
thin = 0 # Number of iterations to thin from the z posterior to speed up computation (0 or more)
nperm = 99 # TODO: 999 (for p value) or 0 permutations
binary_dissimilarity = TRUE # TODO
# INPUT:
path_msom = "data/cache/models/prefinal_msom_jags_nofp_all.rds" # "data/cache/models/msom_nofp_all_2026-02-12_19:37:00.rds"
path_trait_data = "data/cache/2_traits/1_agg_traits/trait_data.csv"
path_occurrence_predictor_plot_data = "data/cache/4_msom/1_assemble_msom_data/V3_occurrence_predictor_plot_data.rds"
path_occurrence_predictor_homerange_data = "data/cache/4_msom/1_assemble_msom_data/V3_occurrence_predictor_homerange_data.rds"
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
stages  = factor(model_data$stages$stratum_4,
                 levels = c("standinit", "compex", "thin", "mature"),
                 labels = c("Stand initiation", "Stem exclusion", "Thinning", "Mature"))

z_raw = model_data$msom$sims.list$z
str(z_raw)
n_iter    = dim(z_raw)[1]
n_sites   = dim(z_raw)[2]
n_seasons = dim(z_raw)[3]
n_species = dim(z_raw)[4]

rm(model_data)
gc()

if (thin > 0) {
  message("Thinning z_raw posterior to every ", thin, " iterations")
  thin_idx = seq(1, n_iter, by = thin)
  z_raw  = z_raw[thin_idx, , , ]
  n_iter   = length(thin_idx)
  str(z_raw)
}

message("Loading species trait data from ", path_trait_data)
species_traits = read_csv(path_trait_data, show_col_types = FALSE) %>% filter(common_name %in% species)

message("Loading occurrence predictor plot scale data from ", path_occurrence_predictor_plot_data)
occurrence_predictor_plot_data = readRDS(path_occurrence_predictor_plot_data)[[1]] %>% filter(site %in% sites)
stopifnot(all(occurrence_predictor_plot_data$site == sites))

if (binary_dissimilarity) {
  message("Summarizing binary occurrence state across years")
  z = array(NA_integer_, dim = c(n_iter, n_sites, n_species))
  for (i in seq_len(n_species)) {
    acc = matrix(0, nrow = n_iter, ncol = n_sites)
    for (t in seq_len(n_seasons)) {
      acc = pmax(acc, z_raw[, , t, i])
    }
    z[, , i] = acc
  }
  str(z)
  
} else {
  # "We retained posterior distributions of latent occurrence states averaged across years to propagate uncertainty through analyses of community composition."
  # Average over years (dim 3) while retaining iterations, sites, species, using accumulation to avoid memory limits
  message("Averaging posterior occurrence states across years as probabilities")
  z = array(0, dim = c(n_iter, n_sites, n_species))
  for (i in seq_len(n_species)) {
    acc = matrix(0, nrow = n_iter, ncol = n_sites)
    for (t in seq_len(n_seasons)) {
      acc = acc + z_raw[, , t, i]
    }
    z[, , i] = acc / n_seasons
  }
  str(z)
}

message("Using dissimilarity index '", ifelse(binary_dissimilarity, "Sorensen", "Bray-Curtis"), "'")

# Posterior estimated richness (averaged across years) ----------------------------------
{
  # Posterior estimated richness
  post_tmean_rich = matrix(0, nrow = n_iter, ncol = n_sites)
  for (i in seq_len(n_species)) {
    post_tmean_rich = post_tmean_rich + z[, , i]
  }
  
  # Posterior mean estimated richness, for summary and visualization
  post_tmean_rich_stats = data.frame(
    site      = sites,
    stage     = stages,
    mean      = apply(post_tmean_rich, 2, mean),
    median    = apply(post_tmean_rich, 2, median),
    bci_2.5   = apply(post_tmean_rich, 2, quantile, probs = 0.025),
    bci_97.5  = apply(post_tmean_rich, 2, quantile, probs = 0.975)
  )
  fig_richness = ggplot(post_tmean_rich_stats, aes(x = stage, y = median, color = stage, fill = stage)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.25, alpha = 0.35) +
    scale_fill_manual(values = alpha(colors_stats, 0.2), name = "Managment stage") +
    scale_color_manual(values = colors_stats, name = "Managment stage") +
    labs(x = "", y = "Posterior median species richness (alpha diversity)") +
    theme(
      panel.grid.major.y = element_line(color = "gray95", linewidth = 0.4)
    ); print(fig_richness)
  
  ggsave("data/cache/figs/fig_richness.pdf", fig_richness, width = 6, height = 6)
  
  post_tmean_rich_summary = post_tmean_rich_stats %>% group_by(stage) %>%
    summarise(
      n        = n(),
      mean     = mean(mean),
      median   = median(median),
      bci_2.5  = mean(bci_2.5),
      bci_97.5 = mean(bci_97.5),
      .groups  = "drop"
    )
  
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
      median    = median(diffs),
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

# Bayesian alternative:
{
  # post_tmean_rich is n_iter × n_sites
  
  # 1. Per-site posterior summaries (for visualization) ---------------------------
  #    Summarize each site's column across iterations — this is unchanged
  post_site_stats = data.frame(
    site     = sites,
    stage    = stages,
    mean     = apply(post_tmean_rich, 2, mean),
    median   = apply(post_tmean_rich, 2, median),
    bci_2.5  = apply(post_tmean_rich, 2, quantile, probs = 0.025),
    bci_97.5 = apply(post_tmean_rich, 2, quantile, probs = 0.975)
  )
  
  fig_richness = ggplot(post_site_stats, aes(x = stage, y = median, color = stage, fill = stage)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.25, alpha = 0.35) +
    scale_fill_manual(values = alpha(colors_stats, 0.2), name = "Managment stage") +
    scale_color_manual(values = colors_stats, name = "Managment stage") +
    labs(x = "", y = "Posterior median species richness (alpha diversity)") +
    theme(
      panel.grid.major.y = element_line(color = "gray95", linewidth = 0.4)
    ); print(fig_richness)
  
  # 2. Per-stage posterior (for reporting) ----------------------------------------
  #    For each iteration, average richness across sites within each stage.
  #    This gives a full n_iter posterior for each stage's mean richness.
  stage_levels = levels(stages)
  
  stage_post = sapply(stage_levels, function(s) {
    rowMeans(post_tmean_rich[, stages == s, drop = FALSE])
  })
  # stage_post: n_iter × n_stages — each column is a posterior distribution
  
  post_stage_stats = data.frame(
    stage    = stage_levels,
    mean     = apply(stage_post, 2, mean),
    median   = apply(stage_post, 2, median),
    bci_2.5  = apply(stage_post, 2, quantile, probs = 0.025),
    bci_97.5 = apply(stage_post, 2, quantile, probs = 0.975)
  )
  
  # 3. Pairwise contrasts ---------------------------------------------------------
  #    Subtract stage posteriors directly — the difference distribution IS the
  #    inference. No test statistic needed.
  stage_pairs = combn(stage_levels, 2, simplify = FALSE)
  
  contrasts = lapply(stage_pairs, function(pair) {
    diffs = stage_post[, pair[1]] - stage_post[, pair[2]]
    data.frame(
      pair     = paste0(pair[1], " - ", pair[2]),
      mean     = mean(diffs),
      median   = median(diffs),
      bci_2.5  = quantile(diffs, 0.025),
      bci_97.5 = quantile(diffs, 0.975),
      prob_gt0 = mean(diffs > 0)   # P(stage B > stage A)
    )
  }) %>% bind_rows() %>% remove_rownames()
  
  print(post_stage_stats, digits = 4)
  print(contrasts, digits = 2)
  
  ggplot(post_stage_stats, aes(x = median, y = stage)) +
    geom_point() + geom_errorbarh(aes(xmin = bci_2.5, xmax = bci_97.5), width = 0.1, position = position_dodge(width = 0.5))
}

## DEBUG: "Calculate the proportion of draws in which the dissimilarities between sites from different stages are larger than dissimilarities between sites from the same stage"
# This is analogous to a PERMANOVA and pairwise PERMANOVAs, but fully Bayesian
{
  n_iter  <- dim(z)[1]
  iter_sample <-  1:n_iter # sample(1:n_iter, 1000)
  
  # ── Index pairs (built once) ──────────────────────────────────────────────────
  
  stage_mat   <- outer(stages, stages, FUN = "==")
  within_idx  <- which(stage_mat  & lower.tri(stage_mat))
  between_idx <- which(!stage_mat & lower.tri(stage_mat))
  
  # All unique stage pairs
  stage_levels <- levels(stages)
  stage_pairs  <- combn(stage_levels, 2, simplify = FALSE)  # 6 pairs
  
  # For each pair, index of between-pair and within-pair (both stages) sites
  pair_idx <- lapply(stage_pairs, function(p) {
    s1 <- stages == p[1]
    s2 <- stages == p[2]
    list(
      between = which(outer(s1, s2) & lower.tri(stage_mat) |
                        outer(s2, s1) & lower.tri(stage_mat)),
      within  = which(stage_mat & lower.tri(stage_mat) &
                        outer(s1 | s2, s1 | s2))
    )
  })
  names(pair_idx) <- sapply(stage_pairs, paste, collapse = " vs ")
  
  # ── Storage ───────────────────────────────────────────────────────────────────
  
  global_diff  <- numeric(length(iter_sample))
  pair_diffs   <- matrix(NA, nrow = length(iter_sample), ncol = length(stage_pairs),
                         dimnames = list(NULL, names(pair_idx)))
  
  within_raw  <- matrix(NA, nrow = length(iter_sample), ncol = length(stage_levels),
                        dimnames = list(NULL, stage_levels))
  between_raw <- matrix(NA, nrow = length(iter_sample), ncol = length(stage_pairs),
                        dimnames = list(NULL, names(pair_idx)))
  
  # ── Main loop ─────────────────────────────────────────────────────────────────
  pb = progress_bar$new(format = progress_bar_format, total = length(iter_sample), clear = FALSE)
  for (i in seq_along(iter_sample)) {
    k    <- iter_sample[i]
    diss <- as.matrix(vegdist(z[k, , ], method = "bray", binary = binary_dissimilarity))
    
    # Global
    global_diff[i] <- mean(diss[between_idx]) - mean(diss[within_idx])
    
    # Pairwise
    for (j in seq_along(pair_idx)) {
      idx <- pair_idx[[j]]
      pair_diffs[i, j]  <- mean(diss[idx$between]) - mean(diss[idx$within])
      between_raw[i, j] <- mean(diss[idx$between])
    }
    
    # Within-group raw dissimilarities
    for (s in stage_levels) {
      idx_s <- which(stages == s)
      d_s   <- diss[idx_s, idx_s]
      within_raw[i, s] <- mean(d_s[lower.tri(d_s)])
    }
    pb$tick()
  }
  
  # ── Summarise ─────────────────────────────────────────────────────────────────
  
  summarise_diff <- function(x, name) {
    data.frame(
      comparison      = name,
      p_greater       = mean(x > 0),
      median_diff     = median(x),
      lower_95        = quantile(x, 0.025),
      upper_95        = quantile(x, 0.975)
    )
  }
  
  summarise_raw <- function(mat) {
    as.data.frame(t(apply(mat, 2, function(x) c(
      median   = median(x),
      lower_95 = quantile(x, 0.025),
      upper_95 = quantile(x, 0.975)
    ))))
  }
  
  results <- do.call(rbind, c(
    list(summarise_diff(global_diff, "Global (all between vs all within)")),
    lapply(names(pair_idx), function(nm) summarise_diff(pair_diffs[, nm], nm))
  ))
  
  rownames(results) <- NULL
  # p_greater - posterior probability that between-group dissimilarity exceeds within-group dissimilarity: "in X% of draws consistent with my data, assemblages from different stages are more dissimilar than assemblages from the same stage."
  # median_diff - median estimate of how much larger between-group dissimilarity is than within-group dissimilarity
  # Answer:
  # 1. Global --> assemblage composition differs between stages, 100% probability
  # 2. Stand initiation is compositionally distinct --> 100% probability, and largest diff values
  # 3. Closed canopy stages --> all differ from each other with high posterior probability, but the effect sizes are much smaller and the credible intervals approach zero in two cases:
  # 3A. compex vs thin - weakest contrast, not distinguishable
  # 3B. compex vs mature - similarly weak, not distinguishable
  # 3C. thin vs mature - strongest contrast, 100% probability; statistically distinguishable but modest differences
  print(results, digits = 3)
  
  cat("\n=== Raw within-group dissimilarities ===\n")
  print(summarise_raw(within_raw), digits = 3)
  
  cat("\n=== Raw between-group dissimilarities ===\n")
  print(summarise_raw(between_raw), digits = 3)
  
  # Answer:
  # Closed-canopy stages are only marginally more different from each other than sites within the same stage are from each other. The assemblages broadly overlap, with stand init as the clear outlier.
}

{ # Q: Are underdev and old sites compositionally dissimlar, or statistically indistinguishable?
  stratum <- occurrence_predictor_plot_data$stratum_5
  
  idx_underdev <- which(stratum == "underdev")
  idx_old      <- which(stratum == "old")
  
  # Within-group indices
  underdev_pairs <- which(outer(stratum == "underdev", stratum == "underdev") & 
                            lower.tri(matrix(0, 224, 224)))
  old_pairs      <- which(outer(stratum == "old", stratum == "old") & 
                            lower.tri(matrix(0, 224, 224)))
  
  # Storage
  diff_vs_within  <- numeric(length(iter_sample))
  between_raw_sub <- numeric(length(iter_sample))
  within_underdev <- numeric(length(iter_sample))
  within_old      <- numeric(length(iter_sample))
  
  for (i in seq_along(iter_sample)) {
    k    <- iter_sample[i]
    diss <- as.matrix(vegdist(z[k, , ], method = "bray", binary = binary_dissimilarity))
    
    between         <- mean(diss[idx_underdev, idx_old])
    w_underdev      <- mean(diss[idx_underdev, idx_underdev][lower.tri(diss[idx_underdev, idx_underdev])])
    w_old           <- mean(diss[idx_old, idx_old][lower.tri(diss[idx_old, idx_old])])
    
    between_raw_sub[i] <- between
    within_underdev[i] <- w_underdev
    within_old[i]      <- w_old
    diff_vs_within[i]  <- between - mean(c(w_underdev, w_old))
  }
  
  cat("P(between > within):", mean(diff_vs_within > 0), "\n")
  
  cat("\nBetween-group (underdev vs old):\n")
  print(quantile(between_raw_sub, c(0.025, 0.5, 0.975)), digits = 3)
  
  cat("\nWithin underdev:\n")
  print(quantile(within_underdev, c(0.025, 0.5, 0.975)), digits = 3)
  
  cat("\nWithin old:\n")
  print(quantile(within_old, c(0.025, 0.5, 0.975)), digits = 3)
  
  cat("\nDifference (between minus mean within) 95% CI:\n")
  print(quantile(diff_vs_within, c(0.025, 0.5, 0.975)), digits = 3)
  
  # Answer:
  # Underdev and old are statistically distinguishable, but their assemblages broadly overlap and the compositional difference is ecologically negligible
  
}
 
{   # Q: Are thinned assemblages more similar to mature than to stem exclusion?
  # Now calculate the posterior probability that Thinning is compositionally closer to Stem exclusion than to Mature, with no baseline correction needed since both comparisons are raw between-group dissimilarities on the same scale.
  # Indices for each stage
  idx_thin <- which(stages == "Thinning")
  idx_stem <- which(stages == "Stem exclusion")
  idx_mat  <- which(stages == "Mature")
  
  # Storage
  thin_vs_stem <- numeric(length(iter_sample))
  thin_vs_mat  <- numeric(length(iter_sample))
  
  for (i in seq_along(iter_sample)) {
    k    <- iter_sample[i]
    diss <- as.matrix(vegdist(z[k, , ], method = "bray", binary = binary_dissimilarity))
    
    thin_vs_stem[i] <- mean(diss[idx_thin, idx_stem])
    thin_vs_mat[i]  <- mean(diss[idx_thin, idx_mat])
  }
  
  # Posterior probability that Thinning is closer to Mature than to Stem exclusion
  mean(thin_vs_mat < thin_vs_stem)
  
  # And the difference (positive = Thinning closer to Mature)
  diff <- thin_vs_stem - thin_vs_mat
  quantile(diff, c(0.025, 0.5, 0.975))
  mean(diff > 0)
  
  # A: We cannot confidently conclude that thinned is more similar to either stem exclusion or mature — they are compositionally roughly equidistant from thinned
  
}

stopifnot(all(occurrence_predictor_plot_data$site == sites))

summarise_cor <- function(x, name) {
  data.frame(
    stage    = name,
    p_pos    = mean(x > 0),
    median_r = median(x),
    lower_95 = quantile(x, 0.025),
    upper_95 = quantile(x, 0.975)
  )
}

# Q: Within each stage, are stands closer in age also more similar in composition?
{
  age <- occurrence_predictor_plot_data$age_mean
  
  age_diss_cor <- matrix(NA, nrow = length(iter_sample), ncol = length(stage_levels),
                         dimnames = list(NULL, stage_levels))
  
  for (i in seq_along(iter_sample)) {
    k    <- iter_sample[i]
    diss <- as.matrix(vegdist(z[k, , ], method = "bray", binary = binary_dissimilarity))
    
    for (s in stage_levels) {
      idx_s    <- which(stages == s & !is.na(age))
      age_s    <- age[idx_s]
      age_diff <- as.vector(dist(age_s))
      diss_s   <- diss[idx_s, idx_s]
      diss_vec <- diss_s[lower.tri(diss_s)]
      
      age_diss_cor[i, s] <- cor(age_diff, diss_vec, method = "spearman")
    }
  }
  
  cor_results <- do.call(rbind, lapply(stage_levels, function(s) {
    summarise_cor(age_diss_cor[, s], s)
  }))
  
  rownames(cor_results) <- NULL
  print(cor_results, digits = 3)
  
  # Answer:
  # Stand age similarity predicted compositional similarity most strongly in stand initiation (median Spearman rs = 0.45, P = 1), and weakly in stem exclusion and mature stages (rs = 0.10–0.11, P > 0.97). However, thinned  and not at all in thinned stands (median ρ = 0.01, P = 0.61), suggesting that management history rather than stand age drives compositional variation in thinned forests.
}

# Q: Does dissimilarity between thinned and mature decrease with thinned stand age?
{
  age <- occurrence_predictor_plot_data$age_mean
  
  idx_thin <- which(stages == "Thinning")
  idx_mat  <- which(stages == "Mature")
  
  age_thin <- age[idx_thin]  # one age per thinned site
  
  # Storage
  slope     <- numeric(length(iter_sample))
  intercept <- numeric(length(iter_sample))
  
  for (i in seq_along(iter_sample)) {
    k    <- iter_sample[i]
    diss <- as.matrix(vegdist(z[k, , ], method = "bray", binary = binary_dissimilarity))
    
    # Mean dissimilarity from each thinned site to all mature sites
    mean_diss_to_mature <- rowMeans(diss[idx_thin, idx_mat])
    
    fit          <- lm(mean_diss_to_mature ~ age_thin)
    intercept[i] <- coef(fit)[1]
    slope[i]     <- coef(fit)[2]
  }
  
  cat("P(slope < 0):", mean(slope < 0), "\n")
  cat("\nSlope 95% CI:\n")
  print(quantile(slope, c(0.025, 0.5, 0.975)), digits = 3)
  
  # Q:
  # - We found no evidence that dissimilarity between thinned and mature assemblages decreased with time since thinning (up to 25 years).
}

# Principal coordinates analysis (PCoA) -----------------------------------------------------------
{
  # Average posterior draws of Bray-curtis dissimilarity matrix
  d_sum = matrix(0, n_sites, n_sites)
  pb = progress_bar$new(format = progress_bar_format, total = n_iter, clear = FALSE)
  for (i in seq_len(n_iter)) {
    d_sum = d_sum + as.matrix(vegdist(z[i, , ], method = "bray", binary = binary_dissimilarity))
    pb$tick()
  }
  d_mean_post = as.dist(d_sum / n_iter)
  
  # PCoA
  pcoa_res = cmdscale(d_mean_post, eig = TRUE, k = 2)
  pct  = round(pcoa_res$eig / sum(pcoa_res$eig[pcoa_res$eig > 0]) * 100, 1)
  pcoa_df = as.data.frame(pcoa_res$points) %>% setNames(c("PCoA1", "PCoA2")) %>% mutate(site = sites, stage = stages)
  
  fig_PCoA_sites_blank = ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2)) +
    geom_vline(xintercept = 0, color = "gray90") + geom_hline(yintercept = 0, color = "gray90") +
    geom_point(size = 2, alpha = 0.8, color = "gray20") +
    labs(x = paste0("PCoA 1 (", pct[1], "%)"), y = paste0("PCoA 2 (", pct[2], "%)")) +
    scale_x_continuous(limits = c(-0.095, 0.155), expand = c(0,0)) +
    scale_y_continuous(limits = c(-0.06, 0.072), expand = c(0,0)) +
    theme(); print(fig_PCoA_sites_blank)
  
  fig_PCoA_sites = ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = stage)) +
    geom_vline(xintercept = 0, color = "gray90") + geom_hline(yintercept = 0, color = "gray90") +
    geom_point(size = 2, alpha = 0.8) +
    stat_ellipse(aes(fill = stage), geom = "polygon", alpha = 0.1, level = 0.95) +
    scale_color_manual(values = colors_stats) +
    scale_fill_manual(values = colors_stats) +
    labs(x = paste0("PCoA 1 (", pct[1], "%)"), y = paste0("PCoA 2 (", pct[2], "%)")) +
    scale_x_continuous(limits = c(-0.095, 0.155), expand = c(0,0)) +
    scale_y_continuous(limits = c(-0.06, 0.072), expand = c(0,0)) +
    theme(legend.position = "none"); print(fig_PCoA_sites)
  
  ggsave("data/cache/figs/fig_PCoA_sites_blank.pdf", fig_PCoA_sites_blank)
  ggsave("data/cache/figs/fig_PCoA_sites.pdf", fig_PCoA_sites)
  
  { # Species ordination correlations
    species_fit = envfit(pcoa_res$points, apply(z, c(2, 3), mean), permutations = nperm)
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
    env_fit = envfit(pcoa_res$points, env_data, permutations = nperm, na.rm = TRUE)
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

# NMDS
{
  # NMDS
  nmds_res = metaMDS(d_mean_post, distance = "bray", k = 2, trymax = 100, trace = FALSE)
  cat("Stress:", nmds_res$stress, "\n")
  
  nmds_k2 = metaMDS(d_mean_post, distance = "bray", trymax = 200, k=2)
  nmds_k2
  nmds_k2$stress
  
  nmds_k3 = metaMDS(d_mean_post, distance = "bray", trymax = 200, k=3)
  nmds_k3
  nmds_k3$stress
  
  nmds_df = as_tibble(scores(nmds_k3, display = "sites")) %>%
    mutate(stage = stages)
  
  fig_nmds_blank = ggplot(nmds_df, aes(x = NMDS1, y = NMDS2)) +
    geom_point(size = 2, color = "gray20") +
    annotate("text", x = Inf, y = -Inf,
             label = paste("Stress =", round(nmds_k3$stress, 2)),
             hjust = 1.1, vjust = -1, size = 3, fontface = "italic") +
    coord_cartesian(xlim = c(-0.095, 0.14), ylim = c(-0.065, 0.075)) +
    theme(legend.position = "none"); print(fig_nmds_blank)
  
  fig_nmds_color = ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = stage)) +
    geom_point(size = 2, alpha = 0.8) +
    stat_ellipse(level = 0.95, aes(fill = stage), geom = "polygon") +
    scale_color_manual(values = colors_stats) +
    scale_fill_manual(values = alpha(colors_stats, 0.1)) +
    labs(color = "Management stage") +
    coord_cartesian(xlim = c(-0.095, 0.14), ylim = c(-0.065, 0.075)) +
    annotate("text", x = Inf, y = -Inf,
             label = paste("Stress =", round(nmds_k3$stress, 3)),
             hjust = 1.1, vjust = -1, size = 3, fontface = "italic") +
    theme(legend.position = "none"); print(fig_nmds_color)
  
  ggsave("data/cache/figs/fig_nmds_blank.pdf", fig_nmds_blank, width = 6, height = 6)
  ggsave("data/cache/figs/fig_nmds_color.pdf", fig_nmds_color, width = 6, height = 6)
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
    d          = vegdist(z[i, , ], method = "bray")
    fit        = adonis2(d ~ stages, permutations = nperm)
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
  # Q: Is the spread of sites around their stage centroid homogeneous across stages?
  dispersion_post = matrix(NA, nrow = n_iter, ncol = nlevels(stages))
  colnames(dispersion_post) = levels(stages)
  dispersion_p <- numeric(n_iter)
  
  pb = progress_bar$new(format = progress_bar_format, total = n_iter, clear = FALSE)
  for (i in seq_len(n_iter)) {
    d  = vegdist(z[i, , ], method = "bray")
    bd = betadisper(d, stages)
    dispersion_post[i, ] = tapply(bd$distances, stages, mean)
    dispersion_p[i]      = permutest(bd, permutations = nperm)$tab["Groups", "Pr(>F)"]
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
  message("Proportion of posterior iterations with significant dispersion: ", round(mean(dispersion_p < 0.05), 2))
  
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
      d     = vegdist(z[i, keep, ], method = "bray")
      fit   = adonis2(d ~ stg_sub, permutations = how(nperm = nperm))
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

# Landscape-level turnover vs nestedness decomposition
{
  beta_multi_post = matrix(NA, nrow = n_iter, ncol = 3)
  colnames(beta_multi_post) = c("beta.SOR", "beta.SIM", "beta.SNE")
  
  pb = progress_bar$new(format = progress_bar_format, total = n_iter, clear = FALSE)
  for (i in seq_len(n_iter)) {
    bp_core            = betapart.core(z[i, , ])
    bm                 = beta.multi(bp_core, index.family = "sorensen")
    beta_multi_post[i, ] = c(bm$beta.SOR, bm$beta.SIM, bm$beta.SNE)
    pb$tick()
  }
  
  beta_multi_df = do.call(rbind, lapply(colnames(beta_multi_post), function(nm) {
    data.frame(
      component = nm,
      mean      = mean(beta_multi_post[, nm]),
      median    = median(beta_multi_post[, nm]),
      bci_2.5   = unname(quantile(beta_multi_post[, nm], 0.025)),
      bci_97.5  = unname(quantile(beta_multi_post[, nm], 0.975))
    )
  }))
  print(beta_multi_df, digits = 3)
  # Interpretation:
  # - beta.SOR = total Sorensen dissimilarity
  # - beta.SIM = turnover component
  # - beta.SNE = nestedness component
  # Findings:
  # - Overall beta diversity (SOR 0.95) is high, so communities are dissimilar across all sites
  # - Nearly all dissimilarity driven by species turnover (SIM 0.9), with nestedness contributing almost nothing (SNE 0.04).
  # - Sites are not simply species-poor subsets of richer sites, but instead genuinely differ in species.
}

# Pairwise turnover vs nestedness decomposition
{
  stage_pairs = combn(levels(stages), 2, simplify = FALSE)
  pair_names  = sapply(stage_pairs, paste, collapse = " vs ")
  components  = c("beta.sor", "beta.sim", "beta.sne")
  
  # Between-stage
  # Q: How compositionally different are pairs of stages from each other?
  between_bp = array(NA, dim = c(n_iter, length(stage_pairs), 3),
                      dimnames = list(NULL, pair_names, components))
  # Within-stage
  # Q: How compositionally homogeneous are sites within the same stage?
  within_bp  = array(NA, dim = c(n_iter, nlevels(stages), 3),
                      dimnames = list(NULL, levels(stages), components))
  
  pb = progress_bar$new(format = progress_bar_format, total = n_iter, clear = FALSE)
  for (i in seq_len(n_iter)) {
    bp_core = betapart.core(z[i, , ])
    bp_pair = beta.pair(bp_core, index.family = "sorensen")
    
    for (comp in components) {
      d_mat = as.matrix(bp_pair[[comp]])
      
      for (j in seq_along(stage_pairs)) {
        pair = stage_pairs[[j]]
        m    = d_mat[stages == pair[1], stages == pair[2]]
        between_bp[i, pair_names[j], comp] = mean(m)
      }
      
      for (s in levels(stages)) {
        m = d_mat[stages == s, stages == s]
        within_bp[i, s, comp] = mean(m[upper.tri(m)])
      }
    }
    pb$tick()
  }
  
  # Summarise
  summarise_array = function(arr) {
    do.call(rbind, lapply(dimnames(arr)[[2]], function(grp) {
      do.call(rbind, lapply(components, function(comp) {
        data.frame(
          group     = grp,
          component = comp,
          mean      = mean(arr[, grp, comp]),
          median    = median(arr[, grp, comp]),
          bci_2.5   = unname(quantile(arr[, grp, comp], 0.025)),
          bci_97.5  = unname(quantile(arr[, grp, comp], 0.975))
        )
      }))
    }))
  }
  
  between_bp_df = summarise_array(between_bp) %>% arrange(component, group)
  within_bp_df  = summarise_array(within_bp)  %>% arrange(component, group)
  
  print(between_bp_df, digits = 2)
  print(within_bp_df,  digits = 2)
  
  format_bp = function(df, group_col) {
    df %>%
      select({{ group_col }} := group, component, median) %>%
      pivot_wider(names_from = component, values_from = median) %>%
      rename(sor = beta.sor, sim = beta.sim, sne = beta.sne) %>%
      mutate(
        pct_turnover  = round(100 * sim / sor, 1),
        pct_nestedness = round(100 * sne / sor, 1)
      )
  }
  
  between_summary = format_bp(between_bp_df, group_col = "pair")
  within_summary  = format_bp(within_bp_df,  group_col = "stage")
  between_summary
  # Interpretation:
  # Higher beta.sor values -> greater compositional turnover between that pair of stages
  within_summary
  # Interpretation:
  # Higher beta.sor values -> more compositionally heterogeneous sites within that stage
  
  # Findings:
  # - Differences among closed-canopy stages are primarily driven by species turnover. Although compex is species-poor relative to thin and mature, those species are not simply a subset -- compex isn't only filtering the mature/thin species, but also selecting for species affiliated with compex specifically.
  # - Increasing nestedness driving differences between standinit and closed canopy indicate that canopy closure filters species out of the standinit pool. Nestedness is strongest in compex stands where stem exclusion environment is most intense.
  # - Canopy closure reduces richness (nestedness), but which species persist under different management regimes still varies substantially (turnover). The practical implication is that thin and mature are not interchangeable from a conservation standpoint despite similar richness — they are likely supporting different species within their shared reduced pool.
  # - Standinit supports a superset community containing many of the species found elsewhere plus a suite of early-successional specialists. The other three strata share a broadly similar forest-interior community that differs among itself mainly through species replacement rather than richness difference.
  # Standinit is compositionally distinct from everything else, and nestedness plays a large role in why. The three standinit contrasts have the highest total dissimilarity (SOR 0.229–0.262), which you already knew from the PERMANOVA. What's new here is that nestedness accounts for 39–49% of that dissimilarity depending on the comparison. This means standinit's community is not simply a reshuffled version of the closed-canopy stages — it retains some species found in other stages while also harboring a distinct set of open-habitat species. The high nestedness component suggests that closed-canopy stages are compositionally nested within standinit to a meaningful degree, which makes ecological sense: early successional sites often support a superset of species that includes both generalists found throughout and specialists tied to open conditions.
  # Among the three closed-canopy stages, dissimilarity is lower and driven almost entirely by turnover. Compex vs mature, compex vs thin, and thin vs mature all have SOR around 0.18 — substantially lower than the standinit contrasts — and 68–75% of that dissimilarity is pure species replacement rather than nestedness. This suggests these stages support broadly similar species pools but with different dominant assemblages, consistent with gradual compositional shifts along a canopy closure or structural complexity gradient rather than wholesale community reorganization.
  # The standinit vs compex contrast stands out with the highest total dissimilarity (SOR 0.262) and the most even turnover/nestedness split (51/49). Compex is the most structurally distinct of the closed-canopy stages — recently transitioned from standinit, with complex vertical structure — which may explain why it shows the strongest nestedness signal relative to standinit compared to thin and mature.
  # nestedness signal would reflect landscape context rather than habitat filtering: The nestedness gradient across standinit contrasts could simply reflect patch permeability — small openings in a closed-canopy matrix are easily traversed by forest interior species, inflating apparent co-occurrence with early successional specialists.
  # alternatively, the nestedness could indicate that several species found in closed-canopy habitat also use open-canopy habitat; The nestedness signal reflects the habitat breadth of species (e.g. generalists, mature species that use early successional for foraging)
  # it may be a mixture of both
}

## RESULTS ----------------------------------------------------------------------------

# Richness
# Q: Do stages differ in alpha diversity (species richness)?
post_tmean_rich_summary
print(contrasts, digits = 2) # Tukey
message("Proportion of posterior iterations with significant ANOVA: ", mean(anova_pvals < 0.05))

# Global PERMANOVA
# Q: What proportion of variance is explained by stage, and are there overall differences among stages (as measured by differences in centroid location)?
print(res_adonis2, digits = 2)
message("Proportion of posterior iterations with significant PERMANOVA: ", mean(p_post < 0.05))
# Within-stage dispersion (mean distance to centroid)
# Q: Are some stages more compositionally variable than others, as measured by the spread of sites around their stage centroid?
print(dispersion_df, digits = 2)
message("Proportion of posterior iterations with significant dispersion: ", mean(dispersion_p < 0.05))
# Pairwise between-stage
# Q: Which stages are different in their composition, and by how much?
print(pairwise_df %>% filter(stat == "R2"), digits = 2)
print(pairwise_df %>% filter(stat == "p"), digits = 2)

# Global beta diversity decomposition
# Q: How compositionally different are sites across the landscape, and are differences driven by turnover or nestedness?
print(beta_multi_df, digits = 2)
# Pairwise beta diversity decomposition
# Q: How compositionally different are pairs of stages from each other, and are differences driven by turnover or nestedness?
print(between_bp_df, digits = 2)
between_summary
# Within-stage beta diversity decomposition
# Q: How compositionally homogeneous are sites within the same stage, and are differences driven by turnover or nestedness?
print(within_bp_df,  digits = 2)
within_summary

stop("DEBUG")

## Combine into tables

# Stage-level metrics (median values)
list(
  post_tmean_rich_summary %>% select(stage, n, median) %>% rename(alpha = median),
  dispersion_df %>% select(stage, median) %>% rename(dispersion = median),
  within_summary %>% select(stage, sor, sim, sne)
) %>% reduce(left_join, by = "stage") %>% mutate(across(where(is.numeric), ~ round(.x, 2)))

# Pairwise metrics
list(
  pairwise_df %>% filter(stat == "R2") %>% select(pair, median) %>% rename(R2 = median),
  pairwise_df %>% filter(stat == "F") %>% select(pair, median) %>% rename(F = median),
  between_summary %>% select(pair, sor, sim, sne)
) %>% reduce(left_join, by = "pair") %>% mutate(across(where(is.numeric), ~ round(.x, 2)))


# TODO: 
# Q: Is proportion of another stage within the median species home range correlated with similarity to that stage?
# canada jay ~ mean home range size
# barred owl ~ 1km
scale = "canada jay"
homerange_data = readRDS(path_occurrence_predictor_homerange_data)[[1]][[scale]]
unique(homerange_data$buffer_radius_m)
{
  # Stage to homerange column mapping
  stage_to_pcnt <- c(
    "Stand initiation" = "pcnt_standinit",
    "Stem exclusion"   = "pcnt_compex",
    "Thinning"         = "pcnt_thin",
    "Mature"           = "pcnt_mature"
  )
  
  # Confirm site order matches
  stopifnot(all(homerange_data$site == sites))
  
  # Storage: one correlation per stage per iteration
  # For each stage s: correlate pcnt_s at non-s sites vs mean dissimilarity to s
  hr_cor <- matrix(NA, nrow = length(iter_sample), ncol = length(stage_levels),
                   dimnames = list(NULL, stage_levels))
  
  for (i in seq_along(iter_sample)) {
    k    <- iter_sample[i]
    diss <- as.matrix(vegdist(z[k, , ], method = "bray", binary = binary_dissimilarity))
    
    for (s in stage_levels) {
      idx_s     <- which(stages == s)   # sites IN stage s (used as reference)
      idx_not_s <- which(stages != s)   # focal sites NOT in stage s
      
      # Mean dissimilarity from each non-s site to all sites in stage s
      mean_diss_to_s <- rowMeans(diss[idx_not_s, idx_s])
      
      # Proportion of stage s in home range for focal sites
      pcnt_s <- homerange_data[[stage_to_pcnt[s]]][idx_not_s]
      
      keep <- !is.na(pcnt_s)
      hr_cor[i, s] <- cor(pcnt_s[keep], mean_diss_to_s[keep], method = "spearman")
    }
  }
  
  # Summarise
  # p_neg = posterior probability that more of stage s nearby -> more similar to stage s
  hr_cor_results <- do.call(rbind, lapply(stage_levels, function(s) {
    x <- hr_cor[, s]
    data.frame(
      stage    = s,
      p_neg    = mean(x < 0),
      median_r = median(x),
      lower_95 = quantile(x, 0.025),
      upper_95 = quantile(x, 0.975)
    )
  }))
  
  rownames(hr_cor_results) <- NULL
  print(hr_cor_results, digits = 3)
  # Answer:
  # 
}


#### ------------------------------------------------------------------------------
# Q: Do smaller patches have higher nestendess?

# Filter to standinit sites
standinit_idx    <- which(stages == "standinit")
closedcanopy_idx <- which(stages != "standinit")

# For each posterior draw, compute mean SNE for each standinit site
# against all closed-canopy sites
site_sne_post <- matrix(NA, nrow = n_iter, ncol = length(standinit_idx))
colnames(site_sne_post) <- sites[standinit_idx]

pb <- progress_bar$new(format = progress_bar_format, total = n_iter, clear = FALSE)
for (i in seq_len(n_iter)) {
  bp_core <- betapart.core(z)
  bp_pair <- beta.pair(bp_core, index.family = "sorensen")
  sne_mat <- as.matrix(bp_pair$beta.sne)
  
  for (k in seq_along(standinit_idx)) {
    site_sne_post[i, k] <- mean(sne_mat[standinit_idx[k], closedcanopy_idx])
  }
  pb$tick()
}

# Posterior mean SNE per standinit site
site_sne_mean <- colMeans(site_sne_post)

# Join with landscape data
spillover_df <- homerange_data %>%
  filter(site %in% sites[standinit_idx]) %>%
  mutate(
    mean_sne          = site_sne_mean,
    pcnt_closedcanopy = pcnt_compex + pcnt_thin + pcnt_mature
  )

# Visualise
ggplot(spillover_df, aes(x = focalpatch_area_homeange_pcnt, y = mean_sne)) +
  geom_point() + geom_smooth(method = "lm") +
  labs(x = "Focal patch area (% of buffer)", y = "Mean nestedness (SNE) vs closed-canopy sites")

ggplot(spillover_df, aes(x = pcnt_closedcanopy, y = mean_sne)) +
  geom_point() + geom_smooth(method = "lm") +
  labs(x = "Surrounding closed-canopy cover (%)", y = "Mean nestedness (SNE) vs closed-canopy sites")

# Regression
lm_spillover <- lm(mean_sne ~ focalpatch_area_homeange_pcnt + pcnt_closedcanopy, data = spillover_df)
summary(lm_spillover)
# A plausible ecological interpretation of the negative relationship is that standinit sites surrounded by more closed-canopy forest are more deeply embedded in a forest-interior landscape context, which may actually exclude early-successional species that would otherwise colonize from the broader landscape. In other words, landscape context filters the open-habitat specialists out rather than letting closed-canopy species in — the opposite of spillover.





# Are thinned sites more similar to mature than to stem exclusion?
{
  ggplot(between_bp_df, aes(x = mean, y = group, color = group)) +
    geom_pointrange(aes(xmin = bci_2.5, xmax = bci_97.5)) +
    facet_wrap(~component) +
    theme_minimal() +
    theme(legend.position = "none")

  diff_samples <- between_bp[, "Thinning vs Mature", "beta.sor"] -
    between_bp[, "Stem exclusion vs Thinning", "beta.sor"]
  
  mean(diff_samples)                    # posterior mean difference
  quantile(diff_samples, c(0.025, 0.975))  # 95% BCI of the difference
  mean(diff_samples < 0)               # P(Thinning more similar to Mature)
  mean(diff_samples > 0)               # P(Thinning more similar to Stem exclusion)
}


{
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  
  # ── Wrangle (leave between_bp_df untouched) ───────────────────────────────────
  stage_order <- c("Stand initiation", "Stem exclusion", "Thinning", "Mature")
  x_stages    <- stage_order[1:3]
  y_stages    <- stage_order[2:4]
  
  plot_df <- between_bp_df %>%
    mutate(
      stage_a = factor(sub(".* vs ", "", group), levels = stage_order),
      stage_b = factor(sub(" vs .*", "", group), levels = stage_order)
    )
  
  sor_vals <- plot_df %>%
    filter(component == "beta.sor") %>%
    select(stage_a, stage_b, sor_mean = mean)
  
  plot_df <- plot_df %>%
    left_join(sor_vals, by = c("stage_a", "stage_b")) %>%
    mutate(
      fill_val = case_when(
        component == "beta.sor" ~ mean,
        TRUE                    ~ mean / sor_mean
      )
    )
  
  # ── Panel labels and colour palettes ─────────────────────────────────────────
  panel_meta <- list(
    beta.sor = list(title = "Total diversity (β sor)", diverging = FALSE, show_legend = TRUE,
                    show_y = TRUE,
                    low = "white",
                    high = "#e66101" #"mediumpurple"
                    ),
    beta.sim = list(title = "Turnover (β sim)",        diverging = TRUE,  show_legend = TRUE,
                    show_y = FALSE),
    beta.sne = list(title = "Nestedness (β sne)",      diverging = TRUE,  show_legend = TRUE,
                    show_y = FALSE)
  )
  
  # ── Helper: build one panel ───────────────────────────────────────────────────
  make_panel <- function(component_name, df, meta) {
    
    d <- filter(df, component == component_name) %>%
      select(stage_a, stage_b, mean, fill_val)
    
    full_grid <- expand.grid(
      stage_a = factor(y_stages, levels = stage_order),
      stage_b = factor(x_stages, levels = stage_order),
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        stage_a = factor(stage_a, levels = stage_order),
        stage_b = factor(stage_b, levels = stage_order)
      ) %>%
      left_join(d, by = c("stage_a", "stage_b")) %>%
      mutate(is_empty = is.na(mean))
    
    p <- ggplot(full_grid, aes(x = stage_b, y = stage_a, fill = fill_val)) +
      geom_tile(
        data  = filter(full_grid, !is_empty),
        color = "white", linewidth = 0.75
      ) +
      geom_tile(
        data  = filter(full_grid, is_empty),
        fill  = "white", color = "white", linewidth = 0.75
      ) +
      geom_text(
        data = filter(full_grid, !is_empty),
        aes(label = sprintf("%.2f", mean)),
        size = 3, color = "grey20"
      ) +
      coord_fixed(ratio = 1) +
      scale_x_discrete(position = "bottom", limits = x_stages) +
      scale_y_discrete(limits = rev(y_stages)) +
      labs(title = meta$title, x = NULL, y = NULL) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title            = element_text(face = "plain", size = 11, hjust = 0.5,
                                             margin = margin(b = 6)),
        axis.text.x           = element_text(angle = 35, hjust = 1, size = 9),
        axis.text.y           = if (TRUE) element_text(size = 9) else element_blank(),
        panel.grid            = element_blank(),
        legend.position       = "bottom",
        legend.key.width      = unit(1, "cm"),
        legend.key.height     = unit(0.35, "cm"),
        # legend.key.width  = unit(0.35, "cm"),   # narrow
        # legend.key.height = unit(1.5,  "cm"),   # tall
        legend.title.position = "top",
        legend.title          = element_text(size = 9, hjust = 0.5),
        legend.text           = element_text(size = 8),
        plot.margin           = margin(8, 4, 8, 4)
      )
    
    if (meta$diverging) {
      p <- p + scale_fill_gradient2(
        low      = "#0571b0", # "#018571",  # "#2166ac",
        mid      = "#f7f7f7",  # "white",
        high     = "#e66101",  # "#a6611a",  # "#b2182b",
        midpoint = 0.5,
        name     = "Percent of βsor",
        limits   = c(0, 1),
        breaks   = c(0, 0.25, 0.50, 0.75, 1.00),
        labels   = c("0%", "25%", "50%", "75%", "100%"),
        na.value = "white"
      )
    } else {
      p <- p + scale_fill_gradient(
        low      = meta$low,
        high     = meta$high,
        name     = "Mean βsor",
        limits   = c(0, 0.16),
        breaks   = c(0, 0.05, 0.10, 0.15),
        labels   = c("0.00", "0.05", "0.10", "0.15"),
        na.value = "grey85"
      )
    }
    
    p
  }
  
  # ── Build panels ──────────────────────────────────────────────────────────────
  p_sor <- make_panel("beta.sor", plot_df, panel_meta[["beta.sor"]])
  p_sim <- make_panel("beta.sim", plot_df, panel_meta[["beta.sim"]])
  p_sne <- make_panel("beta.sne", plot_df, panel_meta[["beta.sne"]])
  
  # ── Combine with patchwork ────────────────────────────────────────────────────
  fig_beta_dcomp <- (p_sor / p_sim / p_sne) +
    plot_layout(guides = "collect")
  # +
  #   plot_annotation(
  #     theme = theme(
  #       plot.title            = element_text(face = "plain", size = 13, hjust = 0.5),
  #       plot.subtitle         = element_text(size = 10, hjust = 0.5, color = "grey40"),
  #       legend.position       = "right",
  #       legend.title.position = "top",
  #       legend.title          = element_text(size = 9, hjust = 0.5),
  #       legend.text           = element_text(size = 8)
  #       ,
  #       legend.key.width      = unit(1, "cm"),
  #       legend.key.height     = unit(0.35, "cm")
  #       # legend.key.width  = unit(0.35, "cm"),   # narrow
  #       # legend.key.height = unit(1.5,  "cm")   # tall
  #     )
  #   )
  print(fig_beta_dcomp)
  
  # ── Save ──────────────────────────────────────────────────────────────────────
  ggsave("data/cache/figs/fig_3C.pdf", fig_beta_dcomp,
         width = 6, height = 8, device = cairo_pdf)
}


# Figures:
fig_3AB = (
  fig_PCoA_sites
) + (
  fig_richness + theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 35, hjust = 1),
  )
) + plot_layout(widths = c(3, 2)) + plot_annotation(tag_levels = "A")
fig_3AB


ggsave("data/cache/figs/fig_3AB.pdf", fig_3AB, device = cairo_pdf,
       width = 8, height = 5.5)




