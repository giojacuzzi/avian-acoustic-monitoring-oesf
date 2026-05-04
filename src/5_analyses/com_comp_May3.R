##############################################################################
# Community composition analyses, propagating posterior uncertainty
#
# CONFIG:
binary_not_psi = TRUE # Summarize latent occurrence states across years as binary (max) or probabilities (psi)
thin = 20 # Number of iterations to thin from the z posterior to speed up computation (0 or more)
nperm = 99 # TODO: 999 (for p value) or 0 permutations
# INPUT:
path_msom                                = "data/cache/models/prefinal_msom_jags_nofp_all.rds"
path_trait_data                          = "data/cache/2_traits/1_agg_traits/trait_data.csv"
path_occurrence_predictor_plot_data      = "data/cache/4_msom/1_assemble_msom_data/V3_occurrence_predictor_plot_data.rds"
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
  message("Thinning z posterior to every ", thin, " iterations")
  thin_idx = seq(1, n_iter, by = thin)
  z_raw = z_raw[thin_idx, , , ]
  n_iter = length(thin_idx)
  str(z_raw)
}

message("Loading species trait data from ", path_trait_data)
species_traits = read_csv(path_trait_data, show_col_types = FALSE) %>% filter(common_name %in% species)

message("Loading occurrence predictor plot scale data from ", path_occurrence_predictor_plot_data)
occurrence_predictor_plot_data = readRDS(path_occurrence_predictor_plot_data)[[1]] %>% filter(site %in% sites)
stopifnot(all(occurrence_predictor_plot_data$site == sites))

# Summarize community posterior ---------------------------------------------------------------

if (binary_not_psi) {
  # Q: "Did the species ever use the site?"
  # Summarize binary occurrence state across years as 1 if a site was
  # occupied in ANY season, and 0 if it was unoccupied in all seasons.
  message("z: Summarizing binary occurrence state across years")
  z = array(dim = c(n_iter, n_sites, n_species))
  for (i in seq_len(n_species)) {
    acc = matrix(0, nrow = n_iter, ncol = n_sites)
    for (t in seq_len(n_seasons)) {
      acc = pmax(acc, z_raw[, , t, i])
    }
    z[, , i] = acc
  }
} else {
  # Q: "How consistently did the species use the site across years?"
  # Average posterior occurrence across years as the mean probability of occurrence across
  # years, propagating posterior uncertainty from seasonal estimates and therefore accounting
  # for seasonal stocasticity.
  message("z: Averaging posterior occurrence states across years as probabilities")
  z = array(dim = c(n_iter, n_sites, n_species))
  for (i in seq_len(n_species)) {
    acc = matrix(0, nrow = n_iter, ncol = n_sites)
    for (t in seq_len(n_seasons)) {
      acc = acc + z_raw[, , t, i]
    }
    z[, , i] = acc / n_seasons
  }
}
str(z)

# TODO:
# Q: "Was the species present at the site more often than not?"
# Alternative -- use a threshold (z_t = ifelse(z > 0.5, 1, 0)), which is compatible with downstream betapart analyses

# stop("READY")

# PCoA and NMDS visualization ----------------------------------------------------------------

# For visualization only, compute dissimilary matricies for each posterior iteration, then take
# the posterior mean dissimlarity matrix. This propagates uncertainty through the dissimilarity
# computation itself, and the mean dissimilarity matrix reflects how sites differ on average
# across the posterior.
{
  dissim_matricies_post = matrix(0, nrow = n_sites, ncol = n_sites)
  pb = progress_bar$new(format = progress_bar_format, total = n_iter, clear = FALSE)
  for (i in seq_len(n_iter)) {
    dissim_matricies_post = dissim_matricies_post + as.matrix(vegdist(z[i, , ], method = "bray", binary = binary_not_psi)) # Sorensen if binary = TRUE and method = "bray"
    pb$tick()
  }
  diss_matrix_mean = dissim_matricies_post / n_iter
}

## PCoA
{
  pcoa = cmdscale(as.dist(diss_matrix_mean), k = 2, eig = TRUE, add = TRUE) # 2 dimensional pcoa with centered (positive) eigenvalues
  pcoa_pcnt = (pcoa$eig / sum(pcoa$eig) * 100) %>% round(1)
  pcoa_df = as.data.frame(pcoa$points) %>% setNames(c("PCoA1", "PCoA2")) %>% mutate(site = sites, stage = stages)
  
  fig_pcoa_sites = ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = stage)) +
    geom_vline(xintercept = 0, color = "gray90") + geom_hline(yintercept = 0, color = "gray90") +
    geom_point(size = 2, alpha = 0.8) + stat_ellipse(aes(fill = stage), geom = "polygon", alpha = 0.1, level = 0.95) +
    scale_color_manual(values = colors_stats) + scale_fill_manual(values = colors_stats) +
    labs(x = paste0("PCoA 1 (", pcoa_pcnt[1], "%)"), y = paste0("PCoA 2 (", pcoa_pcnt[2], "%)")) +
    theme(legend.position = "none") + labs(title = "pcoa")
}
print(fig_pcoa_sites)

## NMDS (note that k = 3)
{
  nmds_k3 = metaMDS(diss_matrix_mean, distance = "bray", trymax = 500, k=3)
  nmds_k3$stress; stressplot(nmds_k3)
  nmds_df = as_tibble(scores(nmds_k3, display = "sites")) %>% mutate(stage = stages, site = sites)
  
  fig_nmds_color = ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = stage)) +
    geom_point(size = 2, alpha = 0.8) + stat_ellipse(level = 0.95, aes(fill = stage), geom = "polygon") +
    # geom_text(aes(label = site), size = 2.5, vjust = -0.5) +
    scale_color_manual(values = colors_stats) + scale_fill_manual(values = alpha(colors_stats, 0.1)) +
    theme(legend.position = "none") + labs(title= "nmds_k3")
}
print(fig_nmds_color)

# Q1: Does assemblage composition vary among stages? ---------------------------------------------------

# Global frequentist PERMANOVA
{
  R2_post = numeric(n_iter) # Proportion of variance explained by stage
  F_post  = numeric(n_iter) # Larger pseudo-F indicates greater difference between groups relative to variation within groups
  pval_post  = numeric(n_iter) # Significant frequentist difference exists
  pb = progress_bar$new(format = progress_bar_format, total = n_iter, clear = FALSE)
  for (i in seq_len(n_iter)) {
    d          = vegdist(z[i, , ], method = "bray", binary = binary_not_psi)
    fit        = adonis2(d ~ stages, permutations = nperm)
    R2_post[i] = fit$R2[1]
    F_post[i]  = fit$F[1]
    pval_post[i]  = fit$`Pr(>F)`[1]
    pb$tick()
  }
  permanova = data.frame(
    stat      = c("R2", "F", "pval"),
    mean      = c(mean(R2_post),             mean(F_post),             mean(pval_post)),
    median    = c(median(R2_post),           median(F_post),           median(pval_post)),
    bci_2.5   = c(quantile(R2_post, 0.025, names = FALSE),  quantile(F_post, 0.025, names = FALSE),  quantile(pval_post, 0.025, names = FALSE)),
    bci_97.5  = c(quantile(R2_post, 0.975, names = FALSE),  quantile(F_post, 0.975, names = FALSE),  quantile(pval_post, 0.975, names = FALSE))
  )
}
print(permanova, digits = 2)

# Pairwise frequentist PERMANOVAs
{
  stage_pairs = combn(levels(stages), 2, simplify = FALSE)
  pair_names  = sapply(stage_pairs, paste, collapse = " vs ")
  
  pairwise_R2_post = matrix(0, nrow = n_iter, ncol = length(stage_pairs), dimnames = list(NULL, pair_names))
  pairwise_F_post  = matrix(0, nrow = n_iter, ncol = length(stage_pairs), dimnames = list(NULL, pair_names))
  pairwise_p_post  = matrix(0, nrow = n_iter, ncol = length(stage_pairs), dimnames = list(NULL, pair_names))
  
  pb = progress_bar$new(format = progress_bar_format, total = n_iter, clear = FALSE)
  for (i in seq_len(n_iter)) {
    d = vegdist(z[i, , ], method = "bray", binary = binary_not_psi)
    for (p in seq_along(stage_pairs)) {
      idx       = which(stages %in% stage_pairs[[p]])
      d_sub     = as.dist(as.matrix(d)[idx, idx])
      fit       = adonis2(d_sub ~ stages[idx], permutations = nperm)
      pairwise_R2_post[i, p] = fit$R2[1]
      pairwise_F_post[i, p]  = fit$F[1]
      pairwise_p_post[i, p]  = fit$`Pr(>F)`[1]
    }
    pb$tick()
  }
  
  pairwise_permanova = bind_rows(lapply(seq_along(stage_pairs), function(p) {
    data.frame(
      comparison = pair_names[p],
      R2_mean    = mean(pairwise_R2_post[, p]),
      R2_bci_2.5 = quantile(pairwise_R2_post[, p], 0.025),
      R2_bci_97.5= quantile(pairwise_R2_post[, p], 0.975),
      F_mean     = mean(pairwise_F_post[, p]),
      F_bci_2.5  = quantile(pairwise_F_post[, p], 0.025),
      F_bci_97.5 = quantile(pairwise_F_post[, p], 0.975),
      pval_mean     = mean(pairwise_p_post[, p]),
      pval_bci_2.5  = quantile(pairwise_p_post[, p], 0.025),
      pval_bci_97.5 = quantile(pairwise_p_post[, p], 0.975)
    )
  }))
}
print(pairwise_permanova, digits = 2)

# Bayesian posterior pairwise_dissim_contrasts: Posterior probability that dissimilarities between sites of
# different stages are greater than dissimilarities between sites of the same stage
{
  iter_sample = 1:n_iter # sample(1:n_iter, 1000)
  
  # For each pair, index of between-pair and within-pair (both stages) sites
  stage_levels = levels(stages)
  stage_pairs  = combn(stage_levels, 2, simplify = FALSE)
  pair_idx = lapply(stage_pairs, function(p) {
    idx1 = which(stages == p[1])
    idx2 = which(stages == p[2])
    list(
      between = expand.grid(i = idx1, j = idx2),
      within  = bind_rows(
        expand.grid(i = idx1, j = idx1) %>% filter(i < j),
        expand.grid(i = idx2, j = idx2) %>% filter(i < j)
      )
    )
  })
  names(pair_idx) = sapply(stage_pairs, paste, collapse = " vs ")
  
  # Global indices
  stage_mat   = outer(stages, stages, FUN = "==")
  within_idx  = which(stage_mat  & lower.tri(stage_mat))
  between_idx = which(!stage_mat & lower.tri(stage_mat))
  
  global_diff = numeric(length(iter_sample))
  pair_diffs  = matrix(NA, nrow = length(iter_sample), ncol = length(stage_pairs),  dimnames = list(NULL, names(pair_idx)))
  within_raw  = matrix(NA, nrow = length(iter_sample), ncol = length(stage_levels), dimnames = list(NULL, stage_levels))
  between_raw = matrix(NA, nrow = length(iter_sample), ncol = length(stage_pairs),  dimnames = list(NULL, names(pair_idx)))
  
  pb = progress_bar$new(format = progress_bar_format, total = length(iter_sample), clear = FALSE)
  for (i in seq_along(iter_sample)) {
    k    = iter_sample[i]
    diss = as.matrix(vegdist(z[k, , ], method = "bray", binary = binary_not_psi))

    global_diff[i] = mean(diss[between_idx]) - mean(diss[within_idx])
    for (j in seq_along(pair_idx)) {
      idx = pair_idx[[j]]
      pair_diffs[i, j]  = mean(diss[cbind(idx$between$i, idx$between$j)]) -
        mean(diss[cbind(idx$within$i,  idx$within$j)])
      between_raw[i, j] = mean(diss[cbind(idx$between$i, idx$between$j)])
    }
    
    # Within-group raw dissimilarities
    for (s in stage_levels) {
      idx_s = which(stages == s)
      d_s   = diss[idx_s, idx_s]
      within_raw[i, s] = mean(d_s[lower.tri(d_s)])
    }
    pb$tick()
  }
  
  summarise_diff = function(x, name) {
    data.frame(
      comparison      = name,
      P_greater       = mean(x > 0),
      median_diff     = median(x),
      lower_95        = unname(quantile(x, 0.025)),
      upper_95        = unname(quantile(x, 0.975))
    )
  }
  
  summarise_raw = function(mat) {
    as.data.frame(t(apply(mat, 2, function(x) c(
      median   = median(x),
      lower_95 = unname(quantile(x, 0.025)),
      upper_95 = unname(quantile(x, 0.975))
    ))))
  }
  
  pairwise_dissim_contrasts = do.call(rbind, c(
    list(summarise_diff(global_diff, "Global (all between vs all within)")),
    lapply(names(pair_idx), function(nm) summarise_diff(pair_diffs[, nm], nm))
  ))
  
  # print(summarise_raw(within_raw), digits = 2) # Note this is the same as within_bp_df %>% filter(component == "beta.sor")
  # ggplot(summarise_raw(within_raw) %>% rownames_to_column("stage"), aes(x = median, y = stage, xmin = lower_95, xmax = upper_95)) + geom_pointrange()
  # print(summarise_raw(between_raw), digits = 2)
}
print(pairwise_dissim_contrasts, digits = 2)

# Frequentist dispersion (multivariate homogeneity of group dispersions)
{
  dispersion_post = matrix(NA, nrow = n_iter, ncol = nlevels(stages))
  colnames(dispersion_post) = levels(stages)
  dispersion_p <- numeric(n_iter)
  
  pb = progress_bar$new(format = progress_bar_format, total = n_iter, clear = FALSE)
  for (i in seq_len(n_iter)) {
    d  = vegdist(z[i, , ], method = "bray", binary = binary_not_psi)
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
  
  pairwise_disp_contrasts = do.call(rbind, lapply(stage_pairs, function(p) {
    diff = dispersion_post[, p[1]] - dispersion_post[, p[2]]
    data.frame(
      stage_1  = p[1],
      stage_2  = p[2],
      mean_diff    = mean(diff),
      bci_2.5      = quantile(diff, 0.025),
      bci_97.5     = quantile(diff, 0.975),
      P_direction  = mean(diff > 0)
    )
  })) %>% remove_rownames()
}
print(dispersion_df, digits = 3)
data.frame(P_pval_lt_0.05 = mean(dispersion_p < 0.05), pval_mean = mean(dispersion_p), bci_2.5 = quantile(dispersion_p, 0.025, names = FALSE), bci_97.5 = quantile(dispersion_p, 0.975, names = FALSE)) %>% print(digits = 2)
print(pairwise_disp_contrasts, digits = 3)

# Q2: Is variation in composition shaped by turnover or nestedness? ------------------------------------

# Global turnover vs nestedness decomposition
# NOTE: The betapart abundance framework does not directly calculate nestedness -- using psi here instead reflects expected richness differences!
# CONSIDER: 
{
  beta_multi_post = matrix(NA, nrow = n_iter, ncol = 3)

  pb = progress_bar$new(format = progress_bar_format, total = n_iter, clear = FALSE)
  if (binary_not_psi) {
    colnames(beta_multi_post) = c("beta.SOR", "beta.SIM", "beta.SNE")
    for (i in seq_len(n_iter)) {
      bp_core            = betapart.core(z[i, , ])
      bm                 = beta.multi(bp_core, index.family = "sorensen")
      beta_multi_post[i, ] = c(bm$beta.SOR, bm$beta.SIM, bm$beta.SNE)
      pb$tick()
    }
  } else {
    colnames(beta_multi_post) = c("beta.BRAY", "beta.BRAY.BAL", "beta.BRAY.GRA")
    for (i in seq_len(n_iter)) {
      bp_core            = betapart.core.abund(z[i, , ])
      bm                 = beta.multi.abund(bp_core, index.family = "bray")
      beta_multi_post[i, ] = c(bm$beta.BRAY, bm$beta.BRAY.BAL, bm$beta.BRAY.GRA)
      pb$tick()
    }
  }
  
  betapart_global = do.call(rbind, lapply(colnames(beta_multi_post), function(nm) {
    data.frame(
      component = nm,
      mean      = mean(beta_multi_post[, nm]),
      median    = median(beta_multi_post[, nm]),
      bci_2.5   = unname(quantile(beta_multi_post[, nm], 0.025)),
      bci_97.5  = unname(quantile(beta_multi_post[, nm], 0.975))
    )
  }))
}
print(betapart_global, digits = 2)

# Pairwise turnover vs nestedness decomposition
# NOTE: The betapart abundance framework does not directly calculate nestedness -- using psi here instead reflects expected richness differences!
{
  stage_pairs = combn(levels(stages), 2, simplify = FALSE)
  pair_names  = sapply(stage_pairs, paste, collapse = " vs ")
  
  pb = progress_bar$new(format = progress_bar_format, total = n_iter, clear = FALSE)
  if (binary_not_psi) {
    
    components  = c("beta.sor", "beta.sim", "beta.sne")
    between_bp = array(NA, dim = c(n_iter, length(stage_pairs), 3), dimnames = list(NULL, pair_names, components))
    within_bp  = array(NA, dim = c(n_iter, nlevels(stages), 3), dimnames = list(NULL, levels(stages), components))
    
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
  } else {
    
    components  = c("beta.bray", "beta.bray.bal", "beta.bray.gra")
    between_bp = array(NA, dim = c(n_iter, length(stage_pairs), 3), dimnames = list(NULL, pair_names, components))
    within_bp  = array(NA, dim = c(n_iter, nlevels(stages), 3), dimnames = list(NULL, levels(stages), components))
    
    for (i in seq_len(n_iter)) {
      bp_core = betapart.core.abund(z[i, , ])
      bp_pair = beta.pair.abund(bp_core, index.family = "bray")
      
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
  }
  
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
  
  if (binary_not_psi) {
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
  } else {
    format_bp = function(df, group_col) {
      df %>%
        select({{ group_col }} := group, component, median) %>%
        pivot_wider(names_from = component, values_from = median) %>%
        rename(bray = beta.bray, bal = beta.bray.bal, gra = beta.bray.gra) %>%
        mutate(
          pct_turnover  = round(100 * bal / bray, 1),
          pct_nestedness = round(100 * gra / bray, 1)
        )
    }
  }
  
  betapart_between = format_bp(between_bp_df, group_col = "pair")
  betapart_within  = format_bp(within_bp_df,  group_col = "stage")
  betapart_between
  betapart_within
}
betapart_between
betapart_within

# Q3: Does richness vary among stages? ------------------------------------

# Posterior estimated richness
post_tmean_rich = matrix(0, nrow = n_iter, ncol = n_sites)
for (i in seq_len(n_species)) {
  post_tmean_rich = post_tmean_rich + z[, , i]
}

# Figure: Posterior mean estimated richness for visualization
{
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
    )
  
  # ggsave("data/cache/figs/fig_richness.pdf", fig_richness, width = 6, height = 6)
}
print(fig_richness)

# Bayesian posterior pairwise_richness_contrasts: Posterior probability that richness is greater in
# one stage than another
{
  # For each iteration, average richness across sites within each stage,
  # producing a full posterior for each stage's mean richness.
  stage_levels = levels(stages)
  stage_post = sapply(stage_levels, function(s) {
    rowMeans(post_tmean_rich[, stages == s, drop = FALSE])
  })
  richness_posterior = data.frame(
    stage    = stage_levels,
    mean     = apply(stage_post, 2, mean),
    median   = apply(stage_post, 2, median),
    bci_2.5  = apply(stage_post, 2, quantile, probs = 0.025),
    bci_97.5 = apply(stage_post, 2, quantile, probs = 0.975)
  )
  
  # Difference in stage-specific posteriors
  stage_pairs = combn(stage_levels, 2, simplify = FALSE)
  
  pairwise_richness_contrasts = lapply(stage_pairs, function(pair) {
    diffs = stage_post[, pair[1]] - stage_post[, pair[2]]
    data.frame(
      pair     = paste0(pair[1], " - ", pair[2]),
      mean     = mean(diffs),
      median   = median(diffs),
      bci_2.5  = quantile(diffs, 0.025),
      bci_97.5 = quantile(diffs, 0.975),
      prob_gt0 = round(mean(diffs > 0), 3)   # P(stage B > stage A)
    )
  }) %>% bind_rows() %>% remove_rownames()
}
print(richness_posterior, digits = 4)
print(pairwise_richness_contrasts, digits = 2)

# Additional post-hoc questions ---------------------------------------------------

summarise_cor <- function(x, name) {
  data.frame(
    stage    = name,
    P_pos    = mean(x > 0),
    median_r = median(x),
    lower_95 = quantile(x, 0.025),
    upper_95 = quantile(x, 0.975)
  )
}

# Q: Within each stage, are stands closer in age also more similar in composition?
{
  iter_sample = 1:n_iter
  
  age <- occurrence_predictor_plot_data$age_mean
  
  age_diss_cor <- matrix(NA, nrow = length(iter_sample), ncol = length(stage_levels),
                         dimnames = list(NULL, stage_levels))
  
  pb = progress_bar$new(format = progress_bar_format, total = n_iter, clear = FALSE)
  for (i in seq_along(iter_sample)) {
    k    <- iter_sample[i]
    diss <- as.matrix(vegdist(z[k, , ], method = "bray", binary = binary_not_psi))
    
    for (s in stage_levels) {
      idx_s    <- which(stages == s & !is.na(age))
      age_s    <- age[idx_s]
      age_diff <- as.vector(dist(age_s))
      diss_s   <- diss[idx_s, idx_s]
      diss_vec <- diss_s[lower.tri(diss_s)]
      
      age_diss_cor[i, s] <- cor(age_diff, diss_vec, method = "spearman")
    }
    pb$tick()
  }
  
  cor_results <- do.call(rbind, lapply(stage_levels, function(s) {
    summarise_cor(age_diss_cor[, s], s)
  }))
  
  rownames(cor_results) <- NULL
}
print(cor_results, digits = 3)

# Q: Are thinned assemblages more similar to mature than to stem exclusion?
{
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
    diss <- as.matrix(vegdist(z[k, , ], method = "bray", binary = binary_not_psi))
    
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
    diss <- as.matrix(vegdist(z[k, , ], method = "bray", binary = binary_not_psi))
    
    # Mean dissimilarity from each thinned site to all mature sites
    mean_diss_to_mature <- rowMeans(diss[idx_thin, idx_mat])
    
    fit          <- lm(mean_diss_to_mature ~ age_thin)
    intercept[i] <- coef(fit)[1]
    slope[i]     <- coef(fit)[2]
  }
  
  cat("P(slope < 0):", mean(slope < 0), "\n")
  cat("\nSlope 95% CI:\n")
  print(quantile(slope, c(0.025, 0.5, 0.975)), digits = 3)
}

# Q: Is proportion of another stage within the mean species home range correlated with similarity to that stage?
{
  # canada jay ~ mean home range size
  # barred owl ~ 1km
  scale = "canada jay"
  homerange_data = readRDS(path_occurrence_predictor_homerange_data)[[1]][[scale]]
  unique(homerange_data$buffer_radius_m)

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
    diss <- as.matrix(vegdist(z[k, , ], method = "bray", binary = binary_not_psi))
    
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
      P_neg    = mean(x < 0),
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


