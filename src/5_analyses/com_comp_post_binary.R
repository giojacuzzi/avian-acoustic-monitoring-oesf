# 4_community_comp.R #############################################################################
# Community composition analyses, propagating posterior uncertainty
#
# CONFIG:
thin = 0 # Number of iterations to thin from the z posterior to speed up computation (0 or more)
nperm = 99 # TODO: 999 (for p value) or 0 permutations
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
stages  = factor(model_data$stages$stratum_4,  levels = c("standinit", "compex", "thin", "mature"))

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

if (TRUE) {
  message("Summarizing occurrence state across years")
  z = array(NA_integer_, dim = c(n_iter, n_sites, n_species))
  for (i in seq_len(n_species)) {
    acc = matrix(0, nrow = n_iter, ncol = n_sites)
    for (t in seq_len(n_seasons)) {
      acc = pmax(acc, z_raw[, , t, i])
    }
    z[, , i] = acc
  }
  str(z)
  
  binary_dissimilarity = TRUE
  
} else {
  # "We retained posterior distributions of latent occurrence states averaged across years to propagate uncertainty through analyses of community composition."
  # Average over years (dim 3) while retaining iterations, sites, species, using accumulation to avoid memory limits
  message("Averaging posterior occurrence states across years")
  z = array(0, dim = c(n_iter, n_sites, n_species))
  for (i in seq_len(n_species)) {
    acc = matrix(0, nrow = n_iter, ncol = n_sites)
    for (t in seq_len(n_seasons)) {
      acc = acc + z_raw[, , t, i]
    }
    z[, , i] = acc / n_seasons
  }
  str(z)
  
  binary_dissimilarity = FALSE
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
  fig_richness = ggplot(post_tmean_rich_stats, aes(x = stage, y = median, color = stage)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.25, alpha = 0.35, aes(color = stage)) +
    scale_color_manual(values = unname(stage_colors)); print(fig_richness)
  
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

# TODO: Posterior estimated richness (each year)

# Principal coordinates analysis (PCoA) -----------------------------------------------------------
{
  # Average posterior draws of Bray-curtis dissimilarity matrix
  d_sum = matrix(0, n_sites, n_sites)
  pb = progress_bar$new(format = progress_bar_format, total = n_iter, clear = FALSE)
  for (i in seq_len(n_iter)) {
    d_sum = d_sum + as.matrix(vegdist(z[i, , ], method = "bray"))
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

#### ------------------------------------------------------------------------------
# TODO: Q: Do assemblages vary between late-successional and old-growth?

#### ------------------------------------------------------------------------------
# Q: Do smaller patches have higher nestendess?

homerange_data = readRDS(path_occurrence_predictor_homerange_data)[[1]][["median"]]

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