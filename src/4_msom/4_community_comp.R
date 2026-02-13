# 4_community_comp.R ####################################################################################
# Community composition analyses
#
# OUTPUT:
out_cache_dir  = "data/cache/4_msom/4_community_comp"
#
# INPUT:
path_msom = "data/cache/models/msom_nofp_all_2026-02-12_19:37:00.rds"
path_trait_data = "data/cache/trait_data/trait_data.csv"
path_unit_key = "data/unit_key.csv"
##################################################################################################################

source("src/global.R")

# if (!dir.exists(dirname(path_out))) dir.create(dirname(path_out), recursive = TRUE)

# Load data for multi-species occupancy model --------------------------------------------------

message("Loading species trait data from ", path_trait_data)
species_traits = read_csv(path_trait_data, show_col_types = FALSE)

message("Loading site unit key from ", path_unit_key)
unit_key = read_csv(path_unit_key, show_col_types = FALSE) %>% rename(site = unit) %>% mutate(site = str_to_lower(site))

message("Loading data for multi-species occupancy model ", path_msom)
model_data = readRDS(path_msom)

msom_summary = model_data$msom_summary
msom = model_data$msom
groups = model_data$groups %>% arrange(common_name)
sites = model_data$sites
species = model_data$species
seasons = c("2020", "2021", "2022", "2023") # TODO: get from model_data

# Indicator species analysis -----------------------------------------------
# https://cran.r-project.org/web/packages/indicspecies/vignettes/IndicatorSpeciesAnalysis.html

message("Conducting indicator species analysis")

## Indicator species analysis
library(indicspecies)
nperm = 10000

# Simplify posterior occurrence site x species matrix as the most probable occupancy state across all draws and seasons
z = msom$sims.list$z
z_binary = apply(z, c(2, 4), function(x) {
  ifelse(mean(x) >= 0.5, 1, 0)
})
dim(z_binary)

occ = as.data.frame(z_binary)
rownames(occ) = sites
colnames(occ) = species
head(occ)

strata = as.factor(unit_key$stratum[ match(sites, unit_key$site) ])

# Single-stratum indicator analysis 
message("Starting single-stratum indicator analysis (", Sys.time(), ")"); t_start = Sys.time()
indval_singlegroup = multipatt(occ, strata, duleg = TRUE, control = how(nperm=nperm))
message("Finished (", Sys.time() - t_start, ")")

# Species with significant p-values are strongly associated with a stratum or combination of strata
# Specificity: "Component ‘A’ is sample estimate of the probability that the surveyed site belongs to the target site group given the fact that the species has been found. This conditional probability is called the specificity or positive predictive value of the species as indicator of the site group." Species with A = 1.0 occur exclusively in sites of a given stratum. In other words, species with high A values are largely restricted to a given stratum. 
# Fidelity: "Component ‘B’ is sample estimate of the probability of finding the species in sites belonging to the site group. This second conditional probability is called the fidelity or sensitivity of the species as indicator of the target site group." Species with B = 1.0 appear in all sites belonging to a given stratum, i.e. it is widespread in the stratum.
# Therefore, white-crowned sparrow is a good indicator of STAND INIT because it almost exclusively occurs in STAND INIT sites (A > 0.95), and most sites belonging to STAND INIT include it (B = 0.75).
# By contrast, rufous hummingbird appears in all STAND INIT sites (B = 1.0), though it is not completely restricted to STAND INIT (A = 0.48).
# Species show significant associations only with STAND INIT and MATURE strata, not COMP EXCL or THINNED.
# Several STAND INIT-associated species occur more exclusively in STAND INIT sites (high A), while most MATURE-associated species do not exclusively occur in MATURE sites (low A), though they appear in nearly all MATURE sites (high B)
summary(indval_singlegroup, indvalcomp=TRUE)

# NOTE: Species like brown creeper, pileated woodpecker, has similar IndVal values across multiple groups (e.g. THINNED and MATURE), causing them to jump between group assignments in single group analyses due to stocasticity between runs. USING A MULTI-GROUP ANALYSIS IS MORE APPROPRIATE HERE

# Show results for all species, regardless of test significance
# Note that of the 72 species in our dataset, only 31 are selected for analysis. The remaining 41 species
# "have their highest IndVal value for the set of all sites. In other words, those species occur in sites
# belonging to all groups. The association with the set of all sites cannot be statistically tested, because
# there is no external group for comparison."
summary(indval_singlegroup, alpha=1)
# The species that are not selected can be seen via $sign with NA p.value
indval_singlegroup$sign

# Multi-stratum indicator analysis (this is less ecologically clear?)
message("Starting multi-stratum indicator analysis (", Sys.time(), ")"); t_start = Sys.time()
indval_multigroup = multipatt(occ, strata, duleg = FALSE, control = how(nperm=nperm))
message("Finished (", Sys.time() - t_start, ")")
# Indicator species associated with a single stratum are specialists of that stratum, while
# those associated with a combination of strata are generalists across those strata.
summary(indval_multigroup, indvalcomp=TRUE)
stop("DEBUG")
indval_multigroup$sign

# Composition PERMANOVA and NMDS ---------------------------------------------------------------------

## Test for significant differences in community composition of different strata
z = msom$sims.list$z # Simplify posterior occurrence site x species matrix as the most probable occupancy state across all draws and seasons
z_binary = apply(z, c(2, 4), function(x) {
  ifelse(mean(x) >= 0.5, 1, 0)
})
dim(z_binary)

# PERMANOVA
# Model R2 is the % of variation in species composition across sites explained by strata as a predictor
# If Pr(>F) < 0.05, differences in composition between strata are statistically significant (at least one stratum differs from the others)
# If not, null hypothesis is true: there is no difference in species composition among strata
res_adonis2 = adonis2(z_binary ~ strata, method="jaccard")

# Pairwise PERMANOVA
strata_levels = unique(strata)
pairwise_combinations = combn(strata_levels, 2, simplify = FALSE)
results = data.frame(Group1=character(),
                      Group2=character(),
                      R2=numeric(),
                      F=numeric(),
                      p.value=numeric(),
                      stringsAsFactors = FALSE)

for (pair in pairwise_combinations) {
  # Subset data for this pair
  idx = strata %in% pair
  z_sub = z_binary[idx, ]
  strata_sub = strata[idx]
  
  # Run adonis2
  ad = adonis2(z_sub ~ strata_sub, method="jaccard", permutations=999)
  
  # Store results
  results = rbind(results, data.frame(
    Group1 = pair[1],
    Group2 = pair[2],
    R2 = ad$R2[1],
    F = ad$F[1],
    p.value = ad$`Pr(>F)`[1]
  ))
}
# Adjust p-values for multiple testing (optional)
results$p.adj = p.adjust(results$p.value, method="fdr")

# Comparisons with small R2 host largely similar assemblages, while those
# with large R2 have very different assemblages.
# Even a small effect (low R2) can be highly significant (i.e. detectable) if groups are consistently different,
# but these differences are minor when compared to differences between strata with higher R2.
results

# NMDS ----------------------------------------------------------------------------------

## Use simplified posterior occurrence site x species matrix
# Simplify posterior occurrence site x species matrix as the most probable occupancy state across all draws and seasons
z = msom$sims.list$z
z_binary = apply(z, c(2, 4), function(x) {
  ifelse(mean(x) >= 0.5, 1, 0)
})
dim(z_binary)

occ = z_binary
rownames(occ) = sites
colnames(occ) = species
head(occ)

# Exclude species with little or no variation in occurrence
sd_threshold = 0.05
species_sd = apply(occ, 2, sd)
species_for_nmds = species_sd > sd_threshold
message("Excluding species with very little variation:")
print(species[species_sd <= sd_threshold])
occ_filtered = occ[, species_for_nmds]

occ_filtered = matrix(as.integer(occ_filtered), nrow = nrow(occ_filtered))
dimnames(occ_filtered) = list(sites, species[species_for_nmds])

nmds_binary_k2 = metaMDS(occ_filtered, distance = "jaccard", trymax = 200, k=2)
nmds_binary_k2 # Stress should be < 0.2 for interpretable results

# Stress plot diagnostic
# Monotone increasing trend → NMDS preserved rank-order distances
# Scatter / residuals → large deviations indicate poor fit
# No systematic curvature → suggests axes adequately represent dissimilarity
stressplot(nmds_binary_k2)

strata = unit_key$stratum[ match(sites, unit_key$site) ]
strata_cols = c(
  "COMP EXCL"  = "forestgreen",
  "MATURE"     = "tan4",
  "STAND INIT" = "orange",
  "THINNED"    = "purple"
)
colors = strata_cols[strata]

# Convex hulls
# Species points are centroids in the ordination space, showing average position of the sites
# where the species occurs. Species near each other are commonly found together, and species
# near a group of sites are typical of that group.
ordiplot(nmds_binary_k2,type="n")
ordihull(nmds_binary_k2, groups = strata, draw = "polygon", col = strata_cols, label = FALSE)
points(nmds_binary_k2, display="sites",col=colors, pch = 16, cex=1)
orditorp(nmds_binary_k2,display="species",col="grey",air=0.01, cex=1)

# 95% SD elipses
ordiplot(nmds_binary_k2,type="n")
ordiellipse(nmds_binary_k2, groups = strata, draw = "polygon", kind = "sd", conf = 0.95, lwd = 2, col = NA, border = strata_cols, label = FALSE)
points(nmds_binary_k2, display="sites",col=colors, pch = 16, cex=0.75)

# Clusters
dist_mat = vegdist(occ_filtered, method = "jaccard")  # Distance matrix
hc = hclust(dist_mat, method = "average")  # Hierarchical clustering with average linkage
ordiplot(nmds_binary_k2,type="n")
ordiellipse(nmds_binary_k2, groups = strata, draw = "polygon", kind = "sd", conf = 0.95, lwd = 2, col = NA, border = strata_cols, label = FALSE)
ordicluster(nmds_binary_k2, hc, col = "gray", lwd = 1)
points(nmds_binary_k2, display="sites",col=colors, pch = 16, cex=1)

## Fit with k = 3
nmds_binary_k3 = metaMDS(occ_filtered, distance = "jaccard", trymax = 200, k=3)
nmds_binary_k3 # Stress should be < 0.2 for interpretable results

# Visualize in 3d
# library(rgl)
# site_scores = scores(nmds_binary_k3, display="sites")
# plot3d(site_scores[,1],
#        site_scores[,2],
#        site_scores[,3],
#        col=colors,
#        size=1,
#        type="s",
#        xlab="NMDS1",
#        ylab="NMDS2",
#        zlab="NMDS3")

# 2D plot showing only axes 1 and 2
ordiplot(nmds_binary_k3,type="n")
ordiellipse(nmds_binary_k3, groups = strata, draw = "polygon", kind = "sd", conf = 0.95, lwd = 2, col = NA, border = strata_cols, label = FALSE)
points(nmds_binary_k3, display="sites",col=colors, pch = 16, cex=0.75)
orditorp(nmds_binary_k3,display="species",col="darkgray",air=0.1, cex=0.75)

## Use posterior mean occurrence probability for NMDS instead (site similarity in expected composition)

# Get posterior mean occurrence probability site x species matrix across all seasons
str(z)
psi_mean = apply(z, c(2,4), mean)
rownames(psi_mean) = sites
colnames(psi_mean) = species

nmds_occprob_k2 = metaMDS(psi_mean, distance = "bray", trymax = 200, k=2)
nmds_occprob_k2

nmds_occprob_k3 = metaMDS(psi_mean, distance = "bray", trymax = 200, k=3)
nmds_occprob_k3

ordiplot(nmds_occprob_k3,type="n")
ordiellipse(nmds_occprob_k3, groups = strata, draw = "polygon", kind = "sd", conf = 0.95, lwd = 2, col = NA, border = strata_cols, label = FALSE)
points(nmds_occprob_k3, display="sites",col=colors, pch = 16, cex=0.75)
orditorp(nmds_occprob_k3,display="species",col="darkgray",air=0.1, cex=0.75)

# Estimate richness -------------------------------------------------------------------------

z = msom$sims.list$z
samples = dim(z)[1]
J = dim(z)[2] # n sites
T = dim(z)[3] # n seasons
I = dim(z)[4] # n species

message("Calculating richness for each grouping")
site_group_richness = tibble()
for (grouping in c("group_all")) {
  
  species_group = tibble(
    common_name = species_traits$common_name,
    group = species_traits[[grouping]],
    group_idx = as.integer(as.factor(group))
  ) %>%
  filter(common_name %in% species) %>%
  slice(match(species, common_name))
  
  G = length(unique(species_group$group))
  
  rich_group = array(NA, c(samples, J, T, G))
  for (g in unique(species_group$group)) {
    species_in_g = which(species %in% (species_group %>% filter(group == g) %>% pull(common_name)))
    gi = (species_group %>% filter(group == g) %>% pull(group_idx))[1]
    rich_group[ , , , gi] = apply(z[ , , , species_in_g, drop = FALSE], c(1,2,3), sum)
  }
  rich_group_mean  = apply(rich_group, c(2,3,4), mean)
  rich_group_lower = apply(rich_group, c(2,3,4), quantile, probs = 0.025)
  rich_group_upper = apply(rich_group, c(2,3,4), quantile, probs = 0.975)
  
  site_richness = apply(rich_group_mean, c(1,3), mean) # mean site richness per site averaged across years
  year_richness = apply(rich_group_mean, c(2,3), mean) # mean site richness per year averaged across sites
  
  d = tibble(
    site = sites,
    rich = as.vector(site_richness)
  ) %>% left_join(unit_key, by = "site")
  
  ggplot(d, aes(x = stratum, y = rich)) +
    geom_boxplot()
}
