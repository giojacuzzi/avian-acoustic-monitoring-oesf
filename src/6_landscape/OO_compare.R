path_original = "data/cache/7_landscape/landscape_data_homerange_scale_2020_clean_strata_4_landscape.rds"
path_optimize = "data/cache/7_landscape/OPT_landscape_data_homerange_scale_2020_clean_strata_4_landscape.rds"

org = read_rds(path_original)
opt = read_rds(path_optimize)

# Check dimensions
identical(lapply(org, dim), lapply(opt, dim))

# Check names
identical(names(org), names(opt))

# Check lengths
all(lengths(org) == lengths(opt))

# Exact match
identical(org, opt)
all(mapply(identical, org, opt))

# With numeric tolerance
all.equal(org, opt)
all(sapply(Map(function(a, b) all.equal(a, b), org, opt), isTRUE))

# If any of the above are false, inspect:
{
  # Exact match, element by element
  mapply(identical, org, opt)
  # With numeric tolerance
  Map(function(a, b) all.equal(a, b), org, opt)
}
