old_plot = read_rds("data/cache/3_gis/3_calc_occurrence_vars/OPT_data_plot_scale_2020_clean_strata_4_sites.rds")
new_plot = read_rds("data/cache/3_gis/3_calc_occurrence_vars/V2_data_plot_scale_2020_clean_strata_4.rds")

stopifnot(sort(names(old_plot)) == sort(names(new_plot)))

new_plot2 <- new_plot[, names(old_plot)]
all.equal(old_plot, new_plot[, names(old_plot)])
identical(old_plot, new_plot2)

old_hr = read_rds("data/cache/3_gis/3_calc_occurrence_vars/V2_data_homerange_scale_2020_clean_strata_4.rds")
new_hr = read_rds("data/cache/3_gis/3_calc_occurrence_vars/OPT_data_homerange_scale_2020_clean_strata_4_sites.rds")

stopifnot(sort(names(old_hr)) == sort(names(new_hr)))

compare_list <- function(old, new) {
  results <- lapply(names(old), function(nm) {
    
    print(nm)
    
    x <- old[[nm]]
    y <- new[[nm]]
    
    # If data.frame → align columns
    if (is.data.frame(x) && is.data.frame(y)) {
      y <- y[, names(x), drop = FALSE]
    }
    
    # Compare
    res <- all.equal(x, y)
    
    list(name = nm, equal = isTRUE(res), msg = res)
  })
  
  results
}
res <- compare_list(old_hr, new_hr)
all(sapply(res, `[[`, "equal"))
Filter(function(x) !x$equal, res)


nm <- "yellow warbler"

x <- old_hr[[nm]]
y <- new_hr[[nm]][, names(x)]

# example column
col <- "focalpatch_area_homeange_pcnt"

summary(x[[col]] - y[[col]])
