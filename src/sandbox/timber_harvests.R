library(tidyverse)
library(sf)
library(terra)
library(mapview)
options(mapview.maxpixels = 9000000)
library(ggrepel)
theme_set(theme_minimal())
library(landscapemetrics)
library(lwgeom)
library(units)

### Study area

crs_m = 32610 # EPSG:32610, UTM Zone 10N, (meters)

# ARU locations
aru_sites = st_read('data/environment/GIS Data/AcousticStations.shp') %>%
  st_drop_geometry() %>% janitor::clean_names() %>% as.data.frame() %>%
  st_as_sf(coords = c("utm_e", "utm_n"), crs = crs_m) %>%
  select(name, ces, treatment, geometry) %>% rename(site = name)
mapview(aru_sites, label = aru_sites$name)

# Timber sale and harvest unit locations

# all sales in the study area that were completed during the period of interest, with a couple extra years on each end. Each record is a timber sale unit polygon, with multiple units making up one timber sale
# TECHNIQUE_ - this is the type of harvest. Most are VRH (regeneration harvest), some are thinnings (these are called “COMMRCL_THIN” or “VARIABL_THIN” or “LATE_RTN_THIN”).
# FMA_DT – this field is the date that the timber sale was certified as completed. The actual harvesting of the trees in a given unit may have occurred a few months up to maybe a year or so earlier, but usually months.
timber_sales = st_read('data/environment/GIS Data/Harvest/PAMarea_TimberSales/PAMarea_TimberSales.shp') %>% janitor::clean_names() %>% st_transform(crs_m)
timber_sales$fma_dt = as.Date(timber_sales$fma_dt)
timber_sales$fma_yr = as.factor(year(timber_sales$fma_dt))

timber_sales$fma_before_pam = timber_sales$fma_dt < as.Date('2020-04-07') # PAM start date
timber_sales$fma_after_pam  = timber_sales$fma_dt > as.Date('2023-08-14') # PAM end date

mapview(timber_sales %>% filter(fma_before_pam == FALSE, fma_after_pam == FALSE), zcol = 'fma_yr') +
  mapview(aru_sites, label = aru_sites$name)

table(timber_sales %>% filter(fma_before_pam == FALSE, fma_after_pam == FALSE) %>% pull(technique))

# TODO: calculate days since start of 2020 monitoring
# TODO: flag for sales that were completed before the start of the 2020 monitoring. These sites should be "safe", i.e. line up with the raster data



mapview(timber_sales, zcol = 'fma_yr') + mapview(aru_sites, label = aru_sites$name)

# more precise dates for when T3 sale units were actually harvested in 2023
harvest_units = st_read('data/environment/GIS Data/Harvest/T3uplandHarvestUnits/HarvestUnitsUpland.shp') %>% janitor::clean_names() %>% st_transform(crs_m)

# “PAMarea_TimberSales.zip” only has the units from sales that have been completely finished, and some of Emily’s records show sales that are partially completed.
