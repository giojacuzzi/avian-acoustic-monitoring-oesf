library(readxl)
library(ggplot2)
library(dplyr)
library(forcats)

data_species_pool = read_excel('data/species_pool/regional_species_pool.xlsx', sheet = 'included')

data_species_pool['season_start'] = as.Date(sapply(strsplit(data_species_pool$season_breeding, "-"), `[`, 1), format = "%d %b")
data_species_pool['season_end'] = as.Date(sapply(strsplit(data_species_pool$season_breeding, "-"), `[`, 2), format = "%d %b")
data_species_pool = data_species_pool %>%
  mutate(season_mid = as.Date((as.numeric(season_start) + as.numeric(season_end)) / 2, origin = "1970-01-01"))
data_species_pool = data_species_pool %>% mutate(status = recode(residency, "Vagrant" = "Migration / Vagrant", "Migration" = "Migration / Vagrant"))
data_species_pool$status = factor(data_species_pool$status, levels = c('Year-round', 'Breeding', 'Migration / Vagrant', 'Nonbreeding'))

data_species_pool <- data_species_pool %>%
  arrange(status, season_start) %>%
  mutate(common_name = factor(common_name, levels = unique(common_name)))

ggplot(data_species_pool) +
  geom_segment(aes(x = season_start, xend = season_end, y = common_name, color = status), size=2) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  geom_vline(xintercept=as.numeric(as.Date("9 Apr", format = "%d %b"))) +
  geom_vline(xintercept=as.numeric(as.Date("26 Jul", format = "%d %b"))) +
  scale_color_manual(values = c('#633372FF', '#D8443CFF', '#F4C40FFF', '#1F6E9CFF')) +
  labs(title = "Breeding periods (regional species pool)", x = "Date", y = "Species", color = "Status") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.text.y = element_text(size = 4)
  )
ggsave('figs/fig_species_breeding_seasons.pdf', width=12, height=8)
