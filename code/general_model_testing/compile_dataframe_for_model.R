# Compile the data for building a model of metabolism
# 2/2022
# A Carter
library(tidyverse)
library(lubridate)
# setwd('C:/Users/alice.carter/git/nhc_50yl/src')

# Load metabolism predictions:
dat <- readRDS("data/metabolism/compiled/met_preds_stream_metabolizer_O2.rds")
sites <- read_csv('data/siteData/NHCsite_metadata.csv') %>%
  slice(c(1:5, 7))
met <- dat$preds

# load additional driver variables:
dd <- tibble()
for(s in sites$sitecode){
  d <- read_csv(paste0('data/metabolism/processed/', s, '.csv'), guess_max = 10000) %>%
  select(DateTime_UTC, site, depth, avg_velocity)
  dd <- bind_rows(dd, d)
}

# load light modeled from StreamLight:
ll <- read_csv('data/daily_modeled_light_all_sites.csv')

# summarize data to daily and combine:
d <- dd %>%
  mutate(DateTime_EST = with_tz(DateTime_UTC, 'EST'),
         date = as.Date(DateTime_EST, tz = 'EST')) %>%
  group_by(site, date) %>%
  summarize(across(any_of(c('depth', 'avg_velocity', 'light')), mean, na.rm = T)) %>%
  ungroup() %>%
  full_join(ll, by = c('site', 'date'))

met <- met %>%
  filter(era == 'now') %>%
  rename(depth_hall = depth) %>%
  left_join( d, by = c('site', 'date')) %>%
  mutate(depth = case_when(is.na(depth) ~ depth_hall,
                           TRUE ~ depth)) %>%
  select(-depth_hall, -method, -era)

write_csv(met, 'data/metabolism/metabolism_and_drivers.csv')

