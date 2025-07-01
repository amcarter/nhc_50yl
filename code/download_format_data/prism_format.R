library(tidyverse)
library(lubridate)

read_csv('~/git/papers/alice_nhc/data/prism_raw/PRISM_PPT_TMEAN_Median_monthly.csv') %>%
  select(date, ppt, tmean) %>%
  group_by(year = year(date)) %>%
  summarize(
    cumulative_precip = sum(ppt, na.rm = TRUE),
    temp_mean = mean(tmean, na.rm = TRUE)
  ) %>%
  write_csv('~/git/papers/alice_nhc/data/watershed/prism_annual.csv')

read_csv('~/git/papers/alice_nhc/data/prism_raw/NLDAS_Hourly_Temp_Precip.csv') %>%
  select(date, precip_mmd = total_precipitation) %>%
  group_by(date) %>%
  summarize(
    precip_mmd = sum(precip_mmd, na.rm = TRUE)
    # temp_mean = mean(tmean, na.rm = TRUE)
  ) %>%
  write_csv('~/git/papers/alice_nhc/data/watershed/nldas2.csv')
