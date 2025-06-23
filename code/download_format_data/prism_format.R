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
