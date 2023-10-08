
library(tidyverse)
library(lubridate)

ysi <- read_csv("data/water_chemistry/all_nhc_ysi_data.csv") %>%
    select(site, DateTime_UTC, DO_mgL, watertemp_C)

sensor <-  read_csv('data/site_data/metabolism.csv') %>%
    mutate(date = as.Date(DateTime_UTC, tz = 'UTC'))

dat <- sensor %>%
    select(DateTime_UTC, site, DO.obs, temp.water) %>%
    left_join(ysi, by = c('DateTime_UTC', 'site')) %>%
    filter(!is.na(DO_mgL)) %>%
    mutate(DO_error = abs(DO.obs - DO_mgL)/DO_mgL,
           temp_error = (temp.water - watertemp_C)/watertemp_C)
summary(dat)

plot(density(dat$DO_error, na.rm = T))
plot(density(dat$temp_error, na.rm = T))

