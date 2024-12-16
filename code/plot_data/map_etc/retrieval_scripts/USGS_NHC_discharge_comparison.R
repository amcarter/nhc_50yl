library(tidyverse)
library(dataRetrieval)
library(lubridate)


sites <- read_csv('data/siteData/NHCsite_metadata.csv')

site_no <- substr(sites$sitecode[12], 6, 13)

q <- read_csv('data/rating_curves/NHC_UNHC_Q.csv', guess_max = 10000) %>%
    mutate(discharge_nhc = case_when(notes == 'nhc modeled' ~ NA_real_,
                                     TRUE ~ discharge_nhc),
           discharge_unhc = case_when(notes == 'unhc modeled' ~ NA_real_,
                                      TRUE ~ discharge_unhc)) %>%
    rename(stage_nhc = NHC, stage_unhc = UNHC)

q <- q %>%
    select(DateTime_UTC, discharge = discharge_nhc)
q$site = 'NHC'


# download raw daily data (discharge and gage height)
dat <- readNWISuv(siteNumbers = site_no, parameterCd = c('00060', '00065'),
                  startDate = '1985-10-01', endDate = '2020-01-01')


dat <- dat %>%
    as_tibble() %>%
    rename(DateTime_UTC = dateTime) %>%
    select(-ends_with('cd')) %>%
    pivot_longer(starts_with('X')) %>%
    mutate(site = paste('USGS', site_no, sep = '-'),
           parm_cd = substr(name, 3, 7),
           variable = case_when(parm_cd == '00060' ~ 'discharge_cfs',
                                parm_cd == '00065' ~ 'stage_f')) %>%
    select(-name, -parm_cd) %>%
    pivot_wider(names_from = 'variable', values_from = 'value') %>%
    mutate(discharge = discharge_cfs / (3.28^3),
           stage = stage_f / 3.28) %>%
    select(DateTime_UTC, site, discharge, stage)

dat <- dat%>%
    mutate(discharge = zoo::na.approx(discharge, na.rm = F),
           stage = zoo::na.approx(stage, na.rm = F))

dat$site = 'BLANDS'
q <- bind_rows(q, dat) %>%
    select(-stage) %>%
    pivot_wider(names_from = 'site',
                values_from = 'discharge') %>%
    arrange(DateTime_UTC) %>%
    mutate(NHC = zoo::na.approx(NHC, x = DateTime_UTC, na.rm = F),
           BLANDS = zoo::na.approx(BLANDS, x = DateTime_UTC, na.rm = F))

ggplot(q, aes(log(BLANDS), log(NHC))) + geom_point()

mod <- lm(log(NHC) ~ log(BLANDS), q)

q$NHC_mod <- exp(-1.2146 + 0.9847 * log(q$BLANDS))



RBIcalc <- function(x){
    sum(x)/length(x)
}

qq <- q %>%
    mutate(date = as.Date(DateTime_UTC)) %>%
    group_by(date) %>%
    summarize(NHC = mean(NHC_mod, na.rm = T),
              BLANDS = mean(BLANDS, na.rm = T)) %>%
    ungroup() %>%
    mutate(year = lubridate::year(date)) %>%
    pivot_longer(cols = c('NHC', 'BLANDS'),
                 names_to = 'site',
                 values_to = 'discharge_m3s')

fall_mean <- qq %>%
    mutate(month = month(date)) %>%
    filter(month %in% c(9, 10, 11)) %>%
    group_by(year, site) %>%
    summarize(fall_mean = mean(discharge_m3s, na.rm = T))

fall_mean %>%
    filter(site == 'NHC') %>%
    write_csv('data/rating_curves/modeled_fall_mean_NHC_flow.csv')

Q_stats <- qq %>%
    filter(site == 'BLANDS') %>%
    group_by(year, site) %>%
    summarize(n = length(discharge_m3s),
              RBI = RBIcalc(discharge_m3s),
              peak_Q = max(discharge_m3s, na.rm = T),
              ar_1 = arima(discharge_m3s, order = c(1,0,0))$coef[1],
              q05 = quantile(discharge_m3s, .05, na.rm = T)) %>%
    ungroup()  %>%
    filter(n >= 250) %>%
    left_join(fall_mean)# Don't keep incomplete years

png(width=7, height=6, units='in', type='cairo', res=300,
    filename='figures/SI/USGS_rbi_trends.png')

Q_stats %>%
    select(-n) %>%
    pivot_longer(cols = c('RBI', 'peak_Q', 'ar_1', 'q05', 'fall_mean'),
                 names_to = 'var', values_to = 'val') %>%
    ggplot(aes(year, val)) +
    geom_point() + geom_line() +
    facet_wrap(.~var, scales = 'free_y')

dev.off()
