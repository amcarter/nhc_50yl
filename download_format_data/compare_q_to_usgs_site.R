# compare median flows at the sites along NHC to watershed areas:
library(tidyverse)
library(dataRetrieval)
library(lubridate)

setwd('C:/Users/Alice Carter/git/nhc_50yl/')
sites <- read_csv('data/siteData/NHCsite_metadata.csv')


# Compare Discharge across sites ####
# read in NHC and UNHC discharge data:

q <- read_csv('data/rating_curves/NHC_UNHC_Q.csv', guess_max = 10000) %>%
  mutate(discharge_nhc = case_when(notes == 'nhc modeled' ~ NA_real_,
                                   TRUE ~ discharge_nhc),
         discharge_unhc = case_when(notes == 'unhc modeled' ~ NA_real_,
                                   TRUE ~ discharge_unhc)) %>%
  rename(stage_nhc = NHC, stage_unhc = UNHC)

q <- q %>%
  pivot_longer(cols =  starts_with(c('discharge', 'stage')),
               names_pattern = '^(discharge|stage)_(.*)$',
               names_to = c('.value','site')) %>% 
  mutate(site = toupper(site))
                                  

# get discharge data for NHC site at blands from NWIS
site_no <- substr(sites$sitecode[12], 6, 13)

# download raw daily data (discharge and gage height)
dat <- readNWISuv(siteNumbers = site_no, parameterCd = c('00060', '00065'),
                  startDate = '2017-03-01', endDate = '2020-04-01')

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
qq <- bind_rows(q, dat) %>%
  arrange(DateTime_UTC) %>%
  select(-notes) %>%
  filter(DateTime_UTC >= ymd_hms('2017-03-01 00:00:00'),
         DateTime_UTC <= ymd_hms('2020-03-21 00:00:00'))

write_csv(qq, 'data/rating_curves/discharge_NHC_UNHC_blands.csv')
# plot the data
qq %>% filter(DateTime_UTC >= ymd_hms('2020-01-05 00:00:00')) %>%
  ggplot(aes(DateTime_UTC, log(discharge), col = site)) +
 geom_line()

meds <- qq %>%
  mutate(logq = log(discharge))%>%
  group_by(site) %>%
  summarize(medq = median(discharge, na.rm = T),
            med_logq = median(logq, na.rm = T)) %>%
  ungroup() %>%
  rename(sitecode = site) %>%
  left_join(sites, by = 'sitecode' )%>%
  select(sitecode, medq, med_logq, ws_area.km2)

# examine how median discharge varies with watershed area:
AQcurve <- data.frame(x = seq(-5,5, by = .1)) 
AQcurve$y = exp(AQcurve$x)*197/1.03
library(ggridges)
png('figures/Q_by_ws_area_NHC.png')
qq %>%
  rename(sitecode = site) %>%
  left_join(sites[,c(2,12)]) %>%
  ggplot(aes(x = log(discharge), y = ws_area.km2,
             col = sitecode)) +
    geom_density_ridges(aes(fill = sitecode), alpha = .3) +
    theme_ridges() + 
    # geom_density(alpha = .3) + 
    geom_point(data = meds, aes(x = med_logq, y = ws_area.km2 ))+
    geom_line(data = AQcurve, aes(x, y), col = 'black')+
    xlim(-5,5)+ ylim(0, 300)+
    ylab('watershed area') + xlab('log(discharge)')
dev.off()

yy <- qq %>% 
  mutate(year = year(DateTime_UTC),
         discharge = case_when(site == 'NHC' & discharge > 6 ~ NA_real_,
                               site == 'UNHC' & discharge > 2.6 ~ NA_real_,
                               TRUE ~ discharge)) %>%
  group_by(site, year) %>% 
  mutate(discharge = zoo::na.approx(discharge, na.rm = F)) %>% 
  summarize(discharge = sum(discharge, na.rm = T)) %>%
  rename(sitecode = site) %>% left_join(meds[,c(1,4)])

ggplot(yy, aes(x = year, y = discharge, fill = sitecode)) +
  geom_bar(stat = 'identity', position = 'dodge')

ggplot(yy, aes(x = ws_area.km2, y = discharge, col = factor(year))) + 
  geom_line()
