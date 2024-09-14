# Build a basic linear mixed effects model of met as f(temp, Q)
library(lme4)
library(tidyverse)
library(lubridate)
setwd('C:/Users/alice.carter/git/nhc_50yl/src')

dat <- readRDS("data/metabolism/compiled/met_preds_stream_metabolizer.rds")
sites <- read_csv('data/siteData/NHCsite_metadata.csv') %>%
  slice(c(1:5, 7))
met <- dat$preds

dd <- tibble()
for(s in sites$sitecode){
  d <- read_csv(paste0('data/metabolism/processed/', s, '.csv'), guess_max = 10000) %>%
  select(DateTime_UTC, site, depth, avg_velocity, light)
  dd <- bind_rows(dd, d)
}

d <- dd %>%
  mutate(DateTime_EST = with_tz(DateTime_UTC, 'EST'),
         date = as.Date(DateTime_EST, tz = 'EST')) %>%
  group_by(site, date) %>%
  summarize(across(any_of(c('depth', 'avg_velocity', 'light')), mean, na.rm = T))
met <- met %>%
  rename(depth_hall = depth) %>%
  left_join( d, by = c('site', 'date')) %>%
  mutate(depth = case_when(is.na(depth) ~ depth_hall,
                           TRUE ~ depth)) %>%
  select(-depth_hall)

# subset out modern dataset to build model
mm <- met %>%
  filter(era == 'now') %>%
  mutate(ER = -ER,
         logQ = log(discharge),
         year = factor(year))
ggplot(mm, aes(date, ER, col = year)) +
  geom_point()

ggplot(mm, aes(site, (K600)))+
  geom_boxplot() #+ facet_wrap(~year)

ggplot(mm, aes(temp.water, log(discharge), col = factor(month)))+
  geom_point() + geom_smooth(method = lm) +
  facet_wrap(~site)

# subset to fall respiration data and add in annual flow values

yy <- mm %>%
  group_by(year, site) %>%
  summarize(logQ_mean = mean(logQ, na.rm = T),
            depth_mean = mean(depth, na.rm = T)) %>%
  ungroup() %>%
  rename(sitecode = site) %>%
  left_join(sites, by = 'sitecode') %>%
  select(site = sitecode, year, logQ_mean, depth_mean, distance_m, slope)

dat <- left_join(mm, yy, by = c('site', 'year')) %>%
  select(-method, -era)
write_csv(dat, 'data/metabolism/metabolism_and_drivers.csv')

fall <- dat %>%
  filter(month %in% c(10)) %>%
  select(date, GPP, ER, logQ, temp.water, site, year, logQ_mean,
         depth_mean, distance_m, slope)

ggplot(fall, aes(temp.water, ER, col = factor(year)))+
  geom_point() + geom_smooth(method = lm) +
  facet_wrap(year~site, scales = 'free_x')
# scale data to model:

sfall <- fall %>%
  mutate(across(.cols = all_of(c('GPP', 'ER', 'temp.water')),
                .fns = scale))
# model for ER ####

mER = lmer(ER ~ temp.water + (temp.water|slope) + (temp.water|logQ_mean),
           data = sfall)
summary(mER)$coeff
confint(mER)
ranef(mER)

ER_mod <- predict(mER)
tmp <- data.frame(index = as.numeric(names(ER_mod)), ER_mod = ER_mod)
sfall <- sfall %>%
  mutate(index = seq(1:nrow(sfall))) %>%
  left_join(tmp) %>%
  select(-index)

ggplot(sfall, aes(ER, ER_mod, col = site))+
  geom_point() + geom_smooth(method = lm) +
  geom_abline(slope = 1, intercept = 0)

# It seems like there is actually no relationship between temperature and ER in  the
#  fall at the three run sites. The model is unable to predict at these sites (unsurprisingly)

# My next step would be to incorporate depth or slope (whichever best differentiates between
# the two sets of sites). I think this type of model will make this a much better paper and is
# worth doing.

# another possibility is to try something like a random forest or other machine
# learning approach, I would just have to see how much of an investment something
# like this is in order to do it from a place of knowing what's up.

# finally, maybe I do want to do a bayesian hierarchical model. It might actually work and
# be cool, it just means that I am still working on my thesis for ages.
