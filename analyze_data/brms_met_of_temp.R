# Basic hierarchical model of metabolism as a function of temperature and Q

# A carter
# 2/2022
library(tidyverse)
library(lme4)

setwd('C:/Users/Alice Carter/git/nhc_50yl/')

dat <- read_csv("data/metabolism/metabolism_and_drivers.csv")
sites <- read_csv('data/siteData/NHCsite_metadata.csv') %>%
  slice(c(1:5, 7))

# Test basic GPP model
P <- dat %>% 
  select(date, site, GPP, GPP.lower, GPP.upper, discharge, 
         temp.water, light, slope) %>%
  group_by(site) %>%
  mutate(GPP_pre = c(NA, GPP[1:(n()-1)]))

P <- filter(P, site == 'CBP') %>%
  mutate(discharge = log(discharge),
         across(.cols = any_of(c('light', 'discharge')), .fns = scale))
mod <- lm(GPP ~ GPP_pre + light + discharge, data = P)
summary(mod)

phi = 0.83
b_light = 0.0364
b_q = .06826
p = numeric()
p[1] = P$GPP[2]
for(i in 2:nrow(P)){
  p[i] = phi * p[i-1] + b_light * P$light[i] + b_q * P$discharge[i]
}
plot(P$date, p, ylim = c(-1,2))
lines(P$date, P$GPP)
plot(P$date, P$light)
m = T))

shapiro.test(dat$GPP)
library(ggpubr)
ggqqplot(dat$ER)
aes(temp.water, GPP, col = factor(month))) +
  geom_point() +
  facet_wrap(~site, scales = 'free_y') +

  geom_smooth(method = lm)

mod <- brms::brm(ER ~ (1 + temp.water| slope +  logQ_mean), data = dat, iter = 1000)
mod
