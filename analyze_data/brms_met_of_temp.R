# build brms models of metabolism as a function of temperature

library(tidyverse)
library(brms)

setwd('C:/Users/Alice Carter/git/nhc_50yl/')

dat <- read_csv("data/metabolism/metabolism_and_drivers.csv")
sites <- read_csv('data/siteData/NHCsite_metadata.csv') %>%
  slice(c(1:5, 7))


ggplot(dat, aes(temp.water, ER, col = month)) +
  geom_point() +
  facet_wrap(~site, scales = 'free_y') 
  # geom_smooth(method = lm)

mod <- brms::brm(ER ~ (1 + temp.water| slope +  logQ_mean), data = dat, iter = 1000)
mod
