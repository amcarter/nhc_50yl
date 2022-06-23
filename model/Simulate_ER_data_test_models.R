# Script for simulating data and running different versions of ER models
# A Carter
# Feb 2022
setwd('C:/Users/alice.carter/git/nhc_50yl/')
library(tidyverse)
library(rstan)

# Load data ####
# Use real driver data from the CBP station

dat <- read_csv('data/metabolism/metabolism_and_drivers.csv')

met <- dat %>%
  arrange(site, date) %>%
  select(site, date, starts_with(c('ER', 'GPP')), temp.water,
         discharge, PAR_surface, depth) %>%
  group_by(site, date) %>%
  summarize(across(everything(), mean, na.rm = T)) %>%
  ungroup()

write_csv(met, 'data/met_and_drivers_for_model.csv')

dat <- read_csv('data/met_and_drivers_for_model.csv')


# define litterfall to be 20 days from Oct 5th - 24th with the total annual litter = 500 gC/m2/y, so daily litterfall = 25 gC/m2/d
cbp <- dat %>%
  filter(site == 'CBP') %>%
  select(-site) %>%
  arrange(date) %>%
  mutate(across(-date, zoo::na.approx, na.rm = F),
         litter = case_when(date >= as.Date('2019-10-05') &
                              date <= as.Date('2019-10-24') ~ 25,
                            TRUE ~ 0),
         logQ = log(discharge),
         tempK = temp.water + 273) %>%
  filter(!is.na(ER))

cbp %>%
  select(date, temp = temp.water, logQ, light = PAR_surface, litter) %>%
  pivot_longer(-date, names_to = 'variable', values_to = 'value') %>%
  ggplot(aes(date, value, col = variable)) +
  geom_line() +
  facet_wrap(.~variable, nrow = 4, scales = 'free_y')

# calculate observation error:
mean(cbp$GPP.upper- cbp$GPP.lower)/4
sd((cbp$GPP.upper- cbp$GPP.lower)/4)

# Simulate data
N <- nrow(cbp)
light = cbp$PAR_surface/max(cbp$PAR_surface)
Q = cbp$discharge/max(cbp$discharge)
temp = cbp$temp.water + 273
litter = cbp$litter

# simulate GPP
P <- numeric()
P[1] = 0.5
for(i in 2:N){
  P[i] = 0.5 * P[i-1] + 0.8 * light[i] + rnorm(1,0, 0.05)
}

# Load Hall driver data for testing hindcast (with simulations)
**** Load Hall data


# Basic ER model without latent state ####
# ER_base_model_1.stan

# define parameters
b1 = .5 * 10e8 # I can only get a measurable amount of heterotrophic respiration with an extremely high b1
b2 = .8
E = 0.63
sigma_proc = 0.1
sigma_obs = 0.05 # adjust this based on the actual measurement error in the data

# Constants
k_b = 8.6173 * 10^-5
AR = 0.5

# Process model
C = numeric()
R = numeric()
R_obs = numeric()
C[1] = 100 # initial accessible carbon storage in stream
R[1] = b1 * C[1] * exp(-E/(k_b * temp[1])) + AR * P[1]

for(i in 2:N){
  C[i] = (C[i-1] + litter[i] - R[i-1] + P[i-1])*(1-b2*Q[i]) + rnorm(1, 0, sigma_proc)
  if(C[i]<0) C[i] = 0
  R[i] = b1 * C[i] * exp(-E/(k_b * temp[i])) + AR * P[i]
}

R_obs = rnorm(N, R, sigma_obs)

data.frame(date = cbp$date, P = P, R = R_obs, C = C) %>%
  pivot_longer(-date, names_to = 'var', values_to = 'value') %>%
  ggplot(aes(date, value, col = var)) +
  geom_line()+
  facet_wrap(.~var, scales = 'free_y', nrow = 3)

# Run Stan model on simulated data:

stan_dat <- list(N = length(R_obs), R_obs = R_obs, P_obs = P,
                 tempK = temp, Q = Q, litter = litter, C0 = C[1])
stan_dat <- list(N = length(R_obs), R_obs = R_obs, P_obs = P,
                 tempK = temp, C = C)

mod <- stan('src/model/stan/ER_base_model_1.stan',
            data = stan_dat,
            chains = 4, init = 0,
            warmup = 500, iter = 1000)


saveRDS(mod, 'src/model/model_runs/simulated_ER_base_mod_1.rds')

# Basic ER model with relative temperature correction ####
# ER_base_model_2.stan

# define parameters
R20 = .01  # Heterotrophic respiration at 20 C
beta = 0.8   # Percent removal of organic carbon from max storm
sigma_proc = 0.1
sigma_obs = 0.08

# Constants
E = 0.63              # activation energy for heterotrophic respiration
k_b = 8.6173 * 10^-5  # Boltzmann's constant in eV/K
AR = 0.5              # the fraction of GPP respired by autotrophs


# Process model
C = numeric()
R = numeric()
R_obs = numeric()
C[1] = 100 # initial accessible carbon storage in stream
R[1] = AR * P[1] + (1 - exp(-R20 * exp(-E/k_b * (1/temp[1] - 1/293)))) * C[1]

for(i in 2:N){
  C[i] = (C[i-1] + litter[i] - R[i-1] + P[i-1])*(1-beta*Q[i]) +
    rnorm(1, 0, sigma_proc)
  if(C[i]<5) C[i] = 5
  R[i] = AR * P[i] + (1 - exp(-R20 * exp(-E/k_b * (1/temp[i] - 1/293))))*C[i]
}

R_obs = rnorm(N, R, sigma_obs)

data.frame(date = cbp$date, P = P, R = R_obs, C = C) %>%
  pivot_longer(-date, names_to = 'var', values_to = 'value') %>%
  ggplot(aes(date, value, col = var)) +
  geom_line()+
  facet_wrap(.~var, scales = 'free_y', nrow = 3)

# Run Stan model on simulated data:

stan_dat <- list(N = length(R_obs), R_obs = R_obs, P_obs = P,
                 tempK = temp, C = C)

mod <- stan('src/model/stan/ER_base_model_2.stan',
            data = stan_dat, init = 0,
            chains = 4, cores = 4,
            warmup = 200, iter = 500)

saveRDS(mod, 'src/model/model_runs/simulated_ER_base_mod_2.rds')

pairs(mod, pars = c('R20', 'sigma_obs'))
traceplot(mod, pars = c('R20', 'sigma_obs'), ncol = 1)

launch_shinystan(mod)
library(shinystan)
ponmlkjihgfedcba:bcdefghijklmnopqrstsrqponmlklmnmlkjijklmnopqrstuvwxyz$yxwvutsrqponmlkjihgfefghijklmnopqrstuvwxyz(yxwvutsrqpopqrstuvwxyz.yxwvutsrqponmlkjijklmnopqrstuvwxyz=zyxwvutsrqponmlkjihgfedcba'bcdefghijklmnopqrstuvwxyz:yxwvutsrqponmlkjihgfedcba/bcdefghijklmnopqrstUtsrstuvwxyz/yxwvutsrqponmlkjihgfedCbabcdefghijklmnopqrqponmlkjihgfefghijklmnopqrstuvwxyz/yxwvutsrqponmlkjihghijklmnopqrstuvwxyz/yxwvutsrqponmlkhihgfedcba1234543210abcdefghijklmnopqrstuvwxyz'yxwvutsrqponmlkjihgfedcba)bcdefghijklmnopqrstuvwxyz`yxwvutsrqponmlkjihgfefghijklmnopqrstuvutsrqponmlkjihgfefghijklmnopqrstuvwxyz:yxWvutsrqponmlkjihgfeedcbcdefghijklkjihgfedefgfedcbabcdefghijklmnonmlkjihgfghijklmnopqrsstsrqponmnopqrqponmlkjihgfefghijklmnopqrstsrqponmlkjijklmnmlkjihgfedcbabcdefghijklmnopqrsrqponmlkjihgfedefghijklmnonmlkjihihgfedcdefghijklkjihgfedefghijklmnopqrqponmlkjihgfefghijklmnopqrsrqponmlkjihghijklmnopqrstutsrqponmlkjihgfedefghijklmnopqrstsrqponmlkjihgfghijklmnopqrqponmlkjihgfedcdefghijklmnopqrstutsrqponmlkjihgfeprint(mod)plot(mod)readRDS('src/model/model_runs/simulated_ss_ER.rds')
