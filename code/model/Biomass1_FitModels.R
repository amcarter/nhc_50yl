## Fitting models to data - JR Blaszczak
## Modified: June 23, 2021
## Heili Lowman

setwd('C:/Users/alice.carter/git/nhc_50yl/src')
# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","here", "viridis"), require, character.only=T)

## Source data
df <- read_csv('data/metabolism/metabolism_and_drivers.csv') %>%
  rename(Q = discharge, light = PAR_surface) %>%
  filter(site == 'CBP') %>%
  data.frame()

####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
## specify number of cores
options(mc.cores=6)

## compile data
stan_data_compile <- function(x){
  data <- list(Ndays=length(x$GPP),
               light = x$light/max(x$light, na.rm = T),
               GPP = x$GPP,
               GPP_sd = (x$GPP.upper - x$GPP.lower)/3.92,
               tQ = x$Q/max(x$Q, na.rm = T))
  return(data)
}

stan_data_l <- stan_data_compile(df)

#########################################
## Run Stan to get parameter estimates - all sites
#########################################

## PM 2 - Latent Biomass (Ricker)
# With Persistence Term (P)
init_Ricker <- function(...) {
  list(c = 0.5, s = 0.5)
}

## export results
PM_outputlist_Ricker <- stan("src/model/Stan_ProductivityModel2_Ricker_fixedinit_obserr.stan",
                            data = stan_data_l,
                            chains = 4,iter = 5000,
                            init = init_Ricker,
                            control = list(max_treedepth = 12))

saveRDS(PM_outputlist_Ricker, "stan_1riv_output_Ricker_2021_06_23.rds")

# End of script.
