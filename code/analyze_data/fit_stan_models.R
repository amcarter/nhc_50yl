# Basic hierarchical model of metabolism as a function of temperature, light, and Q

# A carter
# 2/2022
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# load in the modern dataset:
dat <- read_csv("data/metabolism/compiled/metabolism_and_drivers.csv")
sites <- read_csv('data/siteData/NHCsite_metadata.csv') %>%
    slice(c(1:5, 7))

# load in Hall dataset:

hall <- readRDS("data/metabolism/compiled/met_preds_stream_metabolizer_O2.rds")
hall <- hall$preds %>% filter(era == 'then', site == 'CBP')
hall_QT <- read_csv('data/hall/hall_discharge_temp_daily.csv')


# Test basic GPP model
P <- dat %>%
    select(date, site,ER,  GPP, GPP.lower, GPP.upper, discharge,
           temp.water, light = PAR_surface)%>%#, slope) %>%
    group_by(site) %>%
    left_join(select(sites, site = sitecode, slope = slope_nhd)) %>%
    filter(site %in% c('CBP', 'NHC')) %>%
    mutate(log_Q = log(discharge))

scaling_pars <- P %>%
    ungroup() %>%
    summarize(across(c('log_Q', 'light', 'temp.water'),
                     .fns = c(mean = \(x) mean(x, na.rm = TRUE),
                              sd = \(x) sd(x, na.rm = TRUE))))


min_Q = min(P$discharge[P$site == 'CBP'], na.rm = T)
epsilon = 1e-4

P_scaled <- P %>%
    mutate(across(c('log_Q', 'light', 'temp.water'),
                  .fns = \(x) as.vector(scale(x))),
           sin_time = sin(2*pi *as.numeric(format(date, '%j'))/365),
           cos_time = cos(2*pi *as.numeric(format(date, '%j'))/365),
           month = month(date),
           season = case_when(month %in% c(12, 1, 2) ~ 'Winter',
                              month %in% c(3, 4, 5) ~ 'Spring',
                              month %in% c(6, 7, 8) ~ 'Summer',
                              month %in% c(9, 10, 11) ~ 'Fall'),
           GPP = GPP + epsilon,
           ER = -ER + epsilon,
           GPP_pred = log(GPP),
           site = factor(site, levels = c("CBP", "NHC")))

# add light data to hall:
light <- P %>%
    mutate(doy = format(date, '%j')) %>%
    dplyr::select(doy, light) %>%
    group_by(doy) %>%
    summarize(light = mean(light, na.rm = T))

# add light to Hall daatset and scale
hall_preds <- hall_QT %>%
    mutate(doy = format(date, '%j'),
           # don't allow discharge to be less than 50% of the minimum we observed at the same site
           discharge_m3s = case_when(discharge_m3s < 0.5*min_Q ~ 0.5 * min_Q,
                                     TRUE ~ discharge_m3s)) %>%
    left_join(light, by = 'doy') %>%
    mutate(site = 'CBP', log_Q = log(discharge_m3s)) %>%
    dplyr::select(date, site,  temp.water = water_temp_C, log_Q, light) %>%
    mutate(across(c('log_Q', 'light', 'temp.water'),
                  .fns = \(x) zoo::na.approx(x, na.rm = F)),
           sin_time = sin(2*pi *as.numeric(format(date, '%j'))/365),
           cos_time = cos(2*pi *as.numeric(format(date, '%j'))/365),
           month = month(date),
           season = case_when(month %in% c(12, 1, 2) ~ 'Winter',
                              month %in% c(3, 4, 5) ~ 'Spring',
                              month %in% c(6, 7, 8) ~ 'Summer',
                              month %in% c(9, 10, 11) ~ 'Fall'),
           site = factor(site, levels = c("CBP", "NHC"))) %>%
    slice_head(n = -4)

hall_scaled <- hall_preds %>%
    mutate(log_Q = (log_Q - scaling_pars$log_Q_mean)/scaling_pars$log_Q_sd,
           light = (light - scaling_pars$light_mean)/scaling_pars$light_sd,
           temp.water = (temp.water - scaling_pars$temp.water_mean)/scaling_pars$temp.water_sd)




# remove leading and ending NA's in GPP from each site.
NHC <- P_scaled %>% ungroup() %>%
    slice(-c(1, 362:387, 1506)) %>%
    mutate(across(c('temp.water', 'light', 'log_Q'),
                  \(x) zoo::na.approx(x, na.rm = FALSE)),
           GPP_pred = zoo::na.approx(GPP, na.rm = FALSE))
cbp <- NHC %>% filter(site == 'CBP')

P_obs <- NHC %>% filter(!is.na(GPP))

dat_list <- list(
    N_site = length(unique(NHC$site)),
    N_obs = length(which(!is.na(NHC$GPP))),
    N_mis = length(which(is.na(NHC$GPP))),
    ii_obs = which(!is.na(NHC$GPP)),
    ii_mis = which(is.na(NHC$GPP)),
    P_obs = log(P_obs$GPP),
    light = NHC$light,
    temperature = NHC$temp.water,
    site = as.numeric(NHC$site)
)

mod_stan <- stan("code/analyze_data/ar1_GPP.stan",
                 data = dat_list,
                 chains = 4, iter = 2000)

# summary(mod_stan)
fit <- mod_stan
# Evaluate model and hindcast:
    # Extract posterior estimates
    beta0_post <- rstan::extract(fit, pars = "beta_0")$beta_0
    beta_post <- rstan::extract(fit, pars = "beta")$beta
    phi_post <- rstan::extract(fit, pars = "phi")$phi
    sigma_post <- rstan::extract(fit, pars = "sigma")$sigma

    beta_hat <- data.frame(
        var = c('beta_0', 'beta_light', 'beta_temp'),
        mean = apply(beta_post, 2, mean),
        median = apply(beta_post, 2, median),
        min = apply(beta_post, 2, min),
        max = apply(beta_post, 2, max),
        low = apply(beta_post, 2, quantile, probs = 0.025),
        high = apply(beta_post, 2, quantile, probs = 0.975),
        q0.01 = apply(beta_post, 2, quantile, probs = 0.01),
        q0.05 = apply(beta_post, 2, quantile, probs = 0.05),
        q0.1 = apply(beta_post, 2, quantile, probs = 0.1),
        q0.9 = apply(beta_post, 2, quantile, probs = 0.9),
        q0.95 = apply(beta_post, 2, quantile, probs = 0.95),
        q0.99 = apply(beta_post, 2, quantile, probs = 0.99)
    )

    phi_hat <- data.frame(
        mean = mean(phi_post),
        median = median(phi_post),
        min = min(phi_post),
        max = max(phi_post),
        low = quantile(phi_post, probs = 0.025),
        high = quantile(phi_post, probs = 0.975),
        q0.01 = quantile(phi_post, probs = 0.01),
        q0.05 = quantile(phi_post, probs = 0.05),
        q0.1 = quantile(phi_post, probs = 0.1),
        q0.9 = quantile(phi_post, probs = 0.9),
        q0.95 = quantile(phi_post, probs = 0.95),
        q0.99 = quantile(phi_post, probs = 0.99)
    )

    sigma_hat <- data.frame(
        mean = mean(sigma_post),
        median = median(sigma_post),
        min = min(sigma_post),
        max = max(sigma_post),
        low = quantile(sigma_post, probs = 0.025),
        high = quantile(sigma_post, probs = 0.975),
        q0.01 = quantile(sigma_post, probs = 0.01),
        q0.05 = quantile(sigma_post, probs = 0.05),
        q0.1 = quantile(sigma_post, probs = 0.1),
        q0.9 = quantile(sigma_post, probs = 0.9),
        q0.95 = quantile(sigma_post, probs = 0.95),
        q0.99 = quantile(sigma_post, probs = 0.99)
    )

    par_ests <- list(beta_hat = beta_hat,
                     phi_hat = phi_hat,
                     sigma_hat = sigma_hat)

    # hindcast the historical observations

    draws <- nrow(beta_post)

    # matrix of draws from the posterior-predictive distribution
    post_preds <- matrix(nrow = draws, ncol = nrow(hall_scaled))
    # post_preds2 <- matrix(nrow = draws, ncol = nrow(cbp))

    # fill in first observation based on the mean from the measured historical values at the site
    post_preds[, 1] <- matrix(
        rep(log(mean(hall$GPP, na.rm = T) + epsilon), each = draws),
        nrow = draws, ncol = 1
    )
    # post_preds2[, 1] <- matrix(
    #     rep(log(mean(cbp$GPP, na.rm = T)), each = draws),
    #     nrow = draws, ncol = 1
    # )


    for(i in 1:draws){
        for(t in 2:nrow(hall_scaled)){
            post_preds[i, t] <- beta0_post[i,1] + beta_post[i,2] * hall_scaled$light[t] +
                beta_post[i,3] * hall_scaled$temp.water[t] +
                phi_post[i] * post_preds[i, t-1] #+ rnorm(1, sd = sigma_post[i])
        }
    }

    # for(i in 1:draws){
    #     for(t in 2:nrow(cbp)){
    #         post_preds2[i, t] <- beta0_post[i,1] + beta_post[i,2] * cbp$light[t] +
    #             beta_post[i,3] * cbp$temp.water[t] +
    #             phi_post[i] * post_preds2[i, t-1] #+ rnorm(1, sd = sigma_post[i])
    #     }
    # }

    hindcast <- data.frame(
        date = hall_scaled$date,
        GPP_pred = apply(post_preds, 2, mean),
        GPP_median = apply(post_preds, 2, median),
        GPP_low = apply(post_preds, 2, quantile, probs = 0.025),
        GPP_high = apply(post_preds, 2, quantile, probs = 0.975)
    ) %>%
        mutate(GPP = exp(GPP_pred) - epsilon,
               GPP_low = exp(GPP_low)- epsilon,
               GPP_high = exp(GPP_high - epsilon),
               doy = format(date, '%j')) #%>%
        group_by(doy) %>%
        summarize(across(starts_with('GPP'), \(x) mean(x, na.rm = T))) %>%
        mutate(date = as.Date(paste0('2020-', doy), format = '%Y-%j'))

    todaycast <- data.frame(
        date = cbp$date,
        GPP_pred = apply(post_preds2, 2, mean),
        GPP_median = apply(post_preds2, 2, median),
        GPP_low = apply(post_preds2, 2, quantile, probs = 0.025),
        GPP_high = apply(post_preds2, 2, quantile, probs = 0.975)
    ) %>%
        mutate(GPP = exp(GPP_pred) - epsilon,
               GPP_low = exp(GPP_low)- epsilon,
               GPP_high = exp(GPP_high - epsilon),
               doy = format(date, '%j')) %>%
        group_by(doy) %>%
        summarize(across(starts_with('GPP'), \(x) mean(x, na.rm = T))) %>%
        mutate(date = as.Date(paste0('2020-', doy), format = '%Y-%j'))


    ggplot(hindcast, aes(date, GPP), col = 2) +
        geom_ribbon(aes(ymin = GPP_low, ymax = GPP_high),
                    col = NA, fill = 'grey') +
        geom_point() + geom_line() +
        # geom_point(data = todaycast, col = 1)
        geom_point(data = hall,  col = 2)

    hall_scaled$GPP_pred <- hindcast$GPP_pred


# model for ER:

R_obs <- NHC %>% filter(!is.na(ER))

dat_list <- list(
    N_site = length(unique(P_scaled$site)),
    N_obs = length(which(!is.na(NHC$ER))),
    N_mis = length(which(is.na(NHC$ER))),
    ii_obs = which(!is.na(NHC$ER)),
    ii_mis = which(is.na(NHC$ER)),
    R_obs = log(R_obs$ER),
    log_GPP = log(NHC$GPP_pred),
    temperature = NHC$temp.water,
    log_Q = NHC$log_Q,
    site = as.numeric(NHC$site)
)

mod_stan_ER <- stan("code/analyze_data/ar1_ER.stan",
                 data = dat_list,
                 chains = 4, iter = 2000)

# summary(mod_stan)
fit <- mod_stan_ER
# Evaluate model and hindcast:
    # Extract posterior estimates
    beta0_post <- rstan::extract(fit, pars = "beta_0")$beta_0
    beta_post <- rstan::extract(fit, pars = "beta")$beta
    phi_post <- rstan::extract(fit, pars = "phi")$phi
    sigma_post <- rstan::extract(fit, pars = "sigma")$sigma

    beta_hat <- data.frame(
        var = c('beta_0', 'beta_light', 'beta_temp', 'beta_Q', 'beta_tempQ'),
        mean = apply(beta_post, 2, mean),
        median = apply(beta_post, 2, median),
        min = apply(beta_post, 2, min),
        max = apply(beta_post, 2, max),
        low = apply(beta_post, 2, quantile, probs = 0.025),
        high = apply(beta_post, 2, quantile, probs = 0.975),
        q0.01 = apply(beta_post, 2, quantile, probs = 0.01),
        q0.05 = apply(beta_post, 2, quantile, probs = 0.05),
        q0.1 = apply(beta_post, 2, quantile, probs = 0.1),
        q0.9 = apply(beta_post, 2, quantile, probs = 0.9),
        q0.95 = apply(beta_post, 2, quantile, probs = 0.95),
        q0.99 = apply(beta_post, 2, quantile, probs = 0.99)
    )

    phi_hat <- data.frame(
        mean = mean(phi_post),
        median = median(phi_post),
        min = min(phi_post),
        max = max(phi_post),
        low = quantile(phi_post, probs = 0.025),
        high = quantile(phi_post, probs = 0.975),
        q0.01 = quantile(phi_post, probs = 0.01),
        q0.05 = quantile(phi_post, probs = 0.05),
        q0.1 = quantile(phi_post, probs = 0.1),
        q0.9 = quantile(phi_post, probs = 0.9),
        q0.95 = quantile(phi_post, probs = 0.95),
        q0.99 = quantile(phi_post, probs = 0.99)
    )

    sigma_hat <- data.frame(
        mean = mean(sigma_post),
        median = median(sigma_post),
        min = min(sigma_post),
        max = max(sigma_post),
        low = quantile(sigma_post, probs = 0.025),
        high = quantile(sigma_post, probs = 0.975),
        q0.01 = quantile(sigma_post, probs = 0.01),
        q0.05 = quantile(sigma_post, probs = 0.05),
        q0.1 = quantile(sigma_post, probs = 0.1),
        q0.9 = quantile(sigma_post, probs = 0.9),
        q0.95 = quantile(sigma_post, probs = 0.95),
        q0.99 = quantile(sigma_post, probs = 0.99)
    )

    par_ests <- list(beta_hat = beta_hat,
                     phi_hat = phi_hat,
                     sigma_hat = sigma_hat)

    # hindcast the historical observations

    draws <- nrow(beta_post)

    # matrix of draws from the posterior-predictive distribution
    post_preds <- matrix(nrow = draws, ncol = nrow(hall_scaled))
    # post_preds2 <- matrix(nrow = draws, ncol = nrow(cbp))

    # fill in first observation based on the mean from the measured historical values at the site
    post_preds[, 1] <- matrix(
        rep(log(-mean(hall$ER, na.rm = T) + epsilon), each = draws),
        nrow = draws, ncol = 1
    )
    # post_preds2[, 1] <- matrix(
    #     rep(log(mean(cbp$ER, na.rm = T)), each = draws),
    #     nrow = draws, ncol = 1
    # )


    for(i in 1:draws){
        for(t in 2:nrow(hall_scaled)){
            post_preds[i, t] <- beta0_post[i,1] + beta_post[i,2] * hall_scaled$GPP_pred[t] +
                beta_post[i,3] * hall_scaled$temp.water[t] + beta_post[i,4]*hall_scaled$log_Q[t] +
                beta_post[i,5] * hall_scaled$temp.water[t] * hall_scaled$log_Q[t] +
                phi_post[i] * post_preds[i, t-1] #+ rnorm(1, sd = sigma_post[i])
        }
    }
    # for(i in 1:draws){
    #     for(t in 2:nrow(cbp)){
    #         post_preds2[i, t] <- beta0_post[i,1] + beta_post[i,2] * cbp$light[t] +
    #             beta_post[i,3] * cbp$temp.water[t] + beta_post[i,4]*cbp$log_Q[t] +
    #             phi_post[i] * post_preds2[i, t-1] #+ rnorm(1, sd = sigma_post[i])
    #     }
    # }

    hindcast <- data.frame(
        date = hall_scaled$date,
        ER_pred = apply(post_preds, 2, mean),
        ER_median = apply(post_preds, 2, median),
        ER_low = apply(post_preds, 2, quantile, probs = 0.025),
        ER_high = apply(post_preds, 2, quantile, probs = 0.975)
    ) %>%
        mutate(ER = -exp(ER_pred) + epsilon,
               ER_low = -exp(ER_low) + epsilon,
               ER_high = -exp(ER_high + epsilon),
               doy = format(date, '%j'))# %>%
#         group_by(doy) %>%
#         summarize(across(starts_with('ER'), \(x) mean(x, na.rm = T))) %>%
#         mutate(date = as.Date(paste0('2020-', doy), format = '%Y-%j'))
#
# todaycast <- data.frame(
#         date = cbp$date,
#         ER_pred = apply(post_preds2, 2, mean),
#         ER_median = apply(post_preds2, 2, median),
#         ER_low = apply(post_preds2, 2, quantile, probs = 0.025),
#         ER_high = apply(post_preds2, 2, quantile, probs = 0.975)
#     ) %>%
#         mutate(ER = -exp(ER_pred) + epsilon,
#                ER_low = -exp(ER_low) + epsilon,
#                ER_high = -exp(ER_high + epsilon),
#                doy = format(date, '%j')) %>%
#         group_by(doy) %>%
#         summarize(across(starts_with('ER'), \(x) mean(x, na.rm = T))) %>%
#         mutate(date = as.Date(paste0('2020-', doy), format = '%Y-%j'))


    ggplot(hindcast, aes(date, ER)) +
        geom_ribbon(aes(ymin = ER_low, ymax = ER_high),
                    col = NA, fill = 'grey') +
        geom_point() + geom_line() +
        geom_point(data = hall, col = 2)
