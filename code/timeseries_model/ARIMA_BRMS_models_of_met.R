# Basic hierarchical model of metabolism as a function of temperature and Q

# A carter
# 2/2022
library(tidyverse)
library(brms)

# load in the modern dataset:
dat <- read_csv("data/metabolism/compiled/metabolism_and_drivers.csv")
sites <- read_csv('data/siteData/NHCsite_metadata.csv') %>%
  slice(c(1:5, 7))

# load in Hall dataset:

hall <- readRDS("data/metabolism/compiled/met_preds_stream_metabolizer_O2.rds")
hall <- hall$preds %>% filter(era == 'then', site == 'CBP') %>%
    select(date, site, year, doy, temp.water, starts_with(c('GPP', 'ER'))) %>%
    select(-era) %>%
    group_by(date, site) %>%
    summarize(across(everything(), ~mean(.x, na.rm = T)))
hall_QT <- read_csv('data/hall_data/hall_discharge_temp_daily.csv') %>%
    group_by(date) %>%
    summarize(across(c('water_temp_C', 'discharge_m3s'),
                     ~mean(.x, na.rm = T)))

# historical air temperature trend:
airtemp <- read_csv('data/watershed/noaa_air_temp.csv') %>%
    rename(temp_C = temp_mean)

P <- dat %>%
    select(date, site,ER,  GPP, GPP.lower, GPP.upper, discharge,
           temp.water, light = PAR_surface)%>%#, slope) %>%
    group_by(site) %>%
    left_join(select(sites, site = sitecode, slope = slope_nhd)) %>%
    filter(site %in% c('CBP', 'NHC')) %>%
    mutate(light = case_when(light == 0 ~ NA_real_,
                             TRUE ~ light)) %>%
    mutate(log_light = log(light),
           log_Q = log(discharge),
           diff_Q = c(NA, diff(log_Q)),
           RBI_7 = NA_real_,
           year = lubridate::year(date),
           year = case_when(year == 2020 ~ 2019,
                            year == 2016 ~ 2017,
                            TRUE ~ year))

for(s in c('NHC', 'CBP')){
    P1 <- P[P$site == s,] %>%
        mutate(log_Q = zoo::na.approx(log_Q, na.rm = FALSE))
    for(i in 8:nrow(P1)){
        P1$RBI_7[i] = sum(abs(P1$diff_Q[(i-6):i]))/7
    }
    P$RBI_7[P$site == s] <- P1$RBI_7
}

Q_sum <- P %>%
    group_by(site, year) %>%
    filter(lubridate::month(date) %in% c(9, 10, 11)) %>%
    summarize(med_log_Q = median(log_Q, na.rm = T),
              mean_log_Q = mean(log_Q, na.rm = T))

P <- left_join(P, Q_sum, by = c('site', 'year'))

scaling_pars <- P %>%
    ungroup() %>%
    summarize(across(c('log_Q', 'diff_Q', 'RBI_7', 'light', 'log_light', 'temp.water'),
                     .fns = c(mean = \(x) mean(x, na.rm = TRUE),
                              sd = \(x) sd(x, na.rm = TRUE))))


# select a minimum discharge value to limit the extreme values in the historical (digitized) dataset
min_Q = min(P$discharge[P$site == 'CBP'], na.rm = T)

# select a small value to add to the observations before log transforming:
epsilon = 1e-1

P_scaled <- P %>%
    ungroup() %>%
    mutate(across(c('log_Q', 'diff_Q', 'RBI_7', 'light', 'log_light', 'temp.water'),
                  .fns = \(x) as.vector(scale(x))),
           GPP = GPP + epsilon,
           ER = -ER + epsilon,
           site = factor(site, levels = c("CBP", "NHC")),
           CBP = case_when(site == 'CBP' ~ 1,
                           TRUE ~ 0)) %>%
    slice(-c(1, 362:386, 1506)) # remove leading and ending NA's in metabolism from each site

# P_scaled <- P %>%
#     mutate(across(c('log_Q', 'diff_Q', 'RBI_7', 'light', 'temp.water'),
#                   .fns = \(x) as.vector(scale(x))),
#            sin_time = sin(2*pi *as.numeric(format(date, '%j'))/365),
#            cos_time = cos(2*pi *as.numeric(format(date, '%j'))/365),
#            month = month(date),
#            season = case_when(month %in% c(12, 1, 2) ~ 'Winter',
#                               month %in% c(3, 4, 5) ~ 'Spring',
#                               month %in% c(6, 7, 8) ~ 'Summer',
#                               month %in% c(9, 10, 11) ~ 'Fall'),
#            GPP = GPP + epsilon,
#            ER = -ER + epsilon,
#            site = factor(site, levels = c("CBP", "NHC")))

# plot(density(log(P_scaled$GPP), na.rm = T))
# plot(density(log(-P_scaled$ER), na.rm = T))

P$light_smooth <- (zoo::rollmean(P$light, k = 7, fill = NA))
# add light data to hall:
light <- P %>%
    filter(site == 'NHC') %>%
    mutate(doy = format(date, '%j')) %>%
    dplyr::select(doy, light_smooth) %>%
    group_by(doy) %>%
    summarize(light = mean(light_smooth, na.rm = T))
plot(light)

# add light to Hall daatset and scale
hall_QT <- hall_QT %>%
    mutate(doy = format(date, '%j'),
           # don't allow discharge to be less than 50% of the minimum we observed at the same site
           discharge_m3s = case_when(discharge_m3s < 0.5*min_Q ~ 0.5 * min_Q,
                                     TRUE ~ discharge_m3s)) %>%
    left_join(light, by = 'doy') %>%
    mutate(site = 'CBP',
           log_Q = log(discharge_m3s),
           diff_Q = c(NA, diff(log_Q)),
           RBI_7 = NA_real_)
Q_sum_hall <- hall_QT %>%
    filter(lubridate::month(date) %in% c(9, 10, 11)) %>%
    summarize(med_log_Q = median(log_Q, na.rm = T),
              mean_log_Q = mean(log_Q, na.rm = T))

for(i in 8:nrow(hall_QT)){
    hall_QT$RBI_7[i] = sum(abs(hall_QT$diff_Q[(i-6):i]))/7
}

hall_preds <- hall_QT %>%
    dplyr::select(date, site,  temp.water = water_temp_C, log_Q, diff_Q, RBI_7, light) %>%
    mutate(across(c('log_Q', 'diff_Q', 'RBI_7', 'light', 'temp.water'),
                  .fns = \(x) zoo::na.approx(x, na.rm = F)),
           site = factor(site, levels = c("CBP", "NHC"))) %>%
    slice_head(n = -4)

# hall_preds <- hall_QT %>%
#     dplyr::select(date, site,  temp.water = water_temp_C, log_Q, diff_Q, RBI_7, light) %>%
#     mutate(across(c('log_Q', 'diff_Q', 'RBI_7', 'light', 'temp.water'),
#                   .fns = \(x) zoo::na.approx(x, na.rm = F)),
#            sin_time = sin(2*pi *as.numeric(format(date, '%j'))/365),
#            cos_time = cos(2*pi *as.numeric(format(date, '%j'))/365),
#            month = month(date),
#            season = case_when(month %in% c(12, 1, 2) ~ 'Winter',
#                               month %in% c(3, 4, 5) ~ 'Spring',
#                               month %in% c(6, 7, 8) ~ 'Summer',
#                               month %in% c(9, 10, 11) ~ 'Fall'),
#            site = factor(site, levels = c("CBP", "NHC"))) %>%
#     slice_head(n = -4)

hall_preds$med_log_Q = Q_sum_hall$med_log_Q
hall_preds$mean_log_Q = Q_sum_hall$mean_log_Q

hall_scaled <- hall_preds %>%
    mutate(diff_Q = abs(c(NA, diff(log_Q))),
           log_Q = (log_Q - scaling_pars$log_Q_mean)/scaling_pars$log_Q_sd,
           diff_Q = (diff_Q - scaling_pars$diff_Q_mean)/scaling_pars$diff_Q_sd,
           RBI_7 = (RBI_7 - scaling_pars$RBI_7_mean)/scaling_pars$RBI_7_sd,
           light = (light - scaling_pars$light_mean)/scaling_pars$light_sd,
           temp.water = (temp.water - scaling_pars$temp.water_mean)/scaling_pars$temp.water_sd)


################################################################################

# try with brms:
library(brms)
P_scaled <- P_scaled %>%
    mutate(log_GPP = log(GPP),
           log_ER = log(ER))

bform_GPP <- bf(log_GPP | mi() ~ ar(p = 1) + (1|site) + temp.water + light)
get_prior(bform_GPP, data = P_scaled)
GPP_priors <- c(prior("normal(0,5)", class = "b", coef = "light"),
                prior("normal(0,1)", class = "b", coef = "temp.water"),
                prior("beta(1,1)", class = "ar", lb = 0, ub = 1),
                prior("normal(0,5)", class = "Intercept"),
                prior("cauchy(0,1)", class = "sigma"))

bmod_GPP <- brm(bform_GPP,
                data = P_scaled,
                chains = 4, cores = 4, iter = 4000,
                prior = GPP_priors,
                control = list(adapt_delta = 0.999,
                               stepsize = 0.01,
                               max_treedepth = 14) )

# evaluate model fit:
saveRDS(bmod_GPP, 'data/timeseries_model_fits/brms_GPP_mod.rds')
# bmod_GPP <- readRDS('data/timeseries_model_fits/brms_GPP_mod.rds')

summary(bmod_GPP)

png('figures/SI/BRMS_gpp_posterior_pred_check.png', width = 5, height = 4,
    units = 'in', res = 300)

y_rep <- posterior_predict(bmod_GPP)
n_sims <- nrow(y_rep)
plot(density(P_scaled$log_GPP, na.rm = T), main = 'Posterior predictions of GPP',
     xlab = 'log GPP')
for(s in sample(n_sims, 100)){
    lines(density(y_rep[s,]), col = 'grey', alpha = 0.3)
}
lines(density(P_scaled$log_GPP, na.rm = T))

dev.off()

library(modelr)
library(tidybayes)

P_scaled <- P_scaled %>%
    mutate(siteyear = paste(site, year, sep = " "))

png('figures/SI/BRMS_gpp_model_fit.png', width = 7.5, height = 7,
    units = 'in', res = 300)
P_scaled %>%
    group_by(siteyear) %>%
    add_epred_draws(bmod_GPP, ndraws = 100) %>%
    mutate(.epred = exp(.epred + epsilon)) %>%
    ggplot(aes(date, GPP)) +
    stat_lineribbon(aes(y = .epred)) +
    # geom_line(aes(y = .epred, group = paste(siteyear, .draw)), alpha = .1) +
    geom_point(data = P_scaled, col = '#4F716B') +
    scale_fill_brewer(palette = "Greys") +
    facet_wrap(.~siteyear, scales = 'free', ncol = 1, strip.position = 'right')+
    ylab(expression(paste('GPP (g ', O[2], m^-2, d^-1, ')'))) +
    xlab('Date')+
    theme_bw()

dev.off()

draws_fit <- as_draws_array(bmod_GPP)
GPP_pars <- posterior::summarize_draws(draws_fit) %>%
    filter(variable %in%  c('b_Intercept', 'b_temp.water', 'b_light',
                            'ar[1]', 'sigma', 'r_site[CBP,Intercept]',
                            'r_site[NHC,Intercept]', 'sd_site__Intercept')) %>%
    mutate(model = 'GPP')

pd <- posterior::subset_draws(draws_fit,
                        variable = c('b_Intercept', 'b_temp.water', 'b_light',
                                     'ar[1]', 'sigma', 'r_site[CBP,Intercept]',
                                     'r_site[NHC,Intercept]'),
                        draw = 1:2000)


post_GPP <- data.frame(matrix(pd, nrow = 2000, byrow = FALSE))
colnames(post_GPP) <- c('intercept', 'b_temp.water', 'b_light', 'phi', 'sigma',
                  'cbp_intercept', 'nhc_intercept')

# matrix of draws from the posterior-predictive distribution
draws <- nrow(post_GPP)
post_preds <- matrix(nrow = draws, ncol = nrow(hall_scaled))

# fill in first observation based on the mean from the measured historical values at the site
post_preds[, 1] <- matrix(
    rep(log(mean(hall$GPP, na.rm = T) + epsilon), each = draws),
    nrow = draws, ncol = 1
)

post_preds_err <- post_preds


for(i in 1:draws){
    for(t in 2:nrow(hall_scaled)){
        post_preds[i, t] <- (1 - post_GPP$phi[i]) * (post_GPP$intercept[i] + post_GPP$cbp_intercept[i]) +
            post_GPP$b_light[i] * hall_scaled$light[t] +
            post_GPP$b_temp.water[i] * hall_scaled$temp.water[t] +
            post_GPP$phi[i] * post_preds[i, t-1] -
            post_GPP$phi[i] * (post_GPP$b_light[i] * hall_scaled$light[t-1] +
            post_GPP$b_temp.water[i] * hall_scaled$temp.water[t-1])
    }

    post_preds_err[i,] <- rnorm(nrow(hall_scaled), post_preds[i,], post_GPP$sigma[i])
}



hindcast_GPP <- data.frame(
    date = hall_scaled$date,
    GPP_pred = apply(post_preds, 2, mean),
    GPP_low = apply(post_preds, 2, quantile, probs = 0.025),
    GPP_high = apply(post_preds, 2, quantile, probs = 0.975),
    GPP_err_low = apply(post_preds_err, 2, quantile, probs = 0.025),
    GPP_err_high = apply(post_preds_err, 2, quantile, probs = 0.975)
) %>%
    mutate(GPP = exp(GPP_pred) - epsilon,
           GPP_low = exp(GPP_low)- epsilon,
           GPP_high = exp(GPP_high - epsilon),
           GPP_err_low = exp(GPP_err_low)- epsilon,
           GPP_err_high = exp(GPP_err_high - epsilon))


# ggplot(hindcast_GPP, aes(date, GPP), col = 2) +
#     geom_ribbon(aes(ymin = GPP_err_low, ymax = GPP_err_high),
#                 col = NA, fill = 'grey80') +
#     geom_ribbon(aes(ymin = GPP_low, ymax = GPP_high),
#                 col = NA, fill = 'grey50') +
#     geom_line() +
#     geom_point(data = hall,  col = 2)
#
# plot(hindcast_GPP$date, hindcast_GPP$GPP, type = 'l', ylim = c(0,10))
# for(s in sample(draws, 100)){
#     lines(hindcast_GPP$date, exp(post_preds[s,])-epsilon, col = alpha('grey', 0.3))
# }
#
# lines(hindcast_GPP$date, hindcast_GPP$GPP)
# points(hall$date, hall$GPP, col = 2, pch = 20)


############## model for ER

bform_ER <- bf(log_ER | mi() ~ ar(p = 1) + (1|site) + temp.water*mean_log_Q + light)
# bform_ER <- bf(log_ER | mi() ~ ar(p = 1) + (1|site) + temp.water*med_log_Q + light)
get_prior(bform_ER, data = P_scaled)
ER_priors <- c(prior("normal(0,1)", class = "b"),
               prior("normal(0,1)", class = "b", coef = "temp.water"),
               prior("normal(0,5)", class = "b", coef = "light"),
               prior("normal(0,5)", class = "Intercept"),
               prior("cauchy(0,1)", class = "sigma"),
               prior("beta(1,1)", class = "ar", lb = 0, ub = 1))
bmod_ER <- brm(bform_ER,
               data = P_scaled,
               prior = ER_priors,
               chains = 4, cores = 4, iter = 4000,
               control = list(adapt_delta = 0.999,
                              stepsize = 0.01,
                              max_treedepth = 14))

saveRDS(bmod_ER, 'data/timeseries_model_fits/brms_ER_mod.rds')
# bmod_ER <- readRDS('data/timeseries_model_fits/brms_ER_mod.rds')

summary(bmod_ER)

png('figures/SI/BRMS_er_posterior_pred_check.png', width = 5, height = 4,
    units = 'in', res = 300)

y_rep <- posterior_predict(bmod_ER)
n_sims <- nrow(y_rep)
plot(density(P_scaled$log_ER, na.rm = T), main = 'Posterior predictions of ER',
     xlab = 'log ER')
for(s in sample(n_sims, 100)){
    lines(density(y_rep[s,]), col = 'grey', alpha = 0.3)
}
lines(density(P_scaled$log_ER, na.rm = T))

dev.off()

png('figures/SI/BRMS_er_model_fit.png', width = 7.5, height = 7,
    units = 'in', res = 300)
P_scaled %>%
    group_by(siteyear) %>%
    add_epred_draws(bmod_ER, ndraws = 100) %>%
    mutate(.epred = exp(.epred + epsilon)) %>%
    ggplot(aes(date, ER)) +
    stat_lineribbon(aes(y = .epred)) +
    # geom_line(aes(y = .epred, group = paste(siteyear, .draw)), alpha = .1) +
    geom_point(data = P_scaled, col = '#A2865C') +
    scale_fill_brewer(palette = "Greys") +
    facet_wrap(.~siteyear, scales = 'free', ncol = 1, strip.position = 'right')+
    ylab(expression(paste('ER (g ', O[2], m^-2, d^-1, ')'))) +
    xlab('Date')+
    theme_bw()

dev.off()
# args_y <- list(bmod_ER)
# args_ypred <- list(bmod_ER)
#
#     y <- do_call(get_y, args_y)
#     ypred <- do_call(posterior_epred, args_ypred)
#     if (is_ordinal(family(object, resp = resp[i]))) {
#         ypred <- ordinal_probs_continuous(ypred)
#     }
#     R2[[i]] <- .bayes_R2(y, ypred)
# }
# R2 <- do_call(cbind, R2)
# colnames(R2) <- paste0("R2", resp)
# if (summary) {
#     R2 <- posterior_summary(R2, probs = probs, robust = robust)
# }
# R2
# function (y, ypred, ...)
# {
#     e <- t(t(post_preds_ER) - P_scaled$log_ER)
#     var_ypred <- matrixStats::rowVars(post_preds_ER)
#     var_e <- matrixStats::rowVars(e)
#     as.matrix(var_ypred/(var_ypred + var_e))
# }

draws_fit_ER <- as_draws_array(bmod_ER)
ER_pars <- posterior::summarize_draws(draws_fit_ER) %>%
    filter(variable %in% c('b_Intercept', 'b_temp.water', 'b_light',
                           'b_mean_log_Q', 'b_temp.water:mean_log_Q',
                           'ar[1]', 'sigma',
                           'r_site[CBP,Intercept]', 'r_site[NHC,Intercept]',
                           'sd_site__Intercept') ) %>%
    mutate(model = "ER")

pd <- posterior::subset_draws(draws_fit_ER,
                              variable = c('b_Intercept', 'b_temp.water', 'b_light',
                                           'b_mean_log_Q', 'b_temp.water:mean_log_Q',
                                           'ar[1]', 'sigma',
                                           'r_site[CBP,Intercept]',
                                           'r_site[NHC,Intercept]'),
                              draw = 1:2000)
post_ER <- data.frame(matrix(pd, nrow = 2000, byrow = FALSE))
colnames(post_ER) <- c('intercept', 'b_temp.water', 'b_light', 'b_meanlogQ',
                       'b_temp.water_meanlogQ', 'phi',
                       'sigma', 'cbp_intercept', 'nhc_intercept')

# matrix of draws from the posterior-predictive distribution
draws <- nrow(post_ER)
post_preds_ER <- matrix(nrow = draws, ncol = nrow(hall_scaled))

# fill in first observation based on the mean from the measured historical values at the site
post_preds_ER[, 1] <- matrix(
    rep(log(mean(-hall$ER, na.rm = T) + epsilon), each = draws),
    nrow = draws, ncol = 1
)
post_preds_err_ER <- post_preds_ER

for(i in 1:draws){
    for(t in 2:nrow(hall_scaled)){
        post_preds_ER[i, t] <- (1 - post_ER$phi[i]) * (post_ER$intercept[i] + post_ER$cbp_intercept[i]) +
            post_ER$b_light[i] * hall_scaled$light[t] +
            post_ER$b_temp.water[i] * hall_scaled$temp.water[t] +
            post_ER$b_meanlogQ[i] * hall_scaled$mean_log_Q[t] +
            post_ER$b_temp.water_meanlogQ[i] * hall_scaled$temp.water[t] * hall_scaled$mean_log_Q[t] +
            post_ER$phi[i] * post_preds_ER[i, t-1] -
            post_ER$phi[i] * (post_ER$b_light[i] * hall_scaled$light[t-1] +
                               post_ER$b_temp.water[i] * hall_scaled$temp.water[t-1] +
                               post_ER$b_meanlogQ[i] * hall_scaled$mean_log_Q[t-1] +
                               post_ER$b_temp.water_meanlogQ[i] * hall_scaled$temp.water[t-1] *
                                  hall_scaled$mean_log_Q[t-1])
    }
    post_preds_err_ER[i,] <- rnorm(nrow(hall_scaled), post_preds_ER[i,], post_ER$sigma[i])

}

hindcast_ER <- data.frame(
    date = hall_scaled$date,
    ER_pred = apply(post_preds_ER, 2, mean),
    ER_low = apply(post_preds_ER, 2, quantile, probs = 0.025),
    ER_high = apply(post_preds_ER, 2, quantile, probs = 0.975),
    ER_err_low = apply(post_preds_err_ER, 2, quantile, probs = 0.025),
    ER_err_high = apply(post_preds_err_ER, 2, quantile, probs = 0.975)
) %>%
    mutate(ER = -exp(ER_pred) + epsilon,
           ER_low = -exp(ER_low) + epsilon,
           ER_high = -exp(ER_high) + epsilon,
           ER_err_low = -exp(ER_err_low) + epsilon,
           ER_err_high = -exp(ER_err_high + epsilon))



# plot(hindcast_ER$date, hindcast_ER$ER, type = 'l', ylim = c(-20,0))
# for(s in sample(draws, 100)){
#     lines(hindcast_ER$date, -exp(post_preds_ER[s,]) + epsilon,alpha = 0.1)
# }
#
# lines(hindcast_ER$date, hindcast_ER$ER)
# points(hall$date, hall$ER, col = 2, pch = 20)


# ggplot(hindcast_ER, aes(date, ER), col = 2) +
#     geom_ribbon(aes(ymin = ER_err_low, ymax = ER_err_high),
#                 col = NA, fill = 'grey80') +
#     geom_ribbon(aes(ymin = ER_low, ymax = ER_high),
#                 col = NA, fill = 'grey50') +
#     geom_line() +
#     geom_point(data = hall,  col = 2)

hindcast <- full_join(hindcast_GPP, hindcast_ER, by = 'date')

hindcast <- hall %>%
    select(date, GPP_measured = GPP, ER_measured = ER) %>%
    group_by(date) %>%
    summarize(across(c('GPP_measured', 'ER_measured'), mean)) %>%
    right_join(hindcast, by = 'date')

bind_rows(GPP_pars, ER_pars) %>%
    write_csv('data/timeseries_model_fits/BRMS_hindcast_models_coefficients.csv')

################################################################################
# Make Figures:

hindcast <- left_join(hindcast, select(hall_QT, date, discharge_m3s))

# tiff('figures/BRMS_hindcast_comparison_daily.tiff', width = 7.5, height = 4,
#      units = 'in', res = 300)
# png('figures/BRMS_hindcast_comparison_daily.png', width = 7.5, height = 4,
#      units = 'in', res = 300)

x0date <- as.Date('1968-04-08')
mplot <- ggplot(hindcast) +
    geom_ribbon(aes(x = date, ymin = ER_err_high, ymax = ER_err_low),
                col = NA, fill = 'grey80') +
    geom_ribbon(aes(x = date, ymin = ER_high, ymax = ER_low),
                col = NA, fill = 'grey60') +
    geom_ribbon(aes(x = date, ymin = GPP_err_low, ymax = GPP_err_high),
                col = NA, fill = 'grey80') +
    geom_ribbon(aes(x = date, ymin = GPP_low, ymax = GPP_high),
                col = NA, fill = 'grey60') +
    geom_line(aes(x = date, y = ER, color = "AR1 Hindcast")) +
    geom_line(aes(x = date, y = GPP, color = "AR1 Hindcast")) +
    geom_point(aes(x = date, y = ER_measured, color = "Historical Data")) +
    geom_point(aes(x = date, y = GPP_measured, color = "Historical Data")) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = 'dashed') +
    ylab(expression(paste('Metabolism (g ', O[2], m^-2, d^-1, ')'))) +
    xlab('Date')+
    theme_classic() +
    # Custom legend
    scale_color_manual(
        values = c("AR1 Hindcast" = "black", "Historical Data" = "brown3")
    ) +
    guides(color = guide_legend(override.aes = list(shape = c(NA, 16),
                                                    linetype = c(1, 0)),
                                ncol = 2)) +
    theme(
        legend.background = element_rect(fill = 'transparent'),
        legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.title = element_blank(),
        # legend.box.background = element_rect(color = "black", linewidth = 0.2),
        legend.text = element_text(margin = margin(r = -5, unit = "pt")),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        axis.line.y.left = element_blank(),
        # axis.line.y.left = element_line(color = "black", linewidth = 0.2),
        axis.line.x.bottom = element_blank(),
        plot.margin = margin(b = 0)
    ) +
    scale_y_continuous(limits = c(-4.5, 1.6), expand = c(0.39, 0),
                       oob = scales::squish_infinite) +
    geom_segment(x = x0date, xend = x0date, y = -6.2, yend = 2.8,
                 color = "gray25", linewidth = 0.2) +
    coord_cartesian(clip = "off") +
    scale_x_date(expand = c(0, 0))


qplot <- ggplot(hindcast) +
    # geom_line(aes(x = date, y = discharge_m3s, color = "Discharge")) +
    geom_line(aes(x = date, y = discharge_m3s), color = "darkblue") +
    ylab(expression('Discharge (m'^3 * 's'^-1 * ')')) +
    xlab('Date') +
    theme_classic() +
    theme(axis.line.y.left = element_line(color = "black", linewidth = 0.2),
          axis.line.x.bottom = element_line(color = "black", linewidth = 0.2),
          panel.border = element_blank(),
          plot.margin = margin(t = 0)) +
    scale_x_date(expand = c(0, 0))
    # scale_color_manual(
    #     values = c("Discharge" = "darkblue")
    # )
    # guides(color = guide_legend(override.aes = list(shape = c(NA, 16),
    #                                                 linetype = c(1, 0)),
    #                             ncol = 2)) +
    # theme(
    #     legend.background = element_blank(),
    #     legend.position = c(0.02, 0.95), # Upper-left inside the plot
    #     legend.justification = c(0, 1),
    #     legend.title = element_blank(),
    # )

tiff('figures/BRMS_hindcast_comparison_daily.tiff', width = 7.5, height = 4,
     units = 'in', res = 300)
ggpubr::ggarrange(mplot, qplot, heights = c(3, 1), ncol = 1, align = 'v')
    # annotate("text", x = 0.09, y = 0.85, label = "A") +
    # annotate("text", x = 0.56, y = 0.85, label = "B")
dev.off()

png('figures/BRMS_hindcast_comparison_daily.png', width = 7.5, height = 4,
     units = 'in', res = 300)
ggpubr::ggarrange(mplot, qplot, heights = c(3, 1), ncol = 1, align = 'v') +
    theme(plot.margin = margin(r = 2))
    # annotate("text", x = 0.09, y = 0.85, label = "A") +
    # annotate("text", x = 0.56, y = 0.85, label = "B")
dev.off()


# compare the predictions averaged across months like in Hall 1972
hall_pred_sum <- hindcast %>%
    mutate(month = month(date)) %>%
    # group_by(month) %>%
    summarize(across(starts_with(c('GPP','ER')),
                     .fns = c(mean = \(x) mean(x, na.rm = T),
                              sd = \(x) sd(x, na.rm = T)))) %>%
    ungroup()
hall_sum <- hall %>%
    mutate(month = month(date)) %>%
    group_by(month) %>%
    summarize(across(c('GPP','ER'),
                     .fns = c(mean = \(x) mean(x, na.rm = T),
                              sd = \(x) sd(x, na.rm = T)))) %>%
    ungroup()

hind_sum <- filter(hindcast, year(date) == 1969) %>%
    arrange(date)

# tiff('figures/BRMS_hindcast_comparison_monthly.tiff', width = 7.5, height = 4,
#      units = 'in', res = 300)
png('figures/BRMS_hindcast_comparison_monthly.png', width = 7.5, height = 4,
     units = 'in', res = 300)


ggplot(hall_pred_sum, aes(month, GPP_mean)) +
    # Ribbons for GPP and ER with grey fill (Arima Hindcast)
    geom_ribbon(aes(ymin = GPP_mean - GPP_sd,
                    ymax = GPP_mean + GPP_sd, fill = "AR1 Hindcast"),
                col = NA, fill = 'grey') +
    geom_ribbon(aes(ymin = ER_mean - ER_sd,
                    ymax = ER_mean + ER_sd, fill = "AR1 Hindcast"),
                col = NA, fill = 'grey') +

    # Lines for GPP and ER predictions (black, Arima Hindcast)
    geom_line(aes(y = GPP_mean, color = "AR1 Hindcast")) +
    geom_line(aes(y = ER_mean, color = "AR1 Hindcast")) +
    geom_hline(yintercept = 0) +

    # Lines and error bars for actual data (red, Historical Data)
    geom_line(data = hall_sum, aes(y = GPP_mean, color = "Historical Data")) +
    geom_errorbar(data = hall_sum, aes(ymin = GPP_mean - GPP_sd,
                                       ymax = GPP_mean + GPP_sd, color = "Historical Data"),
                  width = 0.2) +
    geom_line(data = hall_sum, aes(y = ER_mean, color = "Historical Data")) +
    geom_errorbar(data = hall_sum, aes(ymin = ER_mean - ER_sd,
                                       ymax = ER_mean + ER_sd, color = "Historical Data"),
                  width = 0.2) +
    ylab(expression(paste('Metabolism (g ', O[2], m^-2, d^-1, ')'))) +
    xlab('')+

    # Custom legend
    scale_color_manual(
        name = "Legend",
        values = c("AR1 Hindcast" = "black", "Historical Data" = "brown3"),
        labels = c("AR1 Hindcast", "Historical Data")
    ) +
    scale_fill_manual(
        name = "Legend",
        values = c("AR1 Hindcast" = "grey"),
        labels = c("AR1 Hindcast")
    ) +
    scale_x_continuous(
        breaks = 1:12, # Assuming months are represented as 1-12
        labels = month.abb
    ) +
    # Adjust the legend to show color, linetype, and points
    guides(
        color = guide_legend(override.aes = list(linetype = c(1, 1), shape = NA)),
        fill = guide_legend(override.aes = list(linetype = 0))
    ) +

    # Apply theme adjustments
    theme_bw() +
    theme(
        legend.position = c(0.88, 0.88), # Adjust as needed for legend position
        legend.title = element_blank()
    )
dev.off()


################################################################################
# Predict change over time in annual metabolism across 50 years

# average discharge across today and 1970 datasets:
comb_dat <- bind_rows(select(P[P$site == 'CBP',], date, log_Q, temp.water),
          select(hall_preds[hall_preds$site == 'CBP',], date, log_Q, temp.water)) %>%
    mutate(doy = format(date, '%j')) %>%
    left_join(airtemp, by = 'date') %>% ungroup()

LQ <- comb_dat %>%
    group_by(doy) %>%
    summarize(log_Q = mean(log_Q, na.rm = T)) %>%
    left_join(light, by = 'doy')

# interpolate water temperature based on air temperature

temp_mod <- left_join(airtemp,
                     select(comb_dat, date, temp.water),
                     by = 'date')

ggplot(temp_mod, aes(temp_C, temp.water, col = date))+
    geom_point()
ggplot(temp_mod, aes(date, temp_C))+
    geom_line()

t_mod <- lm(temp.water ~ temp_C, temp_mod)


# generate data frame with discharge and light:
met_change <- data.frame(
    date = seq(as.Date('1968-01-01'), as.Date('2019-12-31'), by = 'day')
)

met_change <- met_change %>%
    mutate(doy = format(date, '%j'),
           year = lubridate::year(date)) %>%
    left_join(LQ, by = 'doy') %>%
    left_join(airtemp, by = 'date') %>%
    mutate(light = zoo::na.approx(light, na.rm = F))
met_change$temp.water <- predict(t_mod, newdata = met_change)
Q_sum <- met_change %>%
    group_by(year) %>%
    filter(lubridate::month(date) %in% c(9, 10, 11)) %>%
    summarize(mean_log_Q = mean(log_Q, na.rm = T))

Q_sum2 <- read_csv('data/rating_curves/modeled_fall_mean_NHC_flow.csv') %>%
    mutate(mean_log_Q = log(fall_mean))


met_change_scaled <- met_change %>%
    left_join(Q_sum, by = 'year') %>%
    select(-temp_C) %>%
    mutate(log_Q = (log_Q - scaling_pars$log_Q_mean)/scaling_pars$log_Q_sd,
           light = (light - scaling_pars$light_mean)/scaling_pars$light_sd,
           temp.water = (temp.water - scaling_pars$temp.water_mean)/scaling_pars$temp.water_sd)

met_change_scaled <- met_change %>%
    right_join(Q_sum2, by = 'year') %>%
    select(-temp_C) %>%
    mutate(log_Q = (log_Q - scaling_pars$log_Q_mean)/scaling_pars$log_Q_sd,
           light = (light - scaling_pars$light_mean)/scaling_pars$light_sd,
           temp.water = (temp.water - scaling_pars$temp.water_mean)/scaling_pars$temp.water_sd)

met_change2 <- met_change %>%
    group_by(year) %>%
    summarize(temp.water = mean(temp.water, na.rm = T))%>%
    right_join(Q_sum2, by = 'year')  %>%
    mutate(mean_Q = exp(mean_log_Q)) %>%
    ungroup()

# matrix of draws from the posterior-predictive distribution
post_preds_GPP <- matrix(nrow = draws, ncol = nrow(met_change_scaled))
post_preds_ER <- matrix(nrow = draws, ncol = nrow(met_change_scaled))

# fill in first observation based on the mean from the measured historical values at the site
post_preds_GPP[, 1] <- matrix(
    rep(log(mean(hall$GPP, na.rm = T) + epsilon), each = draws),
    nrow = draws, ncol = 1
)
post_preds_ER[, 1] <- matrix(
    rep(log(mean(-hall$ER, na.rm = T) + epsilon), each = draws),
    nrow = draws, ncol = 1
)
post_preds_err_GPP <- post_preds_GPP
post_preds_err_ER <- post_preds_ER

for(i in 1:draws){
    for(t in 2:nrow(met_change_scaled)){
        post_preds_GPP[i, t] <- (1 - post_GPP$phi[i]) * (post_GPP$intercept[i] + post_GPP$cbp_intercept[i]) +
            post_GPP$b_light[i] * met_change_scaled$light[t] +
            post_GPP$b_temp.water[i] * met_change_scaled$temp.water[t] +
            post_GPP$phi[i] * post_preds_GPP[i, t-1] -
            post_GPP$phi[i] * (post_GPP$b_light[i] * met_change_scaled$light[t-1] +
                                   post_GPP$b_temp.water[i] * met_change_scaled$temp.water[t-1])

        post_preds_ER[i, t] <- (1 - post_ER$phi[i]) * (post_ER$intercept[i] + post_ER$cbp_intercept[i]) +
            post_ER$b_light[i] * met_change_scaled$light[t] +
            post_ER$b_temp.water[i] * met_change_scaled$temp.water[t] +
            post_ER$b_meanlogQ[i] * met_change_scaled$mean_log_Q[t] +
            post_ER$b_temp.water_meanlogQ[i] * met_change_scaled$temp.water[t] *
                met_change_scaled$mean_log_Q[t] +
            post_ER$phi[i] * post_preds_ER[i, t-1] -
            post_ER$phi[i] * (post_ER$b_light[i] * met_change_scaled$light[t-1] +
                                  post_ER$b_temp.water[i] * met_change_scaled$temp.water[t-1] +
                                  post_ER$b_meanlogQ[i] * met_change_scaled$mean_log_Q[t-1] +
                                  post_ER$b_temp.water_meanlogQ[i] * met_change_scaled$temp.water[t-1] *
                                    met_change_scaled$mean_log_Q[t-1])
    }

    post_preds_err_GPP[i,] <- rnorm(nrow(met_change_scaled), post_preds_GPP[i,], post_GPP$sigma[i])
    post_preds_err_ER[i,] <- rnorm(nrow(met_change_scaled), post_preds_ER[i,], post_ER$sigma[i])

}


hindcast <- data.frame(
    date = met_change_scaled$date,
    doy = as.numeric(format(met_change_scaled$date, '%j')),
    year = year(met_change_scaled$date),
    GPP_pred = apply(post_preds_GPP, 2, mean),
    GPP_low = apply(post_preds_GPP, 2, quantile, probs = 0.025),
    GPP_high = apply(post_preds_GPP, 2, quantile, probs = 0.975),
    # GPP_err_low = apply(post_preds_err_GPP, 2, quantile, probs = 0.025),
    # GPP_err_high = apply(post_preds_err_GPP, 2, quantile, probs = 0.975),
    ER_pred = apply(post_preds_ER, 2, mean),
    ER_low = apply(post_preds_ER, 2, quantile, probs = 0.025),
    ER_high = apply(post_preds_ER, 2, quantile, probs = 0.975)
    # ER_err_low = apply(post_preds_err_ER, 2, quantile, probs = 0.025),
    # ER_err_high = apply(post_preds_err_ER, 2, quantile, probs = 0.975)
) %>%
    mutate(GPP = exp(GPP_pred) - epsilon,
           GPP_low = exp(GPP_low)- epsilon,
           GPP_high = exp(GPP_high - epsilon),
           # GPP_err_low = exp(GPP_err_low)- epsilon,
           # GPP_err_high = exp(GPP_err_high - epsilon),
           ER = exp(ER_pred) - epsilon,
           ER_low = exp(ER_low)- epsilon,
           ER_high = exp(ER_high - epsilon)
           # ER_err_low = exp(ER_err_low)- epsilon,
           # ER_err_high = exp(ER_err_high - epsilon)
           )


library(viridis)
met_plot <- hindcast %>%
    mutate(GPP = zoo::rollmean(GPP, k = 7, na.pad = TRUE)) %>%
    mutate(ER = -zoo::rollmean(ER, k = 7, na.pad = TRUE)) %>%
    # group_by(doy) %>%
    # mutate(GPP = zoo::rollmean(GPP, k = 5, na.pad = TRUE)) %>%
    # mutate(ER = zoo::rollmean(ER, k = 5, na.pad = TRUE)) %>%
    pivot_longer(cols = c('GPP', 'ER'),
                 names_to = 'met',
                 values_to = 'gO2') %>%
    filter(year %% 3 == 1) %>%
    mutate(met = factor(met, levels = c('GPP', 'ER'))) %>%
    ggplot(aes(doy, gO2, col = year, group = factor(year))) +
    geom_line(linewidth = 0.7) +
    # geom_line(aes(y = -ER), linewidth = 1) +
    scale_color_viridis(name = 'Year', option = 'C', begin = 0, end = 0.9) +
    facet_wrap(.~met, scales = 'free_y', ncol = 1, strip.position = 'right')+
    ylab(expression(paste('ER (g ', O[2], m^-2, y^-1, ')                    GPP (g ', O[2], m^-2, y^-1, ')'))) +
    xlab('Month') +
    scale_x_continuous(breaks = c(0,31,60,91,121,152,182,213,244,274,305,335),
                       labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          # legend.key.width = unit(0.3, 'in'),
          legend.position = 'inside',
          legend.position.inside = c(.09, .18))
          # legend.direction = 'horizontal')
          # legend.position = 'top')


# tiff('figures/change_in_annual_met_preds_through_time.tiff',
#      width = 6, height = 3.5, res = 300, units = 'in')
# met_plot
# dev.off()

met_change3 <- met_change2 %>%
    rename(`Mean Fall Discharge` = mean_Q, `Mean Annual Temperature` = temp.water) %>%
    pivot_longer(cols = c('Mean Fall Discharge', 'Mean Annual Temperature'),
                 names_to = 'variable', values_to = 'value')
met_dot <- filter(met_change3, year %% 3 == 1)
clim_plot <- met_change3 %>%
    ggplot(aes(year, value, col = year)) +
    geom_line(linewidth = 1) +
    geom_point(data = met_dot, size = 2) +
    geom_point(data = met_dot, col = 'black', pch = 1, size = 2) +
    scale_color_viridis(name = 'Year', option = 'C', begin = 0, end = 0.9) +
    facet_wrap(.~variable, scales = 'free_y', ncol = 1) +
    xlab('year') +
    ylab(expression(paste('Mean Fall Discharge (', m^3, s^-1, ')           Mean Annual Temp. (', degree, 'C)')))+
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          legend.key.width = unit(0.5, 'in'),
          legend.position = 'none')


tiff('figures/Annual_hindcast_trajectory_4panel.tiff', width = 9, height = 5,
     units = 'in', res = 300)
ggpubr::ggarrange(met_plot, clim_plot, widths = c(2,1))#, common.legend = TRUE )
dev.off()
png('figures/Annual_hindcast_trajectory_4panel.png', width = 9, height = 5,
     units = 'in', res = 300)
ggpubr::ggarrange(met_plot, clim_plot, widths = c(2,1))#, common.legend = TRUE )
dev.off()

tiff('figures/Annual_hindcast_trajectory.tiff', width = 4, height = 3,
     units = 'in', res = 300)

hindcast %>%
    group_by(year) %>%
    summarize(GPP_mean = sum(GPP),
              ER_mean = sum(ER),
              GPP_low = sum(GPP_low),
              GPP_high = sum(GPP_high),
              ER_low = sum(ER_low),
              ER_high = sum(ER_high),
              GPP_sd = sd(GPP)*365,
              ER_sd = sd(ER)*365) %>%
    pivot_longer(cols = starts_with(c('GPP', 'ER')),
                 names_to = c('met', 'stat'),
                 names_sep = '_',
                 values_to = 'value') %>%
    pivot_wider(values_from = 'value', names_from = 'stat' ) %>%
    mutate(met = factor(met, levels = c('GPP', 'ER'))) %>%
    ggplot(aes(year, mean))+
    geom_ribbon(aes(ymin = low,
                    ymax = high, fill = met),
                col = NA, alpha = 0.5) +
    geom_line() +
    scale_fill_manual(name = '', values = c("#3b7c70", "#ce9642"))+
    facet_wrap(.~met, ncol = 1, scales = 'free_y', strip.position = 'right')+
    ylab(expression(paste('Annual Metabolism (g ', O[2], m^-2, y^-1, ')'))) +
    xlab('Year') +
    theme_bw() +
    theme(legend.position = 'none')

dev.off()


################################################################################
# Attempt to model using GAMS:(old code)
# Fit the model with Gamma distribution and random effects
model_mgcv <- gam(
    GPP ~ s(temp.water) + s(light) + s(site, bs = "re"),  # Random effect for site
    family = Gamma(link = "log"),                         # Gamma distribution
    data = P_scaled,
    method = "REML",                                       # Use restricted maximum likelihood
    na.action = na.exclude)


summary(model_mgcv)

fitted_values <- mgcv::predict.gam(model_mgcv, type = "response", se.fit = TRUE)
P_scaled$fitted_values <- fitted_values$fit
P_scaled$fit_se <- fitted_values$se.fit

# Plot the original data and the fitted values
ggplot(P_scaled, aes(x = date, y = GPP)) +
    geom_point(color = "blue", alpha = 0.6) +  # Original data points
    geom_line(aes(y = fitted_values), color = "red", linewidth = 1) +  # Fitted values from the model
    geom_ribbon(aes(ymin = fitted_values - 2*fit_se,
                    ymax = fitted_values + 2*fit_se))+
    labs(title = "GAM Fit to Original Data",
         x = "Temperature (Water)",
         y = "GPP") +
    theme_minimal()


GPP_preds <-  predict(model_mgcv, newdata = hall_scaled,
                      type = 'response', se.fit = TRUE)
hall_scaled$GPP_pred <- GPP_preds$fit
hall_scaled$GPP_se <- GPP_preds$se.fit
hall_scaled$GPP <- hall_scaled$GPP_pred - epsilon

ggplot(hall_scaled, aes(date, GPP))+
    geom_ribbon(aes(ymin = GPP - GPP_se*1.96,
                    ymax = GPP + GPP_se*1.96),
                fill = 'grey', col = NA) +
    geom_point() +
    geom_line() +
    geom_point(data = hall, col = 2)

model_mgcv_ER <- gam(
    ER ~ s(temp.water) + s(log_Q) + s(light) + s(site, bs = 're'),
    family = Gamma(link = "log"),
    data = P_scaled,
    method = "REML",
    na.action = na.exclude)
summary(model_mgcv_ER)
plot(model_mgcv_ER)

fitted_values <- mgcv::predict.gam(model_mgcv_ER, type = "response", se.fit = TRUE)
P_scaled$fitted_values <- fitted_values$fit
P_scaled$fit_se <- fitted_values$se.fit

# Plot the original data and the fitted values
ggplot(P_scaled, aes(x = date, y = -ER)) +
    geom_point(color = "blue", alpha = 0.6) +  # Original data points
    geom_line(aes(y = -fitted_values), color = "red", size = 1) +  # Fitted values from the model
    geom_ribbon(aes(ymin = -fitted_values + 2*fit_se,
                    ymax = -fitted_values - 2*fit_se))+
    labs(title = "GAM Fit to Original Data",
         x = "Date",
         y = "ER") +
    theme_minimal()



ER_preds <-  predict(model_mgcv_ER, newdata = hall_scaled,
                     type = 'response', se.fit = TRUE)
hall_scaled$ER_pred <- ER_preds$fit
hall_scaled$ER_se <- ER_preds$se.fit
hall_scaled$ER <- -hall_scaled$ER_pred + epsilon

ggplot(hall_scaled, aes(date, ER))+
    geom_ribbon(aes(ymin = ER - ER_se*1.96,
                    ymax = ER + ER_se*1.96),
                fill = 'grey', col = NA) +
    geom_point() +
    geom_line() +
    geom_point(data = hall, col = 2)
