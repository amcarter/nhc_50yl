# Basic hierarchical model of metabolism as a function of temperature and Q

# A carter
# 2/2022
library(tidyverse)
library(forecast)
library(mgcv)
library(brms)
library(broom)

# load in the modern dataset:
dat <- read_csv("data/metabolism/compiled/metabolism_and_drivers.csv")
sites <- read_csv('data/siteData/NHCsite_metadata.csv') %>%
  slice(c(1:5, 7))

# load in Hall dataset:

hall <- readRDS("data/metabolism/compiled/met_preds_stream_metabolizer_O2.rds")
hall <- hall$preds %>% filter(era == 'then', site == 'CBP')
hall_QT <- read_csv('data/hall_data/hall_discharge_temp_daily.csv')



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


# select a minimum discharge value to limit the extreme values in the historical (digitized) dataset
min_Q = min(P$discharge[P$site == 'CBP'], na.rm = T)

# select a small value to add to the observations before log transforming:
epsilon = 1e-1

P_scaled <- P %>%
    ungroup() %>%
    mutate(across(c('log_Q', 'light', 'temp.water'),
                  .fns = \(x) as.vector(scale(x))),
           GPP = GPP + epsilon,
           ER = -ER + epsilon,
           GPP_pred = log(GPP),
           site = factor(site, levels = c("CBP", "NHC")),
           CBP = case_when(site == 'CBP' ~ 1,
                           TRUE ~ 0)) %>%
    slice(-c(1, 362:386, 1506)) # remove leading and ending NA's in metabolism from each site
plot(density(log(P_scaled$GPP), na.rm = T))
plot(density(log(-P_scaled$ER), na.rm = T))

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
           site = factor(site, levels = c("CBP", "NHC"))) %>%
    slice_head(n = -4)


hall_scaled <- hall_preds %>%
    mutate(log_Q = (log_Q - scaling_pars$log_Q_mean)/scaling_pars$log_Q_sd,
           light = (light - scaling_pars$light_mean)/scaling_pars$light_sd,
           temp.water = (temp.water - scaling_pars$temp.water_mean)/scaling_pars$temp.water_sd,
           CBP = 1)



# fit the model using an AR1 from the forecast package:
for_mod <- forecast::Arima(y = log(P_scaled$GPP),
                           order = c(1,0,0),
                           xreg = matrix(c(P_scaled$temp.water,
                                           P_scaled$light,
                                           P_scaled$CBP),
                                         nrow = nrow(P_scaled), ncol = 3))

# Predict the fitted values
GPP_preds <- predict(for_mod, n_ahead = nrow(hall_scaled),
                     newxreg = matrix(c(hall_scaled$temp.water,
                                        hall_scaled$light,
                                        hall_scaled$CBP),
                                      nrow = nrow(hall_scaled), ncol = 3),
                     se.fit = TRUE)


hall_scaled$GPP_pred <- GPP_preds$pred
hall_scaled$GPP = exp(hall_scaled$GPP_pred) - epsilon

# plot the predictions
ggplot(hall_scaled, aes(date, GPP))+
    geom_point() + geom_line() +
    geom_point(data = hall, col = 2)


# model for ER:

for_mod_ER <- forecast::Arima(y = log(P_scaled$ER),
                              order = c(1,0,0),
                              xreg = matrix(c(P_scaled$light,
                                              P_scaled$temp.water,
                                              P_scaled$log_Q,
                                              P_scaled$CBP),
                                            nrow = nrow(P_scaled),
                                            ncol = 4))
ER_preds <- predict(for_mod_ER, n.ahead = nrow(hall_scaled),
                    newxreg = matrix(c(hall_scaled$light,
                                       hall_scaled$temp.water,
                                       hall_scaled$log_Q,
                                       hall_scaled$CBP),
                                     nrow = nrow(hall_scaled),
                                     ncol = 4))

hall_scaled$ER_pred <- ER_preds$pred
hall_scaled$ER <- -(exp(hall_scaled$ER_pred) - epsilon)


# save the model coefficients:
checkresiduals(for_mod)
summary_fit <- tidy(for_mod) %>%
    bind_rows(data.frame(term = c('sigma2', 'r2'),
                         estimate = c(summary(for_mod)$sigma2,
                                      cor(fitted(for_mod),log(P_scaled$GPP),
                                          use = 'complete.obs')^2))) %>%
    mutate(model = 'GPP',
           term = case_when(term == 'xreg1' ~ 'water_temp',
                            term == 'xreg2' ~ 'light',
                            term == 'xreg3' ~ 'CBP_intercept',
                            TRUE ~ term))
summary_fit_ER <- tidy(for_mod_ER) %>%
    bind_rows(data.frame(term = c('sigma2', 'r2'),
                         estimate = c(summary(for_mod_ER)$sigma2,
                                      cor(fitted(for_mod_ER),log(P_scaled$ER),
                                          use = 'complete.obs')^2))) %>%
    mutate(model = 'ER',
           term = case_when(term == 'xreg1' ~ 'light',
                            term == 'xreg2' ~ 'water_temp',
                            term == 'xreg3' ~ 'log_Q',
                            term == 'xreg4' ~ 'CBP_intercept',
                            TRUE ~ term))
summary_fit <- bind_rows(summary_fit, summary_fit_ER)

write_csv(summary_fit, 'data/Arima_hindcast_models_coefficients.csv')

################################################################################
# Make Figures:

tiff('figures/Arima_hindcast_comparison_daily.tiff', width = 7.5, height = 4,
     units = 'in', res = 300)
    # Q <- ggplot(hall_scaled, aes(date, exp(log_Q)))+
    #     geom_line() +
    #     ylab(expression(paste('Discharge (', m^3, s^-1,')')))+
    #     xlab('Date')+
    #     theme_bw()
    # tempC <- ggplot(hall_scaled, aes(date, temp.water))+
    #     geom_line() +
    #     ylab('Water Temperature (C)')+
    #     theme_bw()+
    #     theme(axis.title.x = element_blank())
    # light <- ggplot(hall_scaled, aes(date, (light - minL)/(maxL-minL)))+
    #     geom_line() +
    #     ylab('Relative Light')+
    #     theme_bw()+
    #     theme(axis.title.x = element_blank())
    ggplot() +
        geom_line(data = hall_scaled, aes(x = date, y = ER, color = "AR1 Hindcast")) +
        geom_line(data = hall_scaled, aes(x = date, y = GPP, color = "AR1 Hindcast")) +
        geom_point(data = hall, aes(x = date, y = ER, color = "Historical Data")) +
        geom_point(data = hall, aes(x = date, y = GPP, color = "Historical Data")) +
        geom_hline(yintercept = 0) +
        ylab(expression(paste('Metabolism (g ', O[2], m^-2, d^-1, ')'))) +
        xlab('Date')+
        theme_bw() +
        # Custom legend
        scale_color_manual(
            values = c("AR1 Hindcast" = "black", "Historical Data" = "brown3")
        ) +
        guides(color = guide_legend(override.aes = list(shape = c(NA, 16), linetype = c(1, 0)))) +
        theme(
            legend.position = c(0.02, 0.97), # Upper-left inside the plot
            legend.justification = c(0, 1),
            legend.title = element_blank(),
            # axis.title.x = element_blank()
        )
    # cowplot::plot_grid(
    #     met, tempC, light, Q,
    #     ncol = 1,             # Define the number of columns
    #     rel_heights = c(3, 1, 1, 1),      # Relative heights for rows
    #     align = 'v',
    #     axis = 'lr'
    # )
dev.off()


# compare the predictions averaged across months like in Hall 1972
hall_pred_sum <- hall_scaled %>%
    select(date, GPP, ER) %>%
    mutate(doy = format(date, '%j'),
           month = month(date)) %>%
    group_by(month) %>%
    summarize(across(c('GPP','ER'),
                     .fns = c(mean = \(x) mean(x, na.rm = T),
                              sd = \(x) sd(x, na.rm = T)))) %>%
    mutate(month_name = month.abb) %>%
    ungroup()

hall_sum <-hall %>%
    select(date, GPP, ER) %>%
    mutate(month = month(date)) %>%
    group_by(month) %>%
    summarize(across(c('GPP','ER'),
                     .fns = c(mean = \(x) mean(x, na.rm = T),
                              sd = \(x) sd(x, na.rm = T)))) %>%
    mutate(month_name = month.abb[month]) %>%
    ungroup()

tiff('figures/Arima_hindcast_comparison_monthly.tiff', width = 7.5, height = 4,
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

    # Horizontal line at y = 0
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

# try with brms:
# regenrate the prediction data frame:

hall_scaled <- hall_preds %>%
    mutate(log_Q = (log_Q - scaling_pars$log_Q_mean)/scaling_pars$log_Q_sd,
           light = (light - scaling_pars$light_mean)/scaling_pars$light_sd,
           temp.water = (temp.water - scaling_pars$temp.water_mean)/scaling_pars$temp.water_sd)

P_scaled <- P_scaled %>%
    mutate(log_GPP = log(GPP),
           log_ER = log(ER))

bform_GPP <- bf(log_GPP | mi() ~ ar(p = 1) + (1|site) + temp.water + light)
get_prior(bform_GPP, data = P_scaled)
GPP_priors <- c(prior("normal(0,5)", class = "b"),
                prior("beta(1,1)", class = "ar", lb = 0, ub = 1),
                prior("normal(0,5)", class = "Intercept"),
                prior("cauchy(0,1)", class = "sigma"))

bmod_GPP <- brm(bform_GPP,
                data = P_scaled,
                chains = 4, cores = 4, iter = 2000,
                prior = GPP_priors,
                control = list(adapt_delta = 0.99,
                               max_treedepth = 12) )

# evaluate model fit:
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


post <- data.frame(matrix(pd, nrow = 2000, byrow = FALSE))
colnames(post) <- c('intercept', 'b_temp.water', 'b_light', 'phi', 'sigma',
                  'cbp_intercept', 'nhc_intercept')

# matrix of draws from the posterior-predictive distribution
draws <- nrow(post)
post_preds <- matrix(nrow = draws, ncol = nrow(hall_scaled))

# fill in first observation based on the mean from the measured historical values at the site
post_preds[, 1] <- matrix(
    rep(log(mean(hall$GPP, na.rm = T) + epsilon), each = draws),
    nrow = draws, ncol = 1
)
post_preds_err <- post_preds

for(i in 1:draws){
    for(t in 2:nrow(hall_scaled)){
        post_preds[i, t] <- (1 - post$phi[i]) * (post$intercept[i] + post$cbp_intercept[i]) +
            post$b_light[i] * hall_scaled$light[t] +
            post$b_temp.water[i] * hall_scaled$temp.water[t] +
            post$phi[i] * post_preds[i, t-1] -
            post$phi[i] * (post$b_light[i] * hall_scaled$light[t-1] +
            post$b_temp.water[i] * hall_scaled$temp.water[t-1])#+ rnorm(1, mean = 0, sd = post$sigma[i])
    }

    post_preds_err[i,] <- rnorm(nrow(hall_scaled), post_preds[i,], post$sigma[i])
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

ggplot(hindcast_GPP, aes(date, GPP), col = 2) +
    geom_ribbon(aes(ymin = GPP_err_low, ymax = GPP_err_high),
                col = NA, fill = 'grey80') +
    geom_ribbon(aes(ymin = GPP_low, ymax = GPP_high),
                col = NA, fill = 'grey50') +
    geom_line() +
    geom_point(data = hall,  col = 2)



plot(hindcast_GPP$date, hindcast_GPP$GPP, type = 'l', ylim = c(0,10))
for(s in sample(draws, 100)){
    lines(hindcast_GPP$date, exp(post_preds[s,])-epsilon, col = alpha('grey', 0.3))
}

lines(hindcast_GPP$date, hindcast_GPP$GPP)
points(hall$date, hall$GPP, col = 2, pch = 20)
# model for ER

bform_ER <- bf(log_ER | mi() ~ ar(p = 1) + (1|site) + temp.water + light + log_Q)
get_prior(bform_ER, data = P_scaled)
ER_priors <- c(prior("normal(0,5)", class = "b"),
               prior("normal(0,5)", class = "Intercept"),
               prior("cauchy(0,1)", class = "sigma"),
               prior("beta(1,1)", class = "ar", lb = 0, ub = 1))
bmod_ER <- brm(bform_ER,
               data = P_scaled,
               prior = ER_priors,
               chains = 4, cores = 4, iter = 2000,
               control = list(adapt_delta = 0.99,
                              max_treedepth = 12))


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


brms::bayes_R2(bmod_ER)

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
                           'b_log_Q', 'ar[1]', 'sigma',
                           'r_site[CBP,Intercept]', 'r_site[NHC,Intercept]',
                           'sd_site__Intercept') ) %>%
    mutate(model = "ER")

pd <- posterior::subset_draws(draws_fit_ER,
                              variable = c('b_Intercept', 'b_temp.water', 'b_light',
                                           'b_log_Q', 'ar[1]', 'sigma',
                                           'r_site[CBP,Intercept]',
                                           'r_site[NHC,Intercept]'),
                              draw = 1:2000)
post <- data.frame(matrix(pd, nrow = 2000, byrow = FALSE))
colnames(post) <- c('intercept', 'b_temp.water', 'b_light', 'b_logQ', 'phi',
                    'sigma', 'cbp_intercept', 'nhc_intercept')

# matrix of draws from the posterior-predictive distribution
draws <- nrow(post)
post_preds_ER <- matrix(nrow = draws, ncol = nrow(hall_scaled))

# fill in first observation based on the mean from the measured historical values at the site
post_preds_ER[, 1] <- matrix(
    rep(log(mean(-hall$ER, na.rm = T) + epsilon), each = draws),
    nrow = draws, ncol = 1
)
post_preds_err_ER <- post_preds_ER

for(i in 1:draws){
    for(t in 2:nrow(hall_scaled)){
        post_preds_ER[i, t] <- (1 - post$phi[i]) * (post$intercept[i] + post$cbp_intercept[i]) +
            post$b_light[i] * hall_scaled$light[t] +
            post$b_temp.water[i] * hall_scaled$temp.water[t] +
            post$b_logQ[i] * hall_scaled$log_Q[t] +
            post$phi[i] * post_preds_ER[i, t-1] -
            post$phi[i] * (post$b_light[i] * hall_scaled$light[t-1] +
                               post$b_temp.water[i] * hall_scaled$temp.water[t-1] +
                               post$b_logQ[i] * hall_scaled$log_Q[t-1]) #+ rnorm(1, sd = post$sigma[i])
    }
    post_preds_err_ER[i,] <- rnorm(nrow(hall_scaled), post_preds_ER[i,], post$sigma[i])

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



plot(hindcast_ER$date, hindcast_ER$ER, type = 'l', ylim = c(-20,0))
for(s in sample(draws, 100)){
    lines(hindcast_ER$date, -exp(post_preds_ER[s,]) + epsilon,alpha = 0.1)
}

lines(hindcast_ER$date, hindcast_ER$ER)
points(hall$date, hall$ER, col = 2, pch = 20)


ggplot(hindcast_ER, aes(date, ER), col = 2) +
    geom_ribbon(aes(ymin = ER_err_low, ymax = ER_err_high),
                col = NA, fill = 'grey80') +
    geom_ribbon(aes(ymin = ER_low, ymax = ER_high),
                col = NA, fill = 'grey50') +
    geom_line() +
    geom_point(data = hall,  col = 2)

hindcast <- full_join(hindcast_GPP, hindcast_ER, by = 'date')

hindcast <- hall %>%
    select(date, GPP_measured = GPP, ER_measured = ER) %>%
    group_by(date) %>%
    summarize(across(c('GPP_measured', 'ER_measured'), mean)) %>%
    right_join(hindcast, by = 'date')

bind_rows(GPP_pars, ER_pars) %>%
    write_csv('data/BRMS_hindcast_models_coefficients.csv')

################################################################################
# Make Figures:

# tiff('figures/BRMS_hindcast_comparison_daily.tiff', width = 7.5, height = 4,
#      units = 'in', res = 300)
png('figures/BRMS_hindcast_comparison_daily.png', width = 7.5, height = 4,
     units = 'in', res = 300)

ggplot(hindcast) +
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
    geom_hline(yintercept = 0) +
    ylab(expression(paste('Metabolism (g ', O[2], m^-2, d^-1, ')'))) +
    xlab('Date')+
    theme_bw() +
    # Custom legend
    scale_color_manual(
        values = c("AR1 Hindcast" = "black", "Historical Data" = "brown3")
    ) +
    guides(color = guide_legend(override.aes = list(shape = c(NA, 16),
                                                    linetype = c(1, 0)),
                                ncol = 2)) +
    theme(
        legend.background = element_rect(fill = 'transparent'),
        legend.position = c(0.02, 0.13), # Upper-left inside the plot
        legend.justification = c(0, 1),
        legend.title = element_blank(),
    )
dev.off()


# compare the predictions averaged across months like in Hall 1972
hall_pred_sum <- hindcast %>%
    mutate(month = month(date)) %>%
    group_by(month) %>%
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
ER_pred <- predict(bmod_ER, newdata = hall_preds)
hall_preds <- cbind(hall_preds, ER_pred) %>%
    rename(ER = Estimate,
           ER_err = Est.Error,
           ER_Q2.5 = Q2.5,
           ER_Q97.5 = Q97.5)

ggplot(hall_preds, aes(date, GPP)) +
    geom_point() +
    geom_point(aes(y = ER), col = 2)

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
    geom_line(aes(y = fitted_values), color = "red", size = 1) +  # Fitted values from the model
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
