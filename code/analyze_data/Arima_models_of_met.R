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
hall_QT <- read_csv('data/hall/hall_discharge_temp_daily.csv')



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
epsilon = 1e-4

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
                           xreg = matrix(c(P_scaled$temp.water, P_scaled$light,
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
        geom_line(data = hall_scaled, aes(x = date, y = ER, color = "Arima Hindcast")) +
        geom_line(data = hall_scaled, aes(x = date, y = GPP, color = "Arima Hindcast")) +
        geom_point(data = hall, aes(x = date, y = ER, color = "Historical Data")) +
        geom_point(data = hall, aes(x = date, y = GPP, color = "Historical Data")) +
        geom_hline(yintercept = 0) +
        ylab(expression(paste('Metabolism (g ', O[2], m^-2, d^-1, ')'))) +
        xlab('Date')+
        theme_bw() +
        # Custom legend
        scale_color_manual(
            values = c("Arima Hindcast" = "black", "Historical Data" = "brown3")
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
                    ymax = GPP_mean + GPP_sd, fill = "Arima Hindcast"),
                col = NA, fill = 'grey') +
    geom_ribbon(aes(ymin = ER_mean - ER_sd,
                    ymax = ER_mean + ER_sd, fill = "Arima Hindcast"),
                col = NA, fill = 'grey') +

    # Lines for GPP and ER predictions (black, Arima Hindcast)
    geom_line(aes(y = GPP_mean, color = "Arima Hindcast")) +
    geom_line(aes(y = ER_mean, color = "Arima Hindcast")) +

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
        values = c("Arima Hindcast" = "black", "Historical Data" = "brown3"),
        labels = c("Arima Hindcast", "Historical Data")
    ) +
    scale_fill_manual(
        name = "Legend",
        values = c("Arima Hindcast" = "grey"),
        labels = c("Arima Hindcast")
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




# try with brms (also old code):
# regenrate the prediction data frame:

hall_scaled <- hall_preds %>%
    mutate(log_Q = (log_Q - scaling_pars$log_Q_mean)/scaling_pars$log_Q_sd,
           light = (light - scaling_pars$light_mean)/scaling_pars$light_sd,
           temp.water = (temp.water - scaling_pars$temp.water_mean)/scaling_pars$temp.water_sd)


# AR1 model with site as a random effect
# determine the parameters for a gamma prior on GPP
a = mean(P_scaled$GPP, na.rm = T)^2/sd(P_scaled$GPP, na.rm = T)
b = mean(P_scaled$GPP, na.rm = T)/sd(P_scaled$GPP, na.rm = T)

bform_GPP <- bf(GPP | mi() ~ ar(p = 1, cov = TRUE) + (1|site) + temp.water + light)
get_prior(bform_GPP, data = P_scaled, family = Gamma(link = 'log'))
GPP_priors <- c(prior("normal(0,5)", class = "b"),
                prior("gamma(1, 1)", class = "shape"))

bmod_GPP <- brm(bform_GPP,
                data = P_scaled,
                family = Gamma(link = "log"),
                chains = 4, cores = 4, iter = 1000,
                prior = GPP_priors,
                control = list(adapt_delta = 0.99) )

summary(bmod_GPP)


GPP_pred <- predict(bmod_GPP, newdata = hall_scaled)
hall_scaled <- cbind(hall_scaled, GPP_pred) %>%
    rename(GPP = Estimate,
           GPP_err = Est.Error,
           GPP_Q2.5 = Q2.5,
           GPP_Q97.5 = Q97.5)
ggplot(hall_preds, aes(date, GPP)) +
    geom_point() + geom_line() +
    geom_ribbon(aes(ymin = GPP_Q2.5,
                    ymax = GPP_Q97.5), fill = 'grey', col = NA)+
    geom_point(data = hall)

bmod_ER <- brm(ER ~ ar(p = 1, gr = site) + (1|site) +
                log_Q * temp.water + GPP,
               data = P_scaled,
               family = Gamma(link = "log"),
               prior = set_prior("normal(0,5)", class = "b"),
               chains = 4, cores = 4, iter = 1000)

summary(bmod_ER)

ER_pred <- predict(bmod_ER, newdata = hall_preds)
hall_preds <- cbind(hall_preds, ER_pred) %>%
    rename(ER = Estimate,
           ER_err = Est.Error,
           ER_Q2.5 = Q2.5,
           ER_Q97.5 = Q97.5)

ggplot(hall_preds, aes(date, GPP)) +
    geom_point() +
    geom_point(aes(y = ER), col = 2)

