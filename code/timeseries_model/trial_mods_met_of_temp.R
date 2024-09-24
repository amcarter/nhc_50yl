# Basic hierarchical model of metabolism as a function of temperature and Q

# A carter
# 2/2022
library(tidyverse)
library(nlme)
library(mgcv)
library(brms)

# load in the modern dataset:
dat <- read_csv("data/metabolism/compiled/metabolism_and_drivers.csv")
sites <- read_csv('data/siteData/NHCsite_metadata.csv') %>%
  slice(c(1:5, 7))

# load in Hall dataset:

hall <- readRDS("data/metabolism/compiled/met_preds_stream_metabolizer_O2.rds")
hall <- hall$preds %>% filter(era == 'then', site == 'CBP')
hall_QT <- read_csv('data/hall_data/hall_discharge_temp_daily.csv')


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


# fit using the nlme package
lme_mod <- lme(fixed = log(GPP) ~ temp.water + light,
              random = ~1|site,
              correlation = corAR1(form = ~1|site),
              data = P_scaled,
              na.action = na.omit)

summary(lme_mod)

# Predict the fitted values
GPP_preds <- predict(lme_mod, newdata = hall_scaled)
hall_scaled$GPP_pred <- GPP_preds
hall_scaled$GPP = exp(hall_scaled$GPP_pred) - epsilon

# plot the predictions
ggplot(hall_scaled, aes(date, GPP))+
    geom_point(col = 'steelblue') + geom_line(col = 'steelblue') +
    geom_point(data = hall, col = 2)


# model for ER:
lme_mod_ER <- lme(fixed = log(ER) ~ temp.water * log_Q + GPP_pred,
                  random = ~1|site,
                  correlation = corAR1(form = ~1|site),
                  data = P_scaled,
                  na.action = na.omit)

summary(lme_mod_ER)
ER_preds <- predict(lme_mod_ER, newdata = hall_scaled)
hall_scaled$ER_pred <- ER_preds
hall_scaled$ER <- -(exp(hall_scaled$ER_pred) - epsilon)

ggplot(hall_scaled, aes(date, ER)) +
    geom_point() + geom_line() +
    geom_point(data = hall, col = 2)


# model using GAMS
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
    # geom_line(aes(y = temp.water/10))
    geom_point(data = hall, col = 2)




# try with brms:
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
site) +
            discharge + temp.water + discharge*temp.water + slope,
            data = P)

m3 <- brm(ER ~ ar(p = 1, gr = site) +
            light + discharge*temp.water + slope,
            data = P)


G1 <- brm(GPP ~ ar(p = 1, gr = site) +
            light + discharge*temp.water + slope,
            data = P)

summary(m2)

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
