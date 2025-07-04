# Compile physical variables from NHC site years
library(streamMetabolizer)
dat <- read_csv("data/metabolism/compiled/metabolism_and_drivers.csv")
dat %>%
    mutate(v_mh = avg_velocity * 60 * 60,
           k600 = K600/depth,
           kO2 = streamMetabolizer::convert_k600_to_kGAS(k600, temp.water,
                                                         gas = "O2"),
           kO2_h = kO2/24,
           reach_length_95 = 3 * v_mh/kO2_h/1000,
           year = case_when(site == "NHC" & year == 2020 ~ 2021,
                            TRUE ~ year)) %>%
    filter(site %in% c("NHC", "CBP"), year != 2021) %>%
    group_by(site) %>%
    summarize(med = median(reach_length_95, na.rm = T),
              mean = mean(reach_length_95, na.rm = T),
              max = max(reach_length_95, na.rm = T),
              min = min(reach_length_95, na.rm = T),
              quant = quantile(reach_length_95, probs = 0.97, na.rm = T))


    select(site, date, avg_velocity, v_mh, temp.water, K600, kO2_h, depth, reach_length_95) %>%
    arrange(desc(reach_length_95))
    ggplot(aes(date, depth)) +
    geom_line() +
    facet_wrap(.~site)

# dat <- readRDS("data/metabolism/compiled/met_preds_stream_metabolizer.rds")
preds <- dat$preds %>%
  filter(era == "now",
         site != "PWC",
         year != 2020)

tmp <- preds %>%
  filter(!is.na(DO.obs)) %>%
  summarize(DO_mean = mean(DO.obs),
            DO_sat_mean = mean(DO.obs/DO.sat),
            n_hypoxic = sum(DO.obs <=2),
            per_hypoxic = n_hypoxic/n()) %>%
  mutate(site = "all",
         time = "year")
tmp <- preds %>%
  filter(!is.na(DO.obs),
         month %in% 10) %>%
  summarize(DO_mean = mean(DO.obs),
            DO_sat_mean = mean(DO.obs/DO.sat),
            n_hypoxic = sum(DO.obs <=2),
            per_hypoxic = n_hypoxic/n()) %>%
  mutate(site = "all",
         time = "Oct") %>%
  bind_rows(tmp)
tmp <- preds %>%
  filter(!is.na(DO.obs)) %>%
  group_by(site, year) %>%
  summarize(DO_mean = mean(DO.obs),
            DO_sat_mean = mean(DO.obs/DO.sat),
            n_hypoxic = sum(DO.obs <=2),
            per_hypoxic = n_hypoxic/n()) %>%
  mutate(time = "year") %>%
  bind_rows(tmp)
DO <- preds %>%
  filter(!is.na(DO.obs),
         month %in% 10) %>%
  group_by(site, year) %>%
  summarize(DO_mean = mean(DO.obs),
            DO_sat_mean = mean(DO.obs/DO.sat),
            n_hypoxic = sum(DO.obs <=2),
            per_hypoxic = n_hypoxic/n()) %>%
  mutate(time = "Oct") %>%
  bind_rows(tmp)

# filter(DO, site %in% c('NHC','UNHC')) %>%
ggplot(DO,aes(year, DO_sat_mean, fill = time)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  facet_wrap(~site)
filter(DO, year == 2019) %>%
ggplot(aes(site, DO_sat_mean, fill = time)) +
  geom_bar(stat = 'identity', position = 'dodge')

write_csv(DO, 'data/metabolism/driver_variables/DO_summary_siteyears.csv')
DO %>%
  filter(time == "year", !is.na(year)) %>%
  summary()
DO %>%
  filter(time == "Oct", !is.na(year)) %>%
  summary()

preds %>%
  select(DO.obs, DO.sat) %>%
  mutate(DO = DO.obs/DO.sat)%>%
  summary()


predmeans <- preds %>%
  mutate(DO.psat = DO.obs/DO.sat) %>%
  select(site, year, GPP, ER, K600,
         discharge, temp.water, DO.obs, DO.psat) %>%
  group_by(site, year) %>%
  summarize(across(everything(), .fns = list(mean = ~mean(., na.rm = T),
                                             min = ~min(., na.rm = T),
                                             max = ~max(., na.rm = T),
                                             cum = ~sum(., na.rm = T) * 365/sum(!is.na(.)) * 60*60*24*365,
                                             cv = ~sd(., na.rm = T)/mean(., na.rm = T)),
                   .names = '{col}_{fn}'))
write_csv(predmeans, 'NHC_2019_metabolism/data/metabolism/driver_variables/physical_vars_summ_siteyears.csv')

summary(predmeans)

# fix mean annual discharge at UNHC
# scale <- predmeans %>%
#   filter(year != 2019) %>%
#   select(discharge_mean, site, year) %>%
#   pivot_wider(names_from = site, values_from = discharge_mean) %>%
#   mutate(scale = NHC/UNHC) %>%
#   summarize(scale = mean(scale))
# scale = scale$scale[1]
# nhcq <- predmeans$discharge_mean[predmeans$site == 'NHC' & predmeans$year == 2019]
# unhcq <- nhcq/scale
# predmeans$discharge_mean[predmeans$site == 'UNHC' & predmeans$year == 2019] <- unhcq
# predmeans$discharge_mean[predmeans$site == 'WB'] <- nhcq/(scale * 58.02/61.64)
# predmeans$discharge_mean[predmeans$site == 'WBP'] <- nhcq/(scale * 58.02/61.38)

png("figures/drivers_Qtemp_across_sites.png", res = 300, units = 'in',
    width = 5, height = 3)
predmeans %>%
  filter(year == 2019) %>%
  select(ends_with("mean"), -starts_with("DO"), -K600_mean) %>%
  mutate(ER_mean = -ER_mean,
         discharge_mean = log(discharge_mean)) %>%
  pivot_longer(cols = c(-site, -GPP_mean, -ER_mean), names_to = "driver",
               values_to = "driver_mean",
               names_pattern = '([a-zA-Z0-9\\.]+)_[a-zA-Z]') %>%
  pivot_longer(c(GPP_mean, ER_mean), names_to = "met", values_to = "value",
               names_prefix = '_mean') %>%
  ggplot(aes(driver_mean, value, col = met)) +
  geom_point() +
  geom_smooth(method = lm) +
  facet_wrap(.~driver, scales = "free") +
  ggtitle("Across Sites in 2019")
dev.off()

png("figures/drivers_Qtemp_across_sites.png", res = 300, units = 'in',
    width = 5, height = 3)
predmeans %>%
  filter(site %in% c('NHC', 'UNHC')) %>%
  select(ends_with("mean"), -starts_with("DO"), -K600_mean) %>%
  mutate(ER_mean = -ER_mean,
         discharge_mean = log(discharge_mean)) %>%
  pivot_longer(cols = c(-site, -GPP_mean, -ER_mean), names_to = "driver",
               values_to = "driver_mean",
               names_pattern = '([a-zA-Z0-9\\.]+)_[a-zA-Z]') %>%
  pivot_longer(c(GPP_mean, ER_mean), names_to = "met", values_to = "value",
               names_prefix = '_mean') %>%
  ggplot(aes(driver_mean, value, col = site)) +
  geom_point() +
  geom_smooth(method = lm) +
  facet_wrap(met~driver, scales = "free") +
  ggtitle("Across years")
dev.off()

predmeans %>%
  filter(site %in% c('NHC', 'UNHC')) %>%
  select(site, ends_with("mean"), -starts_with("DO"), -K600_mean) %>%
  mutate(ER_mean = -ER_mean,
         discharge_mean = log(discharge_mean)) %>%
  ggplot(aes(discharge_mean, temp.water_mean, color = site)) +
  geom_point()


Q_stats <- read_csv('data/rating_curves/annual_Q_stats.csv')
predmeans <- left_join(predmeans, Q_stats[,1:4], by = 'year')
yy <- predmeans %>% filter(site %in% c('NHC', 'UNHC'))
ss <- predmeans %>% filter(year == 2019)

ggplot(yy, aes(discharge_mean, ER_mean, col = site)) +
  geom_point() +
  geom_smooth(method = lm)

ggplot(predmeans, aes(temp.water_mean, GPP_mean)) +
  geom_point() +
  geom_smooth(method = lm)

summary(lm(ER_mean~temp.water_mean, data = predmeans))
