# Plot summarized metabolism data for Hall comparison
dat <- read_csv("data/metabolism/compiled/metabolism_summary_table_gC.csv") %>%
  filter(method == 'uninformed_raymond'| site== "CBP" & is.na(year)) %>%
  mutate(year = ifelse(method == 'uninformed_raymond', year, 1969)) %>%
  select(-c(3:11))

# seasonal distribution of peak metabolism ####

sum <- dat %>% filter(year != 1969, year != 2020, site != 'PWC') %>%
  mutate(across(starts_with('peak'), ~as.numeric(format(., '%j'))+4,
                .names = '{col}_doy'))

png('figures/distribution_peak_gpp_er.png', height = 4, width = 5, res = 300,
    units = 'in', family = 'cairo')
  plot(density(sum$peak_gpp_doy, na.rm = T), xlim = c(0, 365),
       col = 'forestgreen', lwd = 2, xaxt = 'n', xlab = 'Date',
       main = 'Distribution of peak ER and GPP (n = 10)')
  lines(density(sum$peak_er_doy), col = 'sienna', lwd = 2)
  axis(1, at = as.numeric(format(seq(as.Date('2019-01-01'), length.out = 12,
                                    by = 'month'), '%j')), labels = month.abb)
  polygon(c(62, 71, 71, 62), c(0,0,0.02,0.02), col = alpha('sienna', .4),
          border = NA)
  polygon(c(68, 77, 77, 68), c(0,0,0.02,0.02),
          col = alpha('forestgreen', .5), border = NA)
  legend('topright', c('GPP 2017-19', 'ER 2017-19', 'GPP 1969', 'ER 1969'),
         col = c('forestgreen', 'sienna', NA, NA),
         fill = c(NA, NA, alpha('forestgreen', .5), alpha('sienna', .4)),
         lwd = 2, lty = c(1,1,NA, NA), border = NA, bty = 'n',
         x.intersp = c(1,1,0, 0), seg.len = 1.2, inset = .02)
dev.off()

sum %>%
  summarize(across(ends_with('cum'), ~sd(.)/mean(.)))
# add geomorphic variables: ####
geo <- read_csv("data/reach_characterization/nhc_habitat_dimensions_by_reach.csv") %>%
  filter(habitat == "total") %>%
  select(-habitat, -total_length_m)
geo <- sites %>%
  select(site = sitecode, width_mar_m, habitat, slope, distance_m) %>%
  left_join(geo, by = c('site', 'distance_m')) %>%
  right_join(dat)

geosum <- geo %>%
  select(site, width_mar_m, avg_depth_mar, slope, distance_m) %>%
  group_by(site) %>%
  summarize_all(mean, na.rm = T)

summary(geosum)
met_dat <- geo %>%
  select(-starts_with("peak"), -days, -pctcoverage, -daterange, -era) %>%
  pivot_longer(cols = starts_with(c('gpp', 'er')),
               names_to = c('met', "variable"),
               values_to = "value",
               names_sep = '_')  %>%
  filter(year != 1969, year!= 2020, site != 'PWC',
         variable != "min")
  # mutate(gC_m2_time = ifelse(met == 'er', -gC_m2_time, gC_m2_time))
# met_dat %>%
#   filter(year == 2019) %>%
# ggplot( aes(slope, value, col = met)) +
#   # geom_bar(stat = 'identity', position = 'dodge') +
#   geom_line() +
#   facet_wrap(.~variable, scales = 'free_y')

png('figures/drivers_annual_met_by_width_depth_slope.png', width = 5, height = 4,
    res = 300, units = 'in')
  met_dat %>%
    filter(year == 2019,
           variable %in% c('cum', 'median')) %>%
    pivot_longer(cols = any_of(c('avg_width', 'avg_depth_mar', 'slope')),
                 names_to = 'geo_measure',
                 values_to = 'geo_value') %>%
  ggplot( aes(geo_value, value, color = met)) +
    geom_point() +
    geom_smooth(method = lm) +
    facet_grid(variable~geo_measure, scales = 'free')+
    xlab('meters') +
    ylab('gC/m2/time')
dev.off()

# summarize model fits: ####
ss <- met_dat %>%
  filter(year == 2019,
         variable == "mean") %>%
  slice(-c(5,6))%>%
  pivot_wider(names_from = met, values_from = value)


  # ggplot(aes(habitat, value,  color = met) )+
  # geom_boxplot()
summary(lm(gpp~slope, data = ss))
summary(lm(er~slope, data = ss))
summary(lm(gpp~avg_width, data = ss))
summary(lm(er~avg_width, data = ss))
summary(lm(gpp~avg_depth_mar, data = ss))
summary(lm(er~avg_depth_mar, data = ss))
summary(lm(gpp~avg_depth_oct, data = ss))
summary(lm(er~avg_depth_oct, data = ss))

met <- readRDS('data/metabolism/compiled/met_preds_stream_metabolizer.rds')$preds %>%
  filter(era == 'now') %>%
  select(site, year, GPP, ER, discharge, temp.water, DO.obs) %>%
  group_by(site, year) %>%
  summarize(across(everything(), mean, na.rm = T))

summary(lm(ER ~ temp.water, data = met))
summary(lm(GPP ~ temp.water, data = met))
summary(lm(ER+GPP ~ DO.obs, data = met))
summary(lm(GPP ~ DO.obs, data = met))

### Doesn't work below here
# met_dat %>%
#   filter(site %in% c('NHC', 'UNHC')) %>%
# ggplot( aes(site, value, fill = factor(year))) +
#   geom_bar(stat = 'identity', position = 'dodge') +
#   facet_wrap(.~variable, scales = 'free_y')

slist <- sites$sitecode[2:6]
gpp_mean = ss$gpp_mean[ss$site == 'NHC']
er_mean = ss$er_mean[ss$site == 'NHC']
pd <- data.frame()
for(s in slist){
  gpp = abs((ss$gpp_mean[ss$site == s] - gpp_mean)/
              (ss$gpp_mean[ss$site == s] + gpp_mean))
  er = abs((ss$er_mean[ss$site == s] - er_mean)/
              (ss$er_mean[ss$site == s] + er_mean))
  pd = bind_rows(pd, data.frame(gpp = gpp, er = er))
}

all_pd <- pd %>%
  summarize_all(.funs = list(mean = ~mean(.*2),
                             sd = ~sd(.*2))) %>%
  mutate(scale = "sites")
nhc = filter(yy, site == "NHC")
unhc = filter(yy, site == "UNHC")
pd <- data.frame()
for(y in 2018:2019){
  gpp_mean = nhc$gpp_mean[nhc$year == 2017]
  gpp = abs((nhc$gpp_mean[nhc$year == y] - gpp_mean)/
              (nhc$gpp_mean[nhc$year == y] + gpp_mean))
  er_mean = nhc$er_mean[nhc$year == 2017]
  er = abs((nhc$er_mean[nhc$year == y] - er_mean)/
              (nhc$er_mean[nhc$year == y] + er_mean))
  pd = bind_rows(pd, data.frame(gpp = gpp, er = er))

  gpp_mean = unhc$gpp_mean[unhc$year == 2017]
  gpp = abs((unhc$gpp_mean[unhc$year == y] - gpp_mean)/
              (unhc$gpp_mean[unhc$year == y] + gpp_mean))
  er_mean = unhc$er_mean[unhc$year == 2017]
  er = abs((unhc$er_mean[unhc$year == y] - er_mean)/
              (unhc$er_mean[unhc$year == y] + er_mean))
  pd = bind_rows(pd, data.frame(gpp = gpp, er = er))
}

all_pd <- dd %>%
  summarize(gpp_mean = abs(diff(gpp_mean)/sum(gpp_mean))*2,
            er_mean = abs(diff(er_mean)/sum(er_mean))*2) %>%
  mutate(scale = "decades") %>%
  bind_rows(all_pd)
all_pd <- pd %>%
  summarize_all(.funs = list(mean = ~mean(.*2),
                             sd = ~sd(.*2))) %>%
  mutate(scale = "years") %>%
  bind_rows(all_pd)

png('figures/percent_difference_in_cumulative_met_across_scales.png',
    width = 5, height = 4, units = 'in', res = 300)
all_pd %>%
  pivot_longer(cols = -scale, names_to = c('met', 'meas'),
               names_pattern = '([a-z]+)_([a-z]+)',
               values_to = 'percent_difference') %>%
  pivot_wider(names_from = meas, values_from = percent_difference) %>%
ggplot(aes(scale, mean, fill = met)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
             stat = 'identity', position = position_dodge(.9), width = .2) +
  ylab("percent difference between site years")
dev.off()

