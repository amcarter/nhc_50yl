# plot data for NHC watershed
# precip ####

library(tidyverse)
library(zoo)
library(lubridate)

nldas <- read_csv("data/watershed/nldas.csv") %>%
  dplyr::select(datetime = DateTime, value = '1', variable) %>%
  filter(!variable %in% c("wind_speed",
                          "max_relative_humidity",
                          "min_relative_humidity",
                          "surface_downwelling_shortwave_flux_in_air"))
nldas2 <- read_csv("data/watershed/nldas2.csv")

pp_annual = read_csv('data/watershed/prism_raw.csv') %>%
  mutate(year = year(DateTime)) %>%
  dplyr::select(year, ppt_mm = '1') %>%
  filter(year < 1979) %>%
  group_by(year) %>%
  summarize(cumulative_precip = sum(ppt_mm)/1000)

all <- nldas %>%
pivot_wider(names_from = "variable", values_from = "value") %>%
select(datetime, precip_mmd = precipitation_amount)

p90 <- quantile(all$precip_mmd, .9, na.rm = T)
# pp <- all %>%
  # group_by(year = year(datetime)) %>%
pp <- nldas2 %>%
  group_by(year = year(date)) %>%
  summarize(zero_days = length(which(precip_mmd == 0)),
            cumulative_precip = sum(precip_mmd)/1000,
            cum90 = sum(precip_mmd[which(precip_mmd >= p90)]),
            max_precip = max(precip_mmd, na.rm = T),
            percent_extreme = cum90/cumulative_precip/10) %>%
  dplyr::select(-cum90, -max_precip) %>%
  bind_rows(pp_annual)
#   pivot_longer(cols= -year, names_to = "variable", values_to = "value")
# write_csv()

# calculate precip during 2019 drought ####
p19 <- all %>%
  filter(as.numeric(substr(datetime, 1, 4)) == 2019,
         as.numeric(substr(datetime, 6, 7)) %in% 9:10) %>%
  slice(c(6:37)) %>%
  mutate(pre_cum = cumsum(precip_mmd))
plot(p19$datetime, p19$pre_cum, type = "l")
p19


# NLCD Land Cover (stats) [OBSOLETE] ####

dat <- read_csv("data/watershed/nlcd_1992-2016_summary.csv") %>%
  filter(!is.na(category)) %>%
  mutate(category = case_when(id %in% c(22, 23, 24) ~ 'developed',
                              id %in% c(42, 43, 41, 90) ~ 'forested',
                              id %in% c(81, 82) ~ 'agriculture',
                              id %in% c(21, 52, 71) ~ 'grass_shrub',
                              # id %in% c(11) ~ 'water',
                              TRUE ~ 'other')) %>%
  select(-id) %>%
  group_by(category) %>%
  summarize_all(sum, na.rm = T) %>%
  ungroup()

ncells <- sum(dat$CellTally1992, na.rm = T)
nlcd <- dat %>%
  filter(category != 'other') %>%
  mutate(across(starts_with("Cell"), ~ . / ncells)) %>%
  pivot_longer(cols = starts_with('Cell'), names_to = 'year',
               values_to = 'percent') %>%
  mutate(year = as.numeric(substr(year, 10, 13)),
         percent = ifelse(is.na(percent), 0, percent),
         category = factor(category, levels = c('agriculture', 'developed',
                                                'grass_shrub', 'forested')))
nn <- nlcd %>%
  pivot_wider(names_from = category, values_from = percent)

# Riparian Zone ####
NPP <- read_csv("data/watershed/gee_files/npp.csv") %>%
  select(year = date, NPP_mean = annualNPP_mean, NPP_std = annualNPP_std) %>%
  filter(year != 2020)

dat <- read_csv("data/watershed/gee_files/lai.csv") %>%
  filter(!is.na(date)) %>%
  mutate(date = as.Date(date, format = "%Y_%m_%d"),
         year = year(date)) %>%
  select(date, year,lai = Lai_500m_mean)
dates <- data.frame(date = seq(dat$date[1], dat$date[nrow(dat)], by = 'day'))
max <- dat %>%
  group_by(year) %>%
  summarize(lai_max = quantile(lai, 0.975, na.rm = T)) %>%
  ungroup()
dat <- dates %>%
  left_join(dat, by = 'date') %>%
  mutate(year = year(date),
         doy = as.numeric(format(date, '%j')),
         lai = na.approx(lai)) %>%
  left_join(max, by = "year") %>%
  mutate(lai = lai/lai_max) %>%
  as_tibble()

lai <- data.frame()
for(y in unique(dat$year)){
  tmp <- dat %>%
    filter(year == y) %>%
    arrange(date)
  w <- which(tmp$lai >= .5)
  don <- tmp$doy[w[1]]
  doff <- tmp$doy[w[length(w)]]
  tmp <- data.frame(year = y,
                    leaf_on = don,
                    leaf_off = doff)
  lai <- bind_rows(lai, tmp)
}

lai  <- left_join(lai, max, by = "year") %>%
  mutate(leaf_on_days = leaf_off - leaf_on)


riparian <- left_join(NPP, lai, by = "year")

# air temperature ####
air <- read_csv( "data/watershed/noaa_air_temp.csv")

air <-air %>%
  mutate(year = year(date)) %>%
  group_by(year) %>%
  summarize(temp_mean = mean(temp_mean, na.rm = T),
            min_mean = mean(temp_min, na.rm = T)) %>%
  filter(year !=1967)

# Hydrology ####
hydro <- read_csv('data/watershed/nhcblands_usgs_stats.csv')  %>%
  filter(year > 1972)
# join dataframes ####
alldat <- full_join(riparian, nn, by = 'year') %>%
  filter(year !=2020)

ws <- alldat %>%
  select(year, agriculture, developed, forested, grass_shrub,
         lai_max, NPP_mean, leaf_on_days, leaf_on, leaf_off)

cc <- full_join(pp , air, by = 'year')%>%
  full_join(hydro, by = 'year') %>%
  filter(#year !=2020,
         year >1967) %>%
  select(-n,  -peak_Q, -ar_1) %>%
  arrange(year)
summary(lm(cumulative~year, data = cc))

prism_monthly_full <- read_csv('data/watershed/prism_annual.csv')
cc <- cc %>%
    rename(cp = cumulative_precip, tm = temp_mean) %>%
    full_join(prism_monthly_full, by = 'year') %>%
    arrange(year) %>%
    select(-cp, -tm) %>%
    mutate(cumulative_precip = cumulative_precip / 1000) %>%
    filter(between(year, 1967, 2021))

# png("figures/precip.png", width = 7.5, height = 6, units = "in", res = 300)
# cc %>%
#   pivot_longer(-year, names_to = 'variable', values_to = 'value') %>%
# ggplot(aes(x = year, y = value)) +
#   geom_line(lwd = 2) + geom_point() +
#   facet_wrap(.~variable, scales = "free_y", dir = "v", switch = "y") +
#   geom_smooth(method = lm, lwd = 1, col = "black") +
#   xlim(1986,2019)+
#   theme_bw()
# dev.off()



# make multipanel figure (v1, obsolete) ####

dir.create('figures', showWarnings = FALSE)
tiff("figures/Climate_Watershed_Multipanel_figure.tif", width = 6, height = 5.4,
    units = 'in', res = 800, compression = 'lzw')
precip_col ='steelblue'
dat_col = 'grey25'
prod_col = 'forestgreen'

m <- cbind(c(1,1,1,1,2,2,3,3,4,4), c(rep(5,6),6,6,7,7))
layout(m)
# layout.show(8)
par(mar = c(1,2,0,1.5),
    oma = c(3.5, 2, 4, 0),
    adj = 0, ps = 10)
plot(cc$year, cc$temp_mean, type = 'l', lwd = 1.2, xaxt = 'n', yaxt = 'n',
     col = dat_col, ylim = c(8, 17),bty = 'n')
mtext("Climate Change", 3, 1, cex = 0.9)
mm <- lm(temp_mean~year, data = cc)
conf_interval <- data.frame(predict(mm, interval="confidence", level = 0.95))
lines(1968:2019, conf_interval$fit, col = dat_col)
polygon(c(1968:2019, 2019:1968), c(conf_interval$lwr, rev(conf_interval$upr)),
        col = alpha(dat_col, .3), border = NA)
axis(2)
lines(cc$year, cc$min_mean, lwd = 1.2, col = dat_col)
mm <- lm(min_mean~year, data = cc)
conf_interval <- data.frame(predict(mm, interval="confidence", level = 0.95))
lines(1968:2019, conf_interval$fit, col = dat_col)
polygon(c(1968:2019, 2019:1968), c(conf_interval$lwr, rev(conf_interval$upr)),
        col = alpha(dat_col, .3), border = NA)
mtext('Air Temp C', 2, 2.1, cex = .8)
text(1982, 14.2,'Daily mean, slope = 0.4 C/decade', col = dat_col, cex = 1)
text(1982, 8.5,'Daily min, slope = 0.45 C/decade', col = dat_col, cex = 1)
rect(1967, 7.5, 2020, 17.5, xpd = NA)

plot(cc$year, cc$cumulative_precip, type = 'l', lwd = 1.2, col = precip_col, axes = F)
axis(2, at = c(0.9, 1.3, 1.7))
mtext('Annual Precip (m)', 2, 2.1, cex = .8)
plot(cc$year, cc$percent_extreme, type = 'l', lwd = 1.2, col = precip_col, axes = F)
axis(2, at = c(60, 70, 80), labels = c('60%', '70%', '80%'))
mtext('% Extreme', 2, 2.1, cex = .8)
mm <- lm(percent_extreme~year, data = cc)
conf_interval <- data.frame(predict(mm, interval="confidence", level = 0.95))
lines(1979:2019, conf_interval$fit, col = precip_col)
polygon(c(1979:2019, 2019:1979), c(conf_interval$lwr, rev(conf_interval$upr)),
        col = alpha(precip_col, .3), border = NA)
text(1970, 75, 'slope = 2.9%/decade', col = precip_col, cex = 1)
plot(cc$year, cc$zero_days, type = 'l', lwd = 1.2, col = precip_col, axes = F)
axis(2)#, at = c(60, 70, 80), labels = c('60%', '70%', '80%'))
mtext('No Precip Days', 2, 2.1, cex = .8)
mm <- lm(zero_days~year, data = cc)
conf_interval <- data.frame(predict(mm, interval="confidence", level = 0.95))
lines(1979:2019, conf_interval$fit, col = precip_col)
polygon(c(1979:2019, 2019:1979), c(conf_interval$lwr, rev(conf_interval$upr)),
        col = alpha(precip_col, .3), border = NA)
text(1970, 220, 'slope = 12 days/decade', col = precip_col, cex = 1)
rect(1967, 152.5,2020,471,xpd = NA)

# plot(cc$year, cc$cumulative/10^7, type = 'l', lwd = 1.2, col = dat_col, xaxt = 'n',
#      yaxt = 'n', bty = 'n', ylim = c(-5, 12))
# axis(2, at = seq(3, 12, by = 3))
# mtext('Eno River Discharge', 2, 2.1, cex = .8)
# par(new = T)
# plot(cc$year, cc$q05, col = dat_col, type = 'l', lwd = 1.2, axes = F, ylim = c(0, 2.2))
# axis(2, at = c(0, 0.5))
# text(1969, 1.8, 'Cumulative Annual \n(10^7 m3)', col = dat_col, cex = 1.2)
# text(1969, 0.25, '5th quantile \n(m3/s)', col = dat_col, cex = 1.2)
# rect(1967, -0.09, 2020, 2.25, xpd = NA)
axis(1)
# panel 2 Watershed Change
par(mar = c(4,2.5,0,1))
plot(ws$year, ws$forested, ylim = c(0, .90), pch = 19, col = 'forestgreen',
     axt = 'n', yaxt = 'n', bty = 'n', xlim = c(1985.5, 2019), xlab = '')
axis(1)
mtext("Watershed Change", 3, 1, cex = 0.9)
mtext("Full Watershed", 3, 0.1, cex = 0.8)
axis(2, at = seq(0, .8, by = .2), labels = c('0%', '20%', '40%', '60%', '80%'))
text(1987, .73, 'Forested, slope = -2.9%/decade', col = prod_col, cex = 1)
mm = lm(forested~year, data = ws)
conf_interval <- data.frame(predict(mm,interval="confidence", level = 0.95))
years <- c(1992, 2001, 2004, 2006, 2008, 2011, 2014,2016)
lines(years, conf_interval$fit, col = 'forestgreen')
polygon(c(years, rev(years)), c(conf_interval$lwr, rev(conf_interval$upr)),
        col = alpha('forestgreen', .3), border = NA)
points(ws$year, ws$agriculture, pch = 19, col = dat_col)
text(1987, .14, 'Agriculture', col = dat_col, cex = 1)
text(1987, .05, 'Developed, slope = 0.6%/decade', col = 'brown3', cex = 1)
points(ws$year, ws$developed, pch = 19, col = 'brown3')
mm = lm(developed~year, data = ws)
conf_interval <- data.frame(predict(mm,interval="confidence", level = 0.95))
lines(years, conf_interval$fit, col = 'brown3')
polygon(c(years, rev(years)), c(conf_interval$lwr, rev(conf_interval$upr)),
        col = alpha('brown3', .3), border = NA)

mtext('Watershed landcover', 2, 2.1, cex = 0.8)
rect(1985, -.0355, 2020.2, .937, xpd = NA)
par(mar = c(2,2.5,0,1))


plot(ws$year, ws$NPP_mean, type = 'l', lwd = 1.2, col = 'forestgreen',
     xlim = c(1985.5, 2019), axes = F)
mtext("Riparian Area", 3, .2, cex = 0.8)
mtext('NPP kg/m2/y', 2, 2.1, cex = .8)
axis(2)
par(mar = c(1,2.5,0,1))
plot(ws$year, ws$leaf_on, type = 'l', lwd = 1.2, col = prod_col,
     xlim = c(1985.5, 2019), axes = F, ylim = c(92, 170))
axis(2, at = c(92, 123), labels = c('Apr','May'), las = 2)
rect(1985, 89, 2020.2, 270, xpd = NA)
par(new = T)
plot(ws$year, ws$leaf_off,type = 'l', lwd = 1.2, col = 'sienna',
     xlim = c(1985.5, 2019), axes = F, ylim = c(233, 310))

axis(2, at = c(275,306), labels = c('Oct', 'Nov'), las = 2)
# mtext('Leaf Date', 2, 2.1, cex = .8)
mm = lm(leaf_off~year, data = ws)
conf_interval <- data.frame(predict(mm,interval="confidence", level = 0.95))
lines(2000:2019, conf_interval$fit, col = 'sienna')
polygon(c(2000:2019, 2019:2000), c(conf_interval$lwr, rev(conf_interval$upr)),
        col = alpha('sienna', .3), border = NA)

text(1986, 290, 'leaf off date \nslope = -5 days/decade', col = 'sienna', cex = 1)
text(1986, 250, 'leaf on date', col = prod_col, cex = 1)
axis(1)

dev.off()

# make multipanel figure with map (v2) ####


dir.create('figures', showWarnings = FALSE)
png("figures/Climate_Watershed_Multipanel_figure3.png", width = 4, height = 4.2,
    units = 'in', res = 800, type = 'cairo')
# tiff("figures/Climate_Watershed_Multipanel_figure.tif", width = 4, height = 5,
#     units = 'in', res = 800, compression = 'lzw')

precip_col ='steelblue'
dat_col = 'grey25'
prod_col = 'forestgreen'
pt_yrs <- c(1968, 1969, 1970, 2019, 2020)
pt_yrs_nhc <- c(2017, 2018, 2019)


m <- cbind(c(1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 6))
layout(m, heights = c(1,1,1,1,0.1,1,1,1,1,1,1))
library(forecast)
#air
par(mar = c(0,2,0,1),
    oma = c(2, 2, 0.5, 0),
    adj = 0, ps = 10)
plot(cc$year, cc$temp_mean, type = 'l', lwd = 1.2, xaxt = 'n', yaxt = 'n',
     col = dat_col, ylim = c(13.5, 17), bty = 'n')
mm <- lm(temp_mean~year, data = cc)
m_ar1 <- forecast::Arima(cc$temp_mean, order = c(1, 0, 0), xreg = matrix(cc$year, ncol = 1))
# conf_interval <- data.frame(predict(m_ar1, interval="confidence", level = 0.95))
conf_interval <- forecast(m_ar1, xreg = cc$year)
lines(1967:2021, conf_interval$mean, col = dat_col)
polygon(c(1967:2021, 2021:1967), c(conf_interval$lower[,1], rev(conf_interval$upper[,1])),
        col = alpha(dat_col, .3), border = NA)
axis(2, at = 14:17, labels = 14:17, las = 2)
# lines(cc$year, cc$min_mean, lwd = 1.2, col = dat_col)
points(pt_yrs_nhc, cc$temp_mean[cc$year %in% pt_yrs_nhc], pch = 3, col = dat_col, cex = 1.5, xpd = NA)
points(pt_yrs, cc$temp_mean[cc$year %in% pt_yrs], pch = 19, col = dat_col, cex = 1, xpd = NA)
legend(x = 1967, y = 17.5, legend = c("CB", "NHC"), pch = c(19, 3), col = "gray45", bty = 'n', xpd=NA)
# mm <- lm(min_mean~year, data = cc)
# conf_interval <- data.frame(predict(mm, interval="confidence", level = 0.95))
# lines(1968:2019, conf_interval$fit, col = dat_col)
# polygon(c(1968:2019, 2019:1968), c(conf_interval$lwr, rev(conf_interval$upr)),
#         col = alpha(dat_col, .3), border = NA)
mtext('Air Temp C', 2, 2.7, cex = .8)
text(1982, 13.7, 'Daily mean, slope = 0.28 ± 0.06°C/decade', col = dat_col, cex = 1, xpd=NA)

#precip
plot(cc$year, cc$cumulative_precip, type = 'l', lwd = 1.2, col = precip_col,
     axes = F, ylim = c(0.9, 1.7), xpd = NA, ylab = '', xlab = '')
# mm <- lm(cumulative_precip~year, data = cc)
# m_ar1 <- forecast::Arima(cc$cumulative_precip, order = c(1, 0, 0), xreg = matrix(cc$year, ncol = 1))
# conf_interval <- data.frame(predict(mm, interval="confidence", level = 0.95))
points(pt_yrs_nhc, cc$cumulative_precip[cc$year %in% pt_yrs_nhc], pch = 3, col = precip_col, cex = 1.5, xpd = NA)
points(pt_yrs, cc$cumulative_precip[cc$year %in% pt_yrs], pch = 19, col = precip_col, cex = 1, xpd = NA)
axis(2, at = c(0.9, 1.3, 1.7), las = 2)
mtext('Annual Precip (m)', 2, 2.7, cex = .8)

plot(1,1, type = 'n', axes = FALSE, ann = FALSE)

#% extreme
m_ar1 <- forecast::Arima(cc$percent_extreme, order = c(1, 0, 0), xreg = matrix(cc$year, ncol = 1))
plot(cc$year, cc$percent_extreme, type = 'l', lwd = 1.2, col = precip_col, axes = F)
axis(2, at = c(60, 70, 80), labels = c('60%', '70%', '80%'), las = 2)
mtext('% Extreme', 2, 2.7, cex = .8)
mm <- lm(percent_extreme~year, data = cc)
ccc <- filter(cc, ! is.na(percent_extreme))
m_ar1 <- forecast::Arima(ccc$percent_extreme, order = c(1, 0, 0), xreg = matrix(ccc$year, ncol = 1))
# conf_interval <- data.frame(predict(mm, interval="confidence", level = 0.95))
conf_interval <- forecast(m_ar1, xreg = ccc$year)
lines(1979:2021, conf_interval$mean, col = dat_col)
polygon(c(1979:2021, 2021:1979), c(conf_interval$lower[,1], rev(conf_interval$upper[,1])),
        col = alpha(dat_col, .3), border = NA)

# lines(1979:2019, conf_interval$fit, col = precip_col)
points(pt_yrs_nhc, cc$percent_extreme[cc$year %in% pt_yrs_nhc], pch = 3, col = precip_col, cex = 1.5, xpd = NA)
points(pt_yrs, cc$percent_extreme[cc$year %in% pt_yrs], pch = 19, col = precip_col, cex = 1)
# polygon(c(1979:2019, 2019:1979), c(conf_interval$lwr, rev(conf_interval$upr)),
#         col = alpha(precip_col, .3), border = NA)
# text(1968, 62, 'slope = 2.95 ±\n0.68%/decade', col = precip_col, cex = 1)
text(1968, 62, 'slope = 2.16 ±\n0.63%/decade', col = precip_col, cex = 1)

#no precip days
m_ar1 <- forecast::Arima(cc$zero_days, order = c(1, 0, 0), xreg = matrix(cc$year, ncol = 1))
plot(cc$year, cc$zero_days, type = 'l', lwd = 1.2, col = precip_col, axes = F)
# axis(1, padj = -0.5)
axis(2, las = 2)#, at = c(60, 70, 80), labels = c('60%', '70%', '80%'))
mtext('No Precip Days', 2, 2.7, cex = .8)
mm <- lm(zero_days~year, data = cc)
ccc <- filter(cc, ! is.na(zero_days))
m_ar1 <- forecast::Arima(ccc$zero_days, order = c(1, 0, 0), xreg = matrix(ccc$year, ncol = 1))
# conf_interval <- data.frame(predict(mm, interval="confidence", level = 0.95))
conf_interval <- forecast(m_ar1, xreg = ccc$year)
lines(1979:2021, conf_interval$mean, col = dat_col)
polygon(c(1979:2021, 2021:1979), c(conf_interval$lower[,1], rev(conf_interval$upper[,1])),
        col = alpha(dat_col, .3), border = NA)

# m_ar1 <- forecast::Arima(cc$zero_days, order = c(1, 0, 0), xreg = matrix(cc$year, ncol = 1))
# conf_interval <- data.frame(predict(mm, interval="confidence", level = 0.95))
# lines(1979:2019, conf_interval$fit, col = precip_col)
points(pt_yrs_nhc, cc$zero_days[cc$year %in% pt_yrs_nhc], pch = 3, col = precip_col, cex = 1.5, xpd = NA)
points(pt_yrs, cc$zero_days[cc$year %in% pt_yrs], pch = 19, col = precip_col, cex = 1)
# polygon(c(1979:2019, 2019:1979), c(conf_interval$lwr, rev(conf_interval$upr)),
#         col = alpha(precip_col, .3), border = NA)
# text(1968, 131, 'slope = 12.1 ±\n2.4 days/decade', col = precip_col, cex = 1, xpd=NA)
text(1968, 131, 'slope = 9.2 ±\n2.3 days/decade', col = precip_col, cex = 1, xpd=NA)
rect(1966.9, 5.5, 2021.1, 578, xpd = NA)

drought_col <- "firebrick3"
d <- read_csv("data/other_watershed_stuff/PDSI_timeseries.csv")
zz <- d %>%
  select(date, pdsi) %>%
  mutate(pdsi = pdsi * 0.01) %>%
  filter(
    month(date) %in% 9:11,
    year(date) %in% 1967:2021
  ) %>%
  group_by(year = year(date)) %>%
  summarize(pdsi = mean(pdsi))

pdsi <- forecast::Arima(zz$pdsi, order = c(1, 0, 0), xreg = matrix(zz$year, ncol = 1))
plot(zz$year, zz$pdsi, type = "l", lwd = 1.2, col = drought_col, axes = F)
points(pt_yrs_nhc, zz$pdsi[zz$year %in% pt_yrs_nhc], pch = 3, col = drought_col, cex = 1.5, xpd = NA)
points(pt_yrs, zz$pdsi[zz$year %in% pt_yrs], pch = 19, col = drought_col, cex = 1)
axis(1, padj = -0.5)
axis(2, las = 2) # , at = c(60, 70, 80), labels = c('60%', '70%', '80%'))
mtext("PDSI", 2, 2.7, cex = .8)
mm <- lm(pdsi ~ year, data = zz)
# conf_interval <- data.frame(predict(mm, interval = "confidence", level = 0.95))
# lines(1968:2019, conf_interval$fit, col = drought_col)
# polygon(c(1968:2019, 2019:1968), c(conf_interval$lwr, ,rev(conf_interval$upr)),
#   col = alpha(drought_col, .3), border = NA
# )
# text(1970, 220, "slope = 12.1 ± 0.2 days/decade", col = drought_col, cex = 1)
# rect(1967, 152.5, 2020, 537, xpd = NA)

dev.off()


