library(pracma)
library(MetBrewer)
# setwd("C:/Users/alice.carter/git/nhc_50yl/src")
source("code/metabolism/inspect_model_fits.r")

then_col = "brown3"
now_col = "gray"
fall_col = "brown3"

sites <- read_csv("data/siteData/NHCsite_metadata.csv") %>%
    slice(c(1:5,7))

colors = MetBrewer::met.brewer(name="Kandinsky", n=2)

dat <- readRDS("data/metabolism/compiled/met_preds_stream_metabolizer_C.rds")

preds_19 <- dat$preds %>%
  filter(year == 2019)
preds_nhc <- dat$preds %>%
  # filter(site == "NHC")
  filter(site %in% c("NHC", "UNHC"))
preds <- dat$preds %>% filter(date > as.Date('2017-01-01'))
qq <- read_csv('data/rating_curves/interpolatedQ_allsites.csv', guess_max = 10000)
# plot(qq$DateTime_UTC, qq$NHC.Q, log = 'y', type = 'l')
qq <- qq %>%
  mutate(date = as.Date(with_tz(DateTime_UTC, tz = "EST")))%>%
  group_by(date) %>%
  select(-DateTime_UTC, -notes,-notes_rc, -PWC.Q) %>%
  summarize_all(mean, na.rm = T)


# plot annual met with discharge


sites$sscode <- c("8.5 km", "6.9 km", "5 km", "2.5 km", "2.3 km", "0 km")

axis_size = 0.9
Qlim = c(.02, max(qq$NHC.Q, na.rm = T) * 1e7)
ylim = c(-9,4)
florence <- as.Date("2018-09-14")

# png("figures/metQ_across_sites_SM_2019.png",
#     width = 10, height = 6, units = 'in', res = 300)
# tiff("figures/metQ_across_sites_SM_2019.tif",
#     height = 3.6*800, width = 6*800, units = 'px', res = 800,
#     compression = 'lzw')
#   xlim = c(as.Date(c("2019-03-06", "2020-03-20")))
#   par(ps = 10,
#       oma = c(5,4,4,4),
#       mfrow = c(3,2),
#       mar = c(.5,0,0,.5))
#   for(s in c(6,3,5,2,4,1)){
#     ss <- sites$sitecode[s]
#     sn <- sites$sscode[s]
#     tmp <- preds_19 %>%
#       filter(site == ss) %>%
#       arrange(date)
#
#     tmpq <- qq[,c(1, s+1)] %>%
#       filter(date >= xlim[1],
#              date <= xlim[2])
#     colnames(tmpq) <- c('date', 'discharge')
#     plot_metab(tmp, ylim = ylim, xlim = xlim, xaxt = 'n', yaxt = 'n')
#     if(s %in% c(4,5,6)) {
#       axis(2, cex.axis = axis_size, at = c(-6, -3, 0, 3))
#     }
#     if(s %in% c(1,4)){
#       axis(1, at = seq(as.Date('2019-03-01'), by = 'month', length.out = 13),
#            labels = month.abb[c(3:12, 1:3)], cex.axis = axis_size, las = 2)
#       # axis(1, at = seq(as.Date('2019-03-01'), by = '2 months', length.out = 7),
#       #      labels = month.abb[c(3,5,7,9,11,1,3)])
#     }
#     par(new = T)
#     plot(tmpq$date, tmpq$discharge, log = "y", type = "l", ylim = Qlim, xlim = xlim,
#          axes = FALSE, xlab = '', ylab = '')
#     mtext(sn, cex = 0.8, line = -1.2, adj = 0.02)
#     if(s %in% c(1,2,3)){
#       axis(4, cex.axis = axis_size, at = c(0.01, 1, 100),
#            labels = c(0.01, 1, 100), las = 2)
#     }
#   }
#   par(new = T, mfrow = c(1,1), oma = c(1,2,0,2))
#   plot(1,1, axes = F, ann = F, type = 'n')
#   mtext(expression(paste("Metabolism g C/", m^2, "/d")), 2, .9, cex = 0.9)
#   mtext(expression(paste("Discharge (", m^3, "/s)")), 4, 1, cex = 0.9)
#   # mtext("Date", 1, 2.5)
#   mtext("New Hope Creek metabolism across sites in 2019", 3, -2.5)
#   mtext('A', 3, -2.5, adj = 0)
# dev.off()
#
# # png("figures/metQ_across_years_SM_2019.png",
# #     width = 10, height = 6, units = 'in', res = 300)
# tiff("figures/metQ_across_years_SM_2019.tif",
#      height = 3.6*800, width = 6*800, units = 'px', res = 800,
#      compression = 'lzw')
#   par(ps = 10,
#       oma = c(5,4,4,4),
#       mfrow = c(3,2),
#       mar = c(0.5,0,0,0.5))
#   xlim = c(as.Date(c("2017-03-01", "2018-03-01")))
#   for(y in c(2017,2018,2019)){
#     tmp <- preds_nhc %>%
#       filter(site == "UNHC",
#              year == y)
#     tmpq <- qq[,c(1, 7)] %>%
#       filter(date >= xlim[1],
#              date <= xlim[2])
#     colnames(tmpq) <- c('date', 'discharge')
#     plot_metab(tmp, ylim = ylim, xlim = xlim, xaxt = 'n', yaxt = 'n')
#     axis(2, cex.axis = axis_size, at = c(-6, -3, 0, 3))
#     if(y == 2018){
#       abline(v = florence, lty = 2)
#       # text(x = florence, y = -8.5, labels = "Hurricane Florence",
#       #      cex = axis_size, pos = 4)
#     }
#     if(y == 2019){
#       axis(1, at = seq(as.Date('2019-03-01'), by = 'month', length.out = 13),
#            cex.axis = axis_size, labels = month.abb[c(3:12, 1:3)], las = 2)
#       mtext("Upstream (0 km)",1, line = 2.5, cex = 0.9)
#
#     }
#     par(new = T)
#     plot(tmpq$date, tmpq$discharge, log = "y", type = "l", ylim = Qlim, xlim = xlim,
#          axes = FALSE, xlab = '', ylab = '')
#     mtext(y, cex = 0.8, line = -1.2, adj = 0.02)
#
#     tmp <- preds_nhc %>%
#       filter(site == "NHC",
#              year == y)
#     tmpq <- qq[,c(1, 2)] %>%
#       filter(date >= xlim[1],
#              date <= xlim[2])
#     colnames(tmpq) <- c('date', 'discharge')
#     plot_metab(tmp, ylim = ylim, xlim = xlim, xaxt = 'n', yaxt = 'n')
#     if(y == 2018){
#       abline(v = florence, lty = 2)
#       # text(x = florence, y = -8.5, labels = "Hurricane Florence",
#       #      cex = axis_size, pos = 4)
#     }
#     if(y == 2019){
#       axis(1, at = seq(as.Date('2019-03-01'), by = 'month', length.out = 13),
#            cex.axis = axis_size, labels = month.abb[c(3:12, 1:3)], las = 2)
#       mtext("Downstream (8.5 km)",1, line =2.5, cex = 0.9)
#     }
#     par(new = T)
#     plot(tmpq$date, tmpq$discharge, log = "y", type = "l", ylim = Qlim, xlim = xlim,
#          axes = FALSE, xlab = '', ylab = '')
#     mtext(y, cex = 0.8, line = -1.2, adj = 0.02)
#     axis(4, cex.axis = axis_size, at = c(0.01, 1, 100),
#          labels = c(0.01, 1, 100), las = 2)
#     xlim = xlim + 365
#   }
#   par(new = T, mfrow = c(1,1), oma = c(1,2,0,2))
#   plot(1,1, axes = F, ann = F, type = 'n')
#   mtext(expression(paste("Metabolism g C/", m^2, "/d")), 2, .9, cex = .9)
#   mtext(expression(paste("Discharge (", m^3, "/s)")), 4, 1, cex = .9)
#   # mtext("Date", 1, 2.5)
#   mtext("New Hope Creek metabolism across years", 3, -2.5)
#   mtext('B', 3, -2.5, adj = 0)
# dev.off()



# Plot just NHC data across the three years

nhc <- preds_nhc %>%
    mutate(year = year(date)) %>%
    filter(site == 'NHC', year < 2020) %>%
    mutate(hurricane = case_when(date == florence ~ doy,
                                 TRUE ~ NA_real_),
           GPP_fill = na.approx(GPP, na.rm = F),
           ER_fill = na.approx(ER, na.rm = F),
           GPP_high = case_when(GPP_fill > -ER_fill ~ GPP_fill,
                                TRUE ~ NA_real_))

color_pal <- c("GPP" = colors[1], "ER" = colors[2])

    # mutate(date = case_when(year == 2017 ~ date + 365*2,
    #                         year == 2018 ~ date + 365,
    #                         TRUE ~ date)) %>%

met <- ggplot(nhc, aes(doy)) +
    geom_ribbon(aes(ymin = GPP_fill, ymax = -ER_fill),
                alpha = 0.2, color = NA)+
    geom_ribbon(aes(ymax = GPP_high, ymin = -ER_fill),
                fill = 'white',color = NA)+
    geom_line(aes(y = GPP), size = 0.72, col = colors[1]) +
    geom_line(aes(y = -ER), col = colors[2], size = 0.72) +
    geom_vline(aes(xintercept = hurricane), lty = 2) +
    geom_hline(aes(yintercept = 0))+
    facet_grid(year~., scales = 'free_x')+
    ylim(0, 6.5)+
    labs(x = "Day of year",
         y = "Metabolism (gC/m2/d)",
         color = "Legend") +
    theme_classic() +
    scale_color_manual(values = color_pal)


# plot flow duration curves
q_dur <- qq %>%
    select(date, discharge = NHC.Q) %>%
    mutate(year = year(date)) %>%
    filter(year > 2016 & year < 2020)%>%
    arrange(year, discharge) %>%
    mutate(index = NA_real_)

for(y in unique(q_dur$year)){
    n = nrow(filter(q_dur, year == y))
    q_dur$index[q_dur$year == y] <- seq(100, 0, length.out = n)
}

fd <- q_dur %>%
    mutate(Year = factor(year))%>%
ggplot(aes(index, log(discharge))) +
    geom_line(aes(lty = Year)) +
    theme_classic() +
    labs(x = 'Exceedence Frequency',
         y = 'log Discharge (m3/s)')
# fd <-
wt <- nhc %>%
    mutate(Year = factor(year))%>%
ggplot(aes(doy, temp.water)) +
    geom_line(aes(lty = Year)) +
    theme_classic() +
    labs(x = 'Day',
         y = 'Water Temperature (C)')


p2 <- ggpubr::ggarrange(fd, wt, common.legend = TRUE, nrow = 2)
ggpubr::ggarrange(met, p2)




LAI <- read_csv('data/light/drivers/daily_modeled_light_all_sites.csv') %>%
    mutate(year = year(date)) %>%
    filter(site == 'NHC')
lf <- stats::filter(x =diff(LAI$LAI),
                    filter = rep(1,3))
litter <- slice(LAI, -1) %>%
    mutate(lf = case_when(lf > 0 ~ 0,
                          TRUE ~ -lf*100),
           doy = as.numeric(format(date, "%j")),
           Year = factor(year)) %>%
    filter(year ==  2018)

# get usgs discharge
library(dataRetrieval)
siteNumber <- "02085000"
parameterCd <- "00060"

# Raw daily data:
eno <- readNWISdv(
    siteNumber, parameterCd,
    "2018-01-01", "2018-12-31") %>%
    select(date = Date,
           q_cms = X_00060_00003) %>%
    mutate(q_eno = q_cms * 0.3048^3) %>%
    select(-q_cms)


litter <- left_join(litter,
          select(qq, date, discharge = NHC.Q),
          by = 'date') %>%
    left_join(eno)


plot(litter$q_eno, litter$discharge, log = 'xy')
litter <- litter %>%
    mutate(discharge = case_when(discharge > 100 ~ NA_real_,
                                 TRUE ~ discharge))

mod <- lm(log(discharge) ~ log(q_eno), data = litter)
litter$q_pred <- exp(-0.2498 + 0.7797*log(litter$q_eno))

ggplot(litter, aes(date, log(q_pred)))+
    geom_line() +
    geom_line(aes(y = log(discharge)), col = 2)

ggplot(litter, aes(doy, lf))+
    geom_line(col = 'sienna') +
    geom_line(aes(y = log(q_pred)))+
    theme_classic()


