library(pracma)
library(MetBrewer)

# setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism")
# setwd("C:/Users/alice.carter/git/nhc_50yl/src")
setwd("~/git/papers/alice_nhc")
# source("src/metabolism/inspect_model_fits.r")
#source("metabolism/inspect_model_fits.r")
source("code/metabolism/inspect_model_fits.r")

then_col = "brown3"
now_col = "gray"
fall_col = "brown3"

sites <- read_csv("data/siteData/NHCsite_metadata.csv") %>%
    slice(c(1:5,7))

colors = MetBrewer::met.brewer(name="Kandinsky", n=4)

# dat <- readRDS("data/metabolism/compiled/met_preds_stream_metabolizer_C.rds")
dat <- readRDS("data/metabolism/compiled/met_preds_stream_metabolizer_O2.rds")

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

# color_pal <- c("GPP" = colors[1], "ER" = colors[2],
#                "Net Heterotrophy" = alpha('black', 0.2))
gppcol <- colors[1]
ercol <- colors[2]

    # mutate(date = case_when(year == 2017 ~ date + 365*2,
    #                         year == 2018 ~ date + 365,
    #                         TRUE ~ date)) %>%

LAI <- read_csv('data/daily_modeled_light_all_sites.csv') %>%
    mutate(year = year(date)) %>%
    filter(site == 'NHC')
lf <- stats::filter(x =diff(LAI$LAI),
                    filter = rep(1,3))
litter <- slice(LAI, -1) %>%
    mutate(lf = case_when(lf > 0 ~ 0,
                          TRUE ~ -lf*100),
           doy = as.numeric(format(date, "%j")),
           Year = factor(year)) %>%
    filter(year %in% c(2017, 2018, 2019)) %>%
    group_by(year) %>%
    mutate(lf_percent = lf/max(lf),
           lf_type = case_when(lf_percent >= 0.5 ~ 'max',
                               lf_percent >= 0.1 ~ 'med',
                               TRUE ~ NA_character_)) %>%
    left_join(select(nhc, date, discharge, temp.water))


#### Plot components start here
lsum <- litter %>%
    group_by(year, lf_type) %>%
    summarize(start = min(doy),
              end = max(doy)) %>%
    filter(!is.na(lf_type)) %>%
    pivot_longer(cols = c('start', 'end'),
                 names_to = 'name',
                 values_to = 'doy') %>%
    mutate(level = 5,
           date = as.Date(paste(year, doy), format = '%Y %j'))

qsum <- litter %>%
    filter(lf_type == 'max') %>%
    ggplot(aes(factor(year), discharge))+
    geom_boxplot() +
    theme_classic() +
    scale_y_log10()+
    labs(x = '', y = 'Discharge during peak litterfall')
tsum <- litter %>%
    filter(lf_type == 'max') %>%
    ggplot(aes(factor(year), temp.water))+
    geom_boxplot() +
    theme_classic() +
    labs(x = '', y = 'Water temp during peak litterfall')
sumsum <- ggpubr::ggarrange(qsum, tsum)

met <- nhc %>%
    filter(year != 2016) %>%
    ggplot(aes(doy)) +
    geom_point(aes(y = GPP), col = '#4F716B', shape = 20, cex = 0.7) +
    geom_point(aes(y = ER), col = '#A2865C', shape = 20, cex = 0.7) +
    # geom_point(aes(y = GPP), col = 'lemonchiffon4', shape = 95, cex = 2) +
    # geom_point(aes(y = ER), col = 'wheat4', shape = 95, cex = 2) +
    geom_linerange(aes(ymin = GPP.lower, ymax = GPP.upper),
                col = alpha(gppcol, 0.6), linewidth = 0.3) +
    geom_linerange(aes(ymin = ER.lower, ymax = ER.upper),
                col = alpha(ercol, 0.6), linewidth = 0.3) +
    # geom_line(aes(y = GPP), col = 'gray40', linewidth = 0.72) +
    # geom_line(aes(y = ER), col = 'gray40', linewidth = 0.72) +
    # geom_ribbon(aes(ymin = GPP.lower, ymax = GPP.upper),
    #             col = NA, fill = alpha(gppcol, 0.5))+
    # geom_ribbon(aes(ymin = ER.lower, ymax = ER.upper),
    #             col = NA, fill = alpha(ercol, 0.5))+
    geom_segment(aes(x = hurricane, y = -14, xend = hurricane, yend = -7),
                  arrow = arrow(length = unit(0.12, 'inches')),
                  linewidth = 1, col = colors[4])+
    geom_text(aes(x = hurricane, y = -16, label = 'Hurricane Florence'),
              vjust = 1, col = colors[4], size = 4) +
    geom_hline(aes(yintercept = 0), linewidth = 0.3, col = 'gray60')+
    geom_line(data = rename(lsum, Litterfall = lf_type),
              aes(doy, level, lty = Litterfall),
              linewidth = 1, col = colors[3])+
    facet_grid(year~., scales = 'free_x')+
    labs(x = "Month",
         y = expression("Metabolism (g O"[2] * "m"^-2 * "d"^-1 * ")"),
         color = "Metabolism") +
    theme_classic() +
    scale_x_continuous(breaks = c(0,31,60,91,121,152,182,213,244,274,305,335),
                       labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) +
    theme(legend.position = 'top',
          legend.justification = 'center',
          strip.background = element_blank(),
          strip.text = element_blank())
    # theme(legend.position = c(0.22, 0.99),
    #       legend.direction = 'horizontal',
    #       strip.background = element_blank(),
    #       strip.text = element_blank())
          # panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

# cumulative metab plots

cumulplots <- nhc %>%
    filter(year != 2016) %>%
    select(date, starts_with(c('GPP', 'ER'), ignore.case = FALSE), doy, month) %>%
    mutate(across(starts_with("GPP"), ~ if_else(. < 0, 0, .)),
           across(starts_with("ER"), ~ if_else(. > 0, 0, .))) %>%
    mutate(year = as.numeric(strftime(date, format = '%Y'))) %>%
    group_by(year) %>%
    mutate(across(-any_of(c('date', 'year', 'doy', 'month')), ~ na.approx(., na.rm = FALSE, rule = 2)),
           across(-any_of(c('date', 'year', 'doy', 'month')), cumsum, .names = '{col}cumul')) %>%
    ungroup() %>%
    ggplot(aes(x = doy)) +
    geom_line(aes(y = GPPcumul), col = 'gray40', linewidth = 0.72) +
    geom_line(aes(y = ERcumul), col = 'gray40', linewidth = 0.72) +
    geom_ribbon(aes(ymin = GPP.lowercumul, ymax = GPP.uppercumul, fill = 'GPP'),
                col = NA, alpha = 0.5) +
    geom_ribbon(aes(ymin = ER.lowercumul, ymax = ER.uppercumul, fill = 'ER'),
                col = NA, alpha = 0.5) +
    geom_hline(aes(yintercept = 0), linewidth = 0.3, col = 'gray60') +
    facet_grid(year ~ ., scales = 'free_x') +
    labs(x = "Month",
         y = expression("Cumul. Metabolism (g O"[2] * "m"^-2 * ")"),
         fill = "Metabolism") +
    scale_fill_manual(values = c("GPP" = gppcol, "ER" = ercol)) +
    scale_x_continuous(breaks = c(0,31,60,91,121,152,182,213,244,274,305,335),
                       labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) +
    theme_classic() +
    theme(legend.position = 'top',
          legend.justification = 'left',
          legend.box.margin = margin(0, 0, 0, -120))

# library(gtable)
# library(grid)
# g <- ggplotGrob(cumulplots)
# panel_indices <- grep("^panel", g$layout$name)
# panel_layout <- g$layout[panel_indices, ]
# panels_to_box <- g$layout$name %in% c("panel-1-1", "panel-1-2")
# panels_positions <- g$layout[panels_to_box, ]
# top <- min(panels_positions$t)
# left <- min(panels_positions$l)
# bottom <- max(panels_positions$b)
# right <- max(panels_positions$r)
# rect <- rectGrob(
#     gp = gpar(lwd = 2, col = "red", fill = NA)
# )
# g <- gtable_add_grob(
#     g, rect, t = top, l = left, b = bottom, r = right, z = Inf, name = "rect"
# )
# grid.newpage()
# grid.draw(g)


# alice_mostor_plot <- nhc %>%
#     select(date, starts_with(c('GPP', 'ER'), ignore.case = FALSE), doy, hurricane) %>%
#     mutate(across(starts_with("GPP"), ~ if_else(. < 0, 0, .)),
#            across(starts_with("ER"), ~ if_else(. > 0, 0, .))) %>%
#     mutate(year = as.numeric(strftime(date, format = '%Y'))) %>%
#     group_by(year) %>%
#     mutate(GPP.med.mean = GPP, ER.med.mean = ER,
#            GPP.lower.mean = GPP.lower, ER.lower.mean = ER.upper,
#            GPP.upper.mean = GPP.upper, ER.upper.mean = ER.lower,
#            ER = -ER, ER.upper = -ER.upper, ER.lower = -ER.lower) %>%
#     rename(GPP.med = GPP, ER.med = ER) %>%
#     mutate(across(any_of(c('GPP.med', 'ER.med', 'GPP.lower', 'GPP.upper')), ~ na.approx(., na.rm = FALSE, rule = 2)),
#            across(any_of(c('GPP.med', 'ER.med', 'GPP.lower', 'GPP.upper')), cumsum, .names = '{col}.cumul')) %>%
#     ungroup() %>%
#     pivot_longer(cols = ends_with(c('mean', 'cumul')),
#                  names_to = c('met', 'stat', 'plot'),
#                  names_sep = '\\.', values_to = 'value') %>%
#     pivot_wider(values_from = value, names_from = stat) %>%
#     ggplot(aes(x = doy)) +
#     geom_line(aes(y = med, group = met), col = 'gray40', linewidth = 0.72) +
#     geom_ribbon(aes(ymin = lower, ymax = upper, , group = met, fill = met),
#                 col = NA, alpha = 0.5)+
#     geom_segment(aes(x = hurricane, y = -3, xend = hurricane, yend = -1),
#                   arrow = arrow(length = unit(0.12, 'inches')),
#                   linewidth = 1, col = colors[4])+
#     geom_text(aes(x = hurricane, y = -3.5, label = 'Hurricane'),
#               vjust = 1, col = colors[4])+
#     geom_hline(aes(yintercept = 0), linewidth = 0.3, col = 'gray60')+
#     geom_line(data = rename(lsum, Litterfall = lf_type),
#               aes(doy, level, lty = Litterfall),
#               linewidth = 1, col = colors[3])+
#     facet_grid(year~plot, scales = 'free', space = 'free')+
#     labs(x = "Day of year",
#          y = expression("Metabolism (gCm"^-2 * "d"^-1 * ")"),
#          color = "Metabolism") +
#     theme_classic() +
#     theme(legend.position = 'top',
#           legend.justification = c('left', 'top'))
#

    # plot flow duration curves
q_dur <- qq %>%
    select(date, discharge = NHC.Q) %>%
    mutate(year = year(date),
           discharge = case_when(discharge > 100 ~ NA_real_,
                                 TRUE ~ discharge)) %>%
    filter(year > 2016 & year < 2020)%>%
    arrange(year, discharge) %>%
    mutate(index = NA_real_)

for(y in unique(q_dur$year)){
    n = nrow(filter(q_dur, year == y))
    q_dur$index[q_dur$year == y] <- seq(100, 0, length.out = n)
}

t_dur <- nhc %>%
    select(date, temp.water) %>%
    mutate(year = year(date)) %>%
    filter(year > 2016 & year < 2020)%>%
    arrange(year, temp.water) %>%
    mutate(index = NA_real_)

for(y in unique(t_dur$year)){
    n = nrow(filter(t_dur, year == y))
    t_dur$index[t_dur$year == y] <- seq(100, 0, length.out = n)
}

fd <- q_dur %>%
    mutate(Year = factor(year))%>%
ggplot(aes(index, discharge)) +
    geom_line(aes(lty = Year)) +
    theme_classic() +
    scale_y_log10()+
    labs(x = 'Exceedence Frequency',
         y = expression('Discharge (m'^3 * 's'^-1 * ')'))

td <- t_dur %>%
    mutate(Year = factor(year))%>%
ggplot(aes(index, temp.water)) +
    geom_line(aes(lty = Year)) +
    theme_classic() +
    labs(x = 'Exceedence Frequency',
         y = 'Water Temparature (°C)')

wt <- nhc %>%
    mutate(Year = factor(year))%>%
ggplot(aes(doy, temp.water)) +
    geom_line(aes(lty = Year)) +
    scale_x_continuous(breaks = c(0,31,60,91,121,152,182,213,244,274,305,335),
                       labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) +
    theme_classic() +
    labs(x = 'Month',
         y = 'Water Temperature (°C)')


p2 <- ggpubr::ggarrange(fd, wt, common.legend = TRUE, nrow = 2)
# p2 <- ggpubr::ggarrange(fd, td, common.legend = TRUE, nrow = 2)
# p3 <- ggpubr::ggarrange(sumsum, p2, heights = c(1,2), nrow = 2, align = 'h')

png(filename = 'figures/NHC_3year_met.png', width = 10, height = 4.5,
    units = 'in', res = 300)
    ggpubr::ggarrange(met, cumulplots, p2, widths = c(2, 1.2, 1.2), ncol = 3) +
        annotate("text", x = 0.09, y = 0.86, label = "a", fontface = 'bold') +
        annotate("text", x = 0.56, y = 0.86, label = "b", fontface = 'bold') +
        annotate("text", x = 0.97, y = 0.88, label = "c", fontface = 'bold')
dev.off()
tiff(filename = 'figures/NHC_3year_met.tiff', width = 10, height = 4.5,
    units = 'in', res = 300)
    ggpubr::ggarrange(met, cumulplots, p2, widths = c(2, 1.2, 1.2), ncol = 3) +
        annotate("text", x = 0.09, y = 0.86, label = "a", fontface = 'bold') +
        annotate("text", x = 0.56, y = 0.86, label = "b", fontface = 'bold') +
        annotate("text", x = 0.97, y = 0.88, label = "c", fontface = 'bold')
dev.off()

png(filename = 'figures/NHC_qt_during_litterfall.png', width = 6, height = 3,
    res = 300, units = 'in')

    sumsum
dev.off()


ggplot(litter, aes(doy, lf_percent, group = year))+
    geom_line()

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


# concrete bridge metab plot for SI

lsum_cbp <- lsum %>%
    rename(Litterfall = lf_type) %>%
    filter(year == 2019)

cbp <- dat$preds %>%
    filter(site == 'CBP',
           year > 2000) %>%
    ggplot(aes(date)) +
    geom_point(aes(y = GPP), col = '#4F716B', shape = 20, cex = 1) +
    geom_point(aes(y = ER), col = '#A2865C', shape = 20, cex = 1) +
    geom_linerange(aes(ymin = GPP.lower, ymax = GPP.upper),
                   col = alpha(gppcol, 0.6), linewidth = 0.3) +
    geom_linerange(aes(ymin = ER.lower, ymax = ER.upper),
                   col = alpha(ercol, 0.6), linewidth = 0.3) +
    # geom_line(aes(y = GPP), col = 'gray40', linewidth = 0.72) +
    # geom_line(aes(y = ER), col = 'gray40', linewidth = 0.72) +
    # geom_ribbon(aes(ymin = GPP.lower, ymax = GPP.upper),
    #             col = NA, fill = alpha(gppcol, 0.5)) +
    # geom_ribbon(aes(ymin = ER.lower, ymax = ER.upper),
    #             col = NA, fill = alpha(ercol, 0.5)) +
    geom_hline(aes(yintercept = 0), linewidth = 0.3, col = 'gray60')+
    geom_line(data = lsum_cbp,
              aes(date, level, lty = Litterfall),
              linewidth = 1, col = colors[3]) +
    labs(x = "",
         y = expression("Metabolism (g O"[2] * "m"^-2 * "d"^-1 * ")"),
         color = "Metabolism") +
    theme_classic() +
    # scale_x_continuous(breaks = c(0,31,60,91,121,152,182,213,244,274,305,335),
    #                    labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) +
    theme(legend.position = 'top',
          legend.justification = 'left',
          strip.background = element_blank(),
          strip.text = element_blank())

png(filename = 'figures/SI/CBP_met.png', width = 6, height = 4,
    units = 'in', res = 300)
cbp
dev.off()

tiff(filename = 'figures/SI/CBP_met.tiff', width = 6, height = 4,
    units = 'in', res = 300)
cbp
dev.off()
