# Format and compile summary metab data from NHC raymond runs.
# 2021 01 05

library(tidyverse)
library(lubridate)
library(zoo)
setwd('C:/Users/alice.carter/git/nhc_50yl/src')
# load site and met data ####
site_dat <- read_csv("data/siteData/NHCsite_metadata.csv") %>%
  slice(c(1:5,7)) %>%
  select(site = sitecode, distance_m, width_mar_m, slope = slope_nhd)
site_dat <- read_csv("../data/rating_curves/calculated_channel_dimensions_maroct.csv") %>%
  right_join(site_dat, by = "distance_m")

sm_met <- readRDS("data/metabolism/compiled/met_preds_stream_metabolizer_O2.rds")

# Summarize metabolism by month ####
monthly <- sm_met$preds %>%
  filter(year != 2020) %>%
  mutate(year = case_when(era == "then" ~ 1969,
                          TRUE ~ year)) %>%
  # filter(GPP < 5)
  select(-starts_with("K600"), -ends_with(c('lower', 'upper')),
         -level_m, -doy, -method, -date, -era, -depth) %>%
  group_by(site, year, month) %>%
  summarize(across(everything(),
                   .fns = list(mean = ~mean(.,na.rm = T),
                               se = ~sd(.,na.rm = T)/sqrt(n())),
                   .names = "{col}_{fn}")) %>%
  ungroup()

cbp_monthly <- monthly %>%
  filter(site == "CBP")

# plot metabolism by month ####
png("figures/2019monthly_avg_met_all_years_directcalc.png",
    width = 7.5, height = 5, res = 300, units = "in")
monthly %>%
  filter(site %in% c("NHC", "UNHC")) %>%
ggplot(aes(x = month, y = GPP_mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = GPP_mean - GPP_se,
                  ymax = GPP_mean + GPP_se),
              fill = alpha("forestgreen", .3), col = NA) +
  geom_line(aes(y = ER_mean)) +
  geom_ribbon(aes(ymin = ER_mean - ER_se,
                  ymax = ER_mean + ER_se),
              fill = alpha("sienna", .3), col = NA) +
  facet_grid(site ~ year) +
  geom_hline(yintercept = 0, col = "grey") +
  labs(title = "Monthly average metabolism from direct calculation",
       x = "month", y = "gC/m2/d")
dev.off()

now_19 <- monthly %>%
  filter(year == 2019)
now_y <- monthly %>%
  filter(site %in% c("NHC", "UNHC"))
then <- monthly %>%
  filter(year == 1969,
         site == "CBP")
now <- monthly %>%
  filter(year == 2019, site == "CBP")

png("figures/ER_inexcess_GPP_2019cbp.png",
    width = 4.77, height = 4.9, res = 300, units = "in")
  par(mar = c(4,4,1,1), mfrow = c(1,1))
  plot(1, type = 'n', xlim = c(1,12), ylim = c(0,4), ylab = "g O2/m2/d",
       xlab = "month", xaxt = "n",
       main = "Respiration in excess of GPP")
  polygon(c(1:12,12:1),  c(-now$ER_mean, rev(now$GPP_mean)),
          col = alpha("black", .3), border = NA, lwd = 2)
  polygon(c(2.27,3.75,3.75,2.27), c(0,0,3,3), col = 'white', border = NA)
  lines(1:12, -now$ER_mean, lty = 2, lwd = 2)
  points(1:12, -now$ER_mean, pch = 19)
  lines(1:12, now$GPP_mean, lwd = 2)
  points(1:12, now$GPP_mean, pch = 19)
  arrows(x0 = 1:12, y0 = -now$ER_mean,
         x1 = 1:12, y1 = -now$ER_mean - now$ER_se,
         length = 0.05, lwd = 2, angle = 90)
  arrows(x0 = 1:12, y0 = -now$ER_mean,
         x1 = 1:12, y1 = -now$ER_mean + now$ER_se,
         length = 0.05, lwd = 2, angle = 90)
  arrows(x0 = 1:12, y0 = now$GPP_mean,
         x1 = 1:12, y1 = now$GPP_mean + now$GPP_se,
         length = 0.05, lwd = 2, angle = 90)
  arrows(x0 = 1:12, y0 = now$GPP_mean,
         x1 = 1:12, y1 = now$GPP_mean - now$GPP_se,
         length = 0.05, lwd = 2, angle = 90)

  axis(1, at = seq(1, 12, by = 1), labels =
         substr(month.abb[seq(1,12, by = 1)], 1,1))

dev.off()
png("figures/ER_inexcess_GPP_comparison_cbp2.png",
    width = 4.77, height = 4.9, res = 300, units = "in")
  par(mar = c(4,4,1,1), mfrow = c(1,1))
  plot(1, type = 'n', xlim = c(1,12), ylim = c(0,4), ylab = "g O2/m2/d",
       xlab = "month", xaxt = "n",
       main = "Respiration in excess of GPP")
  polygon(c(1:12,12:1),  c(-now$ER_mean, rev(now$GPP_mean)),
          col = alpha("black", .3), border = NA, lwd = 2)
  polygon(c(2.27,3.75,3.75,2.27), c(0,0,3,3), col = 'white', border = NA)
  polygon(c(1:8, 10:12, 12:10, 8:1), c(-then$ER_mean, rev(then$GPP_mean)),
          col = alpha("brown3", .4), border = NA)
  lines(1:12, -now$ER_mean, lty = 2, lwd = 2)
  points(1:12, -now$ER_mean, pch = 19)
  lines(1:12, now$GPP_mean, lwd = 2)
  points(1:12, now$GPP_mean, pch = 19)
  arrows(x0 = 1:12, y0 = -now$ER_mean,
         x1 = 1:12, y1 = -now$ER_mean - now$ER_se,
         length = 0.05, lwd = 2, angle = 90)
  arrows(x0 = 1:12, y0 = -now$ER_mean,
         x1 = 1:12, y1 = -now$ER_mean + now$ER_se,
         length = 0.05, lwd = 2, angle = 90)
  arrows(x0 = 1:12, y0 = now$GPP_mean,
         x1 = 1:12, y1 = now$GPP_mean + now$GPP_se,
         length = 0.05, lwd = 2, angle = 90)
  arrows(x0 = 1:12, y0 = now$GPP_mean,
         x1 = 1:12, y1 = now$GPP_mean - now$GPP_se,
         length = 0.05, lwd = 2, angle = 90)

  axis(1, at = seq(1, 12, by = 1), labels =
         substr(month.abb[seq(1,12, by = 1)], 1,1))

  lines(then$month, -then$ER_mean, lwd = 2, lty = 2, col = 'brown4')
  lines(then$month, then$GPP_mean, lwd = 2, col = 'brown4')
  points(then$month, -then$ER_mean, pch = 19, col = 'brown4')
  points(then$month, then$GPP_mean, pch = 19, col = 'brown4')
  arrows(x0 = then$month, y0 = -then$ER_mean,
         x1 = then$month, y1 = -then$ER_mean + then$ER_se,
         length = 0.05, lwd = 2, angle = 90, col = 'brown4')
  arrows(x0 = then$month, y0 = -then$ER_mean,
         x1 = then$month, y1 = -then$ER_mean - then$ER_se,
         length = 0.05, lwd = 2, angle = 90, col = 'brown4')
  arrows(x0 = then$month, y0 = then$GPP_mean,
         x1 = then$month, y1 = then$GPP_mean - then$GPP_se,
         length = 0.05, lwd = 2, angle = 90, col = 'brown4')
  arrows(x0 = then$month, y0 = then$GPP_mean,
         x1 = then$month, y1 = then$GPP_mean + then$GPP_se,
         length = 0.05, lwd = 2, angle = 90, col = 'brown4')
dev.off()


colors = MetBrewer::met.brewer(name="Kandinsky", n=4)

png("figures/ER_inexcess_GPP_comparison_cbp.png",
    width = 4.5, height = 4.5, res = 300, units = "in")

  par(mar = c(0,4,0,1), oma = c(4,0,1,1), mfrow = c(2,1))
  plot(1, type = 'n', xlim = c(1,12), ylim = c(0,3.5), ylab = "",
       xlab = "month", xaxt = "n", yaxt = 'n', yaxs = 'i')
  axis(2, las = 2, at = 0:4, labels = 0:4)
  mtext('1969', line = -1.5, adj = .05)
  polygon(c(1:8, 10:12, 12:10, 8:1), c(-then$ER_mean, rev(then$GPP_mean)),
          col = alpha("black", .2), border = NA)
  polygon(c(2.7,3.25,3.25,2.7), c(0,0,3,3), col = 'white', border = NA)
  lines(then$month, -then$ER_mean, lwd = 2, col = colors[2])
  lines(then$month, then$GPP_mean, lwd = 2, col = colors[1])
  points(then$month, -then$ER_mean, pch = 19, col =  colors[2])
  points(then$month, then$GPP_mean, pch = 19, col = colors[1])
  arrows(x0 = then$month, y0 = -then$ER_mean,
         x1 = then$month, y1 = -then$ER_mean + then$ER_se,
         length = 0.05, lwd = 2, angle = 90,  col = colors[2])
  arrows(x0 = then$month, y0 = -then$ER_mean,
         x1 = then$month, y1 = -then$ER_mean - then$ER_se,
         length = 0.05, lwd = 2, angle = 90, col = colors[2])
  arrows(x0 = then$month, y0 = then$GPP_mean,
         x1 = then$month, y1 = then$GPP_mean - then$GPP_se,
         length = 0.05, lwd = 2, angle = 90, col = colors[1])
  arrows(x0 = then$month, y0 = then$GPP_mean,
         x1 = then$month, y1 = then$GPP_mean + then$GPP_se,
         length = 0.05, lwd = 2, angle = 90, col = colors[1])
  legend('topright', c('ER', 'GPP'), col = c(colors[2], colors[1]),
         bty = 'n', lty = c(1,1), lwd = 2, cex = 0.7)

  plot(1, type = 'n', xlim = c(1,12), ylim = c(0,3.5), ylab = "",
       xlab = "month", xaxt = "n", yaxt = 'n', yaxs = 'i')
  axis(1, at = seq(1, 12, by = 1),
       labels = substr(month.abb[seq(1,12, by = 1)], 1,1))
  axis(2, las = 2, at = 0:4, labels = 0:4)

  polygon(c(1:12,12:1),  c(-now$ER_mean, rev(now$GPP_mean)),
          col = alpha("black", .2), border = NA, lwd = 2)
  polygon(c(2.27,3.75,3.75,2.27), c(0.1,0.1,3,3), col = 'white', border = NA)
  points(1:12, -now$ER_mean, col = colors[2], pch = 19)
  lines(1:12, -now$ER_mean, col = colors[2], lwd = 2)
  lines(1:12, now$GPP_mean,  col = colors[1],lwd = 2)
  points(1:12, now$GPP_mean, col = colors[1], pch = 19)
  arrows(x0 = 1:12, y0 = -now$ER_mean,
         x1 = 1:12, y1 = -now$ER_mean - now$ER_se,
         length = 0.05, col = colors[2], lwd = 2, angle = 90)
  arrows(x0 = 1:12, y0 = -now$ER_mean,
         x1 = 1:12, y1 = -now$ER_mean + now$ER_se,
         length = 0.05, col = colors[2], lwd = 2, angle = 90)
  arrows(x0 = 1:12, y0 = now$GPP_mean,
         x1 = 1:12, y1 = now$GPP_mean + now$GPP_se,
         length = 0.05, col = colors[1], lwd = 2, angle = 90)
  arrows(x0 = 1:12, y0 = now$GPP_mean,
         x1 = 1:12, y1 = now$GPP_mean - now$GPP_se,
         length = 0.05, col = colors[1], lwd = 2, angle = 90)
  mtext('2019', line = -1.5, adj = .05)
  par(mfrow = c(1,1), new = T)
  mtext('Month', side = 1, line = 2.2)
  mtext(expression('Metabolism (g O'[2] * 'm'^-2 * 'd'^-1 * ')'), side = 2, line = 2.2)
dev.off()


then <- then %>%
  mutate(ER_se = ifelse(is.na(ER_se), 0.3, ER_se),
         GPP_se = ifelse(is.na(GPP_se), 0.2, GPP_se))
png("figures/hall_monthly_met.png",
    width = 6, height = 5, res = 300, units = "in")
plot(1, type = 'n', xlim = c(1,12), ylim = c(-5,4),
     ylab = "Metabolism g O2/m2/d",
     xlab = "month", xaxt = "n")
abline(h = 0)
lines(then$month, then$ER_mean, lwd = 2, col = 'sienna')
lines(then$month, then$GPP_mean, lwd = 2, col = 'forestgreen')
polygon(c(then$month, rev(then$month)), c(then$ER_mean + then$ER_se,
                                          rev(then$ER_mean - then$ER_se)),
        col = alpha('sienna', .5), border = NA)
polygon(c(then$month, rev(then$month)), c(then$GPP_mean + then$GPP_se,
                                          rev(then$GPP_mean - then$GPP_se)),
        col = alpha('forestgreen', .5), border = NA)
axis(1, at = seq(1, 12, by = 1), labels = substr(month.abb[seq(1,12, by = 1)],1,1))
dev.off()


png("figures/hall_monthly_met_flipped3.png",
    width = 6, height = 5, res = 300, units = "in")
plot(1, type = 'n', xlim = c(1,12), ylim = c(-5,4),
     ylab = "Metabolism g O2/m2/d",
     xlab = "month", xaxt = "n")
abline(h = 0)
lines(then$month, then$ER_mean, lwd = 2, col = 'sienna')
lines(then$month, then$GPP_mean, lwd = 2, col = 'forestgreen')
lines(then$month, -then$ER_mean, lwd = 2, col = 'sienna', lty = 2)
axis(1, at = seq(1, 12, by = 1), labels = substr(month.abb[seq(1,12, by = 1)],1,1))
polygon(c(then$month, rev(then$month)),
        c(then$GPP_mean, rev(-then$ER_mean)),
        col = alpha('black', .2), border = NA)
dev.off()

P"],
         x1 = c(1:12)+.03,
         y1 = hmonth_preds$GPP_mean[hmonth_preds$site == "CBP"] -
           hmonth_preds$GPP_sd[hmonth_preds$site == "CBP"],
         length = 0, col = "grey30", lwd = 2, lty = 2)
  arrows(x0 = c(1:8, 10:12) - .03,
         y0 = hall_month$gpp_mean[hall_month$site == "CBP"] +
           hall_month$gpp_sd[hall_month$site == "CBP"],
         x1 = c(1:8, 10:12) -.03,
         y1 = hall_month$gpp_mean[hall_month$site == "CBP"] -
           hall_month$gpp_sd[hall_month$site == "CBP"],
         length = 0, col = "steelblue", lwd = 2, lty = 2)
  arrows(x0 = c(1:8, 10:12) - .01,
         y0 = -hall_month$er_mean[hall_month$site == "CBP"] +
           hall_month$er_sd[hall_month$site == "CBP"],
         x1 = c(1:8, 10:12) -.01,
         y1 = -hall_month$er_mean[hall_month$site == "CBP"] -
           hall_month$er_sd[hall_month$site == "CBP"],
         length = 0, col = "steelblue", lwd = 2)
  legend("topright",
         legend = c("Now", "Then"),
         fill = c(alpha("grey30", .3), alpha("steelblue", .3)),
         bty = "n", cex = 1.3, ncol = 2)
  mtext("CBP site",line =  -1.4, adj = 0.05, cex = 1.2)

  par(mar = c(4,4,0,3))
  plot(1, type = 'n', xlim = c(1,12), ylim = c(0,3), ylab = "gC/m2/d",
       xlab = "month", xaxt = "n")
  polygon(c(1:12,12:1),
            c(-hmonth_preds$ER_mean[hmonth_preds$site == "all"],
              rev(hmonth_preds$GPP_mean[hmonth_preds$site == "all"])),
          col = alpha("grey30", .3), border = "grey30", lwd = 2)
  polygon(c(1:8,10:12, 12:10,8:1),
            c(-hall_month$er_mean[hall_month$site == "all"],
              rev(hall_month$gpp_mean[hall_month$site == "all"])),
          col = alpha("steelblue", .3), border = "steelblue", lwd = 2)
  arrows(x0 = c(1:12) + .01,
         y0 = -hmonth_preds$ER_mean[hmonth_preds$site == "all"] +
           hmonth_preds$ER_sd[hmonth_preds$site == "all"],
         x1 = c(1:12)+.01,
         y1 = -hmonth_preds$ER_mean[hmonth_preds$site == "all"] -
           hmonth_preds$ER_sd[hmonth_preds$site == "all"],
         length = 0, col = "grey30", lwd = 2)
  arrows(x0 = c(1:12) + .03,
         y0 = hmonth_preds$GPP_mean[hmonth_preds$site == "all"] +
           hmonth_preds$GPP_sd[hmonth_preds$site == "all"],
         x1 = c(1:12)+.03,
         y1 = hmonth_preds$GPP_mean[hmonth_preds$site == "all"] -
           hmonth_preds$GPP_sd[hmonth_preds$site == "all"],
         length = 0, col = "grey30", lwd = 2, lty = 2)
  arrows(x0 = c(1:8, 10:12) - .03,
         y0 = hall_month$gpp_mean[hall_month$site == "all"] +
           hall_month$gpp_sd[hall_month$site == "all"],
         x1 = c(1:8, 10:12) -.03,
         y1 = hall_month$gpp_mean[hall_month$site == "all"] -
           hall_month$gpp_sd[hall_month$site == "all"],
         length = 0, col = "steelblue", lwd = 2, lty = 2)
  arrows(x0 = c(1:8, 10:12) - .01,
         y0 = -hall_month$er_mean[hall_month$site == "all"] +
           hall_month$er_sd[hall_month$site == "all"],
         x1 = c(1:8, 10:12) -.01,
         y1 = -hall_month$er_mean[hall_month$site == "all"] -
           hall_month$er_sd[hall_month$site == "all"],
         length = 0, col = "steelblue", lwd = 2)
  mtext("All sites",line =  -1.4, adj = 0.05, cex = 1.2)
dev.off()

png("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/hall_50yl/figures/PR_cbp_nowthen.png",
    width = 6, height = 4, res = 300, units = "in")

plot(c(1:8, 10:12), -hall_month$gpp_mean[hall_month$site =="CBP"]/
       hall_month$er_mean[hall_month$site =="CBP"], type = "l",
     col = "steelblue", lwd = 2, ylim = c(0.1, 2.1), xaxt = "n",
     xlab = "month", ylab = "P:R ratio",
     main = "Productivity to Respiration at CBP")
axis(1, at = seq(1, 12, by = 1), labels = month.abb[seq(1,12, by = 1)])

# arrows(x0 = c(1:8, 10:12, 12:10, 8:1),
#        y0 = -hall_month$gpp_mean[hall_month$site =="CBP"]/
#   hall_month$er_mean[hall_month$site =="CBP"] +
#     sqrt(hall_month$gpp_sd[hall_month$site =="CBP"]^2 +
#     hall_month$er_sd[hall_month$site =="CBP"]^2),
#   c(1:8, 10:12, 12:10, 8:1),
#   -hall_month$gpp_mean[hall_month$site =="CBP"]/
#     hall_month$er_mean[hall_month$site =="CBP"] -
#     sqrt(hall_month$gpp_sd[hall_month$site =="CBP"]^2 +
#            hall_month$er_sd[hall_month$site =="CBP"]^2),
#   col = "steelblue", length = 0)
# lines(c(1:8, 10:12), -hall_month$gpp_mean[hall_month$site =="CBP"]/
#   hall_month$er_mean[hall_month$site =="CBP"],
#   col = "steelblue", lwd = 2, lty = 2)
# lines(1:12, -hmonth_preds$GPP_mean[hmonth_preds$site == "all"]/
#         hmonth_preds$ER_mean[hmonth_preds$site == "all"],
#       lwd = 2, col = "grey30")
lines(1:12, -hmonth_preds$GPP_mean[hmonth_preds$site == "CBP"]/
        hmonth_preds$ER_mean[hmonth_preds$site == "CBP"],
      lwd = 2, col = "grey30")
# arrows(x0 = c(1:8, 10:12, 12:10, 8:1),
#        y0 = -hmonth_preds$GPP_mean[hmonth_preds$site =="CBP"]/
#          hmonth_preds$ER_mean[hmonth_preds$site =="CBP"] +
#          sqrt(hmonth_preds$GPP_sd[hmonth_preds$site =="CBP"]^2 +
#                 hmonth_preds$ER_sd[hmonth_preds$site =="CBP"]^2),
#        c(1:8, 10:12, 12:10, 8:1),
#        -hmonth_preds$GPP_mean[hmonth_preds$site =="CBP"]/
#          hmonth_preds$ER_mean[hmonth_preds$site =="CBP"] -
#          sqrt(hmonth_preds$GPP_sd[hmonth_preds$site =="CBP"]^2 +
#                 hmonth_preds$ER_sd[hmonth_preds$site =="CBP"]^2),
#        col = "grey30", length = 0)

legend("topright",
       legend=c("Now","Then"),
       col = c("grey30", "steelblue"),
       lwd = 2, cex = 1.2, bty = "n")

dev.off()


