# plot monthly metabolism across sites for direct calculation
dat <- readRDS("data/metabolism/compiled/met_preds_stream_metabolizer_O2.rds")

qt <- dat$filled %>%
  filter(!(era == "then" & year %in% c(1968, 1970))) %>%
  mutate(site = case_when(era == "then" ~ "CBP",
                          TRUE ~ site)) %>%
  group_by(site, year, month, era) %>%
  summarize(across(all_of(c("discharge", "temp.water")),
                   .fns = list(mean = ~mean(., na.rm = T),
                               sd = ~sd(., na.rm = T)),
                   .names = '{col}_{fn}')) %>%
  ungroup()

qtt <- dat$filled %>%
  filter(!(era == "then" & year %in% c(1968, 1970))) %>%
  mutate(doy = as.numeric(format(date, '%j')),
         site = case_when(era == "then" ~ "CBP",
                          TRUE ~ site),
         site = factor(site, levels = c('NHC', 'PM', 'WB',
                                        'WBP','UNHC','CBP')))

qtt <- qtt %>%
  group_by(site, year, doy) %>%
  summarize(across(all_of(c("discharge", "temp.water")),
                    ~mean(., na.rm = T))) %>%
  ungroup() %>%
  left_join(select(qtt, site, year, doy, DO.sat_cum))

ctt <- qtt %>%
  filter(site == "CBP") %>%
  mutate(across = "decades")

stt <- qtt %>%
  filter(year == 2019) %>%
  mutate(across = "sites")

ytt <- qtt %>%
  filter(site %in% c('NHC', 'UNHC')) %>%
  mutate(across = "years")

qtt <- bind_rows(ctt, stt, ytt) %>%
  left_join(dat$filled[,c(1,21,22, 24)]) %>%
  mutate(across = factor(across, levels = c("sites", "years", "decades")),
         year = factor(year, levels = c(2019, 2017, 2018, 1969 )))

qqq <- qtt %>%
  group_by(site, year, across) %>%
  arrange(desc(discharge)) %>%
  mutate(index = 1:n())


mon <- dat$preds  %>%
  filter(site != "BLK",
         !(site == "WB" & era == "then")) %>%
  mutate(year = case_when(era == "then" ~ 1969,
                           TRUE ~ year)) %>%
  group_by(site, year, month, era) %>%
  summarize(across(all_of(c("discharge", "temp.water","GPP","ER")),
                   .fns = list(mean = ~mean(., na.rm = T),
                               sd = ~sd(., na.rm = T)),
                   .names = '{col}_{fn}'),
            n = n()) %>%
  ungroup() %>%
  mutate(pr = log10(-GPP_mean/ER_mean),
         month = factor(month.abb[month], levels = month.abb))
         # month = factor(month.abb[month], levels = month.abb[c(3:12, 1:2)]))
# tmp <- data.frame(month = 9,
#                   era = "then",
#                   n = 0,
#                   site = "CBP",
#                   year = 1969)
#
# mon <- bind_rows(mon, tmp)

cbp_mon <- dat$preds %>%
  filter(site == "CBP") %>%
  select(month, era, GPP, ER) %>%
  mutate(ER = -ER) %>%
  pivot_longer(cols = c("GPP", "ER"), names_to = "met", values_to = 'gCm2d') %>%
  group_by(month, era, met) %>%
  summarize(across(all_of("gCm2d"),
                   .fns = list(mean = ~mean(., na.rm = T),
                               sd = ~sd(., na.rm = T)),
                               # lower95 = ~quantile(., 0.025, na.rm = T)),
                   .names = '{col}_{fn}'),
            n = n())

tmp <- data.frame(month = 9,
                  era = "then",
                  met = c("GPP", "ER"),
                  n = 0)

cbp_mon <- bind_rows(cbp_mon, tmp) %>%
  mutate(met = factor(met, levels = c("GPP", "ER")))
         # across(ends_with("er95"), ~case_when(n == 1 ~ NA_real_,
         #                                      TRUE ~ .)))

tiff(filename = 'figures/metabolism_meansd_monthly_SM.tif',
     height = 3.6, width = 6, units = 'in', res = 800)

cbp_mon %>%
    mutate(n = case_when(era == "now" ~ NA_character_,
                         TRUE ~ as.character(n))) %>%
    ggplot( aes(x = month, y = gCm2d_mean, fill = era)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    # geom_text(aes(y = 0.06, label = n),
    #           position = position_dodge(.9)) +
    facet_wrap(.~met, ncol = 1) +
    geom_errorbar(aes(ymin = gCm2d_mean - gCm2d_sd,
                      ymax = gCm2d_mean + gCm2d_sd),
                  width = .2, position = position_dodge(.9)) +
    ylab(expression("Mean Daily Metabolism (g Om"^-2 * "d"^-1 * ")")) +
    xlab("") +
    # ylim(c(0, 2.5)) +
    # labs(title = "Monthly metabolism at Concrete Bridge") +
    scale_fill_manual(values=c(now_col, then_col)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none") +
    scale_x_continuous(breaks = 1:12, labels=month.abb)
dev.off()

png(filename = 'figures/metabolism_meansd_monthly_SM.png',
     height = 3.6, width = 6, units = 'in', res = 800)

cbp_mon %>%
    mutate(n = case_when(era == "now" ~ NA_character_,
                         TRUE ~ as.character(n))) %>%
    ggplot( aes(x = month, y = gCm2d_mean, fill = era)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    # geom_text(aes(y = 0.06, label = n),
    #           position = position_dodge(.9)) +
    facet_wrap(.~met, ncol = 1) +
    geom_errorbar(aes(ymin = gCm2d_mean - gCm2d_sd,
                      ymax = gCm2d_mean + gCm2d_sd),
                  width = .2, position = position_dodge(.9)) +
    ylab(expression("Mean Daily Metabolism (g Om"^-2 * "d"^-1 * ")")) +
    xlab("") +
    # ylim(c(0, 2.5)) +
    # labs(title = "Monthly metabolism at Concrete Bridge") +
    scale_fill_manual(values=c(now_col, then_col)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none") +
    scale_x_continuous(breaks = 1:12, labels=month.abb)
dev.off()

# changes in met across months ####
changes_monthly <- cbp_mon %>%
  pivot_wider(names_from = era, values_from = gCm2d_mean) %>%
  group_by(month, met) %>%
  summarize(now = mean(now, na.rm = T), then = mean(then, na.rm = T)) %>%
  mutate(delta = now - then)# %>%
  select(month, delta, met) %>%
  pivot_wider(names_from = met, values_from = delta)
changes_monthly %>% filter(month %in% c(11,12,1))%>%
  group_by(met) %>% summarize_all(mean)




m_years <- mon %>%
  filter(site %in% c("NHC", "UNHC")) %>%
  mutate(across = "years")
m_sites <- mon %>%
  filter(year == 2019) %>%
  mutate(across = "sites")
c_sites <- cbp_mon %>%
  arrange(era, month) %>%
  mutate(year = case_when(era == "now"~2019,
                          era == "then"~1969),
         gCm2d_mean = case_when(met == "ER"~ -gCm2d_mean,
                                met == "GPP"~gCm2d_mean)) %>%
  select(-era, -gCm2d_sd, -n) %>%
  group_by(met, year) %>%
  mutate(gCm2d_mean = na.approx(gCm2d_mean, x = month)) %>%
  pivot_wider(names_from = met, values_from = gCm2d_mean,
              names_glue = '{met}_mean') %>%
  mutate(pr = log10(-GPP_mean/ER_mean),
         site = "CBP",
         across = "decades",
         month = factor(month.abb[month], levels = month.abb)) %>%
         # month = factor(month.abb[month], levels = month.abb[c(3:12, 1:2)])) %>%
  left_join(mon[,c(1,2,3,7)])
#
# # vertical plot ####
# par(mfrow = c(3, 2), mar = c(1, 2, 1, 2), oma = c(3, 2, 4, 2))
#   ylim = c(min(mon$ER_mean, na.rm = T),
#                 max(mon$GPP_mean, na.rm = T))
#   ylim_p = c(min(mon$pr, na.rm = T),
#                 max(mon$pr, na.rm = T))
#   plot(1:12, type = 'n', xaxt = 'n', ann = F, ylim = ylim)
#   abline(h = 0, lwd = .5)
#   mtext("Daily metabolism", line = 1.2)
#   for(s in unique(sites$sitecode)){
#     col = 'grey50'
#     if(s == "CBP"){ col = "brown3"}
#     tmp <- m_sites[m_sites$site == s,]
#     lines(tmp$month, tmp$GPP_mean, lwd = 1.8, col = col)
#     lines(tmp$month, tmp$ER_mean, lwd = 1.8, col = col)
#   }
#   plot(1:12, type = 'n', xaxt = 'n', ann = F, ylim = ylim_p)
#   abline(h = 1, lwd = .5)
#   mtext("Average P/R", line = 1.2)
#   mtext("across sites", 4, line = 1.2)
#   for(s in unique(sites$sitecode)){
#     col = 'grey50'
#     if(s == "CBP"){ col = "brown3"}
#     tmp <- m_sites[m_sites$site == s,]
#     lines(tmp$month, tmp$pr, lwd = 1.8, col = col)
#   }
#
#   col = "grey50"
#   plot(1:12, type = 'n', xaxt = 'n', ann = F, ylim = ylim)
#   abline(h = 0, lwd = .5)
#   for(s in c("NHC", "UNHC")){
#     tmp <- m_years %>%
#       filter(year == 2019,
#              site == s)
#     lines(tmp$month, tmp$GPP_mean, lwd = 1.8, col = col)
#     lines(tmp$month, tmp$ER_mean, lwd = 1.8, col = col)
#     for(y in c(2017, 2018)){
#       tmp <- m_years %>%
#         filter(year == y,
#                site == s)
#       lines(tmp$month, tmp$GPP_mean, lwd = 1.8, col = col, lty = 2)
#       lines(tmp$month, tmp$ER_mean, lwd = 1.8, col = col, lty = 2)
#     }
#   }
#   plot(1:12, type = 'n', xaxt = 'n', ann = F, ylim = ylim_p)
#   abline(h = 1, lwd = .5)
#   mtext("across years", 4, line = 1.2)
#   for(s in c("NHC", "UNHC")){
#     tmp <- m_years %>%
#       filter(year == 2019,
#              site == s)
#     lines(tmp$month, tmp$pr, lwd = 1.8, col = col)
#     for(y in c(2017, 2018)){
#       tmp <- m_years %>%
#         filter(year == y,
#                site == s)
#       lines(tmp$month, tmp$pr, lwd = 1.8, col = col, lty = 2)
#     }
#   }
#
#   col = "brown3"
#   plot(1:12, type = 'n', xaxt = 'n', ann = F, ylim = ylim)
#   abline(h = 0, lwd = .5)
#   tmp <- mon %>%
#     filter(era == "then",
#            site == "CBP")
#     lines(tmp$month, tmp$GPP_mean, lwd = 1.8, col = col, lty = 2)
#     lines(tmp$month, tmp$ER_mean, lwd = 1.8, col = col, lty = 2)
#   tmp <- mon %>%
#     filter(era == "now",
#            site == "CBP")
#     lines(tmp$month, tmp$GPP_mean, lwd = 1.8, col = col)
#     lines(tmp$month, tmp$ER_mean, lwd = 1.8, col = col)
#
#   plot(1:12, type = 'n', xaxt = 'n', ann = F, ylim = ylim_p)
#   abline(h = 1, lwd = .5)
#   mtext("across decades", 4, line = 1.2)
#   col = "brown3"
#   tmp <- mon %>%
#     filter(era == "then",
#            site == "CBP")
#     lines(tmp$month, tmp$pr, lwd = 1.8, col = col, lty = 2)
#   tmp <- mon %>%
#     filter(era == "now",
#            site == "CBP")
#     lines(tmp$month, tmp$pr, lwd = 1.8, col = col)

# Horizontal plot ####
tiff("figures/multipanel_comps_across_scales_SM.tif",
    height = 8, width = 6, res = 800, units = 'in', compression = 'lzw')
    col_plot = "grey40"
    l_w = 2
    l_w2 = 2
    par(mfrow = c(4, 3),
        mar = c(1, 1, 0, 1),
        oma = c(6, 4, 4, 1),
        ps = 10)
      ylim = c(min(mon$ER_mean, na.rm = T),
                    max(mon$GPP_mean, na.rm = T))
      ylim_p = c(min(mon$pr, na.rm = T),
                    max(mon$pr, na.rm = T))
    # metabolism
      plot(1:12, type = 'n', axes = F, ann = F, ylim = ylim, frame.plot = F)
      abline(h = 0, lwd = .5, col = col_plot)
      mtext("Across sites", line = 1.2)
      mtext(expression(paste("Daily metabolism (g C/", m^2, "/d)")), 2,
            line = 2.4)
      axis(2, col = col_plot, col.axis = col_plot)
      axis(1, col = col_plot, at = c(1,3,5,7,9,11), labels = rep("",6))
      box(which = "plot", bty = "l", col = col_plot)
      for(s in unique(sites$sitecode)){
        col = now_col
        if(s == "CBP"){ col = then_col}
        tmp <- m_sites[m_sites$site == s,]
        lines(tmp$month, tmp$GPP_mean, lwd = l_w, col = col)
        lines(tmp$month, tmp$ER_mean, lwd = l_w, col = col)
      }
      plot(1:12, type = 'n', axes = F, ann = F, ylim = ylim)
      abline(h = 0, lwd = .5, col = col_plot)
      mtext("Across years", line = 1.2)
      col = now_col
      for(s in c("NHC", "UNHC")){
        tmp <- m_years %>%
          filter(year == 2019,
                 site == s)
        lines(tmp$month, tmp$GPP_mean, lwd = l_w, col = col)
        lines(tmp$month, tmp$ER_mean, lwd = l_w, col = col)
        for(y in c(2017, 2018)){
          tmp <- m_years %>%
            filter(year == y,
                   site == s)
          lines(tmp$month, tmp$GPP_mean, lwd = l_w, col = col, lty = 2)
          lines(tmp$month, tmp$ER_mean, lwd = l_w, col = col, lty = 2)
        }
      }
      axis(1, col = col_plot, at = c(1,3,5,7,9,11), labels = rep("",6))
      axis(2, col = col_plot, labels = FALSE)
      box(which = "plot", bty = "l", col = col_plot)
      plot(1:12, type = 'n', axes = F, ann = F, ylim = ylim)
      abline(h = 0, lwd = .5, col = col_plot)
      mtext("Across decades", line = 1.2)
      axis(1, col = col_plot, at = c(1,3,5,7,9,11), labels = rep("",6))
      axis(2, col = col_plot, labels = FALSE)
      box(which = "plot", bty = "l", col = col_plot)
      col = then_col
      tmp <- mon %>%
        filter(era == "then",
               site == "CBP")
        lines(tmp$month, tmp$GPP_mean, lwd = l_w, col = col, lty = 2)
        lines(tmp$month, tmp$ER_mean, lwd = l_w, col = col, lty = 2)
      tmp <- mon %>%
        filter(era == "now",
               site == "CBP")
        lines(tmp$month, tmp$GPP_mean, lwd = l_w, col = col)
        lines(tmp$month, tmp$ER_mean, lwd = l_w, col = col)

    # GPP/ER
      plot(1:12, type = 'n', axes = F, ann = F, ylim = ylim_p)
      abline(h = 0, lwd = .5, col = col_plot)
      mtext("Mean GPP/ER", 2, line = 2.5)
      axis(2, col = col_plot, col.axis = col_plot)
      axis(1, col = col_plot, at = c(1,3,5,7,9,11), labels = FALSE)
      box(which = "plot", bty = "l", col = col_plot)
      for(s in unique(sites$sitecode)){
        col = now_col
        if(s == "CBP"){ col = then_col}
        tmp <- m_sites[m_sites$site == s,]
        lines(tmp$month, tmp$pr, lwd = l_w, col = col)
      }

      plot(1:12, type = 'n',axes = F, ann = F, ylim = ylim_p)
      abline(h = 0, lwd = .5, col = col_plot)
      axis(1, col = col_plot, at = c(1,3,5,7,9,11), labels = rep("",6))
      axis(2, col = col_plot, labels = FALSE)
      box(which = "plot", bty = "l", col = col_plot)
      col = now_col
      for(s in c("NHC", "UNHC")){
        tmp <- m_years %>%
          filter(year == 2019,
                 site == s)
        lines(tmp$month, tmp$pr, lwd = l_w, col = col)
        for(y in c(2017, 2018)){
          tmp <- m_years %>%
            filter(year == y,
                   site == s)
          lines(tmp$month, tmp$pr, lwd = l_w, col = col, lty = 2)
        }
      }
      plot(1:12, type = 'n', axes = F, ann = F, ylim = ylim_p)
      abline(h = 0, lwd = .5, col = col_plot)
      axis(1, col = col_plot, at = c(1,3,5,7,9,11), labels = rep("",6))
      axis(2, col = col_plot, labels = FALSE)
      box(which = "plot", bty = "l", col = col_plot)
      col = then_col
      tmp <- mon %>%
        filter(era == "then",
               site == "CBP")
        lines(tmp$month, tmp$pr, lwd = l_w, col = col, lty = 2)
      tmp <- mon %>%
        filter(era == "now",
               site == "CBP")
        lines(tmp$month, tmp$pr, lwd = l_w, col = col)
    # plot temperature and Q
      ylim = c(min(qt$temp.water_mean, na.rm = T),
               max(qt$temp.water_mean, na.rm = T))
      plot(1:12, type = 'n', axes = F, ann = F, ylim = ylim)
      mtext(expression(paste("Water Temp (", degree, "C)")), 2, line = 2.4)
      axis(2, col = col_plot, col.axis = col_plot)
      axis(1, col = col_plot, at = c(1,3,5,7,9,11),
           labels = month.abb[c(1,3,5,7,9,11)], col.axis = col_plot)
      box(which = "plot", bty = "l", col = col_plot)

      for(s in unique(sites$sitecode)){
        col = now_col
        if(s == "CBP"){ col = then_col}
        tmp <- qt[qt$site == s & qt$year == 2019,] %>%
          arrange(month)
        lines(tmp$month, tmp$temp.water_mean, lwd = l_w2, col = col)
      }

      plot(1:12, type = 'n',axes = F, ann = F, ylim = ylim)
      axis(1, col = col_plot, at = c(1,3,5,7,9,11),
           labels = month.abb[c(1,3,5,7,9,11)], col.axis = col_plot)
      axis(2, col = col_plot, labels = FALSE)
      box(which = "plot", bty = "l", col = col_plot)
      mtext("Month", 1, line = 2.5)

      col = now_col
      for(s in c("NHC", "UNHC")){
        tmp <- qt[qt$site == s & qt$year == 2019,] %>%
          arrange(month)
        lines(tmp$month, tmp$temp.water_mean, lwd = l_w2, col = col)
        for(y in c(2017, 2018)){
          tmp <- qt[qt$site == s & qt$year == y,] %>%
            arrange(month)
          lines(tmp$month, tmp$temp.water_mean, lwd = l_w2, col = col, lty = 2)
        }
      }
      plot(1:12, type = 'n', axes = F, ann = F, ylim = ylim)
      axis(1, col = col_plot, at = c(1,3,5,7,9,11),
           labels = month.abb[c(1,3,5,7,9,11)], col.axis = col_plot)
      axis(2, col = col_plot, labels = FALSE)
      box(which = "plot", bty = "l", col = col_plot)
      col = then_col
      tmp <- qt[qt$site == "CBP" & qt$year == 1969,] %>%
        arrange(month)
      lines(tmp$month, tmp$temp.water_mean, lwd = l_w2, col = col, lty = 2)
      tmp <- qt[qt$site == 'CBP' & qt$year == 2019,] %>%
        arrange(month)
      lines(tmp$month, tmp$temp.water_mean, lwd = l_w2, col = col)

      # discharge
      par(mar = c(0, 1, 3, 1))#, oma = c(6, 4, 3, 1))
      ylim = c(min(qqq$discharge, na.rm = T),
               max(qqq$discharge, na.rm = T))
      plot(1:365, type = 'n', axes = F, ann = F, ylim = ylim,
           log = "y")
      mtext(expression(paste("Discharge (", m^3, "/s)")), 2, line = 2.4)
      axis(2, col = col_plot, col.axis = col_plot)
      axis(1, col = col_plot, col.axis = col_plot,
           at = c(0, 92.25, 182.5, 273.75, 365),
           labels = c("0%", "25%", "50%", "75%", "100%"))
      box(which = "plot", bty = "l", col = col_plot)

      for(s in unique(sites$sitecode)){
        col = now_col
        if(s == "CBP"){ col = then_col}
        tmp <- qqq[qqq$site == s & qqq$year == 2019,] %>%
          arrange(index)
        lines(tmp$index, tmp$discharge, lwd = l_w2, col = col)
      }

      plot(1:365, type = 'n',axes = F, ann = F, ylim = ylim,
           log = "y")
      axis(1, col = col_plot, col.axis = col_plot,
           at = c(0, 92.25, 182.5, 273.75, 365),
           labels = c("0%", "25%", "50%", "75%", "100%"))
      axis(2, col = col_plot, col.axis = col_plot, labels = F)
      mtext("Exceedence frequency", 1, line = 2.5)
      box(which = "plot", bty = "l", col = col_plot)
      col = now_col
      for(s in c("NHC", "UNHC")){
        tmp <- qqq[qqq$site == s & qqq$year == 2019,] %>%
          arrange(index)
        lines(tmp$index, tmp$discharge, lwd = l_w2, col = col)
        for(y in c(2017, 2018)){
          tmp <- qqq[qqq$site == s & qqq$year == y,] %>%
            arrange(index)
          lines(tmp$index, tmp$discharge, lwd = l_w2, col = col, lty = 2)
        }
      }
      plot(1:365, type = 'n', axes = F, ann = F, ylim = ylim,
           log = "y")
      axis(1, col = col_plot, col.axis = col_plot,
           at = c(0, 92.25, 182.5, 273.75, 365),
           labels = c("0%", "25%", "50%", "75%", "100%"))
      axis(2, col = col_plot, col.axis = col_plot, labels = F)
      box(which = "plot", bty = "l", col = col_plot)
      col = then_col
      tmp <- qqq[qqq$site == "CBP" & qqq$year == 1969,] %>%
        arrange(index)
      lines(tmp$index, tmp$discharge, lwd = l_w2, col = col, lty = 2)
      tmp <- qqq[qqq$site == 'CBP' & qqq$year == 2019,] %>%
        arrange(index)
      lines(tmp$index, tmp$discharge, lwd = l_w2, col = col)

      par(new = T, mfrow = c(1,1), oma = c(0,0,0,0))
      plot(1, type = "n", axes = F, ann = F, frame.plot = F)
      legend("bottomleft",
             legend = c("Concrete Bridge", "Other Site"),
             col = c(then_col, rep(now_col, 1)),
             lty = 1, ncol = 2, bty = 'n', lwd = l_w)
      legend("bottomright",
             legend = c("2019", "Other year"),
             col = now_col, lty = c(1,2),
             ncol = 2, bty = 'n', lwd = l_w)
  dev.off()
# # ggplot version####
# allsites <- bind_rows(m_sites, m_years, c_sites ) %>%
#     mutate(site = factor(site, levels = c('NHC', 'PM', 'WB',
#                                           'WBP','UNHC','CBP')),
#            year = factor(year, levels = c(2019, 2017, 2018, 1969 )),
#            across = factor(across, levels = c('sites', 'years', 'decades')))
# png('figures/multipanel_axis.png',
#     height = 3, width = 7.5, res = 300, units = 'in')
#   allsites %>%
#     filter(month %in% month.abb[c(3,5,7,9,11,1)]) %>%
#     ggplot(aes(month, GPP_mean, col = site, lty = year)) +
#     geom_line(lwd = 1.5) +
#     geom_line(aes(y = ER_mean), lwd = 1.5) +
#     geom_hline(yintercept = 0) +
#     facet_wrap(.~across) +
#     scale_color_manual(values = c(rep(now_col, 5),then_col)) +
#     scale_linetype_manual(values = c(1,2,2,2)) +
#     theme_minimal()
# dev.off()
# png('figures/multipanel_metab.png',
#     height = 3, width = 7.5, res = 300, units = 'in')
#   ggplot(allsites, aes(as.numeric(month), GPP_mean, col = site, lty = year)) +
#     geom_line(lwd = 1.5) +
#     geom_line(aes(y = ER_mean), lwd = 1.5) +
#     geom_hline(yintercept = 0) +
#     facet_wrap(.~across) +
#     scale_color_manual(values = c(rep(now_col, 5),then_col)) +
#     scale_linetype_manual(values = c(1,2,2,2)) +
#     theme_minimal()
# dev.off()
# png('figures/multipanel_pr.png',
#     height = 3, width = 7.5, res = 300, units = 'in')
#   ggplot(allsites, aes(as.numeric(month), pr, col = site, lty = year)) +
#     geom_line(lwd = 1.5) +
#     geom_hline(yintercept = 0) +
#     facet_wrap(.~across) +
#     scale_color_manual(values = c(rep(now_col, 5),then_col)) +
#     scale_linetype_manual(values = c(1,2,2,2)) +
#     theme_minimal()
# dev.off()
#
# png('figures/multipanel_temp.png',
#     height = 3, width = 7.5, res = 300, units = 'in')
#
#   ggplot(allsites, aes(as.numeric(month), temp.water_mean, col = site, lty = year)) +
#     geom_line(lwd = 1.5) +
#     facet_wrap(.~across, scales = "free_x") +
#     scale_color_manual(values = c(rep(now_col, 5),then_col)) +
#     scale_linetype_manual(values = c(1,2,2,2)) +
#     theme_minimal()
# dev.off()
#
# png('figures/multipanel_Q.png',
#     height = 3, width = 7.5, res = 300, units = 'in')
#   qqq %>%
#     mutate(site = factor(site, levels = c('NHC', 'PM', 'WB',
#                                           'WBP','UNHC','CBP'))) %>%
#   ggplot(aes(index, log(discharge), col = site, lty = year)) +
#     geom_line(lwd = 1.5) +
#     geom_hline(yintercept = 0) +
#     facet_wrap(.~across) +
#     scale_color_manual(values = c(rep(now_col, 5),then_col)) +
#     scale_linetype_manual(values = c(1,2,2,2)) +
#     theme_minimal()
# dev.off()
