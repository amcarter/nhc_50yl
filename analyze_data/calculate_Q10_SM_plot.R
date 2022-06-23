# calculate Q10 for NHC sites from SM metabolism estimates
dat <- readRDS("data/metabolism/compiled/met_preds_stream_metabolizer_C.rds")

met <- dat$preds %>%
  mutate(NEP = -(GPP + ER))

calc_q10 <- function(met) {
  if(!("r1" %in% colnames(met))){
    met$r1 = met$NEP
  }
  if(!('t1' %in% colnames(met))){
    met$t1 = met$temp.water
  }
  mm <- filter(met, r1 > 0)
  qq <- data.frame()
  for(i in 1:(nrow(mm))){
    resp2 = mm$r1[i]
    temp2 = mm$t1[i]
    q <- mm[-i,] %>%
      mutate(r2 = resp2,
             t2 = temp2)
    qq <- bind_rows(qq, q)
  }
  qqq <- qq %>%
    filter((t2-t1) >= 0.1) %>%
    mutate(deltaT = (t2 - t1)/10,
           rr = r2/r1,
           q10 = rr ^ (1/deltaT))

  m <- lm(log(rr) ~ 0 + deltaT, data = qqq)
  q10_m = m$coefficients[1]
  q10_sd = summary(m)$coefficients[2] * sqrt(nrow(qqq))
  q10_lm = c(exp(c(q10_m - q10_sd, q10_m, q10_m+q10_sd)),
             summary(m)$r.squared)
  names(q10_lm) <- c("lower", "median","upper", "r2")

  p<- ggplot(qqq, aes(deltaT, log(rr))) +
    geom_point() +
    geom_smooth(method = lm)

  tmp <- mm %>%
    mutate(rb = median(r1, na.rm = T),
           tb = median(t1, na.rm = T),
           r_mod = rb * q10_lm[2] ^ ((t1 - tb)/10),
           r_mod_lower = rb * q10_lm[1] ^ ((t1 - tb)/10),
           r_mod_upper = rb * q10_lm[3] ^ ((t1 - tb)/10)) %>%
    select(date, starts_with("r_mod"))
  met <- left_join(met, tmp, by = "date")

  return(list(met = met,
              q10 = q10_lm,
              p = p))

}


# calc Q10 by site by month

met <- dat$preds %>%
  mutate(NEP = -(GPP + ER))
q10_all <- data.frame()
met_rmod <- data.frame()

# pdf("figures/q10_relationships.pdf", onefile = T, height = 11, width = 8.5)
# par(mfrow = c(5,3), mar = c(2,2,3,1))
for(s in unique(sites$sitecode)) {
  smet <- filter(met, site == s)
  for(y in unique(smet$year)){
    ymet <- filter(smet, year == y)
    # mm <- filter(smet, year == y)
    for(m in 1:12){
      mm <- filter(ymet, month == m)
      mm$r1 <- -mm$ER
      n <- length(which(mm$r1 > 0))
      if(n < 5) {next}
      qmod <- calc_q10(mm)
      # print(p + ggtitle(paste(s, y, m, sep = " ")))
      q10 = data.frame(site = s,
                       year = y,
                       month = m,
                       q10_lower = NA_real_,
                       q10 = NA_real_,
                       q10_upper = NA_real_,
                       r2 = NA_real_,
                       n = n)
      q10[,4:7] <- qmod$q10
      q10_all <- bind_rows(q10_all, q10)
      qmod$met <- left_join(qmod$met, q10, by = c("site", "year", "month"))
      met_rmod <- bind_rows(met_rmod, qmod$met)
    }
  }
}
# dev.off()
met <- met_rmod  %>%
  select(-r1, -t1) %>%
  filter(q10 < 8,
         r2 >= 0.1)

write_csv(met, "data/metabolism/Q10vNEP_all_sites_SM.csv")

# png("figures/all_sites_q10_modeled_NEP.png", width = 4, height = 4,
#     res = 300, units = 'in')
  ggplot(met, aes(-ER, r_mod, col = factor(year))) +
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = r_mod_lower, ymax = r_mod_upper), width = 0) +
    geom_abline(intercept = 0, slope = 1, lty = 2, col = "brown3") +
    facet_wrap(~site)+
    xlab("respiration (gC/m2/d)") +
    ylab("Q10 based modeled respiration") +
    labs(title = "All sites based on monthly Q10")

  ggplot(met, aes(month, q10, color = factor(year))) +
    geom_point() +
    facet_wrap(~site)
# dev.off()
  labs(title = paste(met$site[1], met$year[1], month.abb[met$month[1]],
        "  Q10 = ", round(mean(met$q10, na.rm = T), 1),
        "  r2 = ", round(met$r2[1], 3)))

nhc <- met %>%
  filter(site %in% c('NHC', 'UNHC')) %>%
  mutate(year = year(date),
         period = "Year")
nhcf <-  nhc %>%
  filter(month %in% 9:11) %>%
  mutate(period = "Fall (Sep-Nov)") %>%
  bind_rows(nhc)
nhcs <-  nhc %>%
  filter(month %in% 6:8) %>%
  mutate(period = "Summer (Jun-Aug)") %>%
  bind_rows(nhcf)

nhc_seasons <- dat$preds %>%
  filter(site %in% c('NHC', 'UNHC'), year !=2020) %>%
  mutate(NEP = -(GPP + ER),
         season = factor(case_when(month %in% c(1,2,12) ~ "Winter (Dec-Feb)",
                                   month %in% c(3:5) ~ "Spring (Mar-May)",
                                   month %in% c(6:8) ~ "Summer (Jun-Aug)",
                                   month %in% c(9:11) ~ "Fall (Sep-Nov)"),
                         levels = c("Spring (Mar-May)", "Summer (Jun-Aug)",
                                    "Fall (Sep-Nov)", "Winter (Dec-Feb)")))

tiff(filename = 'figures/NEP_seasonal_temp.tif', compression = 'lzw',
     height = 3.4 * 800, width = 9 * 800, units = 'px', res = 800)
# png("figures/NEP_seasonal_temp.png", width = 9, height = 3.4,
#     res = 300, units = 'in')
ggplot(nhc_seasons, aes(temp.water, NEP, col = season, shape = site)) +
  geom_point(size = 1.2) +
  facet_wrap(.~year) +
  scale_shape_manual(values = c(19,21)) +
  xlab('Water Temperature C')+
  # scale_color_manual(values = c(1,1,then_col,1)) +
  # geom_smooth(method = lm, se =F) +
  theme_bw()
dev.off()
# ggplot(nhc_seasons, aes(log(discharge), NEP, col = season)) +
#   geom_point(size = 2) +
#   facet_wrap(.~year) +
#   # scale_shape_manual(values = c(19,21)) +
#   geom_smooth(method = lm, se =F) +
#   theme_bw()

ann_full <- data.frame()
# ss = c('NHC', 'UNHC')
ss = 'NHC'
nhc_fall <- dat$preds %>%
  mutate(NEP = -(GPP + ER),
         season = "Full Year") %>%
  filter(site %in% ss,
         year !=2020,
         !is.na(NEP))
nhc_fall1 <- nhc_fall %>%
  filter(month %in% c(10:11)) %>%
  mutate(season = "Fall (Oct-Nov)")
# nhc_fall <- nhc_fall1 %>%
#   bind_rows(nhc_fall) %>%
#   mutate(season = factor(season, levels = c("Full Year", "Fall (Oct-Nov)")))
ann = data.frame()
for(y in 2017:2019){
  m <- summary(lm(NEP ~ temp.water, data = nhc_fall1[nhc_fall1$year == y,]))
    r = round(m$r.squared, 2)
    p = round(m$coefficients[2,4],2)
    sf = round(m$coefficients[2,1],2)
    mf = paste0("r2 = ", r, ",  p =  ", p)
    # mf = paste0("slope = ",s, ", r2 = ", r, ",  p =  ", p)
  m <- summary(lm(NEP ~ temp.water, data = nhc_fall[nhc_fall$year == y,]))
    r = round(m$r.squared, 2)
    p = round(m$coefficients[2,4],2)
    sy = round(m$coefficients[2,1],2)
    my = paste0("r2 = ", r, ",  p =  ", p)
    # my = paste0("slope = ",s, ", r2 = ", r, ",  p =  ", p)
    ann = bind_rows(ann, data.frame(year = y, mf = mf, my = my,
                                    sf = sf, sy = sy))
}
ann_full <- ann %>%
  mutate(var = 'temp', site = length(ss)) %>%
  bind_rows(ann_full)
# png("figures/NEP_drivers_temp.png", width = 9, height = 3.65,
#     res = 300, units = 'in')
tt <- ggplot(nhc_fall, aes(temp.water, NEP), col = 1) +
  geom_point(size = 2) +
  facet_wrap(.~year) +
  geom_smooth(method = lm, se = F, col = 1) +
  geom_point(data = nhc_fall1, col = fall_col, size = 2) +
  geom_smooth(data = nhc_fall1, method = lm, se =F, col = fall_col) +
  theme_bw() +
  geom_text(data = ann, x = 3, y = 5.4, aes(label = mf),
            col = fall_col, hjust = 0) +
  geom_text(data = ann, x = 3, y = 4.9, aes(label = my),
            col = 1, hjust = 0) +
  ylab("") +
  xlab(expression(paste("Water Temperature (", degree, "C)")))
# dev.off()
# discharge

ann = data.frame()
for(y in 2017:2019){
  m <- summary(lm(NEP ~ log10(discharge), data = nhc_fall1[nhc_fall1$year == y,]))
    r = round(m$adj.r.squared, 2)
    p = round(m$coefficients[2,4],2)
    sf = round(m$coefficients[2,1],2)
    mf = paste0("r2 = ", r, ",  p =  ", p)
    # mf = paste0("slope = ",s, ", r2 = ", r, ",  p =  ", p)
  m <- summary(lm(NEP ~ log10(discharge), data = nhc_fall[nhc_fall$year == y,]))
    r = round(m$adj.r.squared, 2)
    p = round(m$coefficients[2,4],2)
    sy = round(m$coefficients[2,1],2)
    my = paste0("r2 = ", r, ",  p =  ", p)
    # my = paste0("slope = ",s, ", r2 = ", r, ",  p =  ", p)
  ann = bind_rows(ann, data.frame(year = y, mf = mf, my = my,
                                  sf = sf, sy = sy))
}
ann_full <- ann %>%
  mutate(var = 'logQ', site = length(ss)) %>%
  bind_rows(ann_full)
# png("figures/NEP_drivers_Q.png", width = 9, height = 3.65,
#     res = 300, units = 'in')
qq <- ggplot(nhc_fall, aes(log10(discharge), NEP), col = 1) +
  geom_point(size = 2) +
  facet_wrap(.~year) +
  geom_smooth(method = lm, se = F, col = 1) +
  geom_point(data = nhc_fall1, col = fall_col, size = 2) +
  geom_smooth(data = nhc_fall1, method = lm, se =F, col = fall_col) +
  theme_bw() +
  geom_text(data = ann, x = -1, y = 5.4, aes(label = mf),
            col = fall_col, hjust = 0) +
  geom_text(data = ann, x = -1, y = 4.9, aes(label = my),
            col = 1, hjust = 0) +
  ylab(expression(paste("Net Respiration (g C/", m^2, "/d)"))) +
  xlab(expression(paste("Discharge (", m^3, "/s)"))) +
  scale_x_continuous(breaks=c(-2,-1,0),
                     labels=c(0.01, 0.1, 1))
# dev.off()
# DO psat
ann = data.frame()
for(y in 2017:2019){
  m <- summary(lm(NEP ~ DO.obs/DO.sat, data = nhc_fall1[nhc_fall1$year == y,]))
    r = round(m$adj.r.squared, 2)
    p = round(m$coefficients[2,4],2)
    sf = round(m$coefficients[2,1],2)
    mf = paste0("r2 = ", r, ",  p =  ", p)
    # mf = paste0("slope = ",s, ", r2 = ", r, ",  p =  ", p)
  m <- summary(lm(NEP ~ DO.obs/DO.sat, data = nhc_fall[nhc_fall$year == y,]))
    r = round(m$adj.r.squared, 2)
    p = round(m$coefficients[2,4],2)
    sy = round(m$coefficients[2,1],2)
    my = paste0("r2 = ", r, ",  p =  ", p)
    # my = paste0("slope = ",s, ", r2 = ", r, ",  p =  ", p)
  ann = bind_rows(ann, data.frame(year = y, mf = mf, my = my,
                                  sf = sf, sy = sy))
}
ann_full <- ann %>%
  mutate(var = 'DO', site = length(ss)) %>%
  bind_rows(ann_full)
# png("figures/NEP_drivers_DO.png", width = 9, height = 3.65,
#     res = 300, units = 'in')
oo <- ggplot(nhc_fall, aes(DO.obs/DO.sat, NEP), col = 1) +
  geom_point(size = 2) +
  facet_wrap(.~year) +
  geom_smooth(method = lm, se = F, col = 1) +
  geom_point(data = nhc_fall1, col = fall_col, size = 2) +
  geom_smooth(data = nhc_fall1, method = lm, se =F, col = fall_col) +
  theme_bw() +
  geom_text(data = ann, x = 0.5, y = 5.4, aes(label = mf),
            col = fall_col, hjust = 0) +
  geom_text(data = ann, x = 0.5, y = 4.9, aes(label = my),
            col = 1, hjust = 0) +
  ylab("") +
  xlab("DO (% saturation)") +
  scale_x_continuous(breaks=c(0,.25,.5,.75,1),
                   labels=c("0","25","50", "75","100"))

# dev.off()
tiff(filename = 'figures/NEP_drivers_NHC.tif', compression = 'lzw',
     width = 9*800, height = 9.5*800, units = 'px', res = 800)
     ggpubr::ggarrange(tt, qq, oo, ncol = 1, align ='v')
dev.off()

write_csv(ann_full, 'data/NEP_drivers_correlations_NHC_UNHC.csv')
tiff(filename = "figures/NEP_drivers_legend.tif", compression = 'lzw',
    width = 10*800, height = 4.1*800,
    res = 800, units = 'px')
    bind_rows(nhc_fall, nhc_fall1) %>%
    ggplot(aes(log(discharge), NEP, col = season)) +
      geom_point(size = 2) +
      facet_wrap(.~year) +
      scale_color_manual(values = c(then_col,1)) +
      geom_smooth(method = lm, se =F) +
      theme_bw() +
      theme(legend.position = 'top')
dev.off()

m <- summary(lm(NEP ~ log10(discharge), data = nhc_fall1))
r = round(m$adj.r.squared, 2)
p = round(m$coefficients[2,4],2)
s = round(m$coefficients[2,1],1)
mf = paste0("slope = ",s , ", r2 = ", r, ",  p =  ", p)
ann = data.frame(mf = mf)

qq<- ggplot(nhc_fall1, aes(log10(discharge), NEP)) +
  geom_point(size = 1.2) +
  geom_point(aes(col = factor(year)), size = 1.2) +
  # geom_smooth(aes(col = factor(year)),method = lm, se = F) +
  geom_smooth(method = lm, se = F, col = 1) +
  theme_bw()+
  # geom_text(data = ann, x = -0.75, y = 5.4, aes(label = mf),
  #           size = 4, hjust = 0) +
  ylab("") +
  xlab(expression(paste("Discharge (", m^3,"/s)"))) +
  scale_x_continuous(breaks=c(-2,-1,0),
                     labels=c(0.01, 0.1, 1)) +
  # ggtitle("Fall Respiration Across Years") +
  labs(col = "Year") +
  theme(legend.title.align = .6)
tt <- ggplot(nhc_fall1, aes(temp.water, NEP)) +
  geom_point(size = 1.2) +
  geom_point(aes(col = factor(year)), size = 1.2) +
  geom_smooth(aes(col = factor(year)),method = lm, se = F) +
  theme_bw()+
  ylab(expression(paste("Net Respiration (g C/", m^2, "/d)"))) +
  xlab(expression(paste("Water Temperature (",degree,"C)"))) +
  # ggtitle("Autumn Respiration Across Years (Oct - Nov)") +
  labs(col = "Year") +
  theme(legend.title.align = .6)

tiff(filename = "figures/NEP_drivers_autumn_across_years.tif",
     compression = 'lzw', width = 6, height = 5*6/8.4, res = 800, units = 'in')
  ggpubr::ggarrange(tt, qq, common.legend = T, labels = c("A","B"),
                    align = 'h', vjust = 3.5, label.y = 1.08)
dev.off()


nhcall <- dat$preds %>%
  mutate(NEP = -(GPP + ER),
         season = "Full Year") %>%
  filter(year !=2020,
         !is.na(NEP))

ann = data.frame()
for(y in 2017:2019){
  nhc <- filter(nhcall, year == y)
for(s in unique(nhc$site)){
  m <- summary(lm(NEP ~ DO.obs/DO.sat, data = nhc[nhc$site == s,]))
  r = round(m$adj.r.squared, 2)
  p = round(m$coefficients[2,4],2)
  sf = round(m$coefficients[2,1],2)
  mf = paste0("r2 = ", r, ",  p =  ", p)
  ann = bind_rows(ann, data.frame(site = s,year = y, mf = mf, sf = sf))
  }
}
# Q10 plots####
# png("figures/NHC_q10_modeled_NEP.png", width = 7, height = 4,
#     res = 300, units = 'in')
# ggplot(nhcf, aes(NEP, r_mod, col = log(discharge))) +
#   geom_point(size = 1.5) +
#   geom_abline(intercept = 0, slope = 1, lty = 2) +
#   geom_errorbar(aes(ymin = r_mod_lower, ymax = r_mod_upper), width = 0) +
#   xlab("NEP (gC/m2/d)") +
#   ylab("Q10 based modeled NEP") +
#   facet_grid(period~year) +
#   theme_bw() +
#   ylim(0,10)
# dev.off()

# png("figures/NHC_moXyear_q10_modeled_NEP.png", width = 12, height = 6,
#     res = 300, units = 'in')
# ggplot(nhc, aes(NEP, r_mod, col = site)) +
#   geom_point(size = 1.5) +
#   geom_abline(intercept = 0, slope = 1, lty = 2) +
#   geom_errorbar(aes(ymin = r_mod_lower, ymax = r_mod_upper), width = 0) +
#   xlab("NEP (gC/m2/d)") +
#   ylab("Q10 based modeled NEP") +
#   facet_grid(year~month) +
#   theme_bw() +
#   ylim(0,5)
# dev.off()

# png("figures/NHC_discharge_exceedences.png", width = 7, height = 2.5,
#     res = 300, units = 'in')
# nhc %>%
#   arrange(year, discharge) %>%
#   mutate(index = seq(nrow(nhc):1)) %>%
# ggplot(aes(x = index, y = log(discharge), col = log(discharge))) +
#   geom_line(lwd = 2) +
#   facet_grid(.~year, scales = "free") +
#   theme_bw()
# dev.off()

# nhcf <- nhc %>%
#   filter(month %in% 10:11) %>%
#   group_by(year) %>%
#   summarize(across(starts_with('q10'), median, na.rm = T)) %>%
#   mutate(period = "Fall")
# nhcf <- nhc %>%
#   group_by(year) %>%
#   summarize(across(starts_with('q10'), median, na.rm = T)) %>%
#   mutate(period = "Year") %>%
#   bind_rows(nhcf)
# png("figures/NHC_Q10s.png", width = 7, height = 2.6,
#     res = 300, units = 'in')
# ggplot(nhcf, aes(factor(year), y = q10, fill = period)) +
#   geom_bar(stat = 'identity', position = 'dodge',
#            width = .5) +
#   scale_fill_manual(values = c("grey60", "grey15")) +
#   xlab("Year") +
#   ylab("Q10 for NEP") +
#   theme_bw()
# dev.off()

# original NEP driver plots####
# metn <- dat$preds %>%
#   filter(site %in% c("NHC", "UNHC")) %>%
#   mutate(NEP = (GPP + ER),
#          NEP.upper = GPP.upper + ER.upper,
#          NEP.lower = GPP.lower + ER.lower,
#          DO.psat = DO.obs/DO.sat)
# png("figures/NEP_temp_NHC_UNHC.png", width = 7, height = 2.5,
#     res = 300, units = 'in')
# metn %>%
#   # filter(month %in% 7:9) %>%
#   ggplot(aes(temp.water, NEP, shape = site)) +
#   geom_point(size = 2) +
#   facet_grid(.~year) +
#   scale_shape_manual(values = c(19,21)) +
#   theme_bw()
# dev.off()
#
# png("figures/NEP_temp_NHC_UNHC_sum_logQ.png", width = 7, height = 2.5,
#     res = 300, units = 'in')
# metn %>%
#   # filter(month %in% 7:9) %>%
#   mutate(log_Q = log(discharge)) %>%
#   ggplot(aes(temp.water, NEP, col = log_Q)) +
#   geom_point(size = 2) +
#   facet_grid(.~year) +
#   theme_bw()
# dev.off()
#
# png("figures/NEP_temp_NHC_UNHC_sum_DO.png", width = 7, height = 2.5,
#     res = 300, units = 'in')
# metn %>%
#   # filter(month %in% 7:9) %>%
#   ggplot(aes(temp.water, NEP, col = DO.psat)) +
#   geom_point(size = 2) +
#   facet_grid(.~year) +
#   xlab("Water Temp C") +
#   scale_color_gradient(low = "brown3", high = "black")+
#   theme_bw()
# dev.off()

