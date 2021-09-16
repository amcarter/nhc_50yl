# Calculate average width from March survey,
# convert width's and depths to xc area based on riffle vs pool geometry

# library(lubridate)
# library(tidyverse)
# library(ggplot2)
# library(ggpubr)
# 
# source("src/helpers.R")
source("C:/Users/Alice Carter/git/nhc_50yl/src/download_format_data/calc_discharge_from_crosssection.r")

# 1. load raw data files and metadata ####
# only run if updating everything
# longitudinal survey data of widths and habitat types from Feb/Mar 2019
# marlong <- read_csv("data/longitudinal_sampling/nhc_habitat_survey_20190307.csv") %>%
#   mutate(habitat = as.numeric(factor(habitat, levels = c("pool", "run", "riffle"))))
# feblong <- read_csv("data/longitudinal_sampling/nhc_habitat_survey_20190225.csv") %>%
#   mutate(habitat = as.numeric(factor(habitat, levels = c("pool", "run", "riffle"))))
# marsamp <- read_csv("data/longitudinal_sampling/NHCLongitudinalDO_20190308.csv") %>%
#   mutate(habitat = as.numeric(factor(habitat, levels = c("pool", "run", "riffle")))) %>%
#   filter(is.na(trib))
# marsamp <- marsamp[-(1:39),]
# octsamp <- read_csv("data/longitudinal_sampling/NHCLongitudinalDO_20191009.csv")
# 
# depths <- marsamp %>%
#   select(distance_m, mar = depth_cm) %>%
#   full_join(octsamp) %>%
#   select(distance_m, mar, oct = depth_cm) %>%
#   arrange(distance_m)
# 
# # plot(marlong$distance_m, marlong$width_m, type = "l")
# # points(feblong$distance_m, feblong$width_m, col = "brown3", pch = 19)
# # 
# # plot(1, 1, type = "n", xlim = c(0,8500), ylim = c(-12, 12))
# # polygon(c(marlong$distance_m, rev(marlong$distance_m)),
# #         c(marlong$width_m/2, -rev(marlong$width_m)/2),
# #         col = "steelblue", border = NA)
# # points(sites[1:7,]$distance_m, rep(0,7), col = "brown3", pch = 19)
# # points(feblong$distance_m, feblong$width_m/2, pch = 20)
# # 
# # # the february widths look similar, and there are way fewer,
# # # so for consistency, I'll use Mar.
# 
# # group all habitat type observations ####
# long <- feblong %>%
#   select(distance_m, h_feb = habitat) %>%
#   full_join(marlong, by = "distance_m")
# long <- marsamp %>%
#   select(distance_m, h_mar = habitat) %>%
#   full_join(long, by = 'distance_m') %>%
#   arrange(distance_m) %>%
#   select(-date) %>%
#   mutate(hab = NA,
#          n = NA)
# 
# for(i in 1:nrow(long)){
#   long$hab[i] = sum(c(long$h_mar[i], long$h_feb[i], long$habitat[i]), na.rm = T)
#   long$n[i] = sum(!is.na(c(long$h_mar[i], long$h_feb[i], long$habitat[i])))
# }
# long <- long %>%
#   select(-h_mar, -h_feb, -habitat) %>%
#   mutate(tat = case_when(n == 2 ~ floor(hab/n),
#                          n !=2 ~ round(hab/n)))
# 
# for(i in 1:nrow(long)){
#   if(is.nan(long$tat[i])){
#     long$tat[i] <- long$tat[i-1]
#   }
# }
# 
# long <- long %>%
#   mutate(habitat = case_when(tat == 1 ~ "pool",
#                              tat == 2 ~ "run",
#                              tat == 3 ~ "riffle"),
#          width_m = na.approx(width_m, na.rm = F)) %>%
#   select(-tat, -hab, -n)
# 
# long$width_m[1] <- 12
# 
# 
# 
# 
# # join with depth measurements ####
# long <- octsamp %>%
#   mutate(octdepth_m = depth_cm/100) %>%
#   select(distance_m, octdepth_m) %>%
#   right_join(long) %>%
#   arrange(distance_m) %>%
#   mutate(octdepth_m = na.approx(octdepth_m, na.rm = F))
# long <- marsamp %>%
#   mutate(mardepth_m = depth_cm/100) %>%
#   select(distance_m, mardepth_m) %>%
#   right_join(long) %>%
#   arrange(distance_m) %>%
#   mutate(mardepth_m = na.approx(mardepth_m, na.rm = F))
# 
# plot(long$octdepth_m, long$mardepth_m)
# 
# 
# # png("figures/nhc_size_by_habitat_long.png",
# #     width = 8, height = 5, res = 300, units = "in")
# # 
# # long %>%
# #   pivot_longer(cols = c("width_m", "octdepth_m", "mardepth_m"),
# #                names_to = "measure",
# #                values_to = "meters") %>%
# #   ggplot(aes(distance_m, meters)) +
# #     facet_wrap(. ~ measure, nrow = 3, scales = "free_y", ) +
# #     geom_line() +
# #     geom_point(aes(color = habitat)) +
# #     ggtitle("NHC size by habitat type")
# # 
# # dev.off()
# 
# write_csv(long, "data/longitudinal_sampling/nhc_longitudinal_channel_measurements.csv")
# 
# # determine what fraction of the river is what habitat type ####
# calc_habitat_specs <- function(long){
#   pool = run = riffle =
#     wp = wu = wr =
#     dmp = dmu = dmr =
#     dop = dou = dor = 0
# 
#   for(i in 1:(nrow(long)-1)){
#     if(long$habitat[i] == "pool"){
#       p = long$distance_m[i+1] - long$distance_m[i]
#       pool = pool + p
#       wp = wp + (long$width_m[i+1] + long$width_m[i]) * p / 2
#       dmp = dmp + (long$mardepth_m[i+1] + long$mardepth_m[i]) * p / 2
#       dop = dop + (long$octdepth_m[i+1] + long$octdepth_m[i]) * p / 2
#     }
#     if(long$habitat[i] == "run"){
#       u = long$distance_m[i+1] - long$distance_m[i]
#       run = run + u
#       wu = wu + (long$width_m[i+1] + long$width_m[i]) * u / 2
#       dmu = dmu + (long$mardepth_m[i+1] + long$mardepth_m[i]) * u / 2
#       dou = dou + (long$octdepth_m[i+1] + long$octdepth_m[i]) * u / 2
#     }
#     if(long$habitat[i] == "riffle"){
#       r = long$distance_m[i+1] - long$distance_m[i]
#       riffle = riffle + r
#       wr = wr + (long$width_m[i+1] + long$width_m[i]) * r / 2
#       dmr = dmr + (long$mardepth_m[i+1] + long$mardepth_m[i]) * r / 2
#       dor = dor + (long$octdepth_m[i+1] + long$octdepth_m[i]) * r / 2
#     }
#   }
#   len <- sum(pool, riffle, run)
#   out <- data.frame(habitat = c("pool", "run", "riffle", "total"),
#              total_length_m = c(pool, run, riffle, len),
#              avg_width = c(wp/pool, wu/run, wr/riffle,
#                            sum(wp, wr, wu)/len),
#              avg_depth_mar = c(dmp/pool, dmu/run, dmr/riffle,
#                                sum(dmp, dmr, dmu)/len),
#              avg_depth_oct = c(dop/pool, dou/run, dor/riffle,
#                                sum(dmp, dmr, dmu)/len))
# 
#   return(out)
# }
# 
# study_reaches <- data.frame()
# for(r in 1:6){
#   site <- sites$sitecode[r]
#   dist <- sites$distance_m[r]
#   reach <- long %>%
#     filter(distance_m >= dist,
#            distance_m <= dist + 1000)
#   df <- calc_habitat_specs(reach) %>%
#     mutate(site = !!site,
#            distance_m = dist)
#   study_reaches <- bind_rows(study_reaches, df)
# }
# 
# 
# # png("figures/nhc_size_by_habitat_bar.png",
# #     width = 7, height = 6, units = "in", res = 300)
# # study_reaches %>%
# #   filter(habitat != "total") %>%
# #   # mutate(total_length_m = total_length_m / mean(total_length_m),
# #   #        avg_width = avg_width/mean(avg_width),
# #   #        avg_depth_mar = avg_depth_mar/mean(avg_depth_mar),
# #   #        avg_depth_oct = avg_depth_oct/mean(avg_depth_oct)) %>%
# #   pivot_longer(cols = c("total_length_m", "avg_width",
# #                         "avg_depth_mar", "avg_depth_oct"),
# #                names_to = "var",
# #                values_to = "meters") %>%
# # ggplot( aes(x = habitat, y = meters, fill = habitat)) +
# #   geom_bar(stat = "identity") +
# #   facet_grid(var~distance_m, scales = "free_y")
# # dev.off()
# 
# write_csv(study_reaches, "data/longitudinal_sampling/nhc_habitat_dimensions_by_reach.csv")

# Reach Characterization survey data ####

nhc_xc <- read_csv("data/longitudinal_sampling/nhc_channel_crosssections.csv") 

# get discharge for these dates
q <- read_csv("data/rating_curves/NHC_UNHC_Q.csv", guess_max = 10000) %>%
  mutate(date = as.Date(with_tz(DateTime_UTC, tz = "EST"))) %>%
  filter(date %in% unique(nhc_xc$date) | 
           date %in% as.Date(c("2019-03-07", "2019-10-09"))) %>%
  select(-DateTime_UTC, -notes, -notes_rc) %>%
  group_by(date) %>%
  summarize_all(mean, na.rm = T)

# calculate Q at all sites on longitudinal survey days
qall <- sites %>%
  select(site = sitecode, 
         wsarea = ws_area.km2, 
         distance_m) %>%
  mutate(oct9 = c(q$discharge_nhc[4], rep(NA, 5), q$discharge_unhc[4]),
         mar7 = c(q$discharge_nhc[3], rep(NA, 5), q$discharge_unhc[3])) %>%
  mutate(oct9 = na.approx(oct9, x = distance_m), # losing discharge interpolate with distance
         mar7 = na.approx(mar7, x = wsarea))


nhc_reaches <- data.frame()
for(s in c("nhc","unhc")){
  if(s == "nhc"){ Q = q$discharge_nhc[1]}
  if(s == "unhc"){ Q = q$discharge_unhc[2]}
  for(d in seq(0,1000, by = 100)){    
    tmp <- nhc_xc %>%
      filter(distance_m == d, site == s) 
    xc <- calc_xc_discharge(tmp$width_m, tmp$depth_m, rep(NA, nrow(tmp)))
    xc <- tmp %>%
      select(site, distance_m, habitat, baseflow_width) %>%
      filter(!is.na(habitat)) %>%
      bind_cols(xc) %>%
      select(-velocity_avg, -discharge) %>%
      mutate(depth_max = max(tmp$depth_m), 
             discharge_rc = Q)
    nhc_reaches <- bind_rows(nhc_reaches, xc)
  }
}

# a <- ggplot(nhc_reaches, aes(depth_max, depth_avg, color = habitat)) +
#       geom_point() +
#       geom_smooth(method = lm) +
#       xlab("thalweg depth") +
#       ylab("average depth")
# b <- ggplot(nhc_reaches, aes(depth_max, xc_area, color = habitat)) +
#       geom_point() +
#       geom_smooth(method = lm) +
#       xlab("thalweg depth") +
#       ylab("crossectional area")
# c <- ggplot(nhc_reaches, aes(width, depth_avg, color = habitat)) +
#   geom_point() +
#   geom_smooth(method = lm) +
#   xlab("width") +
#   ylab("average depth")
# d <- ggplot(nhc_reaches, aes(width, xc_area, color = habitat)) +
#   geom_point() +
#   geom_smooth(method = lm) +
#   xlab("width") +
#   ylab("crossectional area")
# 
# png("figures/nhc_crossection_relationships.png",
#     width = 7, height = 6, res = 300, units = "in")
# ggpubr::ggarrange(a,b,c,d)
# dev.off()

# save file with reach characterization calculations
nhc_reaches <- select(nhc_reaches, -depth_thalweg)
qqall <- read_csv("data/rating_curves/calculated_discharge_with_levels.csv") %>%
  select(site, depth_avg, depth_max = depth_thalweg, discharge, discharge_rc, 
         velocity_avg, xc_area) %>%
  mutate(habitat = case_when(site %in% c("NHC", "UNHC", "PM", "CBP", "WBP")
                             ~ "pool",
                             site == "WB" ~ "run",
                             TRUE ~ NA_character_)) %>%
  bind_rows(nhc_reaches) %>%
  mutate(site = toupper(site)) %>%
  filter(site != "MC751") 
summary(lm(depth_avg~depth_max, data = qqall))
tiff('figures/SI/thalwegdepth_averagedepth.tif',
    width = 3*800, height = 2.5*800, res = 800, units = 'px')
# png('figures/SI/thalwegdepth_averagedepth.png',
#     width = 5, height = 3.8, res = 300, units = 'in')
  ggplot(qqall, aes(depth_max, depth_avg, color = site)) +
    geom_point() +
    geom_smooth(aes(depth_max, depth_avg), color = 'grey',  method = lm) +
    theme_bw()+
    xlab('Thalweg Depth (m)') +
    ylab('Average Depth (m)') + 
    geom_text(x = 0.38, y = 0.85, label = expression(paste("R"^"2"~" = 0.85")), 
              col = 'black') +
    theme(plot.caption = element_text(family = "Helvetica"),
          axis.text = element_text(size = 7), 
          axis.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9))
dev.off()
write_csv(qqall, 
          "data/longitudinal_sampling/nhc_channel_crosssections_calculated.csv")

# calculate average depths based on thalweg depths ####
long <- read_csv("data/longitudinal_sampling/nhc_longitudinal_channel_measurements.csv")


m <- lm(depth_avg ~ depth_max, qqall)
a <- summary(m)$coefficients[,1]
# pools <- nhc_reaches %>%
#   filter(habitat == "pool")
# rifs <- nhc_reaches %>%
#   filter(habitat != "pool")
mp <- lm(xc_area ~ (depth_max), qqall)
ap <- summary(mp)$coefficients[,1]
# mp1 <- lm(xc_area ~ (I(depth_max ^2) + depth_max), nhc_reaches)
# ap1 <- summary(mp1)$coefficients[,1]
mw <- lm(width_m ~ mardepth_m, long)
aw <- summary(mw)$coefficients[,1]
# mr <- lm(xc_area ~ depth_max, rifs)
# ar <- summary(mr)$coefficients[,1]

cchars <- long %>% 
  mutate(octavg_depth = a[1] + a[2] * octdepth_m,
         maravg_depth = a[1] + a[2] * mardepth_m,
         oct_xcarea = ap[1] + ap[2] * octdepth_m,
         mar_xcarea = ap[1] + ap[2] * mardepth_m,
         oct_width = aw[1] + aw[2] * octdepth_m,
         mar_width = aw[1] + aw[2] * mardepth_m)
         
# plot(cchars$mardepth_m, cchars$maravg_depth, pch = 20)
# points(cchars$octdepth_m, cchars$octavg_depth, pch = 20, col = 3)
# points(long$mardepth_m, long$width_m, pch = 20, col = 5)
# points(nhc_reaches$depth_max, nhc_reaches$depth_avg, pch = 20, col = 5)

write_csv(cchars, "data/rating_curves/calculated_channel_dimensions_maroct.csv")
# the depth to avg depth relationship is great, the others, not so much!
q1 <- q %>%
  pivot_longer(cols = c(2:3), names_to = "site", values_to = "level") %>%
  select(date, site, level)
specs <- q %>% 
  select(date, DateTime_EST, NHC = discharge_nhc, UNHC = discharge_unhc) %>%
  pivot_longer(cols = c(3:4), names_to = "site", values_to = "discharge") %>%
  left_join(q1, by = c("site", "date")) %>%
  mutate(avg_depth = 0,
         avg_xcarea = 0,
         avg_width = 0,
         avg_thwgdepth = 0) %>%
  slice(c(-2,-3))

nhc <- nhc_reaches %>%
  group_by(site) %>%
  select(site,depth_avg,  xc_area, width, depth_max) %>%
  summarize_all(mean)
specs[1:2,6:9] <- nhc[,2:5]

cchars <- data.frame(distance_m = seq(0, max(cchars$distance_m), by = 1)) %>%
  left_join(cchars) %>%
  mutate(mar_avg_depth = na.approx(maravg_depth),
         oct_avg_depth = na.approx(octavg_depth),
         mar_xcarea = na.approx(mar_xcarea),
         oct_xcarea = na.approx(oct_xcarea),
         mar_width = na.approx(width_m),
         oct_width = na.approx(oct_width),
         mar_thwgdepth = na.approx(mardepth_m),
         oct_thwgdepth = na.approx(octdepth_m),
         group = case_when(distance_m %in% seq(0,1000) ~ "NHC",
                           distance_m %in% seq(1570,2570) ~ "PM",
                           distance_m %in% seq(3450,4450) ~ "CBP",
                           distance_m %in% seq(5950,6950) ~ "WB",
                           distance_m %in% seq(6120,7120) ~ "WBP",
                           distance_m %in% seq(6300,7300) ~ "PWC")) %>%
  select(mar_avg_depth, oct_avg_depth, 
         mar_xcarea, oct_xcarea, 
         mar_width, oct_width, 
         mar_thwgdepth, oct_thwgdepth,
         group) %>%
  group_by(group) %>%
  summarize_all(mean)

mchars <- cchars %>% 
  select(site = group, 
         colnames(cchars)[which(grepl("mar", colnames(cchars)))])%>%
  slice(-7) %>%
  left_join(qall[,c(1,5)])
colnames(mchars) <- c("site", "avg_depth","avg_xcarea", 
                      "avg_width", "avg_thwgdepth", "discharge")
mchars$date <- as.Date("2019-03-07")
ochars <- cchars %>% 
  select(site = group, 
         colnames(cchars)[which(grepl("oct", colnames(cchars)))]) %>%
  slice(-7) %>%
  left_join(qall[,c(1,4)])
colnames(ochars) <- c("site", "avg_depth","avg_xcarea", 
                      "avg_width", "avg_thwgdepth", "discharge")
ochars$date <- as.Date("2019-10-09")
specs <- bind_rows(specs[1:2,], mchars, ochars) %>%
  select(-DateTime_EST)
# Assign depth at CBP at spring discharge to be 0.45, which is what Hall measured:
# specs$avg_depth[3] <- 0.45

# Calculate Depth by Q relationships ####
DQ <- data.frame(sitename = c("NHC","PM","CBP","WB","WBP","PWC","UNHC"),
                 c_m = rep(0.409, 7),   # depth at unit discharge
                 f = rep(.294, 7))

qmod = seq(0,1.5, by = .01)
# NHC (3 points)
site <- "NHC"
dat <- specs %>% 
  filter(site == !! site)
c = DQ$c_m[DQ$site == site]
f = DQ$f[DQ$site == site]
# plot(dat$discharge, dat$avg_depth, pch = 19)
# lines(qmod, c * qmod ^ f, lty = 2)
a <- summary(lm(log(avg_depth) ~ log(discharge), dat))$coefficients[,1]
fnhc = a[2]
cnhc = exp(a[1])
# lines(qmod, cnhc * qmod ^ fnhc)
DQ[DQ$site == site, 2:3] <- c(cnhc, fnhc)

# based on depth at unit discharge at NHC, 
# we can calculate depth at unit discharge at the other sites
# ll <- read_csv("rating_curves/all_sites_level_corrected.csv", 
#                guess_max = 100000)
qq <- read_csv("data/rating_curves/interpolatedQ_allsites.csv",
               guess_max = 10000)
for(i in 2:7){
  site = sites$sitecode[i]
  qt <- paste0(site, ".Q")
  a <- summary(lm(I(log(qq[,qt, drop = T])) ~ I(log(qq$NHC.Q))))$coefficients[,1]
  # a2 <- summary(lm(ll[,site,drop = T] ~ ll$NHC))$coefficients[,1]
  qnhc = exp(-a[1]/a[2])
  dnhc <- cnhc*(qnhc ^ fnhc)
  delta <- dnhc - specs$avg_depth[4]
  DQ$c_m[DQ$site == site] <- specs$avg_depth[specs$site == site][1] + delta
  tmp <- data.frame(discharge = 1,
                    avg_depth = c)
  dat <- specs %>% 
    filter(site == !! site) %>%
    bind_rows(tmp)
  c = DQ$c_m[DQ$site == site]
  f = DQ$f[DQ$site == site]
  plot(dat$discharge, dat$avg_depth, pch = 19, main = site)
  lines(qmod, c * qmod ^ f, lty = 2)
  a <- summary(lm(log(avg_depth) ~ log(discharge), dat))$coefficients[,1]
  f1 = a[2]
  c1 = exp(a[1])
  lines(qmod, c1 * qmod ^ f1)
  DQ[DQ$site == site, 2:3] <- c(c1, f1)
}

#not great for  PWC, WB
for(i in c( 4, 6)){
  site = sites$sitecode[i]
  dat <- specs %>% 
    filter(site == !! site)
  # plot(dat$discharge, dat$avg_depth, pch = 19, main = site)
  a <- summary(lm(log(avg_depth) ~ log(discharge), dat))$coefficients[,1]
  f1 = a[2]
  c1 = exp(a[1])
  # lines(qmod, c1 * qmod ^ f1)
  DQ[DQ$site == site, 2:3] <- c(c1, f1)
}

write_csv(specs, "data/rating_curves/depth_discharge_points.csv")
write_csv(DQ, "data/rating_curves/depth_discharge_relationship_LM1953.csv")

# same thing for VQ
qqall <- qqall %>%
  mutate(discharge = ifelse(is.na(discharge), discharge_rc, discharge),
     reach_depth = case_when(site == "NHC" ~ DQ$c_m[1] * discharge ^ DQ$f[1],
                             site == "PM" ~ DQ$c_m[2] * discharge ^ DQ$f[2],
                             site == "CBP" ~ DQ$c_m[3] * discharge ^ DQ$f[3],
                             site == "WB" ~ DQ$c_m[4] * discharge ^ DQ$f[4],
                             site == "WBP" ~ DQ$c_m[5] * discharge ^ DQ$f[5],
                             site == "UNHC" ~ DQ$c_m[7] * discharge ^ DQ$f[7]),
     reach_vel = velocity_avg * depth_avg/reach_depth)

# ggplot(qqall, aes(discharge, reach_vel)) +
#   geom_point(size = 2) +
#   stat_smooth(color = "black") +
#   geom_point(aes( color = site), size = 3)

qplot(discharge, reach_vel, data = qqall) + 
  stat_smooth (aes(outfit = vel_fit <<-..y..))
q_rng <- seq(min(qqall$discharge, na.rm = T), 
             max(qqall$discharge, na.rm = T), length.out = 80)
QV = data.frame(discharge = c(0, q_rng),
                avg_velocity = c(0, vel_fit))
write_csv(QV, "data/rating_curves/velocity_discharge_fit_relationship.csv")
# qql <- read_csv("NHC_2019_metabolism/data/rating_curves/calculated_discharge_with_levels.csv")
# 
# qql <- qql %>%filter(site %in% c("NHC")) %>%
#   arrange(discharge) %>%
#   slice(c(1,2,13,14))
# # ggplot(qql, aes(discharge, velocity_avg, color = site)) +
# #   geom_point()
# 
# # plot(qql$discharge, qql$velocity_avg, pch = 19)
# a <- summary(lm(log(velocity_avg) ~ log(discharge), qql))$coefficients[,1]
# f1 = a[2]
# c1 = exp(a[1])
# # lines(seq(0, 6, by = .1), c1 * seq(0, 6, by = .1) ^ f1)
# DQ = data.frame(sitename = "velocity",
#                 c_m = c1,
#                 f = f1) %>%
#   bind_rows(DQ)
# 
# # qql <- qql %>%
# #   filter(velocity_avg > .42)
# # m <- lm(velocity_avg ~ discharge, data = qql)
# # a <- m$coefficients
# # DQ = data.frame(sitename = "velocity_abv_4",
# #                 c_m = a[1],
# #                 f = a[2]) %>%
# #   bind_rows(DQ)
# 
# # specs <- read_csv("rating_curves/depth_discharge_points.csv")
# # DQ <- read_csv("rating_curves/depth_discharge_relationship_LM1953.csv")
# # 
# # ggplot(specs, aes(avg_depth, avg_width, color = site))+
# #   geom_point()
