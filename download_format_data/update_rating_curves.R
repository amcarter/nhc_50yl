# Update rating curve file with new discharge measurement
source("C:/Users/Alice Carter/git/nhc_50yl/src/download_format_data/calc_discharge_from_crosssection.r")
setwd('C:/Users/Alice Carter/git/nhc_50yl/')

# # 1. initial calcs to create rc_sheet ####
# # only run this section if you have corrected/updated all previous data:
# # load and munge data
# qdat <- read_csv("data/rating_curves/discharge_measurements.csv")
# sp_qdat <- read_csv("data/rating_curves/discharge_measurements_NHC_UNHC.csv")
# metadat <- read_csv("data/siteData/NHCsite_metadata.csv")
# 
# # # # setting negative V equal to zero to be consistent with E Moore's data
# qdat$velocity_ms[qdat$velocity_ms<0] <- 0
# 
# qq <- data.frame()
# for(site in unique(qdat$site)){
#   sdat <- qdat %>%
#   filter(site == !!site)
# 
#   for(date in unique(sdat$date)){
#     dat <- sdat %>%
#       filter(date == !!date)
# 
#     plot_xc(dat)
#     out <- calc_xc_discharge(dat$distance_m, dat$depth_m, dat$velocity_ms)
#     out <- dat %>%
#       slice(1) %>%
#       select(date, time, site, stage_m, notes) %>%
#       bind_cols(out)
#     qq <- bind_rows(qq, out)
#   }
# }
# 
# # the data from this cross sections don't seem right, so I'm leaving it out
# qq <- qq[c(-13),]
# 
# # # 3. run on old SP data for NHC and UNHC ####
# 
# qdat <- sp_qdat
# for(site in unique(qdat$site)){
#   sdat <- qdat %>%
#   filter(site == !!site)
# 
#   for(date in unique(sdat$date)){
#     dat <- sdat %>%
#       filter(date == !!date)
# 
#     plot_xc(dat)
#     out <- calc_xc_discharge(dat$distance_m, dat$depth_m, dat$velocity_ms)
#     out <- dat %>%
#       slice(1) %>%
#       select(date, time, site, stage_m, notes) %>%
#       bind_cols(out)
#     qq <- bind_rows(qq, out)
#   }
# }
# 
# dev.off()
# 
# # pair discharge measurements with metadata
# qq <- metadat %>%
#   select(site = sitecode, location_m = distance_m) %>%
#   right_join(qq, by = "site")
# 
# ggplot(qq, aes(stage_m, discharge)) +
#   geom_point() +
#   facet_wrap(.~site)
# write_csv(qq, "data/rating_curves/calculated_discharge.csv")
# 
# # Add level data from sensor to discharge measurements ####
# qq <- read_csv("data/rating_curves/calculated_discharge.csv")
# qq <- qq %>%
#   mutate(DateTime_EST = round_date(ymd_hms(paste(date, time), tz = "EST"),
#                                    unit = "15 minutes"))
# 
# levels <- read_csv("data/rating_curves/all_sites_level_corrected.csv") %>%
#   mutate(DateTime_EST = with_tz(DateTime_UTC, tz = "EST")) %>%
#   select(-notes, -DateTime_UTC)
#   # pivot_wider(names_from = "site", values_from = "level_m") %>%
#   # select(DateTime_EST, NHC_level = NHC, UNHC_level = UNHC)
# 
# qql <- left_join(qq, levels, by = c('site', 'DateTime_EST'))
# 
# # # inspect data
# qql %>%
#   filter(site %in% c( "NHC")) %>%
#   ggplot(aes(stage_m, discharge)) +
#   geom_point() +
#   facet_wrap(.~site)
# 
# # Assign proper discharges to the pool sections on days when a pool and
# # neighboring riffle were measured
# #
# # PM: that q is just wrong, but we didn't measure a riffle. Will pair with
# # interpolated Q to see what we get
# #
# # CBP
# qql %>%
#   filter(grepl("CBP", site))
# # these two Q's look comprable, I will average them.
# w <- which(grepl("CBP", qql$site) & qql$date == as.Date("2020-06-19"))
# qql$discharge[w] <- mean(qql$discharge[w], na.rm = T)
# w <- which(grepl("CBP", qql$site) & qql$date == as.Date("2020-08-23"))
# qql$discharge[w] <- mean(qql$discharge[w], na.rm = T)
# 
# # WBP
# qql %>%
#   filter(grepl("WBP", site))
# # here I am going to trust the riffle velocity
# w <- which(grepl("WBP", qql$site) & qql$date == as.Date("2020-06-19"))
# qql$discharge[w[1]] <- qql$discharge[w[2]]
# 
# # now correct the average velocites based on new Qs
# # and get rid of the riffle points
# qql$velocity_avg <- qql$discharge/qql$xc_area
# qql <- qql %>%
#   filter(!grepl("riffle", site))
# 
# 
# write_csv(qql, "data/rating_curves/calculated_discharge_with_levels.csv")


# Build Rating Curves ####

qql <- read_csv("data/rating_curves/calculated_discharge_with_levels.csv")
levels <- read_csv("data/rating_curves/all_sites_level_corrected2.csv",
                   guess_max = 1000000)

# unhc <- read_csv("NHC_2019_metabolism/data/metabolism/corrected_level/UNHC_lvl.csv") %>%
#   mutate(level_m = level_m - m_diff)

# write_csv(unhc, "NHC_2019_metabolism/data/metabolism/corrected_level/UNHC_lvl.csv")

levels <- levels %>%  
  filter(site %in% c('NHC', 'UNHC')) %>%
  pivot_wider(names_from = "site", values_from = "level_m") %>%
  mutate(DateTime_EST = with_tz(DateTime_UTC, tz = "EST")) %>%
  arrange(DateTime_EST) 

ggplot(levels, aes(NHC, UNHC, col = DateTime_UTC)) + geom_point()
abline(c(0,1))
# levels %>%
#   select(-DateTime_EST) %>%
#   xts(order.by = levels$DateTime_EST) %>%
#   dygraph() %>%
#   dyRangeSelector()
# Zero flow occurs at a stage of 0.64 at NHC (by observation),

# Build Rating Curves ####
# qql$stage_m[qql$site =="UNHC"][2] <- qql$UNHC[qql$site =="UNHC"][2]
plot(qql$stage_m, qql$level_m)
identify(qql$stage_m, qql$level_m, labels = qql$date)
abline(0,1)
# the stages on 9/22/2016 are much higher than the levels at both sites,
# but the plot below suggests that that is a problem with the level data, 
# not the stage, so I am leaving them in to calculate rating curves.

nhc <- qql %>%
  filter(site =="NHC")

# add a data point based on mannings equation:
nhc <- data.frame(stage_m = 1.6, discharge = 9.81, site = 'NHC',
                  notes = 'mannings eqn') %>% bind_rows(nhc)
# nhc <- data.frame(stage_m = .64, discharge = 0.005) %>% bind_rows(nhc)

unhc <- qql %>%
  filter(site =="UNHC") %>%
  arrange(date) 

plot(unhc$stage_m, unhc$discharge, log = 'xy')
# identify(unhc$stage_m, unhc$discharge, unhc$date)

# the point from 2016-07-15 looks like the stage or discharge is off, 
# I'm going to leave it out

unhc <- unhc[-1,]
plot(unhc$stage_m, unhc$discharge, log = 'xy')

# add a data point based on mannings equation:
unhc <- data.frame(stage_m = 1, discharge = 4.42, site = 'UNHC', 
                   notes = 'mannings eqn') %>% bind_rows(unhc)
# unhc <- data.frame(stage_m = 1.4, discharge = 9.6, site = 'UNHC',
#                    notes = 'mannings eqn') %>% bind_rows(unhc)

png('figures/ratingcurves_mannings.png', width = 480, height = 300)
m <- lm(log(discharge)~log(stage_m), data = nhc)
m_coef <- summary(m)$coefficients[,1]
# 
# What does this curve predict at 0.64 stage?
# nhc_zeroq <- exp(m_coef[1]) * 0.64 ^ m_coef[2]
par(mfrow = c(1, 2))
plot(nhc$stage_m, nhc$discharge, xlim = c(.5, 1.7), ylab = 'discharge', 
     xlab = 'NHC stage', col = 2)
lines(seq(0,2, by = .01), exp(m_coef[1]) * seq(0,2, by = .01) ^ m_coef[2], col = 2)
# remove mannings points
nhc_2 <- nhc[-1,]
m <- lm(log(discharge)~log(stage_m), data = nhc_2)
m_coef <- summary(m)$coefficients[,1]
points(nhc_2$stage_m, nhc_2$discharge)
lines(seq(0,2, by = .01), exp(m_coef[1]) * seq(0,2, by = .01) ^ m_coef[2])
legend('topleft', c('Measured', 'Manning\'s'), lty = c(1,1), col = c(1,2), bty = 'n')
# par(new = T)
# plot(density((levels$NHC), na.rm = T), col = 3, xlim = c(.5, 1.7), axes = F,
#      xlab = '', ylab = '', main = '')

# #calculate expected median UNHC by watershed area
# unhc_medq <- nhc_medq *sites$ws_area.km2[7]/sites$ws_area.km2[1]
# 
m <- lm(log(discharge)~log(stage_m), data = unhc)
m_coef <- summary(m)$coefficients[,1]
plot(unhc$stage_m, (unhc$discharge), xlim = c(.2, 1.5), #ylim = c(0,2), 
     ylab = 'discharge', xlab = 'UNHC stage', col = 2)
lines(seq(0,2, by = .01), exp(m_coef[1]) * seq(0,2, by = .01) ^ m_coef[2], col = 2)
unhc_2 <- unhc[-c(1),]
m <- lm(log(discharge)~log(stage_m), data = unhc_2)
m_coef <- summary(m)$coefficients[,1]

points(unhc_2$stage_m, unhc_2$discharge)
lines(seq(0,2, by = .01), exp(m_coef[1]) * seq(0,2, by = .01) ^ m_coef[2])

# par(new = T)
# plot(density((levels$UNHC), na.rm = T), col = 3, xlim = c(.2, 1.5), axes = F,
#      xlab = '', ylab = '', main = "")
par(new = T, mfrow = c(1,1))
mtext("Rating Curves with bankfull Manning's discharge measurements", 3, 1)
dev.off()

# Read in data from the Wooden Bridge from Tim Covino:
wb <- read_csv('data/rating_curves/covino_wb_ratingcurve.csv') 
wb <- wb %>% mutate(date = as.Date(Datetime),
              stage_m = level_cm/100,
              discharge = discharge_ls/1000,
              site = 'WB') %>%
       select(date, site, stage_m, discharge, velocity_avg = velocity_ms, width = width_m)

comb_rcs <- bind_rows(nhc, unhc, wb)

build_powerlaw_rc <- function(l, q, site){
  m <- lm(log(q)~log(l))
  coefs <- summary(m)$coefficients[,1]
  out<- data.frame(site = site,
                   a = coefs[1, drop = T],
                   b = coefs[2, drop = T],
                   formula = "log(q) = a + b * log(level)",
                   min_Q = min(q, na.rm = T),
                   max_Q = max(q, na.rm = T),
                   min_l = min(l, na.rm = T),
                   max_l = max(l, na.rm = T),
                   n = length(!is.na(q)),
                   row.names = NULL)
  return(out)
}

ZQdat <- data.frame()
ZQdat <- bind_rows(ZQdat, build_powerlaw_rc(nhc$stage_m, nhc$discharge, "NHC"))
ZQdat$min_Q = 0 
ZQdat$min_l = 0.64
ZQdat <- bind_rows(ZQdat, build_powerlaw_rc(unhc$stage_m, unhc$discharge, "UNHC"))
ZQdat <- bind_rows(ZQdat, build_powerlaw_rc(wb$stage_m, wb$discharge, "WB"))

ZQdat <- ZQdat %>% 
  mutate(level_at_0.025 = exp((log(0.025)-a) / b))
comb_rcs %>% tibble() %>%
  mutate(comp_level = case_when(site == 'NHC'~ stage_m - ZQdat$level_at_0.025[1],
                                site == 'UNHC'~ stage_m - ZQdat$level_at_0.025[2],
                                site == 'WB'~ stage_m - ZQdat$level_at_0.025[3])) %>%
  ggplot(aes(comp_level, discharge, col = site)) +
  geom_point() 

plot(seq(0,2, by = .01) - ZQdat$level_at_0.025[1],
     exp(ZQdat$a[1]) * seq(0,2, by = .01) ^ ZQdat$b[1], type = 'l', 
     xlim = c(0,1.1), ylim = c(0,10), xlab = "relative level", ylab = 'Q')
lines(seq(0,2, by = .01) - ZQdat$level_at_0.025[2], 
      exp(ZQdat$a[2]) * seq(0,2, by = .01) ^ ZQdat$b[2], lty = 3)
lines(seq(0,2, by = .01) - ZQdat$level_at_0.025[3],
      exp(ZQdat$a[3]) * seq(0,2, by = .01) ^ ZQdat$b[3], lty = 2)
legend('topleft', c('NHC', 'WB','UNHC'), lty = c(1,2,3), bty = 'n')
mtext('Rating curves adjusted to relative levels')


# Calculate discharge from rating curves ####
# Q = a * level ^ b
# ZQdat_sp <- read_csv(file="data/siteData/NC_streampulseZQ_data.csv")
# ZQdat <- read_csv("NHC_2019_metabolism/data/rating_curves/modified_ZQ_curves.csv")
qq <- levels %>%
  # filter(DateTime_UTC <= ymd_hms("2020-03-02 00:00:00"))%>%
  mutate(discharge_nhc = exp(ZQdat$a[1] + ZQdat$b[1] * log(NHC)),
         # discharge_nhc = case_when(is.na(NHC) ~ NA_real_,
         #                           NHC <= 0.64 ~ 0,
         #                           NHC <= x2 ~ a * (NHC - 0.64),
         #                           TRUE ~ exp(ZQdat$a[1]) + NHC ^ ZQdat$b[1]),
         discharge_unhc = exp(ZQdat$a[2] + ZQdat$b[2] * log(UNHC)),
         notes_rc = case_when(discharge_nhc > ZQdat$max_Q[1] ~
                                "above NHC rc",
                              discharge_unhc > ZQdat$max_Q[2] ~
                                "above UNHC rc"))

# qq %>%
#   # mutate(discharge_nhc = log(discharge_nhc),
#   #        discharge_unhc = log(discharge_unhc)) %>%
#   select(discharge_nhc,  discharge_unhc) %>%
#   xts(order.by = qq$DateTime_UTC) %>%
#   dygraph() %>% dyRangeSelector()
nnhc <- sum(!is.na(qq$NHC))
nunhc <- sum(!is.na(qq$UNHC))

ZQdat <- data.frame(site = c("NHC","UNHC"),
                    above_rc = 
                      c(sum(qq$NHC > ZQdat$max_l[1], na.rm = T)/nnhc,
                        sum(qq$UNHC > ZQdat$max_l[2], na.rm = T)/nunhc),
                    below_rc = 
                      c(sum(qq$NHC < ZQdat$min_l[1], na.rm = T)/nnhc,
                        sum(qq$UNHC < ZQdat$min_l[2], na.rm = T)/nunhc)) %>% 
  left_join(ZQdat)

write_csv(ZQdat, "data/rating_curves/ZQ_curves_with_mannings.csv")

# This chunk was used to pull out good flow data to use with a machine learning
# approach to modeling discharge at these two sites.
# good_flow_for_modeling_Q <-
#   qq %>% mutate(discharge_nhc = case_when(notes == 'nhc modeled' ~ NA_real_,
#                                         notes_rc == 'above NHC rc' ~ NA_real_,
#                                         TRUE ~ discharge_nhc),
#               discharge_unhc = case_when(notes == 'unhc modeled' ~ NA_real_,
#                                         notes_rc == 'above UNHC rc' ~ NA_real_,
#                                         TRUE ~ discharge_unhc)) %>%
#   dplyr::select(DateTime_UTC, discharge_nhc, discharge_unhc) %>%
#   filter(DateTime_UTC >= ymd_hms('2017-01-01 00:00:00'),
#          DateTime_UTC <= ymd_hms('2019-01-01 00:00:00'))
# write_csv(good_flow_for_modeling_Q, 
#           'data/rating_curves/NHC_UNHC_good_flow_years_for_Q_estimation.csv')
par(mfrow = c(1,1))
plot(qq$discharge_unhc, qq$discharge_nhc,log = "xy",
     xlab = "unhc Q", ylab = "nhc Q", pch = 20)
abline(0,1)
qq %>%
  mutate(month = month(DateTime_EST)) %>% 
  filter(DateTime_UTC > ymd_hms('2017-03-01 00:00:00')) %>%
ggplot(aes(log(discharge_unhc), log(discharge_nhc), col = factor(month))) +
  geom_point() +geom_abline(intercept = 0, slope = 1)


# Interpolate discharge ####
write_csv(qq, "data/rating_curves/NHC_UNHC_Q.csv")
# qq <- read_csv("NHC_2019_metabolism/data/rating_curves/NHC_UNHC_Q.csv", guess_max = 10000)
# plot_pres(qq, "discharge_nhc", "discharge_unhc")

sites <- read_csv("data/siteData/NHCsite_metadata.csv") %>%
  slice(1:7)
newQdat <- qq %>%
  mutate(PM.Q = NA_real_,
         CBP.Q = NA_real_,
         WB.Q = NA_real_,
         WBP.Q = NA_real_,
         PWC.Q = NA_real_) %>%
  select(DateTime_UTC, NHC.Q = discharge_nhc,
         PM.Q, CBP.Q, WB.Q, WBP.Q, PWC.Q, 
         UNHC.Q = discharge_unhc, notes) %>%
  # filter(DateTime_UTC >= ymd_hms("2019-03-01 00:00:00")) %>%
  as.data.frame() %>%
  filter(DateTime_UTC <= ymd_hms("2020-04-01 00:00:00"))

w <- which(!is.na(newQdat$NHC.Q))
newQdat <- newQdat[w[1]:w[length(w)],]
for(i in which(!is.na(newQdat$NHC.Q))){
  if(is.na(newQdat$UNHC.Q[i])){ next }
    df <- data.frame(Q = c(newQdat$NHC.Q[i], NA, NA, NA, NA, NA,
                         newQdat$UNHC.Q[i]), 
                   area = c(sites$ws_area.km2[1:7]),
                   distance = c(sites$distance_m[1:7])) 
  if(newQdat$UNHC.Q[i] <= newQdat$NHC.Q[i]){
    df <- df %>% mutate(Q = na.approx(Q, x = area))
  } else {
    df <- df %>% mutate(Q = na.approx(Q, x = distance))
  }

  newQdat[i, 2:8] <- df$Q
  if(i %% 5000 == 0) { print(i/nrow(newQdat))}
}

write_csv(newQdat, "data/rating_curves/interpolatedQ_allsites.csv")
# plot(newQdat$DateTime_UTC, newQdat$NHC.Q, col = "grey80", type = "l", log = "y")
# lines(newQdat$DateTime_UTC, newQdat$NHC.Q, col = "grey80")
# lines(newQdat$DateTime_UTC, newQdat$PM.Q, col = "grey60")
# lines(newQdat$DateTime_UTC, newQdat$CBP.Q, col = "grey50")
# lines(newQdat$DateTime_UTC, newQdat$WB.Q, col = "grey40")
# lines(newQdat$DateTime_UTC, newQdat$WBP.Q, col = "grey35")
# lines(newQdat$DateTime_UTC, newQdat$UNHC.Q, col = "grey20")

# filter high flow days and impossible values #
deltaQ_max = 2
RC <- ZQdat
nhcQ <- read_csv("data/rating_curves/NHC_UNHC_Q.csv", guess_max = 10000) %>%
  mutate(date = as.Date(with_tz(DateTime_UTC, tz = "EST"))) %>%
  group_by(date) %>%
  summarize(nhc_q = mean(discharge_nhc, na.rm = T),
            unhc_q = mean(discharge_unhc, na.rm = T),
            deltaQ = max(discharge_nhc, na.rm = T)/min(discharge_nhc, na.rm = T),
            maxq = max(discharge_nhc, na.rm = T),
            maxqu = max(discharge_unhc, na.rm = T),
            good_flow = ifelse(#nhc_q <= RC$max_Q[1]*1.1 &
                                 # unhc_q <= RC$max_Q[2]*1.1 &
                                 deltaQ < deltaQ_max, 
                               TRUE, FALSE))


flow_dates <- nhcQ %>%
  select(date, nhc_q, deltaQ, good_flow) %>%
  filter(!is.na(good_flow))

ggplot(flow_dates, aes(date, log(nhc_q), col = good_flow))+ geom_point()
sum(flow_dates$good_flow)/nrow(flow_dates)

write_csv(flow_dates, "data/rating_curves/flow_dates_filter.csv")
