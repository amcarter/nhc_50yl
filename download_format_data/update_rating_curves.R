# Update rating curve file with new discharge measurement
source("C:/Users/Alice Carter/git/nhc_50yl/src/download_format_data/calc_discharge_from_crosssection.r")
setwd('C:/Users/Alice Carter/git/nhc_50yl/')

# # 1. initial calcs to create rc_sheet ####
# only run this section if you have corrected/updated all previous data:
# load and munge data
# qdat <- read_csv("data/rating_curves/discharge_measurements.csv")
# sp_qdat <- read_csv("data/rating_curves/discharge_measurements_NHC_UNHC.csv")
# metadat <- read_csv("data/siteData/NHCsite_metadata.csv")
#  
# # # setting negative V equal to zero to be consistent with E Moore's data
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
# qq[,-5]
# 
# # 3. run on old SP data for NHC and UNHC ####
# #qq <- data.frame()
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
# qq[,-5]
# dev.off()
# qq <- metadat %>%
#   select(site = sitecode, location_m = distance_m) %>%
#   right_join(qq, by = "site")
# 
# ggplot(qq, aes(stage_m, discharge)) +
#   geom_point() +
#   facet_wrap(.~site)
# write_csv(qq, "data/rating_curves/calculated_discharge.csv")
# # Add stages ####
# qq <- read_csv("data/rating_curves/calculated_discharge.csv")
# qq <- qq %>%
#   mutate(DateTime_EST = round_date(ymd_hms(paste(date, time), tz = "EST"),
#                                    unit = "15 minutes"))
# 
# levels <- read_csv("data/rating_curves/all_sites_level_corrected.csv") %>%
#   pivot_wider(names_from = "site", values_from = "level_m") %>%
#   mutate(DateTime_EST = with_tz(DateTime_UTC, tz = "EST")) %>%
#   select(-notes)
#   # select(DateTime_EST, NHC_level = NHC, UNHC_level = UNHC)
# 
# qql <- qq %>%
#   left_join(levels)
# 
# # inspect data 
# qql %>%
#   filter(site %in% c( "NHC")) %>%
#   ggplot(aes(stage_m, discharge)) +
#   geom_point() +
#   facet_wrap(.~site)
# 
# qql %>%
#   filter(site == "UNHC") %>%
#   ggplot(aes(stage_m, discharge)) +
#   geom_point(aes(color = UNHC), size = 4)
# 
# # Assign proper discharges to the pool sections on days when a pool and
# # neighboring riffle were measured
# 
# # PM: that q is just wrong, but we didn't measure a riffle. Will pair with
# # interpolated Q to see what we get
# 
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
# qql <- qql %>%
#   mutate(DateTime_UTC = with_tz(DateTime_EST, tz = "UTC")) %>%
#   mutate(level = case_when(site == "NHC" ~ NHC,
#                            site == "PM" ~ PM,
#                            site == "CBP" ~ CBP,
#                            site == "WB" ~ WB,
#                            site == "WBP" ~ WBP,
#                            site == "UNHC" ~ UNHC
#                            )) %>%
#   select(-(15:21))
# 
# # Qdat <- read_csv("data/rating_curves/interpolatedQ_allsites.csv",
# #                  guess_max = 10000) %>%
# #   select(-notes)
# # qql <- qql %>%
# #   left_join(Qdat, by = "DateTime_UTC") %>%
# #   mutate(discharge_rc = case_when(site == "NHC" ~ NHC.Q,
# #                            site == "PM" ~ PM.Q,
# #                            site == "CBP" ~ CBP.Q,
# #                            site == "WB" ~ WB.Q,
# #                            site == "WBP" ~ WBP.Q,
# #                            site == "UNHC" ~ UNHC.Q
# #                            )) %>%
# #   select(-(16:22))
# 
# write_csv(qql, "data/rating_curves/calculated_discharge_with_levels.csv")


# 
# Add a zero flow point for the NHC Q curve ####
# m_diff = 0.04436296 

qql <- read_csv("data/rating_curves/calculated_discharge_with_levels.csv")
levels <- read_csv("data/rating_curves/all_sites_level_corrected.csv",
                   guess_max = 1000000) #%>%
  # mutate(level_m = case_when(site == "UNHC" ~ level_m - m_diff,
  #                            TRUE ~ level_m))
# write_csv(levels, "NHC_2019_metabolism/data/rating_curves/all_sites_level_corrected.csv")

# unhc <- read_csv("NHC_2019_metabolism/data/metabolism/corrected_level/UNHC_lvl.csv") %>%
#   mutate(level_m = level_m - m_diff)

# write_csv(unhc, "NHC_2019_metabolism/data/metabolism/corrected_level/UNHC_lvl.csv")

levels <- levels %>%  
  filter(site %in% c('NHC', 'UNHC')) %>%
  pivot_wider(names_from = "site", values_from = "level_m") %>%
  mutate(DateTime_EST = with_tz(DateTime_UTC, tz = "EST")) %>%
  arrange(DateTime_EST) 

# levels %>%
#   select(-DateTime_EST) %>%
#   xts(order.by = levels$DateTime_EST) %>%
#   dygraph() %>%
#   dyRangeSelector()
# Zero flow occurs at a stage of 0.64 at NHC (by observation),

# Build Rating Curves ####
# qql$stage_m[qql$site =="UNHC"][2] <- qql$UNHC[qql$site =="UNHC"][2]
nhc <- qql %>%
  filter(site =="NHC")

# nhc <- data.frame(stage_m = .64, discharge = 0.005) %>% bind_rows(nhc)
unhc <- qql %>%
  filter(site =="UNHC") %>%
  slice(c( -3, -2))
# unhc$discharge[1]<- c(.001)
# unhc <- data.frame(stage_m = .35, discharge = 0.01) %>% bind_rows(unhc)
# 
# m <- lm(log(discharge)~log(stage_m), data = nhc)
# m_coef <- summary(m)$coefficients[,1]
# 
# # What does this curve predict at 0.62 stage?
# plot(nhc$stage_m, nhc$discharge, xlim = c(.5, 1.5), ylim = c(0,2))
# lines(seq(0,2, by = .01), exp(m_coef[1]) * seq(0,2, by = .01) ^ m_coef[2], col = 2)
# med <- median(levels$NHC, na.rm = T)
# nhc_medq <- exp(m_coef[1]) * med ^ m_coef[2]
# abline(h = nhc_medq, v = med)
# par(new = T)
# plot(density((levels$NHC), na.rm = T), col = 3, xlim = c(.5, 1.5), axes = F)
# 
# #calculate expected median UNHC by watershed area
# unhc_medq <- nhc_medq *sites$ws_area.km2[7]/sites$ws_area.km2[1]
# 
# m <- lm(log(discharge)~log(stage_m), data = unhc)
# m_coef <- summary(m)$coefficients[,1]
# plot(unhc$stage_m, (unhc$discharge), xlim = c(.2, .8), ylim = c(0,2))
# lines(seq(0,2, by = .01), exp(m_coef[1]) * seq(0,2, by = .01) ^ m_coef[2], col = 2)
# u_med <- median(levels$UNHC, na.rm = T)
# med <- (unhc_medq/exp(m_coef[1])) ^ (1/m_coef[2])
# abline(h = unhc_medq, v = med)
# 
# par(new = T)
# plot(density((levels$UNHC), na.rm = T), col = 3, xlim = c(.2, .8), axes = F)
# # med_diff <- u_med-med

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


# Calculate discharge from rating curves ####
# Q = a * level ^ b
# ZQdat_sp <- read_csv(file="siteData/NC_streampulseZQ_data.csv")
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

write_csv(ZQdat, "data/rating_curves/modified_ZQ_curves.csv")

# par(mfrow = c(1,1))
# plot(qq$discharge_unhc, qq$discharge_nhc,log = "xy",
#      xlab = "unhc Q", ylab = "nhc Q", pch = 20)
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
            good_flow = ifelse(nhc_q <= RC$max_Q[1] &
                                 unhc_q <= RC$max_Q[2] &
                                 deltaQ < deltaQ_max, 
                               TRUE, FALSE))


flow_dates <- nhcQ %>%
  select(date, nhc_q, deltaQ, good_flow) %>%
  filter(!is.na(good_flow))

write_csv(flow_dates, "data/rating_curves/flow_dates_filter.csv")
