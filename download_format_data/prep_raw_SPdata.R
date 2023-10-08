# Process raw NHC datafiles:
# files downloaded from SP portal
# prep for running metabolism models
# remotes::install_github("appling/unitted")
# remotes::install_github("USGS-R/streamMetabolizer")
library(streamMetabolizer)
library(LakeMetabolizer)
# library(lubridate)
# library(tidyverse)
# library(dygraphs)
# library(xts)
# library(zoo)
#
# setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data/")
# source("../src/helpers.R")

# 1. load raw data files and metadata ####
# sites <- read_csv("siteData/NHCsite_metadata.csv")
ll <- read_csv("data/rating_curves/all_sites_level_corrected.csv",
               guess_max = 1000000) %>%
  pivot_wider(names_from = site,  values_from = level_m)%>%
  arrange(DateTime_UTC)

Qdat <- read_csv("data/rating_curves/interpolatedQ_allsites.csv",
                 guess_max = 10000)
ZQdat <- read_csv('data/rating_curves/ZQ_curves_with_mannings.csv')
# filter out time points that are above the rating curve:
above_RC <- ll %>%
    mutate(NHC_above = case_when(NHC > ZQdat$max_l[1] ~ TRUE,
                                 TRUE ~ FALSE),
           UNHC_above = case_when(UNHC > ZQdat$max_l[2] ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    select(DateTime_UTC, ends_with('above'))

Qdat_RCcorr <- Qdat %>%
    left_join(above_RC, by = 'DateTime_UTC') %>%
    mutate(NHC.Q = case_when(NHC_above ~ NA_real_,
                             TRUE ~ NHC.Q),
           UNHC.Q = case_when(UNHC_above ~ NA_real_,
                              TRUE ~ UNHC.Q),
           across(any_of(c('PM.Q', 'CBP.Q','WB.Q', 'WBP.Q', 'PWC.Q')),
                  .fns = ~ case_when(NHC_above|UNHC_above ~ NA_real_,
                                   TRUE ~ .x)))
# nhcQ <- read_csv("rating_curves/NHC_UNHC_Q.csv", guess_max = 10000) %>%
#   select(DateTime_UTC, NHC.Q = NHC_Q, UNHC.Q = UNHC_Q, AirPres_kPa) %>%
#   full_join(Qdat) %>%
#   arrange(DateTime_UTC)
filelist <- list.files("data/metabolism/corrected_level")

# get rid of MC751 and Mud for now
# filelist <- filelist[-c(2,3)]

#look at required inputs for a bayesian model in stream Metabolizer
# metab_inputs(type="bayes", input="data")

# 2. Pair with interpolated Q ####
get_Q <- function(dat, Qdat){
  Qname <- paste(dat$site[1], "Q", sep=".")
  Q <- Qdat %>% select (DateTime_UTC,
                        discharge = all_of(Qname)) %>%
    right_join(dat, by="DateTime_UTC") %>%
    arrange(DateTime_UTC)

  return(Q)
}

# 3. Estimate depth based on Leopold and Maddock/Raymond.  ####
# DQ <- data.frame(sitename = c("NHC","PM","CBP","WB","WBP","PWC","UNHC"),
#                  c_m = rep(0.409, 7),   # depth at unit discharge
#                  f = rep(.294, 7))      # exponent in depth discharge relation
# coefficients are default values for Leopold and Maddock 1953 equation D=cQ^f
# defined in Raymond 2012.
# load dataframe with calibrated parameters for each site:

DQ <- read_csv("data/rating_curves/depth_discharge_relationship_LM1953.csv")
VQ <- read_csv("data/rating_curves/velocity_discharge_fit_relationship.csv")

# Calculate velocity from empirical relationship
calc_velocity <- function(dat, VQ){
  tmp <- dat %>%
    bind_rows(VQ) %>%
    arrange(discharge) %>%
    mutate(avg_velocity = na.approx(avg_velocity, x = discharge, na.rm = F)) %>%
    filter(!is.na(DateTime_UTC)) %>%
    arrange(DateTime_UTC)

  return(tmp)
}


# 4. Calculate Level data from water pressure ####
# Depth = pressure_Pa/density/acelleration due to gravity =
#   P (Pa) = ro (kg/m3) * gravity (m/s2) * depth (m)
# density is temperature dependent, this equation in lake metabolizer
# accounts for that based on Martin and McCutcheon 1999
calc_water_level<- function(dat, sites){

  sensor_offset <- sites$sensor_offset_m[sites$sitecode==dat$site[1]]
  ro <- water.density(dat$WaterTemp_C)
  pressure_Pa <- (dat$WaterPres_kPa-dat$AirPres_kPa)*1000
  level_m <- sensor_offset + pressure_Pa/(ro * 9.8)

  return(level_m)
}

# if rerunning for a site already corrected level, use this function
get_l <- function(dat, ll){
  Qname <- dat$site[1]
  L <- ll %>% select (DateTime_UTC, level_m = all_of(Qname)) %>%
    right_join(dat, by="DateTime_UTC") %>%
    arrange(DateTime_UTC)

  return(L)
}

# 5. Function to prepare datafile for drift correction ####

prep_file <- function(filename, sites, Qdat, DQ, VQ){
  dat <- read_csv(paste0("data/metabolism/corrected_level/",
                         filename), guess_max = 10000)
  lat <- sites$latitude[sites$sitecode == dat$site[1]]
  lon <- sites$longitude[sites$sitecode == dat$site[1]]

  # remove leading and ending NAs
  w <- which(!is.na(dat$DO_mgL))
  dat <- dat[min(w):max(w), ]
  dates <- data.frame(DateTime_UTC = seq(min(dat$DateTime_UTC, na.rm = T),
                                         ymd_hms("2020-03-26 20:00:00"),
                                         by = "15 min"))
  dat <- dates %>%
    left_join(dat)
  #Load discharge data
  dat <- get_Q(dat, Qdat)

  # Calculate DO saturation, depth, light
  dat$AirPres_mbar <- dat$air_kPa*10
  dat$DO.sat <- calc_DO_sat(dat$WaterTemp_C, dat$AirPres_mbar,
                            salinity.water = 0, model = "garcia-benson")

  dat$depth <- calc_depth(dat$discharge,
                          c = DQ$c_m[DQ$sitename == dat$site[1]],
                          f = DQ$f[DQ$sitename == dat$site[1]])

  dat <- calc_velocity(dat, VQ)
  # Calculate Level data from water pressure
  # if(corrected_level == T){
  #   dat <- get_l(dat, ll)
  # } else {
  #   dat$level_m <- calc_water_level(dat, sites)
  # }

  dates <- data.frame(DateTime_UTC = seq(min(dat$DateTime_UTC, na.rm = T),
                                         max(dat$DateTime_UTC, na.rm = T),
                                         by = "15 min"))

  dat <- dates %>%
    left_join(dat) %>%
    # select(-site)%>%
    mutate(across(-c("DateTime_UTC", "site"),
                  na.approx, na.rm = F, maxgap = 12)) %>%
    as_tibble()

  #Convert datetime to solar time
  dat$DateTime_EST <- with_tz(dat$DateTime_UTC, tz="EST") # convert to EST timezone
  dat$solar.time <- calc_solar_time(dat$DateTime_EST, longitude=lon)
  dat$light <- calc_light(dat$solar.time, latitude=lat,longitude=lon)

  # rename variables needed for metabolism model
  dat <- dat %>%
    select(-DateTime_EST, -AirPres_mbar, -air_kPa, -air_temp) %>%
    rename(DO.obs = DO_mgL,
           temp.water = WaterTemp_C)

  return(dat)
}


# run with corrected level
dir.create('data/metabolism/processed')
for(i in 1:length(filelist)){
  filename <- filelist[i]
  dat <- prep_file(filename, sites, Qdat_RCcorr, DQ, VQ)
  write_csv(dat, paste0("data/metabolism/processed/", dat$site[1], ".csv"))
}

# # 6. Drift correction functions based on YSI data ####
# # load ysi data
# ysi <- read_csv("data/siteData/all_nhc_ysi_data.csv") %>%
#     mutate(ysi_level = waterdepth_cm/100) %>%
#     select(site, DateTime_UTC, time, ysi_temp = watertemp_C,
#            ysi_DO = DO_mgL, ysi_level) %>%
#     filter(!is.na(site), !is.na(DateTime_UTC))
#
# correct_file <- function(dat, ysi){
#
#     ysid <- ysi %>%
#         filter(site == dat$site[1]) %>%
#         right_join(dat, by = c('site', 'DateTime_UTC')) %>%
#         arrange(DateTime_UTC)
#
#     yd <- ysid %>%
#         select(DateTime_UTC, DO.obs, ysi_DO,
#                temp.water, ysi_temp,
#                level_m, ysi_level) %>%
#         mutate(DO.obs = na.approx(ifelse(!is.na(ysi_DO), NA, DO.obs),
#                               na.rm = F, maxgap = 12),
#                temp.water = na.approx(ifelse(!is.na(ysi_temp), NA, temp.water),
#                                   na.rm = F, maxgap = 12),
#                level_m = na.approx(ifelse(!is.na(ysi_temp), NA, level_m),
#                                na.rm = F, maxgap = 12))
#
#     dat$DO.corr <- yd$DO.corr <- drift_correct(yd, "DO.obs", "ysi_DO")
#     dat$T.corr <- yd$T.corr <- drift_correct(yd, "temp.water", "ysi_temp")
#
#     plot_drift_corr(yd, 'DO.obs', 'ysi_DO')#, 'DO.corr')
#     plot_drift_corr(yd, 'temp.water', 'T.corr', 'ysi_temp')
#
#     return(dat)
# }
#
# # 7. Run for files and check output ####
# # for(i in 1:length(filelist)){
#
# i=2
#
#
# filename <- filelist[i]
# dat <- prep_file(filename, sites, Qdat, DQ)
# dat <- correct_file(dat, ysi)
# ysi %>% filter(site == dat$site[1]) %>%
#   right_join(dat) %>%
#   arrange(DateTime_EST) %>%
#   select(level_m, ysi_level) %>%
#   xts(order.by = dat$DateTime_EST) %>%
#   dygraph() %>%
#   dyRangeSelector()
# check corrections
# cor <- ysi %>%
#   filter(site == dat$site[1]) %>%
#   right_join(dat) %>%
#   arrange(DateTime_EST)
#
# corDO <- cor %>% select(DO.obs, DO.corr, ysi_DO) %>%
#   xts(order.by = dat$DateTime_EST)
# show(dygraph(corDO) %>% dyRangeSelector())
#
# corT <- cor %>%
#   select(temp.water, T.corr, ysi_temp) %>%
#   xts(order.by = dat$DateTime_EST)
# show(dygraph(corT) %>% dyRangeSelector())
#
# df <- dat %>%
#   mutate(DO.obs = DO.corr,
#          temp.water = T.corr) %>%
#   select(-DO.corr, -T.corr)
#
# write_csv(df, paste0("metabolism/processed/drift_corrected/", dat$site[1],"_dfc.csv"))
#
