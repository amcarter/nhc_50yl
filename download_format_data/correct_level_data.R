# Correct level data from NHC sites
# This file is for cleaning and correcting water pressure data from NHC sites
# and converting it into water level. Steps are individualized for each of the
# six sites due to different lower bounds, different timing of moving sensors
# due to drying etc and different data gaps but the process is the same:

# 1) Calculate depth from the water pressure and air pressure data
# 2) remove sensor offsets created by gaps  that are less than 3 hrs
#       by vertically snapping across gaps (usually caused by sensor out of water)
# 3) Correct level data iusing field based depth measurements

# Setup ####
library(lubridate)
library(tidyverse)
library(zoo)
library(xts)
library(dygraphs)
library(LakeMetabolizer)
# library(devtools)
# install_github('streampulse/StreamPULSE', dependencies=TRUE)
library(StreamPULSE)

# uncomment this section if not sourcing file from master script
# setwd('C:/Users/alice.carter/git/nhc_50yl/')
# sites <- read_csv("data/siteData/NHCsite_metadata.csv") %>%
# slice(c(1:7))

plot_pres <- function(NHC, waterpres = "level_m", airpres = "waterdepth_m",
                      extra = NA){
  if(is.na(extra)){
    NHC %>% select(all_of(waterpres), all_of(airpres)) %>%
      xts(order.by = NHC$DateTime_UTC) %>%
      dygraph() %>%
      dyRangeSelector()
  } else {
    NHC %>% select(all_of(waterpres), all_of(airpres), all_of(extra)) %>%
      xts(order.by = NHC$DateTime_UTC) %>%
      dygraph() %>%
      dyRangeSelector()
  }
}

# 2. Load Data ####
# # get airpressure data: This only needs to be rerun if the date window has changed
# NOAA_airpres <- StreamPULSE:::FindandCollect_airpres(sites$latitude[1],
#                                                      sites$longitude[1],
#                                                      ymd_hms("2016-07-14 00:00:00"),
#                                                      ymd_hms("2022-01-01 00:00:00"))
# NHC <- read_csv("data/metabolism/raw/NHC_2020-03-20.csv") %>%
#   select(DateTime_UTC, AirPres_kPa) %>%
#   right_join(NOAA_airpres)
# # there is a systematic offsett from the air pressure measured at our site and the NOAA air pressure,
# # here we correct for this:
# NOAA_airpres <- NHC %>%
#   mutate(air_kPa = air_kPa -
#            mean(NHC$air_kPa - NHC$AirPres_kPa, na.rm = T)) %>%
#   select(-AirPres_kPa)
# write_csv(NOAA_airpres, "data/siteData/NOAA_airpres.csv")
NOAA_airpres <- read_csv("data/siteData/NOAA_airpres.csv")

# load field notes
ysi <- read_csv("data/siteData/all_nhc_ysi_data.csv") %>%
  filter(!is.na(Date)) %>%
  mutate(waterdepth_m = waterdepth_cm/100) %>%
  select(site, DateTime_UTC, waterdepth_m, notes)

# load data files
filelist <- list.files("data/metabolism/raw")
dir.create('data/metabolism/corrected_level')
for(f in 1:nrow(sites)){
  ff  <- filelist[grep(paste0('^', sites$sitecode[f], '_'), filelist)]
  if(length(ff)>1) ff <- ff[grep('_2020-03-20.csv', ff)]
  dat <- read_csv(paste0("data/metabolism/raw/",ff),
                  guess_max = 100000) %>%
    left_join(NOAA_airpres, by = "DateTime_UTC")
  if("AirTemp_C" %in% colnames(dat)){
    dat <- select(dat, -AirTemp_C) }
  if("AirPres_kPa" %in% colnames(dat)){
    dat <- select(dat, -AirPres_kPa) }

# Calculate depth from water pressure and add sensor offset
# Depth = pressure_Pa = kg/ms2/(density_kg/m3*gravity_m/s2)
dat <- dat %>%
  mutate(pressure_Pa = (WaterPres_kPa - air_kPa) * 1000,
         level_m = pressure_Pa/(water.density(WaterTemp_C) * 9.8) +
           sites$sensor_offset_m[f]) %>%
  select(-pressure_Pa) %>%
  left_join(ysi[ysi$site == dat$site[1],])

# plot_pres(dat)

  if(dat$site[1] == "NHC"){
    gaps <- rle_custom(is.na(dat$level_m))
    if(gaps$values[1] == 1) { gaps <- gaps[-1,]}
    if(gaps$values[nrow(gaps)] == 1) { gaps <- gaps[-nrow(gaps),]}
    # snap gaps together that are less than 1 day (usually caused by sensor removal)
    tmp <- gaps %>%
      filter(values == 1) %>%
      mutate(datetime = as.Date(dat$DateTime_UTC[starts]),
             starts = starts - 1,
             stops = stops + 1,
             lvl_start = dat$level_m[starts],
             lvl_stop = dat$level_m[stops],
             jump = lvl_stop - lvl_start) %>%
      filter(abs(jump) >= 0.01,
             lengths < 10 | starts == 121060) %>%
      arrange(starts)
    n <- nrow(dat)
    dat$level_m1 <- dat$level_m
    plot_pres( dat, 'level_m', 'waterdepth_m', 'level_m1')
    for(i in 1:(nrow(tmp))){
      if(i == nrow(tmp)){
        dat$level_m[tmp$stops[i]:n] <- dat$level_m[tmp$stops[i]:n] -
          dat$level_m[tmp$stops[i]] + dat$level_m[tmp$starts[i]]} else {
      dat$level_m[tmp$stops[i]:tmp$starts[i+1]] <-
        dat$level_m[tmp$stops[i]:tmp$starts[i+1]] -
        dat$level_m[tmp$stops[i]] + dat$level_m[tmp$starts[i]]}
    }
    # plot_pres(dat)

    # OPTION A:
    # snap each chunk of data between gaps of >3 days to the average of the measured points
    # tmp <- gaps %>%
    #   filter(values == 1,
    #          lengths > 300)%>%
    #   mutate(datetime = as.Date(dat$DateTime_UTC[starts]))
    # for(i in 1:(nrow(tmp)+1)){
    #   if(i == 1){ start = 1 } else {start = tmp$stops[i - 1]}
    #   if(i == (nrow(tmp)+1)){stop = n} else {stop = tmp$starts[i]}
    #   delta = mean((na.approx(dat$level_m[start:stop], na.rm = F) - dat$waterdepth_m[start:stop]),
    #                na.rm = T)
    #   if(is.na(delta)) { next }
    #   dat$level_m[start:stop] <- dat$level_m[start:stop] - delta
    # }
    # OPTION B: adjust as two chunks, before and after a sensor movement
    shiftdate <- ymd_hms('2017-05-21 00:00:00')
    s <- which(dat$DateTime_UTC == shiftdate)
    delta = mean(na.approx(dat$level_m[1:s], na.rm = F) -
                   dat$waterdepth_m[1:s], na.rm = T)
    dat$level_m[1:s] <- dat$level_m[1:s] - delta
    delta = mean(na.approx(dat$level_m[s:n], na.rm = F) -
                   dat$waterdepth_m[s:n], na.rm = T)
    dat$level_m[s:n] <- dat$level_m[s:n] - delta
  }
  if(dat$site[1] != "NHC"){
    nhc <- read_csv("data/metabolism/corrected_level/NHC_lvl.csv") %>%
        select(DateTime_UTC, level_nhc = level_m)
  }

  if(dat$site[1] == "UNHC"){
    dat <- dat %>%
      left_join(nhc)
    gaps <- rle_custom(is.na(dat$level_m))
    if(gaps$values[1] ==1) { gaps <- gaps[-1,]}
    if(gaps$values[nrow(gaps)] == 1) { gaps <- gaps[-nrow(gaps),]}
    tmp <- gaps %>%
        filter(values == 1,
               lengths < 31) %>%
        mutate(datetime = as.Date(dat$DateTime_UTC[starts]),
               starts = starts - 1,
               stops = stops + 1,
               lvl_start = dat$level_m[starts],
               lvl_stop = dat$level_m[stops],
               jump = lvl_stop - lvl_start) %>%
        filter(abs(jump) >= 0.01) %>%
        slice(-36,-42)%>%
        arrange(starts)
    n <- nrow(dat)
    # dat$level_m1 <- dat$level_m

    for(i in 1:(nrow(tmp))){
      if(i == nrow(tmp)){
        dat$level_m[tmp$stops[i]:n] <- dat$level_m[tmp$stops[i]:n] -
          dat$level_m[tmp$stops[i]] + dat$level_m[tmp$starts[i]]} else {
            dat$level_m[tmp$stops[i]:tmp$starts[i+1]] <-
              dat$level_m[tmp$stops[i]:tmp$starts[i+1]] -
              dat$level_m[tmp$stops[i]] + dat$level_m[tmp$starts[i]]}
    }
    tmp <- gaps %>%
        filter(values == 1,
               # lengths > 300|starts ==27929) %>%
               lengths > 300|starts %in% c(27929, 44759)) %>%
               # lengths > 300|starts %in% c(27929, 32816, 44759)) %>%
        slice(-6)
    for(i in 1:(nrow(tmp)+1)){
        if(i == 1){ start = 1 } else {start = tmp$stops[i - 1]}
        if(i == (nrow(tmp)+1)){stop = n} else {stop = tmp$starts[i]}
        delta = mean((na.approx(dat$level_m[start:stop], na.rm = F) - dat$waterdepth_m[start:stop]),
                     na.rm = T)
        if(is.na(delta)) { next }
        dat$level_m[start:stop] <- dat$level_m[start:stop] - delta
    }
    # plot_pres( dat, 'level_m', 'waterdepth_m', 'level_nhc')
    # plot_pres( dat, 'level_m', 'waterdepth_m', 'level_m1')
    dat <- select(dat, -level_nhc)
  }
  if(dat$site[1] == "PM"){
    dat <- dat %>%
      left_join(nhc)
    # plot_pres(dat, "level_m", "waterdepth_m", "level_nhc")
    gaps <- rle_custom(is.na(dat$level_m))
    if(gaps$values[1] ==1) { gaps <- gaps[-1,]}
    if(gaps$values[nrow(gaps)] == 1) { gaps <- gaps[-nrow(gaps),]}
    tmp <- gaps %>%
      filter(values == 1,
             lengths < 10) %>%
      mutate(datetime = as.Date(dat$DateTime_UTC[starts]),
             starts = starts - 1,
             stops = stops + 1,
             lvl_start = dat$level_m[starts],
             lvl_stop = dat$level_m[stops],
             jump = lvl_stop - lvl_start) %>%
      filter(abs(jump) >= 0.01) %>%
      arrange(starts)
    n <- nrow(dat)
    for(i in 1:(nrow(tmp))){
      dat$level_m[tmp$stops[i]:n] <-  dat$level_m[tmp$stops[i]:n] - tmp$jump[i]
    }
    # snap each chunk of data between gaps of >3 days to the average of the measured points
    tmp <- gaps %>%
      filter(values == 1,
             lengths > 300)
    for(i in 1:(nrow(tmp)+1)){
      if(i == 1){ start = 1 } else {start = tmp$stops[i - 1]}
      if(i == (nrow(tmp)+1)){stop = n} else {stop = tmp$starts[i]}
      delta = mean((na.approx(dat$level_m[start:stop], na.rm = F) - dat$waterdepth_m[start:stop]),
                   na.rm = T)
      if(is.na(delta)) { next }
      dat$level_m[start:stop] <- dat$level_m[start:stop] - delta
    }
    plot_pres(dat, "level_nhc", "waterdepth_m", "level_m")
    dat <- select(dat, -level_nhc)
  }
  if(dat$site[1] == "CBP"){
    dat <- dat %>%
      left_join(nhc)
    # plot_pres(dat, "level_m", "waterdepth_m", "level_nhc")

    gaps <- rle_custom(is.na(dat$level_m))
    if(gaps$values[1] ==1) { gaps <- gaps[-1,]}
    if(gaps$values[nrow(gaps)] == 1) { gaps <- gaps[-nrow(gaps),]}
    tmp <- gaps %>%
      filter(values == 1,
             lengths < 10|
                 starts == which(dat$DateTime_UTC ==
                                     ymd_hms("2019-10-05 01:00:00"))) %>%
      mutate(datetime = as.Date(dat$DateTime_UTC[starts]),
             starts = starts - 1,
             stops = stops + 1,
             lvl_start = dat$level_m[starts],
             lvl_stop = dat$level_m[stops],
             jump = lvl_stop - lvl_start) %>%
      filter(abs(jump) >= 0.01) %>%
      arrange(starts)
    n <- nrow(dat)
    for(i in 1:(nrow(tmp))){
      dat$level_m[tmp$stops[i]:n] <-  dat$level_m[tmp$stops[i]:n] - tmp$jump[i]
    }
    # plot_pres(dat, "level_m", "waterdepth_m", "level_nhc")

    # snap each chunk of data between gaps of >3 days to the average of the measured points
    tmp <- gaps %>%
      filter(values == 1,
             lengths > 300)
    for(i in 1:(nrow(tmp)+1)){
      if(i == 1){ start = 1 } else {start = tmp$stops[i - 1]}
      if(i == (nrow(tmp)+1)){stop = n} else {stop = tmp$starts[i]}
      delta = mean((na.approx(dat$level_m[start:stop], na.rm = F) - dat$waterdepth_m[start:stop]),
                   na.rm = T)
      if(is.na(delta)) { next }
      dat$level_m[start:stop] <- dat$level_m[start:stop] - delta
    }
    dat <- select(dat, -level_nhc)
  }
  if(dat$site[1] == "WB"){
    dat <- dat %>%
      left_join(nhc)
    # plot_pres(dat, "level_m", "waterdepth_m", "level_nhc")
    gaps <- rle_custom(is.na(dat$level_m))
    if(gaps$values[1] ==1) { gaps <- gaps[-1,]}
    if(gaps$values[nrow(gaps)] == 1) { gaps <- gaps[-nrow(gaps),]}
    tmp <- gaps %>%
      filter(values == 1,
             lengths < 10|
                 starts == which(dat$DateTime_UTC ==
                                     ymd_hms("2019-04-15 18:30:00"))|
                 starts == which(dat$DateTime_UTC ==
                                     ymd_hms("2019-05-30 21:45:00"))|
                 starts == which(dat$DateTime_UTC ==
                           ymd_hms("2020-02-27 17:45:00"))) %>%
      mutate(datetime = as.Date(dat$DateTime_UTC[starts]),
             starts = starts - 1,
             stops = stops + 1,
             lvl_start = dat$level_m[starts],
             lvl_stop = dat$level_m[stops],
             jump = lvl_stop - lvl_start) %>%
      filter(abs(jump) >=0.01) %>%
      arrange(starts)
    n <- nrow(dat)
    for(i in 1:(nrow(tmp))){
      dat$level_m[tmp$stops[i]:n] <-  dat$level_m[tmp$stops[i]:n] - tmp$jump[i]
    }
    # snap each chunk of data between gaps of >3 days to the average of the measured points
    delta = mean((na.approx(dat$level_m, na.rm = F) - dat$waterdepth_m),
                 na.rm = T)
    dat$level_m <- dat$level_m - delta
    dat <- select(dat, -level_nhc)
  }
  if(dat$site[1] == "WBP"){
    dat <- dat %>%
      left_join(nhc)
    # plot_pres(dat, "level_m", "waterdepth_m", "level_nhc")
    gaps <- rle_custom(is.na(dat$level_m))
    if(gaps$values[1] ==1) { gaps <- gaps[-1,]}
    if(gaps$values[nrow(gaps)] == 1) { gaps <- gaps[-nrow(gaps),]}
    tmp <- gaps %>%
      filter(values == 1,
             lengths < 10| starts == 5953) %>%
      mutate(datetime = as.Date(dat$DateTime_UTC[starts]),
             starts = starts - 1,
             stops = stops + 1,
             lvl_start = dat$level_m[starts],
             lvl_stop = dat$level_m[stops],
             jump = lvl_stop - lvl_start) %>%
      filter(abs(jump) >= 0.01) %>%
      slice(-10,-11) %>%
      arrange(starts)
    n <- nrow(dat)
    for(i in 1:(nrow(tmp))){
      dat$level_m[tmp$stops[i]:n] <-  dat$level_m[tmp$stops[i]:n] - tmp$jump[i]
    }

    # plot_pres(dat, "level_nhc", "waterdepth_m", "level_m")
    # snap each chunk of data between gaps of >3 days to the average of the measured points
    tmp <- gaps %>%
      filter(values == 1,
             lengths > 100)
    for(i in 1:(nrow(tmp)+1)){
      if(i == 1){ start = 1 } else {start = tmp$stops[i - 1]}
      if(i == (nrow(tmp)+1)){stop = n} else {stop = tmp$starts[i]}
      delta = mean((na.approx(dat$level_m[start:stop], na.rm = F) - dat$waterdepth_m[start:stop]),
                   na.rm = T)
      if(is.na(delta)) { next }
      dat$level_m[start:stop] <- dat$level_m[start:stop] - delta
    }
    dat <- select(dat, -level_nhc)
  }
  if(dat$site[1] == "PWC"){
    dat <- dat %>%
      left_join(nhc)
    plot_pres(dat, "level_m", "waterdepth_m", "level_nhc")
    dat$level_m[dat$level_m < 0.69] <- NA
    tmp <- data.frame(starts = c(which(dat$DateTime_UTC ==
                                         ymd_hms("2019-05-30 21:00:00")),
                                 which(dat$DateTime_UTC ==
                                         ymd_hms("2019-07-19 17:00:00"))
    ))
    tmp$stops <- tmp$starts + c( 1,2)

    gaps <- rle_custom(is.na(dat$level_m))
    if(gaps$values[1] ==1) { gaps <- gaps[-1,]}
    if(gaps$values[nrow(gaps)] == 1) { gaps <- gaps[-nrow(gaps),]}
    tmp <- gaps %>%
      filter(values == 1) %>%
      slice(1,3,5) %>%
      bind_rows(tmp) %>%
      mutate(datetime = as.Date(dat$DateTime_UTC[starts]),
             starts = starts - 1,
             stops = stops + 1,
             lvl_start = dat$level_m[starts],
             lvl_stop = dat$level_m[stops],
             jump = lvl_stop - lvl_start) %>%
      arrange(starts)
    n <- nrow(dat)
    # dat$level_m1 -> dat$level_m
    dat$level_m1 <- dat$level_m
    for(i in 1:(nrow(tmp))){
      dat$level_m[tmp$stops[i]:n] <-  dat$level_m[tmp$stops[i]:n] - tmp$jump[i]
    }

    plot_pres(dat, "level_nhc", "waterdepth_m", "level_m")

    delta = mean((na.approx(dat$level_m, na.rm = F) - dat$waterdepth_m),
                 na.rm = T)
    dat$level_m <- dat$level_m - delta

    dat <- select(dat, -level_nhc)
  }

  # plot_pres(dat, "level_m", 'waterdepth_m', "level_m1")
  dat <- dat %>%
    mutate(level_m = na.approx(level_m, maxgap = 96*3, na.rm = F),
           site = dat$site[1]) %>%
    select(-WaterPres_kPa, -waterdepth_m, -notes) %>%
    group_by(DateTime_UTC, site) %>%
    summarize_all(mean, na.rm = T) %>%
    ungroup()

write_csv(dat, paste0("data/metabolism/corrected_level/",
                      dat$site[1], "_lvl.csv"))

}

# 3. compile all levels ####
filelist <- list.files("data/metabolism/corrected_level/")
dd <- data.frame()
for(f in 1:length(filelist)){
  d <- read_csv(paste0("data/metabolism/corrected_level/",
                       filelist[f]), guess_max = 10000) %>%
  select(DateTime_UTC, level_m, site)
  dd <- bind_rows(dd, d)
}

write_csv(dd, "data/rating_curves/all_sites_level_corrected.csv")
