
library(tidyverse)
library(googledrive)
library(rgee)
library(sf)
library(jsonlite)
library(data.table)
library(lubridate)
library(mapview)
library(ggplot2)
library(magick)

setwd('~/git/papers/alice_nhc/data/gee/')

# setup ####

conf <- jsonlite::fromJSON('~/git/macrosheds/data_acquisition/config.json',
                           simplifyDataFrame = FALSE)
googledrive::drive_auth(email = conf$gee_login_mike)
rgee::ee_Initialize(user = conf$gee_login_ms, drive = TRUE)


# prep ####

#load shapes
nhc_wbd <- st_read('../watershed_boundary/NHC.shp', quiet = TRUE)

#create new asset folder
user_info <- rgee::ee_user_info(quiet = TRUE)
asset_folder <- file.path(user_info$asset_home, 'NHC')
rgee::ee_manage_create(asset_folder)

#upload shapes
sf_as_ee(nhc_wbd,
         via = 'getInfo_to_asset',
         assetId = file.path(asset_folder, 'nhc_wbd'),
         overwrite = FALSE,
         quiet = TRUE)

asset_path <- rgee::ee_manage_assetlist(asset_folder)
asset <- ee$FeatureCollection(asset_path$ID[1])

# compute ####

gee_id = 'ECMWF/ERA5_LAND/DAILY_AGGR'
band = 'temperature_2m'
res = 11132

imgcol <- ee$ImageCollection(gee_id)$select(band)

flat_img <- imgcol$map(function(image){
    image$reduceRegions(
        collection = asset,
        reducer = ee$Reducer$median(),
        scale = res
    )
})$flatten()

ee_task <- ee$batch$Export$table$toDrive(
    collection = flat_img,
    description = 'air_temp_reanalysis',
    fileFormat = 'CSV',
    folder = 'GEE',
    fileNamePrefix = 'NHC_air_temp_reanalysis'
)

ee_task$start()
ee_monitoring(ee_task, max_attempts = Inf, quiet = TRUE)

# retrieve, munge ####

gee_dl <- '~/ssd2/NHC_air_temp.csv'
googledrive::drive_download(file = 'GEE/NHC_air_temp_reanalysis.csv', gee_dl)

d <- read_csv(gee_dl, show_col_types = FALSE) %>%
    mutate(date = ymd(substr(`system:index`, 1, 8))) %>%
    select(FID, median, date)

d <- d %>%
    # mutate(median = median * 0.0001) %>%  #scale factor of gee product
    mutate(median = median -273.15) %>%
    # group_by(FID, date) %>%
    # summarize(median = sum(median, na.rm = TRUE),
    #           count = n(),
    #           .groups = 'drop') %>%
    # mutate(median = (median / (count * 16)) * 365) %>% #if n-daily aggregating to annual
    rename(temp_C = median) %>%
    select(date, temp_C)

write_csv(d, 'historical_air_temperature.csv')
