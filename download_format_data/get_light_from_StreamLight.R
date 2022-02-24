# Get light data using the StreamLight package from P Savoy

# setup ####
#Use the devtools packge to install StreamLightUtils
# devtools::install_github("psavoy/StreamLightUtils")
# devtools::install_github("psavoy/StreamLight")

library(StreamLightUtils)
library(StreamLight)
library(lubridate)
library(tidyverse)
setwd('C:/Users/Alice Carter/git/nhc_50yl')
source('data/light/modified_streamLight_functions.R')
sitedat <- read_csv('data/siteData/NHCsite_metadata.csv') %>%
  rename(Site_ID = sitecode, Lat = latitude, Lon = longitude,
         startDate = startdate.UTC)
# Download and Process NLDAS data for incoming radiation ####
# Set the download location (add your own directory)
working_dir <- "C:/Users/Alice Carter/git/nhc_50yl/data/light/NLDAS"

# # for one site:
# #Download NLDAS data at NC_NHC
# NLDAS_DL(
#   save_dir = working_dir,
#   Site_ID = "NC_NHC",
#   Lat = 35.9925, 
#   Lon = -79.0460, 
#   startDate = "2017-01-01"
# )
# #Process the downloaded data
# NLDAS_processed <- NLDAS_proc(
#   read_dir = working_dir, 
#   Site_IDs = "NC_NHC"
# )

# for multiple sites
#Read in a table with initial site information
sites <- sitedat %>%
  slice(c(1:5,7)) %>%
  select(Site_ID, Lat, Lon, startDate) %>%
  mutate(startDate = as.character(date(startDate)),
         epsg_crs = rep(4326, 6)) %>%
  data.frame()

#Download NLDAS data at NC_NHC
NLDAS_DL_bulk(
  save_dir = working_dir,
  site_locs = sites
)

#List of successfully downloaded sites
NLDAS_list <- stringr::str_sub(list.files(working_dir), 1, -11)

#Processing the downloaded NLDAS data
NLDAS_processed <- StreamLightUtils::NLDAS_proc(read_dir = working_dir, NLDAS_list)

# Download and process MODIS LAI ####
#Make a table for the MODIS request 
request_sites <- sites[, c("Site_ID", "Lat", "Lon")] 

#Export your sites as a .csv for the AppEEARS request  
write.table(
  request_sites, 
  "data/light/NC_sites.csv", 
  sep = ",", 
  row.names = FALSE,
  quote = FALSE, 
  col.names = FALSE
)

# Save the zip file downloaded from AppEEARS
# Unpack the data after downloading
working_dir <- "C:/Users/Alice Carter/git/nhc_50yl/data/light/MODIS"

MOD_unpack <- AppEEARS_unpack_QC(
  zip_file = "nhc-sites.zip", 
  zip_dir = working_dir, 
  request_sites[, "Site_ID"]
)


# Process the downloaded data
MOD_processed <- AppEEARS_proc2(
  unpacked_LAI = MOD_unpack,  
  fit_method = "Gu", 
  plot = TRUE
)


# Make an input datafile for streamlight ####
working_dir <- "C:/Users/Alice Carter/git/nhc_50yl/data/light/drivers"
make_driver(sites, NLDAS_processed, MOD_processed, 
            TRUE, working_dir)

# add parameters to site file ####
site_parm <- sitedat %>%
  select(Site_ID, Lat, Lon, CRS, Width = width_mar_m)

# add summary data from daily site data
daily <- read_csv('C:/Users/Alice Carter/git/nhc_50yl/data/metabolism/metabolism_and_drivers.csv') 

sum <- daily %>%
  filter(date >= date('2019-03-06') & date < date('2020-03-06')) %>%
  select(Site_ID = site, depth) %>%
  group_by(Site_ID) %>%
  summarize(WL = median(depth, na.rm = T))

site_parm <- right_join(site_parm, sum, by = 'Site_ID')

# Calculate site azimuth using Google Maps. The angle I am using is the 
# line that connects the sensor location with the point 1000 m upstream.

site_parm <- site_parm %>%
  mutate(epsg_code = rep(4326, 6),
         Azimuth = c(320, 250, 280, 300, 234, 221), 
         BH = rep(0.1, 6),
         BS = rep(100, 6))

# get canopy data:
# extract_height(Site_ID = site_parm[,'Site_ID'], Lat = site_parm[,'Lat'], 
#                Lon = site_parm[,'Lon'], site_crs = site_parm[,'epsg_code'])
# that didn't work
# from streampulse field measurements:D
site_parm <- site_parm %>%
  mutate(TH = rep(20, 6),
         overhang = rep(7, 6),
         overhang_height = rep(13, 6),
         x = rep(1, 6))

# running streamlight ####
#Function for batching over multiple sites
batch_model <- function(Site, params, read_dir, save_dir){
  #Get the model driver
  driver_file <- readRDS(paste(read_dir, "/", Site, "_driver.rds", sep = ""))
  
  #Get model parameters for the site
  site_p <- params[params[, "Site_ID"] == Site, ]
  
  #Run the model
  modeled <- stream_light(
    driver_file, 
    Lat = site_p[, "Lat"], 
    Lon = site_p[, "Lon"],
    channel_azimuth = site_p[, "Azimuth"], 
    bottom_width = site_p[, "Width"], 
    BH = site_p[, "BH"],
    BS = site_p[, "BS"], 
    WL = site_p[, "WL"], 
    TH = site_p[, "TH"], 
    overhang = site_p[, "overhang"],
    overhang_height = site_p[, "overhang_height"], 
    x_LAD = site_p[, "x"]
  )
  
  #Save the output
  saveRDS(modeled, paste(save_dir, "/", Site, "_predicted.rds", sep = ""))
  
} #End batch_model 

#Applying the model to all sites
model_rd <- working_dir
model_sd <- working_dir

#Running the model
lapply(
  site_parm[, "Site_ID"], 
  FUN = batch_model, 
  params = site_parm,
  read_dir = model_rd,
  save_dir = model_sd
) 

#Take a look at the output
CBP_predicted <- readRDS(paste(working_dir, '/CBP_predicted.rds', sep = ''))
CBP_predicted[1:2, ]

CBP_predicted %>%
  mutate(date = as.Date(local_time)) %>%
  group_by(date) %>%
  summarize(par_surface = sum(PAR_surface),
            par_inc = sum(PAR_inc)) %>%
  ungroup() %>%
  ggplot(aes(date, par_inc)) +
    geom_point(col = 'grey') +
    geom_point(aes(y = par_surface), col = 'gold') +
    theme_minimal()

# Combine the outputs into one file
working_dir <- 'C:/Users/Alice Carter/git/nhc_50yl/data/light/drivers/'
dat <- data.frame()
for(site in data.frame(site_parm)[,'Site_ID']){
  pred <- readRDS(paste(working_dir, site, '_predicted.rds', sep = ''))
  ss <- pred %>% 
  select(local_time, LAI, PAR_inc, PAR_surface) %>%
  mutate(date = as.Date(local_time, tz = 'EST'),
         site = site) %>%
  group_by(site, date) %>%
  summarize(LAI = mean(LAI, na.rm = T),
            PAR_inc = sum(PAR_inc),
            PAR_surface = sum(PAR_surface)) %>%
  ungroup()
  
  dat <- bind_rows(dat, ss)
}
dat %>% filter(date > as.Date('2019-03-01')) %>%
ggplot(aes(date, PAR_surface, col = site)) +
  geom_point()

write_csv(dat, paste0(working_dir, 'daily_modeled_light_all_sites.csv'))

