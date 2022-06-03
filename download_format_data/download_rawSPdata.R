# Pull raw data from streamPulse website
# All data used are publicly available and can be accessed using this script by
# any users.

# library(devtools)
# install_github('streampulse/StreamPULSE', dependencies=TRUE)

# # Load packages.
# library(StreamPULSE)
# library(lubridate)
# library(tidyverse)

sites <- read_csv("data/siteData/NHCsite_metadata.csv") %>% slice(1:7)
dir.create('data/metabolism')
dir.create('data/metabolism/raw')

# download all the data for each at site:
for(i in 1:nrow(sites)){
  vars <- as.character(unlist(
    query_available_data("NC",sites[i,]$sitecode)$variables))
  dat <- request_data(sites[i,]$siteID,
                      startdate = as.Date(sites[i,]$startdate.UTC),
                      enddate = as.Date("2020-03-21"),
                      variables=vars)
  dd <- dat$data %>%
    mutate(value = ifelse(flagtype %in% c("Bad Data", "Questionable"),
                          NA, value)) %>%
    select(DateTime_UTC, site, value, variable)  %>%
    arrange(DateTime_UTC)
  w <- range(which(!is.na(dd$value)))
  dd <- dd[w[1]:w[2],]
  dd <- dd %>%
    pivot_wider(names_from = variable, values_from = value)
  datetimes <- data.frame(DateTime_UTC = seq(dd$DateTime_UTC[1],
                                             dd$DateTime_UTC[nrow(dd)],
                                             by = "15 min"))
  dd <- left_join(datetimes, dd, by = "DateTime_UTC")
  write_csv(dd, paste0("C:/Users/alice.carter/git/nhc_50yl/data/metabolism/raw/",
                       dd$site[1], "_", as.Date(dat$specs$enddate), ".csv"))

}
