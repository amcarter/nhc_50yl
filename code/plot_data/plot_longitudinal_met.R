# longitudinal exploration of 2019 NHC metabolism
# plots of metabolism vs distance binned by month and season 
# 3d plots of metabolic surface through space and time
library(tidyverse)
library(ggpubr)

# setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data")
source("src/metabolism/inspect_model_fits.r")

# it would be nice to recompile the metabolism file so it has the PWC site
dat <- readRDS("data/metabolism/compiled/met_preds_stream_metabolizer.rds")
sites <- read_csv("data/siteData/NHCsite_metadata.csv") %>%
  slice(c(1:7))

preds <- dat$preds %>% 
  filter(year == 2019)

dists <- sites %>% 
  select(site = sitecode, distance_m) %>%
  mutate(distance_m = 8450 - distance_m)

preds <- left_join(preds, dists, by = 'site') %>%
  mutate(month = as.factor(month))

seas <- preds %>%
  mutate(season = case_when(month %in% c(12,1,2) ~ 'winter',
                            month %in% c(3,4,5) ~ 'spring',
                            month %in% c(6,7,8) ~ 'summer',
                            month %in% c(9,10,11) ~'autumn'))

ggplot(preds, aes(factor(distance_m), GPP,  fill = month)) +
  geom_boxplot() 
  geom_smooth(se = FALSE, method = 'lm')
er <- ggplot(preds, aes(distance_m, ER, group = month, color = month)) +
  geom_point() +
  geom_smooth(se = FALSE)
gpp <- ggplot(preds, aes(distance_m, GPP, group = month, color = month)) +
  geom_point() +
  geom_smooth(se = FALSE)

ggpubr::ggarrange(gpp, er, ncol = 1, common.legend = TRUE)

# 3D metabolism plot 

library(plotly)
even_dists <- data.frame(distance_m = seq(0, 8450, by = 10))
GPP_dat <- preds %>% select(date, GPP, distance_m) %>%
  group_by(distance_m) %>%
  mutate(GPP = zoo::na.approx(GPP, na.rm = F),
         GPP = zoo::rollmean(GPP, 7, fill = NA)) %>%
  filter(date >= as.Date('2019-03-10'),
         date <= as.Date('2020-03-03')) %>%
  group_by(date, distance_m) %>%
  summarize(GPP = mean(GPP, na.rm = T)) %>%
  pivot_wider(names_from = 'date', values_from = 'GPP')

GPP_filled <- left_join(even_dists, GPP_dat, by = 'distance_m') %>%
  mutate(across(.fns = zoo::na.approx, na.rm = F)) %>%
  filter(distance_m %% 50 == 0) %>%
  select(-distance_m) %>%
  as.matrix


ER_dat <- preds %>% select(date, ER, distance_m) %>%
  group_by(distance_m) %>%
  mutate(ER = zoo::na.approx(ER, na.rm = F),
         ER = zoo::rollmean(ER, 7, fill = NA)) %>%
  filter(date >= as.Date('2019-03-10'),
         date <= as.Date('2020-03-03')) %>%
  group_by(date, distance_m) %>%
  summarize(ER = mean(ER, na.rm = T)) %>%
  pivot_wider(names_from = 'date', values_from = 'ER')

ER_filled <- left_join(even_dists, ER_dat, by = 'distance_m') %>%
  mutate(across(.fns = zoo::na.approx, na.rm = F)) %>%
  filter(distance_m %% 50 == 0) %>%
  select(-distance_m) %>%
  as.matrix

axx <- list(
  ticketmode = 'array',
  ticktext = c("Apr-19", "Jun-19", "Aug-19", "Oct-19", "Dec-19", "Feb-20"),
  tickvals = c(23, 84, 145, 206, 267, 329),
  range = c(1,360),
  title = "Date"
)

axy <- list(
  ticketmode = 'array',
  ticktext = c("0", "1", "2", "3", "4", "5", "6", "7", "8"),
  tickvals = c(1, 21, 41, 61, 81, 101, 121, 141, 161),
  range = c(1, 170),
  title = 'Distance downstream (km)'
)

axz_GPP <- list(
  nticks = 4,
  range = c(0, 1.5),
  title = "GPP gO2/m2/d"
)
axz_ER <- list(
  nticks = 5,
  range = c(-5, 0),
  title = "ER gO2/m2/d"
)

gpp<- plot_ly(z = GPP_filled, type = 'surface') %>% 
  layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz_GPP))

er<- plot_ly(z = ER_filled, type = 'surface') %>%
  layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz_ER))

# save the widget
library(htmlwidgets)
saveWidget(gpp, file= "figures/surface_plots/3dgpp.html")
saveWidget(er, file= "figures/surface_plots/3der.html")
