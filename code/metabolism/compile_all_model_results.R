## convert model results to Carbon and compile all methods/years
source("code/metabolism/inspect_model_fits.r")
# Compile model outputs ####
convert = FALSE
# convert = TRUE

if(convert == TRUE){
  O2toC <- 12.0107/(2*15.999)
  }
if(convert == FALSE){
  O2toC <- 1
}

sm_fit <- readRDS("data/metabolism/compiled/raymond_met.rds")
# dir_fit <- readRDS("NHC_2019_metabolism/data/metabolism/hall/hall_met_60min_2021_01.rds")

sm_preds <- sm_fit$preds %>%
  as_tibble() %>%
  select(-ends_with("Rhat"), -errors, -good_flow) %>%
  rename(discharge = discharge.daily) %>%
  mutate(site = toupper(site),
         doy = as.numeric(format(date, "%j")),
         month = as.numeric(format(date, "%m")),
         across(starts_with(c('GPP', 'ER')), ~ . * O2toC),
         era = "now") %>%
  filter(date <= as.Date("2020-03-25"))
sm_preds_filled <- sm_fit$cumulative %>%
  as_tibble() %>%
  rename(discharge = discharge.daily) %>%
  select(-discharge.daily_cum, -temp.min_cum, -temp.water_cum) %>%
  mutate(site = toupper(site),
         doy = as.numeric(format(date, "%j")),
         month = as.numeric(format(date, "%m")),
         across(starts_with(c("ER", "GPP")), ~ .* O2toC))

# dir_preds <- dir_fit$preds %>%
#   as_tibble() %>%
#   select(-good_flow) %>%
#   mutate(doy = as.numeric(format(date, '%j')),
#          month = as.numeric(format(date, '%m')),
#          across(c('GPP', 'ER'), ~ . * O2toC),
#          era = 'now')
# dir_preds_filled <- dir_fit$cumulative %>%
#   as_tibble() %>%
#   rename(gpp = GPP, er = ER, gpp_cum = GPP_cum, er_cum = ER_cum) %>%
#   mutate(doy = as.numeric(format(date, '%j')),
#          month = as.numeric(format(date, '%m')),
#          across(starts_with(c('gpp', 'er')), ~ . * O2toC),
#          era = 'now') %>%
#   select(-temp.water_cum) %>%
#   left_join(dir_preds[,c(1,2,8)], by = c("site", "date")) %>%
#   arrange("site", "date") %>%
#   mutate(discharge = na.approx(discharge, na.rm = F))

# dir_fit$summary$method <- "direct_calculation"

sm_met_sum <- sm_fit$summary %>%
  as_tibble() %>%
  mutate(site = toupper(site)) %>%
  # bind_rows(dir_fit$summary) %>%
  rename(peak_gpp = gpp_max10d, peak_er = er_max10d) %>%
  mutate(days = round(total_days * pctcoverage, 0),
         across((starts_with(c('gpp', 'er')) & !ends_with('cv')), ~ . * O2toC))

hall_preds <- read_csv('data/hall_data/hall_table_15.csv') %>%
  select(date, site, depth = depth_m, GPP = GPP_gO2m2d, ER = ER_gO2m2d) %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y"),
         site = case_when(site == "Concrete" ~ "CBP",
                          site == "Blackwood" ~ "BLK",
                          site == "Wood Bridge" ~ "WB"),
         across(c('GPP', 'ER'), ~ . * O2toC)) %>%
  group_by(site, date) %>%
  summarize_all(mean, na.rm = T) %>%
  ungroup()# %>%
    # filter(site == 'CBP')
hall_qt <- read_csv("data/hall_data/hall_discharge_temp_daily.csv") %>%
  rename(temp.water = water_temp_C, discharge = discharge_m3s) %>%
    group_by(date) %>%
    select(-notes) %>%
    summarize(across(everything(), ~mean(.x, na.rm = T))) %>%
    mutate(site = 'CBP') %>%
    ungroup()
hall_preds <- hall_preds %>%
  left_join(hall_qt, by = c("date", "site")) %>%
  mutate(ER = -ER,
         era = "then",
         doy = as.numeric(format(date, '%j')),
         month = as.numeric(format(date, '%m')),
         year = case_when(date <= as.Date("1969-04-13") ~ 1968,
                          date <= as.Date("1970-04-13") ~ 1969,
                          TRUE ~ 1970)) %>%
  filter(GPP < 3)


# dir_preds <- bind_rows(dir_preds, hall_preds)
sm_preds <- bind_rows(sm_preds, hall_preds)
# calculate cumulative data for Hall ####
hall_cbp <- hall_preds %>%
  rename(gpp = GPP, er = ER) %>%
  mutate(er = -er) %>%
  filter(site == "CBP") %>%
  arrange(date)
met <- hall_cbp %>%
  group_by(year) %>%
  summarize(across(c('gpp','er'), .fns = list(mean = ~mean(., na.rm = T),
                                              median = ~median(., na.rm = T),
                                              max = ~max(., na.rm = T),
                                              min = ~min(., na.rm = T),
                                              cv = ~sd(., na.rm = T)/mean(., na.rm = T)),
                   .names = '{col}_{fn}'),
            days = length(gpp)) %>%
  ungroup() %>%
  mutate(site = "CBP") %>%
  slice(-3)
met <- hall_preds %>%
  mutate(gpp = GPP, er = -ER) %>%
  group_by(site) %>%
  summarize(across(c('gpp','er'), .fns = list(mean = ~mean(., na.rm = T),
                                              median = ~median(., na.rm = T),
                                              max = ~max(., na.rm = T),
                                              min = ~min(., na.rm = T),
                                              cv = ~sd(., na.rm = T)/mean(., na.rm = T)),
                   .names = '{col}_{fn}'),
            days = length(gpp)) %>%
  bind_rows(met) %>%
  mutate(method = "direct_calculation",
         era = "then",
         bad_flow = c(0,1,0,0,1))

cum <- data.frame(date = seq(hall_cbp$date[1],
                             hall_cbp$date[nrow(hall_cbp)],
                             by = "day")) %>%
  as_tibble() %>%
  left_join(hall_cbp) %>%
  select(date, depth, gpp, er, doy) %>%
  mutate(across(-date, na.approx, na.rm = F),
         year = case_when(date <= as.Date("1969-04-13") ~ 1968,
                          date <= as.Date("1970-04-13") ~ 1969,
                          TRUE ~ 1970),
         era = "then",
         month = as.numeric(format(date, '%m'))) %>%
  left_join(hall_qt[,c(1,2,4)], by = "date")
cum$gpp_cum <- c(cumsum(cum$gpp[cum$year == 1968]),
             cumsum(cum$gpp[cum$year == 1969]),
             cumsum(cum$gpp[cum$year == 1970]))
cum$er_cum <- c(cumsum(cum$er[cum$year == 1968]),
            cumsum(cum$er[cum$year == 1969]),
            cumsum(cum$er[cum$year == 1970]))

# dir_preds_filled <- bind_rows(dir_preds_filled, cum)
sm_preds_filled <- bind_rows(sm_preds_filled, cum)



met$gpp_cum = c(NA, NA, NA, sum(cum$gpp[cum$year == 1968]),
                sum(cum$gpp[cum$year == 1969]))
met$er_cum = c(NA, NA, NA, sum(cum$er[cum$year == 1968]),
                sum(cum$er[cum$year == 1969]))

tmp <- hall_cbp %>%
  group_by(doy) %>%
  summarise(across(c('gpp', 'er'), ~mean(., na.rm = T))) %>%
  ungroup()

doy <- data.frame(doy = c(347:365, 1:365, 1:17)) %>%
  left_join(tmp) %>%
  as_tibble() %>%
  mutate(across(-doy, ~na.approx(.))) %>%
  slice(-c(1:19, 385:401)) %>%
  mutate(across(-doy, ~cumsum(.), .names = "{col}_cum"))
met$gpp_cum[2] <- doy$gpp_cum[365]
met$er_cum[2] <- doy$er_cum[365]
met <- mutate(met,
              total_days = 365,
              pctcoverage = days/total_days)

# find peak weeks here
n <- nrow(cum)
l = (n-9)
weekly <- tibble(date = cum$date[1:l],
                 GPP_week = rep(NA_real_, l),
                 ER_week = rep(NA_real_, l))
for(i in 1:l){
  weekly$GPP_week[i] <- sum(cum$gpp[i:(i+9)])
  weekly$ER_week[i] <- sum(cum$er[i:(i+9)])
}
peakweek <- weekly %>% left_join(cum[,c(1,6)], by = 'date') %>%
  group_by(year) %>%
  summarize(peak_gpp = date[which.max(GPP_week)],
            peak_er = date[which.max(ER_week)]) %>%
  mutate(site = 'CBP') %>% slice(-3)
weekly <- tibble(doy = 1:356,
                 GPP_week = rep(NA_real_, 356),
                 ER_week = rep(NA_real_, 356))
for(i in 1:356){
  weekly$GPP_week[i] <- sum(doy$gpp[i:(i+9)])
  weekly$ER_week[i] <- sum(doy$er[i:(i+9)])
}
week <- weekly %>%
  mutate(date = seq(as.Date('1969-01-01'), by = 'day', length.out = 356)) %>%
  summarize(peak_gpp = date[which.max(GPP_week)],
            peak_er = date[which.max(ER_week)]) %>%
  mutate(site = 'CBP')
tmp <- data.frame(site = c('BLK', 'CBP', 'WB'))
met <- right_join(week, tmp) %>%
  bind_rows(peakweek) %>%
  full_join(met, by = c('site', 'year'))
summary <- bind_rows(sm_met_sum, met)

if(convert == TRUE){
  write_csv(summary, "data/metabolism/compiled/metabolism_summary_table_gC.csv")
  saveRDS(list(preds = sm_preds,
             filled = sm_preds_filled),
        "data/metabolism/compiled/met_preds_stream_metabolizer_C.rds")
}
if(convert == FALSE){
  write_csv(summary, "data/metabolism/compiled/metabolism_summary_table_gO2.csv")
  saveRDS(list(preds = sm_preds,
             filled = sm_preds_filled),
        "data/metabolism/compiled/met_preds_stream_metabolizer_O2.rds")
}

# saveRDS(list(preds = dir_preds,
#              filled = dir_preds_filled),
#         "NHC_2019_metabolism/data/metabolism/compiled/met_preds_direct_calc.rds")


# Numbers for results section ####
# 3.1 Met patterns and drivers:
summary <- summary %>%
  filter(!(era == 'now' & method == 'direct_calculation'))
total_ests <- summary %>%
  filter(method == "uninformed_raymond") %>%
  summarize(across(where(is.numeric), sum)) %>%
  select(c(2:9, 21)) %>%
  mutate(bad_flow = lost_GPP - missing_data - bad_Rhat - neg_GPP,
         site = 'total')
total_ests <- summary %>%
  filter(method == "uninformed_raymond" | site== "CBP" & is.na(year)) %>%
  select(c(1:2,4:11,25, 28)) %>%
  mutate(bad_flow = lost_GPP - missing_data - bad_Rhat - neg_GPP) %>%
  bind_rows(total_ests)

write_csv(total_ests,
          'data/metabolism/compiled/coverage_of_met_estimates.csv')

yy <- summary %>%
  filter(method == "uninformed_raymond", site %in% c('NHC', 'UNHC')) %>%
  mutate(scale = "years")
ss <- summary %>%
  filter(method == "uninformed_raymond", year == 2019) %>%
  mutate(scale = "sites")
dd <- summary %>%
  filter(site == 'CBP', !(year %in% 1968:1969)) %>%
  mutate(year = ifelse(method == 'uninformed_raymond', 2019, 1969),
         scale = "decades")
scales <- bind_rows(yy, ss, dd) %>%
  select(site, year, scale, days, pctcoverage, starts_with(c('gpp', 'er_')))

annual_summary <- scales %>%
  select(scale, ends_with(c('cum', "cv"))) %>%
  group_by(scale) %>%
  summarize(across(ends_with(c('cum', 'cv')), .fns = list(mean = ~mean(.),
                                                 min = ~min(.),
                                                 max = ~max(.),
                                                 cv = ~sd(.)/mean(.)),
                                     .names = '{col}_{fn}')) %>%
  pivot_longer(starts_with(c('gpp', 'er')), names_to = c('met','var', 'stat'),
               values_to = 'value',
               names_pattern = '([a-z]+)_([a-z]+)_([a-z]+)') %>%
  pivot_wider(names_from = var, values_from = value)

cc <- ggplot(annual_summary, aes(scale, cum, fill = met)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  facet_wrap(.~stat, scale = 'free_y', nrow = 1) +
  scale_fill_manual(values = c('sienna', 'forestgreen'))+
  ylab("cumulative annual metabolism")
vv <- ggplot(annual_summary, aes(scale, cv, fill = met)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  facet_wrap(.~stat, scale = 'free_y', nrow = 1) +
  scale_fill_manual(values = c('sienna', 'forestgreen'))+
  ylab("Within year cv of metabolism")

png('figures/comparison_cumulativemet_withinyearCV_across_scales.png',
    width = 7, height = 5, res = 300, units = 'in')
    ggpubr::ggarrange(cc, vv, common.legend = T, ncol = 1)
dev.off()

write_csv(annual_summary,
          'data/metabolism/compiled/across_scales_cum_comp.csv')


scales %>% filter(site %in% c('NHC', 'UNHC')) %>%
  group_by(site) %>%
  summarize(gpp_cv = sd(gpp_cum)/mean(gpp_cum),
            er_cv = sd(er_cum)/mean(gpp_cum),
            gpp_mean = mean(gpp_cum),
            er_mean = mean(er_cum),
            gpp_sd = sd(gpp_cum),
            er_sd = sd(er_cum))
scales %>% filter(year == 2019) %>%
  summarize(gpp_cv = sd(gpp_cum)/mean(gpp_cum),
            er_cv = sd(er_cum)/mean(gpp_cum))
# n <- nrow(cum)
# l = (n-9)
# weekly <- tibble(date = cum$date[1:l],
#                  GPP_week = rep(NA_real_, l),
#                  ER_week = rep(NA_real_, l))
# for(i in 1:l){
#   weekly$GPP_week[i] <- sum(cum$GPP[i:(i+9)])
#   weekly$ER_week[i] <- sum(cum$ER[i:(i+9)])
# }
# met$gpp_max10d <- weekly$date[which.max(weekly$GPP_week)]
# met$er_max10d <- weekly$date[which.max(weekly$ER_week)]
# se <- sum(is.na(cum$GPP))
# met$gpp_cum <- cum$GPP_cum[n-se]*365/(n-se)
# se <- sum(is.na(cum$ER))
# met$er_cum <- cum$ER_cum[n-se]*365/(n-se)
# met$daterange <- as.character(paste(cum$date[1], "-", cum$date[nrow(cum)]))
# met$total_days <- sum(!is.na(hall_preds$GPP))
# met$pctcoverage <- sum(!is.na(hall_preds$GPP))/nrow(cum)
# met$site = "all"
#
# met_sum <- met
#
# for(site in c("CBP", "WB", "BLK")){
#   preds <- hall_preds %>%
#     filter(site == !! site)
#   met <- preds %>%
#     summarize(gpp_mean = mean(GPP, na.rm = T),
#               gpp_median = median(GPP, na.rm = T),
#               gpp_max = max(GPP, na.rm = T),
#               er_mean = -mean(ER, na.rm = T),
#               er_median = -median(ER, na.rm = T),
#               er_max = -min(ER, na.rm = T))
#   cum <- data.frame(date = seq(preds$date[1],
#                              preds$date[nrow(preds)],
#                              by = "day")) %>%
#     as_tibble() %>%
#     left_join(preds) %>%
#     select(date, depth, GPP, ER) %>%
#     mutate(across(-date, na.approx, na.rm = F)) %>%
#     mutate(across(-date, cumsum, .names = "{col}_cum"))
#
#   n <- nrow(cum)
#   l = (n-9)
#   weekly <- tibble(date = cum$date[1:l],
#                    GPP_week = rep(NA_real_, l),
#                    ER_week = rep(NA_real_, l))
#   for(i in 1:l){
#     weekly$GPP_week[i] <- sum(cum$GPP[i:(i+9)])
#     weekly$ER_week[i] <- sum(cum$ER[i:(i+9)])
#   }
#   met$gpp_max10d <- weekly$date[which.max(weekly$GPP_week)]
#   met$er_max10d <- weekly$date[which.max(weekly$ER_week)]
#   se <- sum(is.na(cum$GPP))
#   met$gpp_cum <- cum$GPP_cum[n-se]*365/(n-se)
#   se <- sum(is.na(cum$ER))
#   met$er_cum <- cum$ER_cum[n-se]*365/(n-se)
#   met$daterange <- as.character(paste(cum$date[1], "-", cum$date[nrow(cum)]))
#   met$total_days <- sum(!is.na(preds$GPP))
#   met$pctcoverage <- sum(!is.na(preds$GPP))/nrow(cum)
#   met$site = site
#
#   met_sum <- bind_rows(met_sum, met)
# }
#
# # Group summary data####
#
# sm_fit$summary <- sm_fit$summary %>% mutate(site = toupper(site))
# summary <- bind_rows(sm_fit$summary, dir_fit$summary, met_sum)
#
# write_csv(summary, "NHC_2019_metabolism/data/metabolism/compiled/metabolism_summary_table_2021_01.csv")

# ggplot(summary, aes(x = site, y = gpp_median))+
#   geom_bar()
# png("figures/tmp/er_mintemp_by_month.png", width = 7, height = 5,
#     res = 300, units = "in")
# all_preds %>%
#   filter(GPP < 3) %>%
#   ggplot(aes(temp.min, ER, color = era)) +
#   geom_point() +
#   theme_bw() +
#   facet_wrap(.~month, scales = "free")
# dev.off()
#
# all_preds %>%
#   filter(GPP < 3) %>%
#   ggplot(aes(temp.water, GPP, color = era)) +
#   geom_point() +
#   theme_bw() +
#   facet_wrap(.~month, scales = "free")

# Compile Hall Data ####
# hallf <- data.frame(doy = 366,
#                    gpp_gcm2d = 0.258,
#                    er_gcm2d = 0.394) %>%
#   bind_rows(hall)
#
# hallf <- data.frame(doy = seq(1:366)) %>%
#   full_join(hallf) %>%
#   arrange(doy) %>%
#   mutate(across(-doy, na.approx, na.rm = F)) %>%
#   filter(doy != 366)
# ncon <- sum(hall$site =="Concrete")
# nblk <- sum(hall$site =="Blackwood")
# nwb <- sum(hall$site =="Wood Bridge")
#
# dd <- as.Date("1968-04-14")
# h68 <- hall %>%
#   filter(date >= dd,
#          date < dd + 365,
#          site == "Concrete")
# h69 <- hall %>%
#   filter(date >= dd+ 365,
#          date < dd + 2*365,
#          site == "Concrete")
# h_all <- hall %>%
#   mutate(doy = format(date, "%j")) %>%
#   group_by(doy) %>%
#   summarize(gpp_gcm2d = mean(gpp_gcm2d, na.rm = T),
#             er_gcm2d = mean(er_gcm2d, na.rm = T))
# rd <- range(hall$date[hall$site =="Wood Bridge"])
# dates <- data.frame(date = seq(rd[1], rd[2], by = "day"))
# nrow(dates)
# h_con <- hall %>%
#   filter(site == "Concrete") %>%
#   select(-site) %>%
#   full_join(dates) %>%
#   arrange(date) %>%
#   mutate(across(-date, na.approx, na.rm = F))
#
# h_c <- hall %>%
#   filter(site == "Concrete")
#
# h_wb <- hall %>%
#   filter(site == "Wood Bridge") %>%
#   mutate(doy = format(date, "%j")) %>%
#   group_by(doy) %>%
#   summarize(gpp_gcm2d = mean(gpp_gcm2d, na.rm = T),
#             er_gcm2d = mean(er_gcm2d, na.rm = T))
#
# h68 <- h_con %>%
#   filter(date >= dd ,
#          date < dd + 365)
# h68 <- h_con %>%
#   filter(date >= dd +365,
#          date < dd + 365*2)
#
# h68 <- hall
# met <- h68 %>%
#     summarize(gpp_mean = mean(gpp_gcm2d, na.rm = T),
#               gpp_median = median(gpp_gcm2d, na.rm = T),
#               gpp_max = max(gpp_gcm2d, na.rm = T),
#               er_mean = -mean(er_gcm2d, na.rm = T),
#               er_median = -median(er_gcm2d, na.rm = T),
#               er_max = -max(er_gcm2d, na.rm = T))
# met$site = "all"
# met$year = NA_real_
#
# cum <- h68 %>%
#  # select(-site) %>%
#  # mutate(across(-date, na.approx, na.rm = F)) %>%
#  mutate(across(-doy, cumsum, .names = "{col}_cum"))
#
# n <- nrow(cum)
# l = (n-9)
# weekly <- tibble(date = cum$date[1:l],
#                  GPP_week = rep(NA_real_, l),
#                  ER_week = rep(NA_real_, l))
# for(i in 1:l){
#   weekly$GPP_week[i] <- sum(cum$gpp_gcm2d[i:(i+9)])
#   weekly$ER_week[i] <- sum(cum$er_gcm2d[i:(i+9)])
# }
# met$gpp_max10d <- weekly$date[which.max(weekly$GPP_week)]
# met$er_max10d <- weekly$date[which.max(weekly$ER_week)]
# met$gpp_cum <- cum$gpp_gcm2d_cum[n]*365/(n)
# met$er_cum <- cum$er_gcm2d_cum[n]*365/(n)
# met$daterange <- as.character(paste(cum$date[1], "-", cum$date[nrow(cum)]))
# met$pctcoverage <- nrow(h68)/n
#
#  # met_hall <- data.frame()
# met_hall <- bind_rows(met_hall, met)

####


# png("../figures/metabolism_contours_K_estimates_v2.png", height = 5, width = 5,
#     units = "in", res = 300)
#
#   plot_kde_metab(hall_preds,  col = "steelblue", lim = 3.5)
#   par(new = T)
#   plot_kde_metab(dir_preds, col = "darkred", lim = 3.5)
#   legend("topright", cex = 1.2,
#          c("Then", "now"),
#          fill = c(alpha("steelblue", .75), alpha("darkred", .75)),
#          border = NA, bty = "n")
#   mtext("NHC Metabolism Estimates (n = 1473)", cex = 1.2)
# dev.off()
