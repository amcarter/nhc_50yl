
# Compile model outputs ####

source("code/metabolism/inspect_model_fits.r")
flow_dates <- read_csv("data/rating_curves/flow_dates_filter.csv")
bad_DO_fits <- data.frame(date = c(as.Date(c('2016-12-29', '2017-03-31', '2017-04-19',
                                           '2019-04-23', '2017-07-06', '2017-07-07',
                                           '2017-07-08', '2017-07-10', '2017-07-11',
                                           '2017-07-18', '2017-07-20', '2017-07-30',
                                           '2017-07-31', '2017-08-24', '2017-08-25',
                                           '2017-09-01', '2017-09-03', '2017-10-30',
                                           '2017-10-31', '2017-11-13', '2017-11-30',
                                           '2017-12-02', '2017-12-03', '2017-12-06',
                                           '2017-12-18', '2017-12-24', '2017-12-30',
                                           '2018-01-17', '2018-01-20', '2018-01-21',
                                           '2018-07-28', '2018-07-29', '2019-06-06',
                                           '2019-06-27', '2019-07-01', '2019-08-19',
                                           '2019-08-22', '2019-08-23', '2019-09-10',
                                           '2019-10-22', '2019-10-23', '2019-10-24',
                                           '2019-10-26', '2019-10-31', '2019-11-06')),
                                   seq(as.Date('2017-01-05'), as.Date('2017-01-19'), by = 'day'),
                                   seq(as.Date('2017-07-14'), as.Date('2017-07-16'), by = 'day'),
                                   seq(as.Date('2017-07-20'), as.Date('2017-07-25'), by = 'day'),
                                   seq(as.Date('2017-08-10'), as.Date('2017-08-13'), by = 'day'),
                                   seq(as.Date('2017-09-06'), as.Date('2017-09-11'), by = 'day'),
                                   seq(as.Date('2017-09-14'), as.Date('2017-09-25'), by = 'day'),
                                   seq(as.Date('2017-10-02'), as.Date('2017-10-04'), by = 'day'),
                                   seq(as.Date('2017-11-03'), as.Date('2017-11-05'), by = 'day'),
                                   seq(as.Date('2017-11-25'), as.Date('2017-11-28'), by = 'day'),
                                   seq(as.Date('2018-06-18'), as.Date('2018-06-26'), by = 'day'),
                                   seq(as.Date('2018-06-30'), as.Date('2018-07-04'), by = 'day'),
                                   seq(as.Date('2018-07-09'), as.Date('2018-07-14'), by = 'day'),
                                   seq(as.Date('2018-07-17'), as.Date('2018-07-22'), by = 'day'),
                                   seq(as.Date('2018-08-27'), as.Date('2018-08-30'), by = 'day'),
                                   seq(as.Date('2018-09-02'), as.Date('2018-09-08'), by = 'day'),
                                   seq(as.Date('2019-07-03'), as.Date('2019-07-07'), by = 'day'),
                                   seq(as.Date('2019-07-11'), as.Date('2019-07-21'), by = 'day'),
                                   seq(as.Date('2019-07-27'), as.Date('2019-08-05'), by = 'day'),
                                   seq(as.Date('2019-08-09'), as.Date('2019-08-14'), by = 'day'),
                                   seq(as.Date('2019-08-28'), as.Date('2019-09-07'), by = 'day'),
                                   seq(as.Date('2019-09-13'), as.Date('2019-09-17'), by = 'day'),
                                   seq(as.Date('2019-09-22'), as.Date('2019-09-26'), by = 'day'),
                                   seq(as.Date('2019-11-17'), as.Date('2019-11-22'), by = 'day'))) %>%
    mutate(nhc = 'bad_DO')
bad_DO_fits <- data.frame(date = as.Date(c('2019-04-02', '2019-04-07', '2019-04-10',
                                           '2019-04-11', '2019-05-10', '2019-06-07',
                                           '2019-06-20', '2019-08-29', '2019-09-15',
                                           '2019-09-22', '2019-09-23', '2019-09-24',
                                           '2019-09-29', '2019-10-03', '2019-10-12',
                                           '2019-10-13', '2019-10-15', '2019-10-27',
                                           '2019-10-28', '2019-11-06', '2019-11-08',
                                           '2019-11-12', '2019-11-17', '2019-11-19',
                                           '2019-11-21', '2019-11-28', '2019-11-29',
                                           '2019-12-11'))) %>%
    mutate(cbp = 'bad_DO') %>%
    full_join(bad_DO_fits, by = 'date')

flow_dates <- left_join(flow_dates, bad_DO_fits, by = 'date') %>%
    mutate(pm = NA_character_,
           wb = NA_character_,
           wbp = NA_character_,
           unhc = NA_character_)

filelist <- list.files("data/metabolism/modeled/finalQ")#[c(-5, -7)]
GPP_min = 0
ER_max = 0

met_summary <- data.frame()
all_preds <- data.frame()
all_filled_preds <- data.frame()

sn <- data.frame(site = c('nhc', 'pm', 'cbp','wb','wbp','unhc'),
                 sn = c('NHC_8.5', 'NHC_6.9', 'NHC_5',
                        'NHC_2.5', 'NHC_2.3', 'NHC_0'))
# pdf("figures/model_diagnostics_raymond.pdf", width = 9, height = 6)
obserr <- data.frame()
for(file in filelist) {
    # uninformed raymond ests
    fit <- readRDS(paste0("data/metabolism/modeled/finalQ/", file))
    tmp <- str_match(string = file,
                     pattern = '^[a-z]+_([a-z]+)_([a-z]+_[a-z]+)')
    site <- tmp[2]
    method <- tmp[3]

    # Incorporate something here to look at obs and proc error
    obs <- get_fit(fit)$overall %>%
        select(err_obs_iid_sigma_50pct, err_obs_iid_sigma_Rhat,
               err_proc_iid_sigma_50pct, err_proc_iid_sigma_Rhat) %>%
        mutate(site = site)
    obserr <- bind_rows(obserr, obs)

    if(site == "wbp"){
        fit@data <- fit@data %>% filter(date <= as.Date("2020-03-25"))
        fit@fit$daily <- fit@fit$daily %>% filter(date <= as.Date("2020-03-25"))
        fit@data_daily <- fit@data_daily %>% filter(date <= as.Date("2020-03-25"))
    }

    group_year = FALSE
    if(site %in% c('nhc', 'unhc')){group_year = TRUE}

    # plot_zoom(fit@data)
    out <- filter_model(fit, flow_dates, group_year = group_year)
    preds <- out[[1]]
    coverage <- out[[2]] %>%
        mutate(site = site,
              method = method) %>%
        select(-any_of("year"))

    out <- fill_summarize_met(preds, group_year = group_year)
    cum <- out[[1]] %>%
        mutate(site = site,
               method = method)
    preds <- preds %>%
        mutate(site = site,
               method = method)

    # bad <- unique(preds$date[is.na(preds$GPP) & is.na(preds$ER)])
    # dat <- fit@data %>%
    #     filter(!(date %in% bad))
    # plot_zoom(dat)
    # mcmc <- get_mcmc(fit)
    # rstan::traceplot(mcmc, pars = c("K600_daily[11]",
    #                                 "K600_daily[12]",
    #                                 "K600_daily[13]",
    #                                 "K600_daily[14]",
    #                                 "K600_daily[15]",
    #                                 "K600_daily[16]",
    #                                 "K600_daily[17]",
    #                                 "K600_daily[18]",
    #                                 "K600_daily[19]",
    #                                 "K600_daily[20]"), nrow = 5)
    # plot_diagnostics(fit, preds, paste(site, year, method),
    #                  ylim = c(-15, 7), lim = 7)
    tiff(paste0('figures/SI/model_fit_', site, '.tiff'),
         width = 6*800, height = 3*800, res = 800, units = 'px',
         compression = 'lzw')
    par(ps = 10,
        mfrow = c(1,2),
        mar = c(3, 4, 2, 0),
        oma = c(0, 0, 0, 2))
    plot_binning2(fit, preds)
    sitename <- toupper(site)
    mtext(sitename, 2, line = 3)
    plot_KvER(preds)

    dev.off()
    met_sum <- bind_cols(coverage, out[[2]])

    met_summary <- bind_rows(met_summary, met_sum)
    all_preds <- bind_rows(all_preds, preds)
    all_filled_preds <- bind_rows(all_filled_preds, cum)

}
# dev.off()

all_preds <- mutate(all_preds,
                    year = case_when(is.na(year) ~ 2019,
                                     TRUE ~ year))
met_summary <- met_summary %>%
    filter(year %in% c(2017, 2018, 2019))

saveRDS(list(preds = all_preds,
             summary = met_summary,
             cumulative = all_filled_preds),
        "data/metabolism/compiled/raymond_met.rds")


all_preds %>%
    filter(year %in% c(2017, 2018, 2019),
           site %in% c('cbp', 'nhc')) %>%
    summary()

preds <- all_preds %>%
    filter(year %in% c(2017, 2018, 2019),
           site %in% c('cbp', 'nhc')) %>%
    mutate(NEP = GPP+ER,
           PR = -GPP/ER,
           PR = case_when(is.infinite(PR) ~ NA,
                          TRUE ~ PR))

summary(preds)
length(which(preds$NEP < 0))/length(which(!is.na(preds$NEP)))

preds %>%
    group_by(year, site) %>%
    summarize(across(c('GPP','ER'),
                     .fns = c(mean = \(x) mean(x, na.rm = T),
                              max = \(x) max(x, na.rm = T),
                              min = \(x) min(x, na.rm = T),
                              cv = \(x) sd(x, na.rm = T)/mean(x, na.rm = T))))
met_summary %>% filter(site %in% c('cbp', 'nhc'))

# png('figures/SI/KQ_KER_plots_SMfits.png', width = 7.5, height = 7,
#     res = 300, units = 'in')
# par(mfrow = c(4, 2),
#     mar = c(2.5, 3, 0, 0),
#     oma = c(0,3,4,2))
# for(file in filelist) {
#   # uninformed raymond ests
#   tmp <- str_match(string = file,
#                    pattern = '^[a-z]+_([a-z]+)_?([0-9]+)?_([a-z]+_[a-z]+)')
#   site <- tmp[2]
#   method <- tmp[4]
#   year = 2019
#
#   if(site %in% c("nhc", "unhc")) next
#
#   fit <- readRDS(paste0("data/metabolism/modeled/finalQ/", file))
#   if(site == "wbp"){
#     fit@data <- fit@data %>% filter(date <= as.Date("2020-03-20"))
#     fit@fit$daily <- fit@fit$daily %>% filter(date <= as.Date("2020-03-20"))
#     fit@data_daily <- fit@data_daily %>% filter(date <= as.Date("2020-03-20"))
#   }
#   out <- filter_model(fit, flow_dates)
#   preds <- out[[1]]
#   coverage <- data.frame(site = site,
#                          year = year,
#                          method = method)
#   coverage <- bind_cols(coverage, out[[2]])
#
#   out <- fill_summarize_met(preds)
#   cum <- out[[1]] %>%
#     mutate(site = site,
#            year = year,
#            method = method)
#   preds <- preds %>%
#     mutate(site = site,
#            year = year,
#            method = method)
#
#   # plot_diagnostics(fit, preds, paste(site, year, method),
#   #                  ylim = c(-15, 7), lim = 7)
#   plot_binning2(fit, preds)
#   sitename <- sn$sn[sn$site == site]
#   mtext(sitename, 2, line = 4.5, cex = 1.3)
#   mtext('K600', 2, line = 2.2, cex = 0.8)
#   plot_KvER(preds)
#   met_sum <- bind_cols(coverage, out[[2]])
#
#   met_summary <- bind_rows(met_summary, met_sum)
#   all_preds <- bind_rows(all_preds, preds)
#   all_filled_preds <- bind_rows(all_filled_preds, cum)
#
# }
# dev.off()
# png('figures/SI/KQ_KER_plots_SMfits2.png', width = 7.5, height = 10,
#     res = 300, units = 'in')
# for(file in filelist) {
#   # uninformed raymond ests
#   tmp <- str_match(string = file,
#                    pattern = '^[a-z]+_([a-z]+)_?([0-9]+)?_([a-z]+_[a-z]+)')
#   site <- tmp[2]
#   method <- tmp[4]
#   if(!(site %in% c("nhc", "unhc"))) next
#     year = as.numeric(tmp[3])
#
#   fit <- readRDS(paste0("data/metabolism/modeled/finalQ/", file))
#   if(site == "wbp"){
#     fit@data <- fit@data %>% filter(date <= as.Date("2020-03-20"))
#     fit@fit$daily <- fit@fit$daily %>% filter(date <= as.Date("2020-03-20"))
#     fit@data_daily <- fit@data_daily %>% filter(date <= as.Date("2020-03-20"))
#   }
#   out <- filter_model(fit, flow_dates)
#   preds <- out[[1]]
#   coverage <- data.frame(site = site,
#                          year = year,
#                          method = method)
#   coverage <- bind_cols(coverage, out[[2]])
#
#   out <- fill_summarize_met(preds)
#   cum <- out[[1]] %>%
#     mutate(site = site,
#            year = year,
#            method = method)
#   preds <- preds %>%
#     mutate(site = site,
#            year = year,
#            method = method)
#
#   # plot_diagnostics(fit, preds, paste(site, year, method),
#   #                  ylim = c(-15, 7), lim = 7)
#   plot_binning2(fit, preds)
#   sitename <- sn$sn[sn$site == site]
#   mtext(paste(sitename, year), 2, line = 4.5, cex = 1.3)
#   mtext('K600', 2, line = 2.2, cex = 0.8)
#   plot_KvER(preds)
#   met_sum <- bind_cols(coverage, out[[2]])
#
#   met_summary <- bind_rows(met_summary, met_sum)
#   all_preds <- bind_rows(all_preds, preds)
#   all_filled_preds <- bind_rows(all_filled_preds, cum)
#
# }
# dev.off()
#
# # #inspect fits ####
# # ray_met <- all_preds %>%
# #   mutate(doy = format(date, "%j")) %>%
# #   group_by(doy) %>%
# #   summarize(GPP.upper = quantile(GPP, .975, na.rm = T),
# #             GPP.lower = quantile(GPP, .025, na.rm = T),
# #             ER.upper = quantile(ER, .975, na.rm = T),
# #             ER.lower = quantile(ER, .025, na.rm = T),
# #             GPP = median(GPP, na.rm = T),
# #               ER = median(ER, na.rm = T))
