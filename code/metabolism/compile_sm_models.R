## Model Metabolism #
# adapted from JRB script
# This version runs metabolism on NHC sites using the K600 values from Hall 1972
setwd('C:/Users/alice.carter/git/nhc_50yl/')
source("src/metabolism/inspect_model_fits.r")
# Compile model outputs ####

flow_dates <- read_csv("data/rating_curves/flow_dates_filter.csv")

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
                     pattern = '^[a-z]+_([a-z]+)_?([0-9]+)?_([a-z]+_[a-z]+)')
    site <- tmp[2]
    method <- tmp[4]

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
    out <- filter_model(fit, flow_dates)
    preds <- out[[1]]
    coverage <- data.frame(site = site,
                           method = method)
    coverage <- bind_cols(coverage, out[[2]])

    out <- fill_summarize_met(preds)
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
        sitename <- sn$sn[sn$site == site]
        mtext(sitename, 2, line = 3)
        plot_KvER(preds)

    dev.off()
    met_sum <- bind_cols(coverage, out[[2]])

    met_summary <- bind_rows(met_summary, met_sum)
    all_preds <- bind_rows(all_preds, preds)
    all_filled_preds <- bind_rows(all_filled_preds, cum)

}
# dev.off()

saveRDS(list(preds = all_preds,
             summary = met_summary,
             cumulative = all_filled_preds),
        "data/metabolism/compiled/raymond_met.rds")




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
