# Bootstrap Metabolism comparison NHC SP data and Hall data #####
library(readr)
library(viridis)
library(beanplot)
library(scales)
library(dplyr)
library(ggplot2)
library(tidyr)

#setup ####
#historic data
dat <- readRDS("data/metabolism/compiled/met_preds_stream_metabolizer_O2.rds")
# dat <- readRDS("NHC_2019_metabolism/data/metabolism/compiled/met_preds_direct_calc.rds")
nhc_68_70 = dat$preds %>%
    filter(era == "then",
           site == "CBP") %>%
    select(-era)

nhc_68 <- nhc_68_70 %>%
    filter(year == 1968)
nhc_69 <- nhc_68_70 %>%
    filter(year == 1969)

gpp_68_70 = nhc_68_70$GPP
gpp_68 = nhc_68$GPP
gpp_69 = nhc_69$GPP

er_68_70 = nhc_68_70$ER
er_68 = nhc_68$ER
er_69 = nhc_69$ER

nep_68_70 = gpp_68_70 + er_68_70
nep_68 = gpp_68 + er_68
nep_69 = gpp_69 + er_69

# contemporary data
nhc_new <- dat$preds %>%
    filter(era == "now",
           site == "CBP")

nhc_site <- dat$preds %>%
    filter(era == "now",
           site == "NHC")

gpp_new = nhc_new$GPP
er_new = nhc_new$ER
nep_new = gpp_new + er_new

dates_new = nhc_new$date

#distribution plots ####
# png(width=9, height=6, units='in', type='cairo', res=300,
#     filename='figures/metab_distributions_v2.png')

# defpar = par(mfrow=c(2,3))
#
# #plot GPP dists, then and now
# plot(density(gpp_68_70, na.rm=TRUE), bty='l', col='sienna3',
#      main='GPP 1968-70 vs. 2019', xlab='GPP')
# lines(density(gpp_new, na.rm=TRUE), col='blue')
# legend('topright',
#        legend=c('68-70; n=76',
#                 paste0('19; n=', length(which(!is.na(gpp_new))))),
#        col = c('sienna3','blue'), lty = 1, bty = 'n',
#        seg.len = 1, cex = 0.9, lwd = 2)
#
# #plot ER dists, then and now
# plot(density(er_68_70, na.rm=TRUE), bty='l', col='sienna3',
#      main='ER 1968-70 vs. 2019', xlab='ER')
# lines(density(er_new, na.rm=TRUE), col='blue')
# legend('topleft',
#        legend=c('68-70; n=76',
#                 paste0('19; n=', length(which(!is.na(er_new))))),
#        col = c('sienna3','blue'),
#        lty = 1, bty = 'n', seg.len = 1, cex = 0.9, lwd = 2)
#
# #plot NEP dists, then and now
# plot(density(nep_68_70, na.rm=TRUE), bty='l', col='sienna3',
#     main='NEP 1968-70 vs. 2019', xlab='NEP', ylim = c(0,4.5))
# lines(density(nep_new, na.rm=TRUE), col='blue')
# legend('topleft',
#        legend=c('68-70; n=76',
#                 paste0('19; n=', length(which(!is.na(nep_new))))),
#        col = c('sienna3','blue'),
#        lty = 1, bty = 'n', seg.len = 1, cex = 0.9, lwd = 2)
#
# #plot GPP dists by year
# # cols = viridis(6)
# cols = c(rep('sienna3', 2), 'blue')
# plot(density(gpp_68, na.rm=TRUE),bty='l', col=cols[1],
#     main='GPP by year', xlab='GPP', ylim = c(0,3), xlim = c(-.1, 1))
# lines(density(gpp_69, na.rm=TRUE), col=cols[2])
# lines(density(gpp_new, na.rm=TRUE), col=cols[3])
# legend('topright',
#     legend=c(paste0('68; n=', length(gpp_68)),
#              paste0('69; n=', length(gpp_69)),
#              paste0('19; n=', length(which(!is.na(gpp_new))))),
#     col = cols, lty = 1, bty = 'n',
#     seg.len = 1, cex = 0.9, lwd = 2)
#
# #plot ER dists by year
# plot(density(er_68, na.rm=TRUE), bty='l', col=cols[1],
#     main='ER by year', xlab='ER', xlim = c(-1.2, .2))
# lines(density(er_69, na.rm=TRUE), col=cols[2])
# lines(density(er_new, na.rm=TRUE), col=cols[3])
# legend('topleft',
#     legend=c(paste0('68; n=', length(gpp_68)),
#              paste0('69; n=', length(gpp_69)),
#              paste0('19; n=', length(which(!is.na(er_new))))),
#     col = cols, lty = 1, bty = 'n',
#     seg.len = 1, cex = 0.9, lwd = 2)
#
# #plot NEP dists by year
# plot(density(nep_68, na.rm=TRUE), bty='l', col=cols[1],
#     main='NEP by year', xlab='NEP', ylim = c(0, 4.5), xlim = c(-.8, .2))
# lines(density(nep_69, na.rm=TRUE), col=cols[2])
# lines(density(nep_new, na.rm=TRUE), col=cols[3])
# legend('topleft',
#     legend = c(paste0('68; n=', length(gpp_68)),
#                paste0('69; n=', length(gpp_69)),
#                paste0('19; n=', length(which(!is.na(nep_new))))),
#     col=cols, lty=1, bty='n', seg.len=1, cex=0.9, lwd=2)
#
# dev.off()

#plot temporal coverage for historic data ####

historic_dates = nhc_68_70$date[! is.na(nhc_68_70$GPP)]
historic_year_agg = as.character(historic_dates)
substr(historic_year_agg, 1, 4) = '2020'
historic_year_agg = as.Date(historic_year_agg)
hy_num = as.numeric(historic_year_agg)

new_dates = nhc_new$date[!is.na(nhc_new$GPP)]
new_year_agg = as.character(new_dates)
substr(new_year_agg, 1, 4) = '2020'
new_year_agg = sort(as.Date(new_year_agg))
ny_num = as.numeric(new_year_agg)

tiff(filename = 'figures/seasonalcoverage.tif',
     width = 6, height = 2.5, res = 300, units = 'in')
# width = 6 * 800, height = 2.05 * 800, res = 800, units = 'px')

lims = c(min(ny_num), max(ny_num))

par(mfrow=c(2,1), mar=c(0,0,0,0), oma=c(3,4,2,3))

beanplot(ny_num, horizontal = TRUE, col = now_col, xaxt = 'n',
         frame.plot = FALSE, ylim = lims)
mtext('2019-20', 2, line = 1, adj = .7)
mtext(paste('n =', length(! is.na(ny_num))), 2, cex = .7)
mtext('Annual Sampling Coverage', 3, line = .1, cex = .8)

beanplot(hy_num, horizontal = TRUE, col=then_col, xaxt = 'n',
         frame.plot = FALSE, ylim = lims)
# mtext(paste('n =', length(! is.na(hy_num))), 2, cex = .7)
mtext(paste('n =', 57), 2, cex = .7)
mtext('1968-70', 2, line = 1, adj = 1)
axis(1, at=seq(as.Date('2020-01-01'), as.Date('2021-01-01'),
               length.out=13)[1:13], labels = FALSE)
axis(1, at=seq(as.Date('2020-01-01'), as.Date('2021-01-01'), length.out=13)[1:12],
     labels=month.abb, cex.axis = .8)

dev.off()

png(width = 6, height=2.5, units='in', type='cairo', res=300,
    filename='figures/seasonalcoverage.png')

lims = c(min(ny_num), max(ny_num))

par(mfrow=c(2,1), mar=c(0,0,0,0), oma=c(3,4,2,3))

beanplot(ny_num, horizontal = TRUE, col = now_col, xaxt = 'n',
         frame.plot = FALSE, ylim = lims)
mtext('2019-20', 2, line = 1, adj = .7)
mtext(paste('n =', length(! is.na(ny_num))), 2, cex = .7)
mtext('Annual Sampling Coverage', 3, line = .1, cex = .8)

beanplot(hy_num, horizontal = TRUE, col=then_col, xaxt = 'n',
         frame.plot = FALSE, ylim = lims)
# mtext(paste('n =', length(! is.na(hy_num))), 2, cex = .7)
mtext(paste('n =', 57), 2, cex = .7)
mtext('1968-70', 2, line = 1, adj = 1)
axis(1, at=seq(as.Date('2020-01-01'), as.Date('2021-01-01'),
               length.out=13)[1:13], labels = FALSE)
axis(1, at=seq(as.Date('2020-01-01'), as.Date('2021-01-01'), length.out=13)[1:12],
     labels=month.abb, cex.axis = .8)

dev.off()

# bootstrap some confidence bounds ####
# bootstrapped proportional to year, excluding september

n_68 = length(gpp_68_70[!is.na(gpp_68_70)])
n_new = length(er_new[!is.na(er_new)])
month_props <- c(31, 28, 31, 30, 31, 30, 31, 31, 0, 31, 30, 31)
names(month_props) <- c('01', '02', '03', '04', '05', '06',
                        '07', '08', '09', '10', '11', '12')
gpp_68_70_bymo = split(gpp_68_70[!is.na(gpp_68_70)],
                       factor(substr(nhc_68_70$date, 6, 7)))
er_68_70_bymo = split(er_68_70[!is.na(er_68_70)],
                      factor(substr(nhc_68_70$date, 6, 7)))
gpp_new_bymo = split(gpp_new[!is.na(gpp_new)],
                     factor(substr(dates_new[!is.na(gpp_new)], 6, 7)))
er_new_bymo = split(er_new[!is.na(er_new)],
                    factor(substr(dates_new[!is.na(er_new)], 6, 7)))

nhc_68_70$doy = as.numeric(format(nhc_68_70$date, "%j"))
nhc_new$doy = as.numeric(format(nhc_new$date, "%j"))

nsamp = 10000


summarize_bootstrap_met <- function(met_dat){
    # make an empty matrix for samples
    samps_mat <- matrix(NA_real_, ncol = nsamp, nrow = 365)
    row.names(samps_mat) <- as.character(1:365)
    samp_mat_GPP = samp_mat_ER = samps_mat

    # iterate taking samples for each day
    for(d in 1:365){
        dat_tmp <- met_dat %>%
            select(date, doy, GPP, ER) %>%
            mutate(day_diff = case_when(doy - d < 183 & doy - d >= 0 ~ doy - d,
                                        d - doy < 183 & d - doy > 0 ~ d - doy,
                                        365 - doy + d < 183 ~ 365 - doy + d,
                                        365 - d + doy < 183 ~ 365 - d + doy,
                                        TRUE ~ NA_real_),
                   prob = 1/(day_diff + 1)^2,
                   prob_norm = prob/sum(prob)) %>% arrange(doy) #%>% print(n = 100)
        samp_mat_GPP[d,] <- sample(x = dat_tmp$GPP, size = nsamp,
                                    replace = TRUE, prob = dat_tmp$prob_norm)
        samp_mat_ER[d,] <- sample(x = dat_tmp$ER, size = nsamp,
                                   replace = TRUE, prob = dat_tmp$prob_norm)
    }
    # summarize data for each day
    sum_dat <- data.frame(date = seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "day"),
               doy = as.numeric(row.names(samp_mat_GPP)),
               GPP_med = apply(samp_mat_GPP, 1, median),
               GPP_mean = apply(samp_mat_GPP, 1, mean),
               GPP_upper = apply(samp_mat_GPP, 1, function(x) quantile(x, probs = 0.975)),
               GPP_lower = apply(samp_mat_GPP, 1, function(x) quantile(x, probs = 0.025)),
               ER_med = apply(samp_mat_ER, 1, median),
               ER_mean = apply(samp_mat_ER, 1, mean),
               ER_upper = apply(samp_mat_ER, 1, function(x) quantile(x, probs = 0.975)),
               ER_lower = apply(samp_mat_ER, 1, function(x) quantile(x, probs = 0.025)))

    return(sum_dat)
}

historic_sum_dat <- summarize_bootstrap_met(nhc_68_70) %>%
    mutate(data = "CB_68_70")
new_sum_dat <- summarize_bootstrap_met(filter(nhc_new, !is.na(GPP) & !is.na(ER))) %>%
    mutate(data = "CB_2019")
nhc_17 <- summarize_bootstrap_met(filter(nhc_site, year == 2017 & !is.na(GPP) & !is.na(ER))) %>%
    left_join(select(filter(nhc_site, year == 2017), doy, GPP, ER)) %>% mutate(data = "NHC_2017")
nhc_18 <- summarize_bootstrap_met(filter(nhc_site, year == 2018 & !is.na(GPP) & !is.na(ER))) %>%
    left_join(select(filter(nhc_site, year == 2018), doy, GPP, ER)) %>% mutate(data = "NHC_2018")
nhc_19 <- summarize_bootstrap_met(filter(nhc_site, year == 2019 & !is.na(GPP) & !is.na(ER))) %>%
    left_join(select(filter(nhc_site, year == 2019), doy, GPP, ER)) %>% mutate(data = "NHC_2019")

bootstrapped_dat <- bind_rows(full_join(historic_sum_dat, select(nhc_68_70, doy, GPP, ER)),
          full_join(new_sum_dat, select(nhc_new, doy, GPP, ER)),
          nhc_17, nhc_18, nhc_19)

write_csv(bootstrapped_dat,
          "data/metabolism/compiled/bootstrapped/SM_met_means_bootstrapped_seasonally_proportions.csv")

bootstrapped_dat %>%
    mutate(month = month(date)) %>%
    group_by(data, month) %>%
    summarize(
        GPP_sd = sd(GPP_mean),
        GPP_mean = mean(GPP_mean),
        ER_sd = sd(ER_mean),
        ER_mean = mean(ER_mean)) %>%
    ggplot(aes(month, GPP_mean))+
    geom_line() +
    geom_line(aes(y = ER_mean)) +
    geom_ribbon(aes(ymin = GPP_mean - GPP_sd*1.96,
                    ymax = GPP_mean + GPP_sd*1.96), fill = "forestgreen", alpha = 0.3)+
    geom_ribbon(aes(ymin = ER_mean - ER_sd*1.96,
                    ymax = ER_mean + ER_sd*1.96), fill = "sienna1", alpha = 0.3) +
    facet_wrap(.~data)

bootstrapped_dat %>%
    group_by(data) %>% summarize(GPP_max = max(GPP, na.rm = T),
                                 GPP_min = min(GPP, na.rm = T),
                                 GPP_mean = mean(GPP, na.rm = T),
                                 ER_max = min(ER, na.rm = T),
                                 ER_min = max(ER, na.rm = T),
                                 ER_mean = mean(ER, na.rm = T))

bootstrapped_dat %>%
    ggplot(aes(doy, GPP_mean)) +
        geom_point(pch = 20) +
        geom_ribbon(aes(ymin = GPP_lower, ymax = GPP_upper), fill = "lightgreen", alpha = 0.5) +
        geom_point(aes(y = ER_med), pch = 20) +
        geom_ribbon(aes(ymin = ER_lower, ymax = ER_upper), fill = "sienna1", alpha = 0.5) +
        geom_point(aes(y = GPP), pch = 1) +
        geom_point(aes(y = ER), pch = 1) +
        facet_wrap(.~data) +
    theme_bw()

ggsave('figures/seasonally_bootstrapped_metabolism.png', width = 10, height = 6)


ggplot(historic_sum_dat, aes(doy, GPP_med)) +
    geom_point() +
    geom_ribbon(aes(ymin = GPP_lower, ymax = GPP_upper), fill = "lightgreen", alpha = 0.5) +
    geom_point(aes(y = ER_med)) +
    geom_ribbon(aes(ymin = ER_lower, ymax = ER_upper), fill = "sienna1", alpha = 0.5) +
    geom_point(data = nhc_68_70, aes(doy, GPP), pch = 1)+
    geom_point(data = nhc_68_70, aes(doy, ER), pch = 1)
ggplot(new_sum_dat, aes(doy, GPP_med)) +
    geom_point() +
    geom_ribbon(aes(ymin = GPP_lower, ymax = GPP_upper), fill = "lightgreen", alpha = 0.5) +
    geom_point(aes(y = ER_med)) +
    geom_ribbon(aes(ymin = ER_lower, ymax = ER_upper), fill = "sienna1", alpha = 0.5) +
    geom_point(data = nhc_new, aes(doy, GPP), pch = 1)+
    geom_point(data = nhc_new, aes(doy, ER), pch = 1)

# calculate cumulative metabolism:
bind_rows(historic_sum_dat, new_sum_dat, nhc_17, nhc_18, nhc_19) %>%
    group_by(data) %>%
    summarize(GPP_cum = sum(GPP_med),
              ER_cum = sum(ER_med))

bind_rows(historic_sum_dat, new_sum_dat, nhc_17, nhc_18, nhc_19) %>%
    group_by(data)

dd <- historic_sum_dat
get_max_windows <- function(dd){
    dd$rollER = dd$rollGPP = NA_real_

    for(i in 1:(nrow(dd)-9)){
        dd$rollER[i] <- sum(dd$ER_med[i:(i + 9)])
        dd$rollGPP[i] <- sum(dd$GPP_med[i:(i + 9)])
    }

    roll_dates <- data.frame(GPP = dd$date[which.max(dd$rollGPP)],
               ER = dd$date[which.min(dd$rollER)])

    return(roll_dates)
}
# make rolling 10 day window calcs.

get_max_windows(historic_sum_dat)
get_max_windows(new_sum_dat)
get_max_windows(nhc_17)
get_max_windows(nhc_18)
get_max_windows(nhc_19)


# bootstrapped method by year:

mean_vect_er_68_70 = mean_vect_er_new = mean_vect_gpp_68_70 =
    mean_vect_gpp_new = c()
turboset <- list(c(1:12), c(10,11), c(3,4), c(7,8), c(1,2))
names <- c("year", "oct_nov", "mar_apr", "jul_aug", "jan_feb")
CI_prop <- data.frame()



for(ss in 1:5){
    set <- turboset[[ss]]
    for(i in 1:nsamp){
        samp_68_70_er = samp_new_er =
            samp_68_70_gpp = samp_new_gpp = c()

        for(m in names(month_props[set])){
            nn <- month_props[m]/sum(month_props[set])
            t_er_68_70 <- er_68_70_bymo[[m]]
            t_gpp_68_70 <- gpp_68_70_bymo[[m]]
            t_er_new <- er_new_bymo[[m]]
            t_gpp_new <- gpp_new_bymo[[m]]

            samp_68_70_er = c(samp_68_70_er,
                              sample(t_er_68_70,
                                     size = round(nn * n_68, 0),
                                     replace=TRUE))
            samp_new_er = c(samp_new_er,
                            sample(t_er_new,
                                   size = round(nn * n_new, 0),
                                   replace=TRUE))
            samp_68_70_gpp = c(samp_68_70_gpp,
                               sample(t_gpp_68_70,
                                      size = round(nn * n_68, 0),
                                      replace=TRUE))
            samp_new_gpp = c(samp_new_gpp,
                                 sample(t_gpp_new,
                                        size = round(nn * n_new, 0),
                                        replace=TRUE))
            }
        mean_vect_er_68_70[i] = mean(samp_68_70_er)
        mean_vect_er_new[i] = mean(samp_new_er)
        mean_vect_gpp_68_70[i] = mean(samp_68_70_gpp)
        mean_vect_gpp_new[i] = mean(samp_new_gpp)
        if(i %% 1000 == 0){ print(i) }
    }

    CI_p = data.frame('CI95_lower'=numeric(4), 'median'=numeric(4),
        'CI95_upper'=numeric(4), 'met' = c(rep("GPP", 2), rep("ER", 2)),
        'era' = rep(c('then', 'now'), 2), 'prop' = names[ss])
    CI_p[1,1:3] = quantile(sort(mean_vect_gpp_68_70),
                           probs=c(0.025, 0.5, 0.975))
    CI_p[2,1:3] = quantile(sort(mean_vect_gpp_new),
                           probs=c(0.025, 0.5, 0.975))
    CI_p[3,1:3] = -quantile(sort(mean_vect_er_68_70),
                            probs=c(0.025, 0.5, 0.975))
    CI_p[4,1:3] = -quantile(sort(mean_vect_er_new),
                            probs=c(0.025, 0.5, 0.975))

    CI_prop <- bind_rows(CI_prop, CI_p)
    if(ss == 1){
        CI = data.frame(met = rep(c(rep("GPP", 3), rep("ER", 3)),2),
                        era = c(rep("now", 6), rep("then", 6)),
                        val = c(quantile(mean_vect_gpp_new,
                                         probs = c(.025, .5, .975)),
                                -quantile(mean_vect_er_new,
                                          probs = c(.025, .5, .975)),
                                quantile(mean_vect_gpp_68_70,
                                         probs = c(.025, .5, .975)),
                                -quantile(mean_vect_er_68_70,
                                          probs = c(.025, .5, .975))))
        }
}

write.csv(CI_prop, 'data/metabolism/compiled/bootstrapped/SM_met_means_bootstrapped_by_month_proportions.csv')

CI <- CI %>%
    mutate(era = factor(era, levels = c("then", "now")),
           met = factor(met, levels = c("GPP", "ER")),
           b = "Bootstrapped 95% CI")

tiff(filename = 'figures/bootstrapped_CI_daily_mean_SM.tif',
     height = 5, width = 2.2, res = 800, units = 'in')

    ggplot(CI, aes(met, val, fill = era)) +
        geom_boxplot(coef = 3)+
        # stat_boxplot(geom ='errorbar', width = .1,
        #              position = position_dodge(.75)) +
        facet_wrap(.~b) +
        scale_fill_manual(values = c(then_col, now_col)) +
        theme_bw() +
        xlab("") +
        ylim(0, 2.2) +
        ylab(expression("Mean Daily Metabolism (g O"[2] * "m"^-2 * "d"^-1 * ")")) +
        theme(legend.position = "none",
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
dev.off()

png(filename = 'figures/bootstrapped_CI_daily_mean_SM.png',
     height = 5, width = 2.2, res = 800, units = 'in')

    ggplot(CI, aes(met, val, fill = era)) +
        geom_boxplot(coef = 3)+
        facet_wrap(.~b) +
        scale_fill_manual(values = c(then_col, now_col)) +
        theme_bw() +
        xlab("") +
        ylim(0, 2.2) +
        ylab(expression("Mean Daily Metabolism (g O"[2] * "m"^-2 * "d"^-1 * ")")) +
        theme(legend.position = "none",
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
dev.off()

# numbers for results section
cbp_mon <- dat$preds %>%
    filter(site == "CBP") %>%
    select(month, era, GPP, ER) %>%
    mutate(ER = -ER) %>%
    pivot_longer(cols = c("GPP", "ER"), names_to = "met", values_to = 'gCm2d') %>%
    group_by(month, era, met) %>%
    summarize(across(all_of("gCm2d"),
                     .fns = list(mean = ~mean(., na.rm = T),
                                 sd = ~sd(., na.rm = T)),
                     # lower95 = ~quantile(., 0.025, na.rm = T)),
                     .names = '{col}_{fn}'),
              n = n())

mm = c(1,11,12)
mm = c(4,5)
mm = c(6,7)
mm = 3

cbp_mon %>% filter(month %in% mm) %>%
    group_by(era, met) %>%
    summarize(mean = mean(gCm2d_mean, na.rm = T),
              cv = sd(gCm2d_mean, na.rm = T)/mean)

dat$preds %>%
    filter(site == 'CBP',
           month %in% mm) %>%
    pivot_longer(cols = any_of(c('GPP', 'ER')),
                 names_to = 'met',
                 values_to = 'gCm2d_mean') %>%
    group_by(era, met) %>%
    summarize(mean = mean(gCm2d_mean, na.rm = T),
              cv = sd(gCm2d_mean, na.rm = T)/mean)

# png('figures/bootstrapped_CI_daily_mean_streamMetabolizer.png',
#     height = 6, width = 3, type = cairo,  res = 300, units = 'in')
#
#     par(mfrow = c(1,2), mar = c(3,2,4,1), oma = c(0,3,1,0))
#     tmp <- CI_prop %>%
#         filter(prop == "by hall sampling")
#     boxplot(t(tmp[,1:3]), col = c(alpha(then_col,.8),now_col ),
#             ylim = c(0.2, .9),
#             xaxt = "n", xlab = "")
#     axis(1, at=c(1.5, 3.5), labels=c("GPP", "ER"), line=0)
#     mtext("Proportional to Hall 1972 sampling")
#     mtext("CI around mean (g C/m2/d)", side = 2, line = 3)
#
#     tmp <- CI_prop %>%
#         filter(prop == "year")
#     boxplot(t(tmp[,1:3]), col = c(alpha(then_col,.8),now_col ),
#             ylim = c(0.2, 0.9),
#             xaxt = "n", xlab = "")
#     axis(1, at=c(1.5, 3.5), labels=c("GPP", "ER"), line=0)
#     mtext("Proportional to days per month")
#     legend("bottomright",cex=1, bty = "n",
#            c("then", "now"),
#            fill = c(alpha(then_col,.75),now_col ))
#     par(new = T, mfrow = c(1,1))
#     mtext('Bootstrapped daily metabolism estimates at Concrete Bridge',
#           line = 2, cex = 1.1)
# dev.off()

# png('figures/bootstrapped_CI_daily_mean_streamMetabolizer_byseason.png',
#     height = 4, width = 8, res = 300, units = 'in')
#     par(mfrow = c(1,4), mar = c(3,2,4,1), oma = c(0,3,1,0))
#     tmp <- CI_prop %>%
#         filter(prop == "jan_feb")
#     boxplot(t(tmp[,1:3]), col = c(alpha(then_col,.8),now_col ),
#             ylim = c(0, 1),
#             xaxt = "n", xlab = "")
#     axis(1, at=c(1.5, 3.5), labels=c("GPP", "ER"), line=0)
#     mtext("Jan - Feb")
#     mtext("CI around mean (g C/m2/d)", side = 2, line = 3)
#     # legend("bottom",cex=1, bty = "n",
#     #        c("then", "now"),
#     #        fill = c(alpha(then_col,.75),now_col ))
#
#     tmp <- CI_prop %>%
#         filter(prop == "mar_apr")
#     boxplot(t(tmp[,1:3]), col = c(alpha(then_col,.8),now_col ),
#             ylim = c(0, 1),
#             xaxt = "n", xlab = "")
#     axis(1, at=c(1.5, 3.5), labels=c("GPP", "ER"), line=0)
#     mtext("Mar - Apr")
#
#     tmp <- CI_prop %>%
#         filter(prop == "jul_aug")
#     boxplot(t(tmp[,1:3]), col = c(alpha(then_col,.8),now_col ),
#             ylim = c(0, 1),
#             xaxt = "n", xlab = "")
#     axis(1, at=c(1.5, 3.5), labels=c("GPP", "ER"), line=0)
#     mtext("Jul - Aug")
#
#     tmp <- CI_prop %>%
#         filter(prop == "oct_nov")
#     boxplot(t(tmp[,1:3]), col = c(alpha(then_col,.8),now_col ),
#             ylim = c(0, 1),
#             xaxt = "n", xlab = "")
#     axis(1, at=c(1.5, 3.5), labels=c("GPP", "ER"), line=0)
#     mtext("Oct - Nov")
#
#     par(new = T, mfrow = c(1,1))
#     mtext('Bootstrapped daily metabolism estimates at Concrete Bridge',
#           line = 3.4, cex = 1.1)
#
# dev.off()

# CI by months

