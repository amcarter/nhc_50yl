# assess the fit and K values for SM metabolism model runs
# 11/17/2020
library(ks)
library(zoo)
library(tidyverse)
library(streamMetabolizer)
#setwd('C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data')

# Cleaning and summarizing ####
filter_model <- function(fit, flow_dates, GPP_min = 0, ER_max = 0){
  q <- fit@data[,c(1,3,4,6,8)] %>%
    group_by(date) %>%
    summarize(discharge.daily = mean(discharge, na.rm = T),
              temp.min = min(temp.water, na.rm = T),
              temp.water = mean(temp.water, na.rm = T),
              DO.obs = mean(DO.obs, na.rm = T),
              DO.sat = mean(DO.sat, na.rm = T)) %>%
    ungroup() %>%
    left_join(flow_dates[,c(1,4)]) 
  preds <- get_fit(fit)$daily  %>%
    select(date, GPP = GPP_daily_50pct, 
           GPP.lower = GPP_daily_2.5pct, GPP.upper = GPP_daily_97.5pct,
           ER = ER_daily_50pct, 
           ER.lower = ER_daily_2.5pct, ER.upper = ER_daily_97.5pct,
           K600 = K600_daily_50pct, 
           K600_2.5 = K600_daily_2.5pct, K600_97.5 = K600_daily_97.5pct,
           GPP_Rhat = GPP_daily_Rhat, ER_Rhat = ER_daily_Rhat, 
           K600_Rhat = K600_daily_Rhat, errors) %>%
    left_join(q) %>%
    mutate(badGPP = case_when(good_flow == FALSE ~ 1,
                              GPP.upper < GPP_min ~ 1,
                              GPP_Rhat > 1.05 ~ 1,
                              K600_Rhat > 1.05 ~ 1,
                              is.na(GPP) ~ 1,
                              TRUE ~ 0),
           badER = case_when(good_flow == FALSE ~ 1,
                             ER.lower > ER_max ~ 1,
                             ER_Rhat > 1.05 ~ 1,
                             K600_Rhat > 1.05 ~ 1,
                             is.na(ER) ~ 1,
                             TRUE ~ 0))
  
  
  coverage <- data.frame(missing_data = sum(preds$errors != ""),
                         bad_flow = sum(preds$good_flow == FALSE, na.rm = T),
                         neg_GPP = sum(preds$GPP.upper < GPP_min, na.rm = T),
                         pos_ER = sum(preds$ER.lower > ER_max, na.rm = T),
                         bad_Rhat = sum(preds$K600_Rhat > 1.05 |
                                          preds$ER_Rhat > 1.05 |
                                          preds$GPP_Rhat > 1.05, na.rm = T),
                         lost_GPP = sum(preds$badGPP),
                         lost_ER = sum(preds$badER),
                         total_days = nrow(preds))
  preds <- preds %>% 
    mutate(across(starts_with("ER", ignore.case = F), 
                  ~ case_when(badER == 1 ~ NA_real_,
                              TRUE ~ .)),
           across(starts_with("GPP"), ~ case_when(badGPP == 1 ~ NA_real_,
                                                 TRUE ~ .)),
           across(starts_with("K600"), ~ case_when(badER == 1 ~ NA_real_,
                                                  TRUE ~ .)),
           GPP = case_when(GPP < 0 ~ 0,
                           TRUE ~ GPP),
           ER = case_when(ER > 0 ~ 0,
                          TRUE ~ ER)) %>%
    select(-badER, -badGPP)
        
  return(list(preds, coverage))

}

fill_summarize_met <- function(preds){
  w <- range(c(which(!is.na(preds$GPP)), which(!is.na(preds$ER))))
  preds <- preds[w[1]:w[2],] 
  met <- preds %>%
    summarize(gpp_mean = mean(GPP, na.rm = T),
              gpp_median = median(GPP, na.rm = T),
              gpp_max = max(GPP, na.rm = T),
              gpp_cv = sd(GPP, na.rm = T)/mean(GPP, na.rm = T),
              er_mean = -mean(ER, na.rm = T),
              er_median = -median(ER, na.rm = T),
              er_max = -min(ER, na.rm = T),
              er_cv = -sd(ER, na.rm = T)/mean(ER, na.rm = T))
  
  cum <- data.frame(date = seq(preds$date[1], 
                               preds$date[nrow(preds)], 
                               by = "day")) %>%
    as_tibble() %>%
    left_join(preds) %>%
    select(-ends_with("Rhat"), -starts_with("K600"), -good_flow, -errors) %>%
    mutate(across(-date, na.approx, na.rm = F)) %>%
    mutate(across(starts_with("GPP"), ~ case_when(. < 0 ~ 0,
                                                  TRUE ~ .))) %>%
    mutate(across(starts_with("ER"), ~ 
                    case_when(. > 0 ~ 0,
                              TRUE ~ -.))) %>%
    mutate(across(-date, cumsum, .names = "{col}_cum")) 
  n <- nrow(cum)
  l = (n-9)
  weekly <- tibble(date = cum$date[1:l],
                   GPP_week = rep(NA_real_, l),
                   ER_week = rep(NA_real_, l))
  for(i in 1:l){
    weekly$GPP_week[i] <- sum(cum$GPP[i:(i+9)]) 
    weekly$ER_week[i] <- sum(cum$ER[i:(i+9)]) 
  }
  met$gpp_max10d <- weekly$date[which.max(weekly$GPP_week)]
  met$er_max10d <- weekly$date[which.max(weekly$ER_week)]
  se <- sum(is.na(cum$GPP))
  met$gpp_cum <- cum$GPP_cum[n-se]*365/(n-se)
  se <- sum(is.na(cum$ER))
  met$er_cum <- cum$ER_cum[n-se]*365/(n-se)
  met$daterange <- as.character(paste(cum$date[1], "-", cum$date[nrow(cum)]))
  met$pctcoverage <- sum(!is.na(preds$GPP))/nrow(cum)
  
  return(list(cum = cum, met = met))
}


# plots####
plot_zoom <- function(dat, vars = c("DO.obs", "DO.mod"),
                      ordby = "solar.time"){
  ob <- dat[,ordby, drop = T]
  tt <- dat %>%
    select(any_of(vars)) %>%
    xts(order.by = ob)
  if("DO.mod" %in% vars){
    tt %>%
      dygraph() %>%
      dySeries("DO.mod", drawPoints = TRUE, strokeWidth = 0) %>%
      dyRangeSelector()
  } else {
    tt %>%
      dygraph() %>%
      dyRangeSelector()
  }
}
plot_rhats <- function(preds){
  rh <- preds %>%
    select(date, ends_with('Rhat'))
  ylim = range(c(rh$K600_Rhat, rh$GPP_Rhat, rh$ER_Rhat, 1.05), 
               na.rm = T)
  plot(x = rh$date, y = rh$K600_Rhat, 
       ylab = "rhat (convergence metric)", xlab = "date",
       type = "l", lwd = 1.5, ylim = ylim)        
  lines(rh$date, rh$GPP_Rhat, col = "forestgreen", lwd = 1.5)
  lines(rh$date, rh$ER_Rhat, col = "sienna", lwd = 1.5)
  abline(h = 1.05, lty = 2, col = "red")
  mtext("Rhat below 1.05 is good", 3, 0, adj = 0, cex = .8)
  legend("topleft", 
         legend = c("K600", "GPP", "ER"),
         col = c(1, "forestgreen", "sienna"),
         lty = 1, bty = "n", lwd = 1.5)
}
plot_binning <- function(fit, preds){
  mm_fit <- get_fit(fit)
  
  SM_output <- mm_fit$daily
  SM_KQbin <-  mm_fit$KQ_binned
  SM_day <- get_data_daily(fit)
  SM_specs <- get_specs(fit)
  
  day <- data.frame(SM_day$discharge.daily, 
                    SM_output$K600_daily_50pct, 
                    SM_output$GPP_50pct,
                    SM_output$K600_daily_Rhat,
                    rep('daily', dim(SM_output)[1]))
  
  colnames(day)<-c('Q', 'K600', 'GPP','Rhat', 'Group')
  
  # gg<-ggplot(day, aes(x=log(Q), y = GPP, col=Rhat))+
  #   geom_point() +
  #   geom_hline(yintercept = 0.2)
  
  nodes<-data.frame(exp(SM_specs$K600_lnQ_nodes_centers), 
                    exp(SM_KQbin$lnK600_lnQ_nodes_50pct),
                    exp(SM_KQbin$lnK600_lnQ_nodes_2.5pct),
                    exp(SM_KQbin$lnK600_lnQ_nodes_97.5pct),
                    exp(SM_specs$K600_lnQ_nodes_meanlog)) 
  prior_sd <- 1.96*exp(SM_specs$K600_lnQ_nodes_sdlog[1])
  colnames(nodes)<-c('Q', 'K600','K600_2.5', 'K600_97.5',  'K600_prior')
  pm <- par()$mar
  plot(log(preds$discharge.daily), preds$K600,
       xlab = "logQ", ylab = "K600", ylim = c(0, 25),
       col = "grey", pch = 19)
  polygon(log(c(nodes$Q, rev(nodes$Q))), 
          c(nodes$K600_prior + prior_sd, rev(nodes$K600_prior - prior_sd)),
          col = alpha('brown3', 0.2), border=NA)
  points(log(preds$discharge.daily), preds$K600,
         col = "grey25", pch = 19)
  points(log(nodes$Q), nodes$K600_prior, col = "brown3", cex = 1.5)
  points(log(nodes$Q), nodes$K600, col = "brown3", pch = 19, cex = 1.5)
  arrows(log(nodes$Q), nodes$K600_2.5, log(nodes$Q), nodes$K600_97.5, 
         length = 0, col = 'brown3', lwd = 2)
  par(new = T, mar = c(0,0,0,0))
  plot(1,1, type = 'n', axes = FALSE, xlab = '', ylab = '')
  legend("top", legend = c("data", "prior", "posterior"),
         col = c("grey25", "brown3", "brown3"),
         pch = c(19, 1, 19), cex = 1.2,  bty = 'n', ncol = 3)
  par(mar = pm)
}
plot_binning2 <- function(fit, preds){
  mm_fit <- get_fit(fit)
  
  SM_output <- mm_fit$daily
  SM_KQbin <-  mm_fit$KQ_binned
  SM_day <- get_data_daily(fit)
  SM_specs <- get_specs(fit)
  
  day <- data.frame(SM_day$discharge.daily, 
                    SM_output$K600_daily_50pct, 
                    SM_output$GPP_50pct,
                    SM_output$K600_daily_Rhat,
                    rep('daily', dim(SM_output)[1]))
  
  colnames(day)<-c('Q', 'K600', 'GPP','Rhat', 'Group')
  
  # gg<-ggplot(day, aes(x=log(Q), y = GPP, col=Rhat))+
  #   geom_point() +
  #   geom_hline(yintercept = 0.2)
  
  nodes<-data.frame(exp(SM_specs$K600_lnQ_nodes_centers), 
                    exp(SM_KQbin$lnK600_lnQ_nodes_50pct),
                    exp(SM_KQbin$lnK600_lnQ_nodes_2.5pct),
                    exp(SM_KQbin$lnK600_lnQ_nodes_97.5pct),
                    exp(SM_specs$K600_lnQ_nodes_meanlog)) 
  prior_sd <- 1.96*exp(SM_specs$K600_lnQ_nodes_sdlog[1])
  colnames(nodes)<-c('Q', 'K600','K600_2.5', 'K600_97.5',  'K600_prior')
  plot(log(preds$discharge.daily), preds$K600,
       xlab = "", ylab = "", ylim = c(0, 25),
       col = "grey", pch = 20, cex = 0.5, cex.axis = 0.6)

  mtext(expression(paste("K600 (d"^"-1"*")")), 2, line = 2, cex = 0.8)
  mtext(expression(paste("log discharge (m"^"3"~"s"^"-1"*")")), 1, cex = 0.8, 
        line = 2)
  legend(-3, 31, legend = c("prior", "posterior", "data"),
         col = c( "brown3", "brown3", "grey25"), xpd = T,
         pch = c(1, 19, 20), cex = 0.7,  bty = 'n', ncol = 3)
  
  polygon(log(c(nodes$Q, rev(nodes$Q))), 
          c(nodes$K600_prior + prior_sd, rev(nodes$K600_prior - prior_sd)),
          col = alpha('brown3', 0.2), border=NA)
  points(log(preds$discharge.daily), preds$K600,
         col = "grey25", pch = 20, cex = 0.6)
  points(log(nodes$Q), nodes$K600_prior, col = "brown3", cex = 0.8)
  points(log(nodes$Q), nodes$K600, col = "brown3", pch = 19, cex = 0.8)
  arrows(log(nodes$Q), nodes$K600_2.5, log(nodes$Q), nodes$K600_97.5, 
         length = 0, col = 'brown3', lwd = 2)
}
plot_KvER <- function(preds){
  pcor <- round(cor(preds$ER, preds$K600, 
                    use = "na.or.complete"), 2)
  plot( preds$K600,-preds$ER, xlab = "", ylab = "",
        pch = 20, cex = 0.8, cex.axis = 0.6,  col = 'grey25')
  mtext(expression(paste("K600 (d"^"-1"*")")), 1, line = 2, cex = 0.8)
  mtext(expression(paste("Ecosystem Respiration (g"~O[2]~"m"^"-2"~"d"^"-1"*")")), 
        2, line = 1.9, cex = 0.8)

  mtext(paste0("pearson's correlation = ", pcor),
        side = 3, line = 0, adj = 1, cex = .8)
}


plot_metab <- function(met, ylim = NULL, xlim = NULL, doy = F, main = "", error = T, 
                       xaxt = NULL, yaxt = NULL){

  if(error){
    yrange = range(c(met$GPP.upper, met$ER.lower), na.rm = T)
  } else {
    yrange = range(c(met$GPP, met$ER), na.rm = T)
  }
  if(!is.null(ylim)){yrange = ylim}
  if(!is.null(xlim)){
    xlim = range(met$date)
  }
  if(doy == T){
    met <- met %>%
      mutate(date = doy) %>%
      arrange(date)
    }
  plot(met$date, met$GPP, main = main,
       type = "l", lwd = 2, col = "forestgreen", xaxt = xaxt, yaxt = yaxt,
       ylim = yrange, xlim = xlim, xlab = "date", ylab = "gO2/m2/d",
       frame.plot = T)  
  lines(met$date, met$ER, 
        lwd = 2, col = "sienna")
  if(error){
    polygon(na.approx(c(met$date, rev(met$date)), na.rm = F), 
            na.approx(c(met$GPP.lower, rev(met$GPP.upper)), na.rm = F),
            col = alpha("forestgreen", 0.4), border = NA)
    polygon(na.approx(c(met$date, rev(met$date)), na.rm = F), 
            na.approx(c(met$ER.lower, rev(met$ER.upper)), na.rm = F),
            col = alpha("sienna", 0.4), border = NA)
  }
  abline(h = 0)
  
}
plot_hall_metab <- function(met, ylim = NULL, 
                            site = c("CBP", "WB", "UNHC"), doy = F, error = T) {

  ss <- data.frame(hall_site = c("Concrete", "Wood Bridge", "Blackwood"),
                   site = c("CBP", "WB", "UNHC"))
  if(!site %in% ss$site){
    print(paste(site, "not in Hall data"))
    break
  }
  
  hall <- read_csv("C:/Users/Alice Carter/git/nhc_50yl/data/hall/hall_table_15.csv") %>%
    rename(hall_site = site) %>%
    mutate(ER_gO2m2d = -ER_gO2m2d) %>%
    left_join(ss) %>%
    select(-hall_site, -newdate) %>%
    filter(site %in% !!site) %>%
    mutate(doy = format(as.Date(date, format = "%m/%d/%Y"), "%j")) %>%
    select(-date)
  
  if(!"doy" %in% colnames(met)){
    met <- met %>%
      mutate(doy = format(date, "%j"))
  }
  
  met <- met %>%
    full_join(hall, by = "doy")
  
  if(is.null(ylim)){
    if(error){
      ylim <- range(c(met$GPP.upper, met$ER.lower, 
                    met$GPP_gO2m2d, met$ER_gO2m2d), na.rm = T)
    } else{
      ylim <- range(c(met$GPP, met$ER, 
                    met$GPP_gO2m2d, met$ER_gO2m2d), na.rm = T)
    }
  }
  
  plot_metab(met, ylim, doy = doy, error = error)
  
  if(doy == T){
    met <- met %>%
      mutate(date = doy) %>%
      arrange(date)
  }
  
  points(met$date, met$GPP_gO2m2d, pch = 19, col = "forestgreen")
  points(met$date, met$ER_gO2m2d, pch = 19, col = "sienna")
}

K600toO2<-function(temp, K600) {
  ((600/(1800.6-(120.1*temp)+(3.7818*temp^2)-(0.047608*temp^3)))^0.5)*K600
  }

plot_k <- function(preds, xlim = NULL){
  hall <- read_csv("C:/Users/Alice Carter/git/nhc_50yl/data/hall/hall_tableA2_k_morphology_extended.csv") 
  k = preds$K600
  if(is.null(xlim)){
    xlim <- range(c(k, hall$K600), na.rm = T)
  } else {xlim = c(0, xlim)}
  plot(density(k, na.rm = T), xlab = "K600 (day-1)", 
       main = "K values", lwd = 2,  xlim = c(xlim[1]-1, xlim[2]+1))
  par(new = T)
  plot(density(hall$K600, na.rm = T), xlab = "",xaxt = "n", yaxt = "n", 
       main = "", lwd = 2,col = "brown3", xlim = c(xlim[1]-1, xlim[2]+1))
  legend("topright", bty = "n", lty = 1, 
         c("now","hall"), col = c(1, "brown3"))
}

plot_kde_metab <- function(met, lim = NULL, col = "grey25"){
  
  kernel <- kde(na.omit(met[,c("GPP","ER")]))
  if(is.null(lim)) {
    lim <- quantile(c(met$GPP, -met$ER), .99, na.rm = T) 
  } 
  
  plot(kernel, xlab = "GPP (gO2m2d)", ylab = "ER (gO2m2d)", 
       ylim = c(-lim, 0), xlim = c(0, lim), 
       display = "filled.contour",
       cont=c(30,60,90), #lwd = 1,
       col = c(NA, 
               alpha(col, .25), 
               alpha(col, .5), 
               alpha(col, .75)))
       
  abline(0,-1)
}
plot_kde_hall_metab <- function(met, lim = NULL, 
                                site = c("CBP", "WB", "UNHC"), legend = T){
  
  ss <- data.frame(hall_site = c("Concrete", "Wood Bridge", "Blackwood"),
                   site = c("CBP", "WB", "UNHC"))
  if(!site %in% ss$site){
    print(paste(site, "not in Hall data"))
    break
  }
  
  hall <- read_csv("C:/Users/Alice Carter/git/nhc_50yl/data/hall/hall_table_15.csv") %>%
    rename(hall_site = site) %>%
    mutate(ER_gO2m2d = -ER_gO2m2d) %>%
    left_join(ss) %>%
    select(-hall_site, -newdate) %>%
    filter(site %in% !!site)
    
  kernel_hall <- kde(na.omit(hall[,c("GPP_gO2m2d","ER_gO2m2d")]))
  
  if(is.null(lim)) {
    lim <- quantile(c(hall$GPP_gO2m2d, -hall$ER_gO2m2d, 
                   met$GPP, -met$ER), .99, na.rm = T) 
  }
  
  plot_kde_metab(met, lim)
  par(new=T)
  plot(kernel_hall, 
       display = "filled.contour",
       xlab = "", ylab = "", 
       ylim = c(-lim, 0), xlim = c(0, lim),
       cont = c(30,60,90), 
       #lwd = 1,
       col = c(NA, 
               alpha("darkred", .25),                 
               alpha("darkred", .5), 
               alpha("darkred", .75)))
  abline(0,-1)
  if(legend){
  legend("topright", cex = 1.4,
         c(paste0("2019  n = ",nrow(na.omit(met))),
           paste0("1969  n = ", nrow(na.omit(hall)))),
         fill = c(alpha("grey25", .75), alpha("darkred", .75)), 
         border = NA, bty = "n")
  }
   
}
plot_diagnostics <- function(fit, preds, site, ylim = NULL, lim = NULL){
  m <- rbind(c(1,1,1,1,2,2),
             c(3,3,4,4,5,5))
  layout(m)
  par(mar = c(4,4,2,2))

  plot_hall_metab(preds, ylim = ylim)
  plot_kde_hall_metab(preds, lim = lim)
  plot_rhats( preds)
  plot_binning(fit, preds)
  plot_KvER(preds)
  mtext(site, outer = T, line = -1.5)
}

#site <- "CBP"

# Look at the metabolism predictions and fits for each year of data
# fit <- readRDS(paste0("metabolism/modeled/", site, "_nreg_v1.rds"))
# plot_metab_preds(predict_metab(fit))
# plot_binning(fit)
# plot_rhats(fit)
# plot_KvER(fit)
# 
# plot_DO_preds(fit, style = "dygraphs", y_var = "conc")
# 
# get_fit(fit)$overall %>%
#   select(ends_with('Rhat'))


# old plots ####
# plot_rhats <- function(fit){ 
#   rh <- get_fit(fit)$daily %>%
#     select(date, ends_with('Rhat'))
#   ylim = range(c(rh$K600_daily_Rhat, rh$GPP_daily_Rhat, rh$ER_daily_Rhat, 1.05), 
#                na.rm = T)
#   plot(x = rh$date, y = rh$K600_daily_Rhat, 
#        ylab = "rhat (convergence metric)", xlab = "date",
#        type = "l", lwd = 1.5, ylim = ylim)        
#   lines(rh$date, rh$GPP_daily_Rhat, col = "forestgreen", lwd = 1.5)
#   lines(rh$date, rh$ER_daily_Rhat, col = "sienna", lwd = 1.5)
#   abline(h = 1.05, lty = 2, col = "red")
#   mtext("Rhat below 1.05 is good", 3, 0, adj = 0, cex = .8)
#   legend("topleft", 
#          legend = c("K600", "GPP", "ER"),
#          col = c(1, "forestgreen", "sienna"),
#          lty = 1, bty = "n", lwd = 1.5)
# }
# plot_binning <- function(fit){
#   mm_fit <- get_fit(fit)
#   
#   SM_output <- mm_fit$daily
#   SM_KQbin <-  mm_fit$KQ_binned
#   SM_day <- get_data_daily(fit)
#   SM_specs <- get_specs(fit)
#   
#   day <- data.frame(SM_day$discharge.daily, 
#                     SM_output$K600_daily_50pct, 
#                     SM_output$GPP_50pct,
#                     SM_output$K600_daily_Rhat,
#                     rep('daily', dim(SM_output)[1]))
#   
#   colnames(day)<-c('Q', 'K600', 'GPP','Rhat', 'Group')
#   
#   # gg<-ggplot(day, aes(x=log(Q), y = GPP, col=Rhat))+
#   #   geom_point() +
#   #   geom_hline(yintercept = 0.2)
#   
#   nodes<-data.frame(exp(SM_specs$K600_lnQ_nodes_centers), 
#                     exp(SM_KQbin$lnK600_lnQ_nodes_50pct),
#                     exp(SM_specs$K600_lnQ_nodes_meanlog))
#   colnames(nodes)<-c('Q', 'K600', 'K600_prior')
#   pm <- par()$mar
#   plot(log(day$Q), day$K600,
#        xlab = "logQ", ylab = "K600", 
#        col = "grey", pch = 19, cex = 1.5)
#   points(log(nodes$Q), nodes$K600, col = "brown3", pch = 19, cex = 1.5)
#   points(log(nodes$Q), nodes$K600_prior, col = "brown3", cex = 1.5)
#   par(new = T, mar = c(0,0,0,0))
#   plot(1,1, type = 'n', axes = FALSE, xlab = '', ylab = '')
#   legend("top", legend = c("daily", "prior", "posterior"),
#          col = c("grey", "brown3", "brown3"),
#          pch = c(19, 1, 19), cex = 1,  bty = 'n', ncol = 3)
#   par(mar = pm)
# }
# plot_KvER <- function(fit){
#   KvER <- get_fit(fit)
#   pcor <- round(cor(KvER$daily$K600_daily_mean, 
#                     KvER$daily$ER_daily_mean, 
#                     use = "na.or.complete"),2)
#   plot(-KvER$daily$ER_daily_mean, KvER$daily$K600_daily_mean,
#        xlab = "ER (gO2/m2/d)", ylab = "K600 (/d)")
#   mtext(paste0("pearson's correlation = ", pcor),
#         side = 3, line = 0, adj = 1, cex = .8)
# }