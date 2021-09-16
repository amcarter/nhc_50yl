
# R helpers for NHC_2019_metabolism project

# run length encoder: ####
# for finding start and end points of NA sequences for interpolating discharge

rle_custom = function(x){
  r = rle(x)
  ends = cumsum(r$lengths)
  r = as.data.frame(cbind(values=r$values,
                         starts=c(1, ends[-length(ends)] + 1),
                         stops=ends, lengths=r$lengths, deparse.level=1))
  return(r)
}



# Drift corrector ####
# yd is a dataframe with a column of the data you want corrected and 
# a column of the point measurements for correction.
drift_correct <-  function(yd,  colname1, colname2){ 
  gap <- rle_custom(is.na(yd[,colname1, drop = T])) %>%
    filter(values == 0)
  ends <- data.frame(pt = c(gap$starts, gap$stops),
                     dat = FALSE) 
  #zero <- rle_custom(yd[,colname1] <= 0)
  ysi_pts <- data.frame(pt = which(!is.na(yd[, colname2])),
                        dat = TRUE) %>%
    bind_rows(ends) %>%
    arrange(pt)
  n <- nrow(yd)
  
  for(i in 1:(nrow(ysi_pts)-1)){
    a = ysi_pts$pt[i]
    b = ysi_pts$pt[i+1]
    if(a == b) next
    meas_p1 <- yd[a, colname2, drop = TRUE]
    sens_p1 <- yd[a, colname1, drop = TRUE]
    meas_p2 <- yd[b, colname2, drop = TRUE]
    sens_p2 <- yd[b, colname1, drop = TRUE]
    if(is.na(meas_p1) & is.na(meas_p2)) next
    
    if(ysi_pts$dat[i]){ # the first point is a measurement
      if(is.na(sens_p1)) next
      if(ysi_pts$dat[i+1]){ # the second point is a measurement
        if(is.na(sens_p2)) next
        #interpolate between the two, shift the rest of the ts to match
        delt <- data.frame(pt = a:b, 
                           dd = NA)
        delt$dd[1] <- meas_p1 - sens_p1
        delt$dd[nrow(delt)] <- meas_p2 - sens_p2
        delt$dd <- na.approx(delt$dd)      
        
        yd[a:b, colname1] <- yd[a:b, colname1] + delt$dd
        
        yd[(b+1):n, colname1] <- 
          yd[(b+1):n, colname1] + delt$dd[nrow(delt)]
      } else{ # the second point is a sensor break
        
        yd[a:n, colname1] <- yd[a:n, colname1] + meas_p1 - sens_p1
        
      }
    } else { # the first point is a sensor break
      if(ysi_pts$dat[i+1]){ # the second point is data
        yd[a:n, colname1] <- yd[a:n, colname1] + meas_p2 - sens_p2
        
      } #don't need an else, if the 2nd point is also a sensor break, do nothing
    }
  }
  return(yd[, colname1, drop = T])
}

#################
# #PlotDOPreds
# Al_plot_DO_preds <- function (DO_preds, y_var = c("conc", "pctsat", "ddodt"), 
#           style = c("ggplot2", "dygraphs"), y_lim = list(conc = c(NA, NA), pctsat = c(NA, NA), ddodt = c(NA, NA)), 
#           date_start = NA, date_end = NA, use_saved = TRUE) 
# {
#   if (is(DO_preds, "metab_model")) {
#     DO_preds <- predict_DO(DO_preds, date_start = date_start, 
#                            date_end = date_end, use_saved = use_saved)
#   }
#   # style <- match.arg(style)
#   # y_var <- match.arg(y_var, several.ok = TRUE)
#   params <- list(xlab = "Local time", ylab = "Predictions (lines) and observations (points)", 
#                  colors = list(conc = c("#CE9C59", "#A64B00","#FF7400"), 
#                                pctsat = c("#7CA586", "#007929","#23BC47"), 
#                                ddodt = c("#4A5869", "#05326D","#4282D3")))
#   DO.obs <- DO.pure <- DO.mod <- DO.sat <- ".dplyr.var"
#   DO_preds_conc <- mutate(DO_preds, as = "conc", var = "DO (mg/L)", 
#                           lab = "DO (mg/L)", col.pure = params$colors$conc[1], 
#                           col.mod = params$colors$conc[2], col.obs = params$colors$conc[3], 
#                           pure = if (exists("DO.pure", DO_preds)) 
#                             DO.pure
#                           else NA, mod = DO.mod, obs = DO.obs)
#   DO_preds_pctsat <- mutate(DO_preds, as = "pctsat", 
#                             var = "DO (% sat)", lab = "DO (% sat)", col.pure = params$colors$pctsat[1], 
#                             col.mod = params$colors$pctsat[2], col.obs = params$colors$pctsat[3], 
#                             pure = if (exists("DO.pure", DO_preds)) 
#                               100 * DO.pure/DO.sat
#                             else NA, mod = 100 * DO.mod/DO.sat, obs = 100 * DO.obs/DO.sat)
#   DO_preds_ddodt <- mutate(DO_preds[-1, ], as = "ddodt", 
#                            var = "dDO/dt (mg/L/d)", lab = "dDO/dt~(mg~L^-1~d^-1)", 
#                            col.pure = params$colors$ddodt[1], col.mod = params$colors$ddodt[2], 
#                            col.obs = params$colors$ddodt[3], pure = if (exists("DO.pure", 
#                                                                                DO_preds)) 
#                              diff(DO_preds$DO.pure)/as.numeric(diff(DO_preds$solar.time), 
#                                                                units = "days")
#                            else NA, mod = diff(DO_preds$DO.mod)/as.numeric(diff(DO_preds$solar.time), 
#                                                                            units = "days"), obs = diff(DO_preds$DO.obs)/as.numeric(diff(DO_preds$solar.time), 
#                                                                                                                                    units = "days")) %>% mutate(pure = ifelse(diff(DO_preds$date) == 
#                                                                                                                                                                                0, pure, NA), mod = ifelse(diff(DO_preds$date) == 0, 
#                                                                                                                                                                                                           mod, NA), obs = ifelse(diff(DO_preds$date) == 0, obs, 
#                                                                                                                                                                                                                                  NA))
#   var <- ".dplyr.var"
#   DO_preds_all <- bind_rows(DO_preds_conc, DO_preds_pctsat, 
#                             DO_preds_ddodt) %>% mutate(var = ordered(var, c(conc = "DO (mg/L)", 
#                                                                             pctsat = "DO (% sat)", ddodt = "dDO/dt (mg/L/d)")[y_var]))
#   plot_out <- switch(style, ggplot2 = {
#     if (!requireNamespace("ggplot2", quietly = TRUE)) stop("call install.packages('ggplot2') before plotting with style='ggplot2'")
#     . <- solar.time <- pure <- mod <- date <- col.pure <- col.mod <- col.obs <- obs <- ".ggplot.var"
#     preds_ggplot <- v(DO_preds_all) %>% filter(as %in% y_var)
#     if ("conc" %in% names(y_lim)) {
#       lim <- y_lim[["conc"]][1]
#       if (!is.na(lim)) preds_ggplot <- filter(preds_ggplot, 
#                                               as != "conc" | (pure >= lim & mod >= lim & 
#                                                                 obs >= lim))
#       lim <- y_lim[["conc"]][2]
#       if (!is.na(lim)) preds_ggplot <- filter(preds_ggplot, 
#                                               as != "conc" | (pure <= lim & mod <= lim & 
#                                                                 obs <= lim))
#     }
#     if ("pctsat" %in% names(y_lim)) {
#       lim <- y_lim[["pctsat"]][1]
#       if (!is.na(lim)) preds_ggplot <- filter(preds_ggplot, 
#                                               as != "pctsat" | (pure >= lim & mod >= 
#                                                                   lim & obs >= lim))
#       lim <- y_lim[["pctsat"]][2]
#       if (!is.na(lim)) preds_ggplot <- filter(preds_ggplot, 
#                                               as != "pctsat" | (pure <= lim & mod <= 
#                                                                   lim & obs <= lim))
#     }
#     if ("ddodt" %in% names(y_lim)) {
#       lim <- y_lim[["ddodt"]][1]
#       if (!is.na(lim)) preds_ggplot <- filter(preds_ggplot, 
#                                               as != "ddodt" | (pure >= lim & mod >= lim & 
#                                                                  obs >= lim))
#       lim <- y_lim[["ddodt"]][2]
#       if (!is.na(lim)) preds_ggplot <- filter(preds_ggplot, 
#                                               as != "ddodt" | (pure <= lim & mod <= lim & 
#                                                                  obs <= lim))
#     }
#     g <- ggplot2::ggplot(preds_ggplot, ggplot2::aes(x = solar.time, 
#                                                     group = date))
#     if (any(!is.na(preds_ggplot$pure))) g <- g + ggplot2::geom_line(ggplot2::aes(y = pure, 
#                                                                                  color = col.pure), size = 0.8, na.rm = TRUE)
#     g + ggplot2::geom_point(ggplot2::aes(y = obs, color = col.obs), 
#                             alpha = 0.6, na.rm = TRUE) + ggplot2::geom_line(ggplot2::aes(y = mod, 
#                                                                                          color = col.mod), size = 0.8, na.rm = TRUE) + ggplot2::scale_color_identity(guide = "none") + 
#       ggplot2::theme_bw() + ggplot2::facet_grid(var ~ ., 
#                                                 scales = "free_y") + ggplot2::xlab(params$xlab) + 
#       ggplot2::ylab(params$ylab)
#   }, dygraphs = {
#     if (!requireNamespace("dygraphs", quietly = TRUE)) stop("call install.packages('dygraphs') before plotting with style='dygraphs'")
#     if (!requireNamespace("xts", quietly = TRUE)) stop("call install.packages('xts') before plotting with style='dygraphs'")
#     . <- ".dplyr.var"
#     preds_xts <- DO_preds_all %>% filter(as %in% y_var) %>% 
#       arrange(solar.time) %>% group_by(date) %>% do(., 
#                                                     {
#                                                       out <- .[c(seq_len(nrow(.)), nrow(.)), ]
#                                                       out[nrow(.) + 1, c("pure", "mod", 
#                                                                          "obs")] <- NA
#                                                       out
#                                                     }) %>% ungroup()
#     prep_dygraph <- function(y_var) {
#       . <- solar.time <- pure <- mod <- obs <- ".dplyr.var"
#       prepped <- preds_xts %>% filter(as == y_var) %>% 
#         select(pure, mod, obs, solar.time) %>% 
#         mutate(solar.time = lubridate::force_tz(solar.time,Sys.getenv("TZ"))) %>% 
#         xts::xts(x = select(.,tzone = Sys.getenv("TZ"))
#       if (all(is.na(prepped[, "pure"]))) prepped <- prepped[, c("mod", "obs")]
#       prepped
#       }
#     if (length(y_var) > 1) {
#       y_var <- y_var[1]
#       warning("can only plot one dygraph y_var at a time for now; plotting ", 
#               y_var)
#     }
#     y_var_long <- preds_xts %>% filter(as == y_var) %>% slice(1) %>% 
#       .[["var"]] %>% as.character()
#     y_var_col <- params$colors[[y_var]]
#     dat <- xts::xts(preds_xts%>%select(DO.obs, DO.mod), order.by=preds_xts$solar.time)
#     ymin <- max(c(min(c(unclass(dat)), na.rm = TRUE), y_lim[[y_var]][1]), 
#                 na.rm = TRUE)
#     ymax <- min(c(max(c(unclass(dat)), na.rm = TRUE), y_lim[[y_var]][2]), 
#                 na.rm = TRUE)
#     d <- dygraphs::dygraph(dat, xlab = params$xlab, ylab = y_var_long, 
#                            group = "plot_DO_preds")
#     if (ncol(dat) == 3) d <- d %>% dygraphs::dySeries("pure", 
#                                                       drawPoints = FALSE, label = paste0("Pure ", 
#                                                                                          y_var_long), color = y_var_col[1])
#     d %>% dygraphs::dySeries("DO.mod", drawPoints = FALSE, 
#                              label = paste0("Modeled ", y_var_long), color = y_var_col[2]) %>% 
#       dygraphs::dySeries("DO.obs", drawPoints = TRUE, 
#                          strokeWidth = 0, label = paste0("Observed ", 
#                                                          y_var_long), color = y_var_col[3]) %>% dygraphs::dyAxis("y", 
#                                                                                                                  valueRange = (c(ymin, ymax) + (ymax - ymin) * c(-0.05, 
#                                                                                                                                                                  0.15))) %>% dygraphs::dyOptions(colorSaturation = 1) %>% 
#       dygraphs::dyLegend(labelsSeparateLines = TRUE, width = 300) %>% 
#       dygraphs::dyRangeSelector(height = 20)
#   })
#   plot_out
# }
# 
# 
# # compress streammetabolizer output into smaller format