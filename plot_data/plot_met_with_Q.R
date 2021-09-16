library(pracma)
# setwd("C:/Users/Alice Carter/Dropbox (Duke Bio_Ea)/projects/NHC_2019_metabolism/data")
source("src/metabolism/inspect_model_fits.r")


dat <- readRDS("data/metabolism/compiled/met_preds_stream_metabolizer.rds")

preds_19 <- dat$preds %>% 
  filter(year == 2019)
preds_nhc <- dat$preds %>%
  filter(site %in% c("NHC", "UNHC"))
preds <- dat$preds %>% filter(date > as.Date('2017-01-01'))
qq <- read_csv('data/rating_curves/interpolatedQ_allsites.csv', guess_max = 10000)
plot(qq$DateTime_UTC, qq$NHC.Q, log = 'y', type = 'l')
qq <- qq %>%
  mutate(date = as.Date(with_tz(DateTime_UTC, tz = "EST")))%>%
  group_by(date) %>%
  select(-DateTime_UTC, -notes, -PWC.Q) %>%
  summarize_all(mean, na.rm = T)

# plot annual met with discharge


sites$sscode <- c("8.5 km", "6.9 km", "5 km", "2.5 km", "2.3 km", "0 km")

axis_size = 0.9
Qlim = c(.02, max(qq$NHC.Q, na.rm = T) * 1e7)
ylim = c(-9,4)

# png("figures/metQ_across_sites_SM_2019.png",
#     width = 10, height = 6, units = 'in', res = 300)
tiff("figures/metQ_across_sites_SM_2019.tif",
    height = 5.25*800, width = 8.75*800, units = 'px', res = 800,
    compression = 'lzw')
  xlim = c(as.Date(c("2019-03-06", "2020-03-20")))
  par(ps = 10,
      oma = c(5,5,4,6),
      mfrow = c(3,2),
      mar = c(.5,0,0,.5))
  for(s in c(6,3,5,2,4,1)){
    ss <- sites$sitecode[s]
    sn <- sites$sscode[s]
    tmp <- preds_19 %>%
      filter(site == ss)
    tmpq <- qq[,c(1, s+1)] %>%
      filter(date >= xlim[1], 
             date <= xlim[2])
    colnames(tmpq) <- c('date', 'discharge')
    plot_metab(tmp, ylim = ylim, xlim = xlim, xaxt = 'n', yaxt = 'n')
    if(s %in% c(4,5,6)) {
      axis(2, cex.axis = axis_size, at = c(-6, -3, 0, 3))
    } 
    if(s %in% c(1,4)){
      axis(1, at = seq(as.Date('2019-03-01'), by = 'month', length.out = 13),
           labels = month.abb[c(3:12, 1:3)], cex.axis = axis_size, las = 2)
      # axis(1, at = seq(as.Date('2019-03-01'), by = '2 months', length.out = 7),
      #      labels = month.abb[c(3,5,7,9,11,1,3)])
    } 
    par(new = T)
    plot(tmpq$date, tmpq$discharge, log = "y", type = "l", ylim = Qlim, xlim = xlim,
         axes = FALSE, xlab = '', ylab = '')
    mtext(sn, cex = 0.9, line = -1.5, adj = 0.02)
    if(s %in% c(1,2,3)){
      axis(4, cex.axis = axis_size, at = c(0.01, 1, 100), 
           labels = c(0.01, 1, 100), las = 2)
    } 
  }
  par(new = T, mfrow = c(1,1), oma = c(4,4,0,4))
  plot(1,1, axes = F, ann = F, type = 'n')
  mtext("Metabolism gC/m2/d", 2, 2.8)
  mtext("Discharge (m3/s)", 4, 2.8)
  # mtext("Date", 1, 2.5)
  mtext("New Hope Creek metabolism across sites in 2019", 3, -2.5)
dev.off()

# png("figures/metQ_across_years_SM_2019.png",
#     width = 10, height = 6, units = 'in', res = 300)
tiff("figures/metQ_across_years_SM_2019.tif",
     height = 5.25*800, width = 8.75*800, units = 'px', res = 800,
     compression = 'lzw')
  par(ps = 10,
      oma = c(5,5,4,6),
      mfrow = c(3,2),
      mar = c(0.5,0,0,0.5))
  florence <- as.Date("2018-09-14")
  xlim = c(as.Date(c("2017-03-01", "2018-03-01")))
  for(y in c(2017,2018,2019)){
    tmp <- preds_nhc %>%
      filter(site == "UNHC",
             year == y)
    tmpq <- qq[,c(1, 7)] %>%
      filter(date >= xlim[1], 
             date <= xlim[2])
    colnames(tmpq) <- c('date', 'discharge')
    plot_metab(tmp, ylim = ylim, xlim = xlim, xaxt = 'n', yaxt = 'n')
    axis(2, cex.axis = axis_size, at = c(-6, -3, 0, 3))
    if(y == 2018){
      abline(v = florence, lty = 2)
      # text(x = florence, y = -8.5, labels = "Hurricane Florence", 
      #      cex = axis_size, pos = 4)
    }
    if(y == 2019){
      axis(1, at = seq(as.Date('2019-03-01'), by = 'month', length.out = 13),
           cex.axis = axis_size, labels = month.abb[c(3:12, 1:3)], las = 2)
      mtext("Upstream (0 km)",1, line = 4) 
      
    } 
    par(new = T)
    plot(tmpq$date, tmpq$discharge, log = "y", type = "l", ylim = Qlim, xlim = xlim,
         axes = FALSE, xlab = '', ylab = '')
    mtext(y, cex = 1, line = -1.5, adj = 0.02)

    tmp <- preds_nhc %>%
      filter(site == "NHC",
             year == y)
    tmpq <- qq[,c(1, 2)] %>%
      filter(date >= xlim[1], 
             date <= xlim[2])
    colnames(tmpq) <- c('date', 'discharge')
    plot_metab(tmp, ylim = ylim, xlim = xlim, xaxt = 'n', yaxt = 'n')
    if(y == 2018){
      abline(v = florence, lty = 2)
      # text(x = florence, y = -8.5, labels = "Hurricane Florence",
      #      cex = axis_size, pos = 4)
    }
    if(y == 2019){
      axis(1, at = seq(as.Date('2019-03-01'), by = 'month', length.out = 13),
           cex.axis = axis_size, labels = month.abb[c(3:12, 1:3)], las = 2)
      mtext("Downstream (8.5 km)",1, line =4) 
    } 
    par(new = T)
    plot(tmpq$date, tmpq$discharge, log = "y", type = "l", ylim = Qlim, xlim = xlim,
         axes = FALSE, xlab = '', ylab = '')
    mtext(y, cex = 0.9, line = -1.5, adj = 0.02)
    axis(4, cex.axis = axis_size, at = c(0.01, 1, 100), 
         labels = c(0.01, 1, 100), las = 2)
    xlim = xlim + 365
  }
  par(new = T, mfrow = c(1,1), oma = c(4,4,0,4))
  plot(1,1, axes = F, ann = F, type = 'n')
  mtext("Metabolism gC/m2/d", 2, 2.8)
  mtext("Discharge (m3/s)", 4, 2.8)
  # mtext("Date", 1, 2.5)
  mtext("New Hope Creek metabolism across years", 3, -2.5)
dev.off()

