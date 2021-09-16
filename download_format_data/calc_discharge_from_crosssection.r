# Calculate Discharge from cross sectional velocity measurements
# 1. fns to calc Q and plot cross section ####

plot_xc <- function(dat){
  par(mfrow = c(2,1), mar = c(0,4,4,1))
  if(is.na(dat$velocity_ms[1])){
    plot(1,1, axes = F, type = "n")
  } else {
  barplot(dat$velocity_ms, 
          width = c(diff(dat$distance_m), diff(dat$distance_m)[nrow(dat)-1]),
          main = paste(dat$site[1], dat$date[1]), 
          ylab = "velocity(m/s)")
  }
  par(mar = c(4,4,0,1))
  plot(dat$distance_m, -dat$depth_m, type = "l", 
       ylab = "depth (m)", xlab = "width (m)")
}

calc_xc_discharge <- function(dists, depths, vels){
  d_mids <- data.frame(dists = dists[1:(length(dists)-1)] +
                         diff(dists)/2)
  # insert midpoints so that trapezoids can be calculated
  # on either side of the velocity measurements
  dd <- data.frame(dists = dists,
                   depths = depths, 
                   v = vels) %>%
    full_join(d_mids) %>%
    arrange(dists) %>%
    mutate(depths = na.approx(depths),
           a = NA)
  # rep the vel measurements so they apply to both sides of 
  # each trapezoid
  for(i in 1:nrow(dd)){
    dd$v[i] <- ifelse(is.na(dd$v[i]), dd$v[i+1], dd$v[i])
  }

  for(i in 1:(length(dd$dists)-1)){
    dd$a[i] <- (diff(dd$dists)[i]) * 
    (abs(diff(dd$depths)[i])/2 + min(dd$depths[i:(i+1)]))
  }
  
  dd <- dd[1:(nrow(dd)-1),]    
  out <- data.frame(
    width = max(dists),
    xc_area = sum(dd$a),
    depth_avg = sum(dd$a)/max(dists),
    velocity_avg = sum(dd$a * dd$v)/sum(dd$a),
    discharge = sum(dd$a * dd$v), 
    depth_thalweg = max(dd$depths))
    
  
  return(out)
}

add_point_to_rc <- function(out){
  qq <- read_csv("data/rating_curves/calculated_discharge.csv")
  if(row.names(out) %in% row.names(qq)){ 
    next } else {
      stop(paste("this dataframe has columns that it shouldn't. 
                  Allowable columns include:", row.names(qq)))
    }
  qq <- bind_rows(qq, out)
  write_csv(qq, "data/rating_curves/calculated_discharge.csv")
}
