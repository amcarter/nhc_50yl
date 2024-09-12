library(nhdplusTools)
library(raster)
library(tidyverse)
library(sf)
library(tmap)
library(rgdal)
library(maps)

tmap_options(check.and.fix = TRUE)

setwd('~/git/papers/alice_nhc')

sites <- read.csv("data/map_files/NHCsite_50yl_coordinates.csv",
                  header=TRUE, stringsAsFactors=FALSE)

sites_sf <- st_as_sf(sites, coords=c("longitude", "latitude"), remove=FALSE, crs=4326)
carter_sites <- filter(sites_sf,
                       type == 'now',
                       sitename %in% c('NewHopeCreek', 'ConcreteBridgePool')) %>%
    mutate(shortname = c('NHC', 'CB'))
hall_sites <- filter(sites_sf, type == 'hall')
stream_line <- st_read("data/map_files/stream_lines.shp")
# riparian_boundary <- st_read("data/other_watershed_stuff/riparian.shp")
study_reaches_line <- st_read("data/map_files/study_reaches.shp")
duke_forest_boundary <- st_read("data/map_files/2019_boundary.shp")
korstian_div <- filter(duke_forest_boundary, DIVISION == 'Korstian')
watershed_boundary <- st_read("data/watershed_boundary/nhc_wb_streamstats.shp")
# nlcd <- raster('data/map_files/NLCD_2019_Land_Cover_L48_20210604_PjC3eOp816qGLgHGHoRx.tiff')
# nlcd <- raster('data/nlcd/NLCD2016.tif')


#remove the piece of the streamlines that's downstream of the lowest site
stream_line_shortened <- stream_line
ll = as(stream_line_shortened$geometry[[1]], 'list')
ll[1][[1]] = ll[1][[1]][1:21, ]
stream_line_shortened$geometry[[1]] = st_multilinestring(ll)

study_reaches_line_shortened <- study_reaches_line
ll = as(study_reaches_line_shortened$geometry[[1]], 'list')
ll[1][[1]] = ll[1][[1]][1:21, ]
study_reaches_line_shortened$geometry[[1]] = st_multilinestring(ll)

#make new riparian boundary
PROJ4 = paste0('+proj=laea +lat_0=', round(mean(sites$latitude), 4), ' +lon_0=',
               round(mean(sites$longitude), 4),
             ' +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
riparian_boundary = study_reaches_line_shortened %>%
    st_transform(crs = PROJ4) %>%
    st_buffer(dist = 250) %>%
    st_transform(crs = 4326)

# tmap_mode("view")
tmap_mode("plot")
# par(bg=NA)

# map_with_colon = tm_shape(watershed_boundary) + tm_polygons(alpha=0, border.col="black", lwd=1) +
#     tm_shape(korstian_div) + tm_polygons(alpha=0.3, col = 'springgreen3',
#                                          border.col="transparent", lwd=.5) +
#     tm_shape(riparian_boundary) + tm_polygons(alpha=0, col="black", lwd=1.5,
#                                               border.col='steelblue3', border.alpha=0.8) +
#     tm_shape(study_reaches_line_shortened) + tm_lines(col='steelblue3', lwd=2.5) +
#     tm_shape(stream_line_shortened) + tm_lines(col='black', alpha=0.5, lwd=0.5) +
#     tm_shape(carter_sites) + tm_symbols(shape=1, col="red2", size=0.6, border.lwd=2) +
#     tm_shape(hall_sites) + tm_symbols(shape=3, col="black", size=0.6, border.lwd=2) +
#     tm_scale_bar(text.size = 1, position=c(0, 0)) +
#     tm_compass(type="arrow", position=c(0.875, 0.05), show.labels=1,
#                size=3, text.size=1) +
#     tm_style(style='white') +
#     tm_layout(frame=TRUE, bg.color="white") +
#     tm_add_legend(type='symbol', labels = '  Study sites', col = 'red2', size = 0.7,
#                   shape=1) +
#     tm_add_legend(type='symbol', labels = '  Hall 1972 sites', col = 'black',
#                   size=0.5, shape=3, border.lwd=2) +
#     tm_add_legend(type='line', labels = 'Study reach', col = 'steelblue3', lwd = 2.5) +
#     tm_add_legend(type='line', labels = 'Riparian zone', col = 'steelblue3', lwd = 1) +
#     tm_add_legend(type='fill', labels = 'Duke Forest', col = 'springgreen3', alpha=0.3,
#                   border.col='transparent') +
#     tm_legend(show=TRUE, position=c('left', 'top'), outside=FALSE, bg.color='gray97',
#               frame=TRUE, text.size=0.8)

qq = sf::st_bbox(watershed_boundary)

locs <- data.frame(x = qq[c(1,3)], y = qq[c(2,4)]) %>%
    st_as_sf(coords = c('x', 'y'), crs = 4326)

dem = elevatr::get_elev_raster(locs, z = 12)

## crop and mask
dem2 <- raster::crop(dem, watershed_boundary)
dem3 <- raster::mask(dem2, watershed_boundary)


qq = sf::st_bbox(watershed_boundary)
qq[1] = qq[1] - 0.01
# qq[4] = qq[4] + 0.005
map_without_colon_legend = tm_shape(watershed_boundary, bbox = qq) +
    tm_polygons(alpha=0, border.col="black", lwd=1) +
    tm_shape(dem3) +
    tm_raster(palette = terrain.colors(20), style = "cont",
              title = "Elevation (m)")+
    tm_shape(korstian_div) + tm_polygons(alpha=0.5, col = 'lightblue',
                                         border.col="grey50", lwd=.5) +
    # tm_shape(study_reaches_line_shortened) + tm_lines(col='steelblue3', lwd=2.5) +
    tm_shape(stream_line_shortened) + tm_lines(col='black', alpha=0.5, lwd=0.5) +
    tm_shape(carter_sites) + tm_symbols(shape=20, col="red2", size=0.1, border.lwd=1.6) +
    tm_text(text = 'shortname', ymod = 0.6, size = 0.8) +
    tm_scale_bar(text.size = 0.6, breaks = c(0, 1, 2, 3),
                 position=c(0.69, 0)) +
    tm_compass(type="arrow", position=c("RIGHT", "bottom", show.labels=3),
               size=1.6, text.size=0.8) +
    tm_style(style='white') +
    tm_layout(frame=TRUE, bg.color="white") +
    # tm_add_legend(type='symbol', labels = '  Study sites          ', col = 'red2', #size = 0.4,
    #               shape=20) +
    # tm_add_legend(type='line', labels = 'Study reach', col = 'steelblue3', lwd = 2.5) +
    tm_add_legend(type='fill', labels = 'Duke Forest', col = 'lightblue', alpha=0.5,
                  border.col='grey50') +
    tm_legend(show=TRUE, position=c('left', 'bottom'), outside=FALSE, bg.color='gray97',
              frame=TRUE, text.size=.7)+
    tm_layout(legend.title.size = 0.8,
              legend.text.size = 0.7)

# tmap_save(map_without_colon_legend, filename="figures/map_without_colon_legend.tiff",
tmap_save(map_without_colon_legend, filename="figures/map_without_colon_legend.png",
          bg="white", dpi = 800, height = 2, width = 4, units = 'in')

# map_without_colon = tm_shape(watershed_boundary, bbox = qq) +
#     tm_polygons(alpha=0, border.col="black", lwd=1) +
#     tm_shape(dem3) +
#     tm_raster(palette = terrain.colors(20), style = "cont",
#               title = "Elevation (m)", legend.show = FALSE)+
#     tm_shape(korstian_div) + tm_polygons(alpha=0.5, col = 'lightblue',
#                                          border.col="grey50", lwd=.5) +
#     # tm_shape(study_reaches_line_shortened) + tm_lines(col='steelblue3', lwd=2.5) +
#     tm_shape(stream_line_shortened) + tm_lines(col='black', alpha=0.5, lwd=0.5) +
#     tm_shape(carter_sites) + tm_symbols(shape=20, col="red2", size=0.1, border.lwd=1.6) +
#     tm_text(text = 'shortname', ymod = 0.6) +
#     tm_scale_bar(text.size = 0.6, breaks = c(0, 1, 2, 3),
#                  position=c(0.75, 0)) +
#     tm_compass(type="arrow", position=c("RIGHT", "bottom", show.labels=3),
#                size=1.6, text.size=0.8) +
#     tm_style(style='white') +
#     tm_layout(frame=TRUE, bg.color="white")

# tmap_save(map_without_colon, filename="figs/map_without_colon.tiff",
#           bg="white", dpi = 800, height = 2, width = 3.5, units = 'in')
