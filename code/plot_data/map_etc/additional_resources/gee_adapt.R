
setwd('~/git/papers/alice_nhc/')

site_codes <- c('NHC') #make sure these are in same order as the rows in the sf object
collection_name <- 'NHC'

# Load landcover defs
color_key <- read_csv('data/nlcd/pixel_color_key.csv',
                      show_col_types = FALSE,
                      col_names = FALSE) %>%
    rename(class_code = X1, color = X2, desc = X3)

nlcd_summary <- color_key %>%
    rename(id = class_code) %>%
    mutate(id = as.character(id))

boundaries <- st_read("data/watershed_boundary/nhc_wb_streamstats.shp") %>%
    mutate(site_code = site_codes) %>%
    select(site_code) %>%
    st_transform(crs = 4326) %>%
    st_make_valid()

# Upload watersheds to GEE
user_info <- rgee::ee_user_info(quiet = TRUE)
asset_folder <- glue('{a}/macrosheds_ws_boundaries/{d}',
                     a = user_info$asset_home,
                     d = collection_name)

gee_file_exist <- try(rgee::ee_manage_assetlist(asset_folder), silent = TRUE)
if(inherits(gee_file_exist, 'try-error')){
    sm(rgee::ee_manage_create(asset_folder))
}

for(i in 1:nrow(boundaries)){
    one_boundary <- boundaries[i, ]
    asset_path <- file.path(asset_folder, site_codes[i])
    sf_as_ee(one_boundary,
             via = 'getInfo_to_asset',
             assetId = asset_path,
             overwrite = TRUE,
             quiet = TRUE)
}

all_ee_task <- c()
needed_files <- c()
nlcd_all <- tibble()
nlcd_epochs <- c(2016)
for(s in 1:length(site_codes)){

    contents <- ee$data$listImages('USGS/NLCD_RELEASES')
    avail_releases <- sapply(contents$images, function(x) x$id)

    asset_path <- rgee::ee_manage_assetlist(asset_folder)

    # if(nrow(asset_path) == 1){
        ws_boundary_asset <- ee$FeatureCollection(asset_path$ID)
        filter_ee <- ee$Filter$inList('site_code', c(site_codes[s], site_codes[s]))
        site_ws_asset <- ws_boundary_asset$filter(filter_ee);
    # } else {
        # ws_boundary_asset <- str_split_fixed(asset_path$ID, '/', n = Inf)
        # ws_boundary_asset <- ws_boundary_asset[,ncol(ws_boundary_asset)]
        #
        # this_asset <- asset_path[grep(site_codes[s], ws_boundary_asset),]
        # site_ws_asset <- ee$FeatureCollection(this_asset$ID)
    # }

    for(e in nlcd_epochs){

        imgcol <- ifelse(as.numeric(str_extract(e, '[0-9]+')) <= 2016,
                         'USGS/NLCD_RELEASES/2016_REL',
                         glue('USGS/NLCD_RELEASES/{e}_REL/NLCD'))

        img <- ee$ImageCollection(imgcol)$
            select('landcover')$
            filter(ee$Filter$eq('system:index', e))$
            first()$
            # clip(ws_boundary_asset)
            clip(site_ws_asset)

        ee_description <- glue('{d}_{s}_{e}',
                               d = collection_name,
                               s = site_codes[s],
                               e = e)

        file_name <- paste0('nlcdX_X', e, 'X_X', site_codes[s])
        ee_task <- ee$batch$Export$image$toDrive(image = img,
                                                 description = ee_description,
                                                 folder = 'GEE',
                                                 # crs = 'ESPG:4326',
                                                 region = site_ws_asset$geometry(),
                                                 # region = ws_boundary_asset$geometry(),
                                                 fileNamePrefix = file_name,
                                                 maxPixels = 105921861)

        needed_files <- c(needed_files, file_name)
        all_ee_task <- c(all_ee_task, ee_description)

        suppressMessages(googledrive::drive_rm(paste0('GEE/', file_name, '.tif')))

        start_mess <- ee_task$start()
        # ee_monitoring(ee_task)
    }
}

needed_files <- paste0(needed_files, '.tif')

task_running <- rgee::ee_manage_task()

task_running <- task_running %>%
    filter(DestinationPath %in% all_ee_task)

while(any(task_running$State %in% c('RUNNING', 'READY'))){
    task_running <- rgee::ee_manage_task()

    task_running <- task_running %>%
        filter(DestinationPath %in% all_ee_task)

    Sys.sleep(10)
}

Sys.sleep(60)
task_running <- rgee::ee_manage_task()
if(any(task_running$State %in% c('RUNNING', 'READY'))){
    stop('weird. somehow the while loop above is exiting early')
}

temp_rgee <- tempfile(fileext = '.tif')

for(i in 1:length(needed_files)){

    rel_file <- needed_files[i]

    string <- str_match(rel_file, '(.+?)X_X(.+?)X_X(.+?)\\.tif')[2:4]
    year <- string[2]
    site <- string[3]

    file_there <- googledrive::drive_get(paste0('GEE/', rel_file))

    if(nrow(file_there) == 0) next

    expo_backoff(
        expr = {
            googledrive::drive_download(file = paste0('GEE/', rel_file),
                                        temp_rgee,
                                        overwrite = TRUE)
        },
        max_attempts = 5
    ) %>% invisible()

    nlcd_rast <- terra::rast(temp_rgee)
    nlcd_rast[as.vector(terra::values(nlcd_rast)) == 0] <- NA

    googledrive::drive_rm(rel_file)

    tabulated_values = terra::values(nlcd_rast) %>%
        table() %>%
        as_tibble() %>%
        rename(id = '.',
               CellTally = 'n')

    if(length(str_split_fixed(year, '_', n = Inf)[1, ]) == 2){
        year <- as.numeric(str_split_fixed(year, pattern = '_', n = Inf)[1, 1])
    }

    if(year == 1992){

        nlcd_e = full_join(nlcd_summary_1992,
                           tabulated_values,
                           by = 'id') %>%
            mutate(sum = sum(CellTally, na.rm = TRUE))

        nlcd_e_1992names <- nlcd_e %>%
            mutate(percent = round((CellTally / sum) * 100, 1)) %>%
            mutate(percent = ifelse(is.na(percent), 0, percent)) %>%
            select(var = macrosheds_1992_code, val = percent) %>%
            mutate(year = !!year)

        nlcd_e_norm_names <- nlcd_e %>%
            group_by(macrosheds_code) %>%
            summarize(CellTally1992 = sum(CellTally, na.rm = TRUE),
                      .groups = 'drop') %>%
            mutate(sum = sum(CellTally1992, na.rm = TRUE)) %>%
            mutate(percent = round((CellTally1992 / sum) * 100, 1)) %>%
            mutate(percent = ifelse(is.na(percent), 0, percent)) %>%
            select(var = macrosheds_code, val = percent) %>%
            mutate(year = !!year)

        nlcd_e <- rbind(nlcd_e_1992names, nlcd_e_norm_names) %>%
            mutate(site_code = !!site)

    } else {

        nlcd_e <- full_join(nlcd_summary,
                            tabulated_values,
                            by = 'id')

        nlcd_e <- nlcd_e %>%
            mutate(percent = round((CellTally * 100) / sum(CellTally, na.rm = TRUE), 1)) %>%
            mutate(percent = ifelse(is.na(percent), 0, percent)) %>%
            select(var = macrosheds_code, val = percent) %>%
            mutate(year = !!year) %>%
            mutate(site_code = !!site)
    }

    nlcd_all <- rbind(nlcd_all, nlcd_e)
}

nlcd_final <- nlcd_all %>%
    mutate(year = as.numeric(year)) %>%
    select(year, site_code, var, val)

nlcd_final <- append_unprod_prefix(nlcd_final, prodname_ms)

nlcd_final <- bind_older_ws_traits(nlcd_final)

save_general_files(final_file = nlcd_final,
                   domain_dir = nlcd_dir)
