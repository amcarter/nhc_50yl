chili <- hindcast %>%
    mutate(GPP = zoo::rollmean(GPP, k = 7, na.pad = TRUE)) %>%
    mutate(ER = -zoo::rollmean(ER, k = 7, na.pad = TRUE)) %>%
    # group_by(doy) %>%
    # mutate(GPP = zoo::rollmean(GPP, k = 5, na.pad = TRUE)) %>%
    # mutate(ER = zoo::rollmean(ER, k = 5, na.pad = TRUE)) %>%
    pivot_longer(cols = c('GPP', 'ER'),
                 names_to = 'met',
                 values_to = 'gO2') %>%
    filter(year %% 3 == 1) %>%
    mutate(met = factor(met, levels = c('GPP', 'ER')))

met_change3 <- met_change2 %>%
    rename(`Mean Fall Discharge` = mean_Q, `Mean Annual Temperature` = temp.water) %>%
    pivot_longer(cols = c('Mean Fall Discharge', 'Mean Annual Temperature'),
                 names_to = 'variable', values_to = 'value')
met_dot <- filter(met_change3, year %% 3 == 1)

line_yrs = sort(unique(chili$year))
yrs = sort(unique(met_change3$year))
cols = viridis::plasma(11, begin = 0, end = 0.9)
ylim_gpp = range(chili$gO2[chili$met == 'GPP'])
ylim_er = range(chili$gO2[chili$met == 'ER'])
xlim_qt = range(met_change3$year)
ylim_q = range(met_change3$value[met_change3$variable == 'Mean Fall Discharge'])
ylim_t = range(met_change3$value[met_change3$variable != 'Mean Fall Discharge'])

gppi <- chili %>%
    filter(year <= line_yrs[1], met == 'GPP') %>%
    ggplot(aes(doy, gO2, col = factor(year), group = factor(year))) +
    ylim(ylim_gpp) +
    ylab(expression(paste('GPP (g ', O[2], m^-2, y^-1, ')'))) +
    xlab('Month') +
    scale_x_continuous(breaks = c(0,31,60,91,121,152,182,213,244,274,305,335),
                       labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          legend.position = 'none')

eri = chili %>%
    filter(year <= line_yrs[1], met == 'ER') %>%
    ggplot(aes(doy, gO2, col = factor(year), group = factor(year))) +
    ylim(ylim_er) +
    ylab(expression(paste('ER (g ', O[2], m^-2, y^-1, ')'))) +
    xlab('Month') +
    scale_x_continuous(breaks = c(0,31,60,91,121,152,182,213,244,274,305,335),
                       labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          legend.position = 'none')

met_plt = ggpubr::ggarrange(gppi, eri, nrow = 2)

cntr = 0
for(i in seq_along(yrs)){

    if(yrs[i] %in% line_yrs){

        cntr = cntr + 1

        gppi <- chili %>%
            filter(year <= line_yrs[cntr], met == 'GPP') %>%
            ggplot(aes(doy, gO2, col = factor(year), group = factor(year))) +
            # ggplot(aes(doy, gO2, col = year, group = factor(year))) +
            geom_line(linewidth = 0.7) +
            ylim(ylim_gpp) +
            # scale_color_viridis(name = 'Year', option = 'C', begin = 0, end = 0.9) +
            # facet_wrap(.~met, scales = 'free_y', ncol = 1, strip.position = 'right')+
            ylab(expression(paste('GPP (g ', O[2], m^-2, y^-1, ')'))) +
            xlab('Month') +
            scale_x_continuous(breaks = c(0,31,60,91,121,152,182,213,244,274,305,335),
                               labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) +
            scale_color_manual(values = cols[1:cntr]) +
            theme_bw() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  strip.text = element_blank(),
                  legend.position = 'none')
                  # legend.position.inside = c(.09, .18))

        eri = chili %>%
            filter(year <= line_yrs[cntr], met == 'ER') %>%
            ggplot(aes(doy, gO2, col = factor(year), group = factor(year))) +
            geom_line(linewidth = 0.7) +
            ylim(ylim_er) +
            ylab(expression(paste('ER (g ', O[2], m^-2, y^-1, ')'))) +
            xlab('Month') +
            scale_x_continuous(breaks = c(0,31,60,91,121,152,182,213,244,274,305,335),
                               labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) +
            scale_color_manual(values = cols[1:cntr]) +
            theme_bw() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  strip.text = element_blank(),
                  legend.position = 'none')

        met_plt = ggpubr::ggarrange(gppi, eri, nrow = 2)
    }

    qi <- met_change3 %>%
        filter(year <= yrs[i], variable == 'Mean Fall Discharge') %>%
        ggplot(aes(year, value, col = year)) +
        geom_line(linewidth = 1) +
        ylim(ylim_q) +
        xlim(xlim_qt) +
        geom_point(data = met_dot[met_dot$year <= yrs[i], ], size = 2) +
        geom_point(data = met_dot[met_dot$year <= yrs[i], ], col = 'black', pch = 1, size = 2) +
        # scale_color_viridis(name = 'Year', option = 'C', begin = 0, end = 0.9) +
        # facet_wrap(.~variable, scales = 'free_y', ncol = 1) +
        xlab('year') +
        ylab(expression(paste('Mean Annual Temp. (', degree, 'C)')))+
        scale_color_gradientn(colors = cols[1:cntr]) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background = element_blank(),
              strip.text = element_blank(),
              legend.position = 'none')

    ti <- met_change3 %>%
        filter(year <= yrs[i], variable != 'Mean Fall Discharge') %>%
        ggplot(aes(year, value, col = year)) +
        geom_line(linewidth = 1) +
        ylim(ylim_t) +
        xlim(xlim_qt) +
        geom_point(data = met_dot[met_dot$year <= yrs[i], ], size = 2) +
        geom_point(data = met_dot[met_dot$year <= yrs[i], ], col = 'black', pch = 1, size = 2) +
        xlab('year') +
        ylab(expression(paste('Mean Fall Discharge (', m^3, s^-1, ')')))+
        scale_color_gradientn(colors = cols[1:cntr]) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background = element_blank(),
              strip.text = element_blank(),
              legend.position = 'none')

    clim_plt = ggpubr::ggarrange(qi, ti, nrow = 2)

    ggpubr::ggarrange(met_plt, clim_plt, ncol = 2, widths = c(2,1))
    ggsave(paste0("figures/gifs/raw/", sprintf(i, fmt = '%02d'), ".png"), width = 9, height = 5, dpi = 100)
}

# system(paste('convert -loop 1 -dispose None -delay 200 figures/gifs/raw/1.png',
system(paste('convert -loop 0 -dispose None',
             '-delay 30 figures/gifs/raw/*.png',
             '-delay 150 figures/gifs/raw/35.png',
             'figures/gifs/metab.gif'))
             # '-delay 50 figures/gifs/blank.png figures/gifs/metab.gif'))

system('ffmpeg -i figures/gifs/metab.gif -movflags faststart -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" figures/gifs/metab.mp4')
