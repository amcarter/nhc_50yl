library(ggplot2)
library(lubridate)
library(ggthemes)
dat %>%
    filter(site == 'UNHC',
           between(date, as.Date('2019-03-01'), as.Date('2020-02-28'))) %>%
    # filter(PAR_surface == max(PAR_surface)) %>% select(site)
    ggplot(aes(x = date, y = PAR_surface)) +
    geom_line(color = 'goldenrod1', linewidth = 1.5) +
    theme_few()

ggsave('/tmp/light_plot_oneoff.png', width = 12, height = 8)
