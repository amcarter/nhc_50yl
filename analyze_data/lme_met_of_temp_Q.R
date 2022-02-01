# Build a basic linear mixed effects model of met as f(temp, Q)
library(lme4)
dat <- readRDS("data/metabolism/compiled/met_preds_stream_metabolizer.rds")

met <- dat$preds 


# subset out modern dataset to build model
mm <- met %>%
  filter(era == 'now') %>%
  mutate(ER = -ER,
         logQ = log(discharge),
         year = factor(year))
ggplot(mm, aes(date, ER, col = year)) +
  geom_point()

ggplot(mm, aes(temp.water, ER, col = factor(month)))+
  geom_point() + geom_smooth(method = lm) +
  facet_wrap(~site)

# subset to fall respiration data and add in annual flow values

yy <- mm %>%
  group_by(year) %>%
  summarize(logQ_mean = mean(logQ, na.rm = T)) %>%
  ungroup() 

fall <- mm %>%
  filter(month %in% c(10,11)) %>%
  select(date, GPP, ER, logQ, temp.water, site, year) %>%
  left_join(yy, by = 'year')

ggplot(fall, aes(date, ER, col = factor(year)))+
  geom_point() + geom_smooth(method = lm) +
  facet_wrap(year~site, scales = 'free_x')
# scale data to model:

sfall <- fall %>% 
  mutate(across(.cols = all_of(c('GPP', 'ER', 'temp.water')),
                .fns = scale))
# model for ER ####

mER = lmer(ER ~ temp.water + (temp.water|logQ_mean) + (temp.water|site), 
           data = sfall)
summary(mER)$coeff
confint(mER)

ranef(mER)
ER_mod <- predict(mER) 
tmp <- data.frame(index = as.numeric(names(ER_mod)), ER_mod = ER_mod)
sfall <- sfall %>%
  mutate(index = seq(1:nrow(sfall))) %>%
  left_join(tmp) %>%
  select(-index)

ggplot(sfall, aes(ER, ER_mod, col = site))+
  geom_point() + geom_smooth(method = lm) +
  geom_abline(slope = 1, intercept = 0)

# It seems like there is actually no relationship between temperature and ER in  the 
#  fall at the three run sites. The model is unable to predict at these sites (unsurprisingly)

# My next step would be to incorporate depth or slope (whichever best differentiates between
# the two sets of sites). I think this type of model will make this a much better paper and is
# worth doing. 

# another possibility is to try something like a random forest or other machine
# learning approach, I would just have to see how much of an investment something
# like this is in order to do it from a place of knowing what's up.

# finally, maybe I do want to do a bayesian hierarchical model. It might actually work and
# be cool, it just means that I am still working on my thesis for ages.