## a) Data Manipulation and Descriptive Stats Packages
# load tidyverse packages
library('tidyverse')

# load here packages
library('here')

## b) Graph Plotting and Visualization Packages
# load ggplot2 packages
library ('ggplot2')

# load gridExtra packages
library('gridExtra')

## c) Modeling Packages
# load lme4 packages
library ('lme4')

# load broom packages
library('broom')
library('broom.mixed')


image.png

### 5.3 Graph of estimates of associations between maternal prenatal social
# stressors and Pan Tissue EAA in mid-childhood

## a) Min temp before model: extract estimates and tidy the data frame
mass_min_before_thermo_adj_lmer_est <- tidy(mass_min_before_thermo_adj_lmer,
                               conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(thermo_bef_min_temp)')

## b) Rename variables for estimates of interest
mass_min_before_thermo_adj_lmer_est['term'][mass_min_before_thermo_adj_lmer_est['term'] ==
                                 'scale(thermo_bef_min_temp)'] <- 'Before day 6'


## c) Min temp before model: Label the estimates in data frame
mass_min_before_thermo_adj_lmer_est$model <- c('Before day 6 model')

## d) Min temp after model: extract estimates and tidy the data frame
mass_min_after_thermo_adj_lmer_est <- tidy(mass_min_after_thermo_adj_lmer,
                                          conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(thermo_aft_min_temp)')


## e) Rename variables for estimates of interest
mass_min_after_thermo_adj_lmer_est['term'][mass_min_after_thermo_adj_lmer_est['term'] ==
                                              'scale(thermo_aft_min_temp)'] <- 'After day 6'


## f) Min temp before model: Label the estimates in data frame
mass_min_after_thermo_adj_lmer_est$model <- c('After day 6 model')


## d) Min temp both model: extract estimates and tidy the data frame
mass_min_both_thermo_adj_lmer_est <- tidy(mass_min_both_thermo_adj_lmer,
                                           conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(thermo_bef_min_temp)' | term == 'scale(thermo_aft_min_temp)')


## e) Rename variables for estimates of interest
mass_min_both_thermo_adj_lmer_est['term'][mass_min_both_thermo_adj_lmer_est['term'] ==
                                              'scale(thermo_bef_min_temp)'] <- 'Before day 6'


mass_min_both_thermo_adj_lmer_est['term'][mass_min_both_thermo_adj_lmer_est['term'] ==
                                             'scale(thermo_aft_min_temp)'] <- 'After day 6'


## f) Min temp before model: Label the estimates in data frame
mass_min_both_thermo_adj_lmer_est$model <- c('Both time periods model')

## j) Combine regression estimates into a tidy table
mass_min_both_thermo_adj_lmer_tbl <- bind_rows(mass_min_before_thermo_adj_lmer_est,
                                               mass_min_after_thermo_adj_lmer_est,
                                               mass_min_both_thermo_adj_lmer_est)

## k) Re-code *nominal* factor (with ordered levels)  
# Set levels (odering) of 'model' variable
mass_min_both_thermo_adj_lmer_tbl <-
  transform(mass_min_both_thermo_adj_lmer_tbl,
            model = factor(model,
                           levels = c('Before day 6 model',
                                      'After day 6 model',
                                      'Both time periods model')))

mass_min_both_thermo_adj_lmer_tbl <-
  transform(mass_min_both_thermo_adj_lmer_tbl,
            term = factor(term,
                          levels = c('Before day 6', 
                                     'After day 6')))

## l) Graph of estimates of associations between maternal prenatal social
# stressors and Pan Tissue EAA in mid-childhood
mass_min_both_thermo_adj_lmer_plot <-
  ggplot(mass_min_both_thermo_adj_lmer_tbl, aes(x = term, y = estimate,
                                     color = term)) +
  theme_classic() + 
  geom_hline(yintercept = 0, color = 'red',
             linetype = 2) + # line at null behind coefs
  geom_point(size = 6,
             position=position_dodge(width = 0.0)) +
  geom_errorbar(aes(ymin=(conf.low),
                    ymax=(conf.high)), width=.1,
                position=position_dodge(0.0), size = 1) +
  scale_color_manual(values=c('#440154', '#7ad151','#440154','#7ad151')) +
  facet_grid(cols = vars(model), scales = 'free', space = 'free') +
  theme(strip.background =element_rect(fill= 'white'))+
  theme(strip.text = element_text(colour = 'black')) +
  #coord_flip() + # flip x and y axes
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18)) +
  theme(text = element_text(size=18)) +
  # bold and size title and axes labels
  theme(legend.position = 'none') +
  theme(axis.ticks = element_blank()) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(size = 16, angle=45,
                                 margin = margin(t = 50, r = 0,
                                                 b = 0, l = 0), face = 'bold'), 
        axis.text.y = element_text(size = 16, angle=0,
                                   margin = margin(t = 25, r = 0,
                                                   b = 0, l = 10)),
        axis.title = element_text(size = 18), legend.title=element_blank(),
        legend.text=element_text(size=18),
        legend.position = 'none', #c(0.91, 0.94),
        legend.key = element_blank()) +
  labs(colour = "Parent feeding level") +
  xlab(expression(italic('Standardized minimum temperature (SD)'))) +
  ylab(expression
       (atop(bold('Beta estimate and 95% CI'),
             paste(italic('Nestling mass (g)')))))
  # remove axis ticks

print(mass_min_both_thermo_adj_lmer_plot)

## m) Save Plot
# use ggsave to save the plot
ggsave('mass_min_both_thermo_adj_lmer_plot.png', plot = mass_min_both_thermo_adj_lmer_plot, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 9, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 



## a) Max temp before model: extract estimates and tidy the data frame
mass_max_before_thermo_adj_lmer_est <- tidy(mass_max_before_thermo_adj_lmer,
                                            conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(thermo_bef_max_temp)')

## b) Rename variables for estimates of interest
mass_max_before_thermo_adj_lmer_est['term'][mass_max_before_thermo_adj_lmer_est['term'] ==
                                              'scale(thermo_bef_max_temp)'] <- 'Before day 6'


## c) Min temp before model: Label the estimates in data frame
mass_max_before_thermo_adj_lmer_est$model <- c('Before day 6 model')

## d) Min temp after model: extract estimates and tidy the data frame
mass_max_after_thermo_adj_lmer_est <- tidy(mass_max_after_thermo_adj_lmer,
                                           conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(thermo_aft_max_temp)')


## e) Rename variables for estimates of interest
mass_max_after_thermo_adj_lmer_est['term'][mass_max_after_thermo_adj_lmer_est['term'] ==
                                             'scale(thermo_aft_max_temp)'] <- 'After day 6'


## f) Min temp before model: Label the estimates in data frame
mass_max_after_thermo_adj_lmer_est$model <- c('After day 6 model')


## d) Min temp both model: extract estimates and tidy the data frame
mass_max_both_thermo_adj_lmer_est <- tidy(mass_max_both_thermo_adj_lmer,
                                          conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(thermo_bef_max_temp)' | term == 'scale(thermo_aft_max_temp)')


## e) Rename variables for estimates of interest
mass_max_both_thermo_adj_lmer_est['term'][mass_max_both_thermo_adj_lmer_est['term'] ==
                                            'scale(thermo_bef_max_temp)'] <- 'Before day 6'


mass_max_both_thermo_adj_lmer_est['term'][mass_max_both_thermo_adj_lmer_est['term'] ==
                                            'scale(thermo_aft_max_temp)'] <- 'After day 6'


## f) Min temp before model: Label the estimates in data frame
mass_max_both_thermo_adj_lmer_est$model <- c('Both time periods model')

## j) Combine regression estimates into a tidy table
mass_max_both_thermo_adj_lmer_tbl <- bind_rows(mass_max_before_thermo_adj_lmer_est,
                                               mass_max_after_thermo_adj_lmer_est,
                                               mass_max_both_thermo_adj_lmer_est)

## k) Re-code *nominal* factor (with ordered levels)  
# Set levels (odering) of 'model' variable
mass_max_both_thermo_adj_lmer_tbl <-
  transform(mass_max_both_thermo_adj_lmer_tbl,
            model = factor(model,
                           levels = c('Before day 6 model',
                                      'After day 6 model',
                                      'Both time periods model')))

mass_max_both_thermo_adj_lmer_tbl <-
  transform(mass_max_both_thermo_adj_lmer_tbl,
            term = factor(term,
                          levels = c('Before day 6', 
                                     'After day 6')))

## l) Graph of estimates of associations between maternal prenatal social
# stressors and Pan Tissue EAA in mid-childhood
mass_max_both_thermo_adj_lmer_plot <-
  ggplot(mass_max_both_thermo_adj_lmer_tbl, aes(x = term, y = estimate,
                                                color = term)) +
  theme_classic() + 
  geom_hline(yintercept = 0, color = 'red',
             linetype = 2) + # line at null behind coefs
  geom_point(size = 6,
             position=position_dodge(width = 0.0)) +
  geom_errorbar(aes(ymin=(conf.low),
                    ymax=(conf.high)), width=.1,
                position=position_dodge(0.0), size = 1) +
  scale_color_manual(values=c('#440154', '#7ad151','#440154','#7ad151')) +
  facet_grid(cols = vars(model), scales = 'free', space = 'free') +
  theme(strip.background =element_rect(fill= 'white'))+
  theme(strip.text = element_text(colour = 'black')) +
  #coord_flip() + # flip x and y axes
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18)) +
  theme(text = element_text(size=18)) +
  # bold and size title and axes labels
  theme(legend.position = 'none') +
  theme(axis.ticks = element_blank()) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(size = 16, angle=45,
                                   margin = margin(t = 50, r = 0,
                                                   b = 0, l = 0), face = 'bold'), 
        axis.text.y = element_text(size = 16, angle=0,
                                   margin = margin(t = 25, r = 0,
                                                   b = 0, l = 10)),
        axis.title = element_text(size = 18), legend.title=element_blank(),
        legend.text=element_text(size=18),
        legend.position = 'none', #c(0.91, 0.94),
        legend.key = element_blank()) +
  labs(colour = "Parent feeding level") +
  xlab(expression(italic('Standardized maximum temperature (SD)'))) +
  ylab(expression
       (atop(bold('Beta estimate and 95% CI'),
             paste(italic('Nestling mass (g)')))))


print(mass_max_both_thermo_adj_lmer_plot)

## m) Save Plot
# use ggsave to save the plot
ggsave('mass_max_both_thermo_adj_lmer_plot.png', plot = mass_max_both_thermo_adj_lmer_plot, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 9, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 


## a) Min temp before model: extract estimates and tidy the data frame
mass_iqr_before_thermo_adj_lmer_est <- tidy(mass_iqr_before_thermo_adj_lmer,
                                            conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(thermo_bef_iqr_temp)')

## b) Rename variables for estimates of interest
mass_iqr_before_thermo_adj_lmer_est['term'][mass_iqr_before_thermo_adj_lmer_est['term'] ==
                                              'scale(thermo_bef_iqr_temp)'] <- 'Before day 6'


## c) Min temp before model: Label the estimates in data frame
mass_iqr_before_thermo_adj_lmer_est$model <- c('Before day 6 model')

## d) Min temp after model: extract estimates and tidy the data frame
mass_iqr_after_thermo_adj_lmer_est <- tidy(mass_iqr_after_thermo_adj_lmer,
                                           conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(thermo_aft_iqr_temp)')


## e) Rename variables for estimates of interest
mass_iqr_after_thermo_adj_lmer_est['term'][mass_iqr_after_thermo_adj_lmer_est['term'] ==
                                             'scale(thermo_aft_iqr_temp)'] <- 'After day 6'


## f) Min temp before model: Label the estimates in data frame
mass_iqr_after_thermo_adj_lmer_est$model <- c('After day 6 model')


## d) Min temp both model: extract estimates and tidy the data frame
mass_iqr_both_thermo_adj_lmer_est <- tidy(mass_iqr_both_thermo_adj_lmer,
                                          conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(thermo_bef_iqr_temp)' | term == 'scale(thermo_aft_iqr_temp)')


## e) Rename variables for estimates of interest
mass_iqr_both_thermo_adj_lmer_est['term'][mass_iqr_both_thermo_adj_lmer_est['term'] ==
                                            'scale(thermo_bef_iqr_temp)'] <- 'Before day 6'


mass_iqr_both_thermo_adj_lmer_est['term'][mass_iqr_both_thermo_adj_lmer_est['term'] ==
                                            'scale(thermo_aft_iqr_temp)'] <- 'After day 6'


## f) Min temp before model: Label the estimates in data frame
mass_iqr_both_thermo_adj_lmer_est$model <- c('Both time periods model')

## j) Combine regression estimates into a tidy table
mass_iqr_both_thermo_adj_lmer_tbl <- bind_rows(mass_iqr_before_thermo_adj_lmer_est,
                                               mass_iqr_after_thermo_adj_lmer_est,
                                               mass_iqr_both_thermo_adj_lmer_est)

## k) Re-code *nominal* factor (with ordered levels)  
# Set levels (odering) of 'model' variable
mass_iqr_both_thermo_adj_lmer_tbl <-
  transform(mass_iqr_both_thermo_adj_lmer_tbl,
            model = factor(model,
                           levels = c('Before day 6 model',
                                      'After day 6 model',
                                      'Both time periods model')))

mass_iqr_both_thermo_adj_lmer_tbl <-
  transform(mass_iqr_both_thermo_adj_lmer_tbl,
            term = factor(term,
                          levels = c('Before day 6', 
                                     'After day 6')))

## l) Graph of estimates of associations between maternal prenatal social
# stressors and Pan Tissue EAA in mid-childhood
mass_iqr_both_thermo_adj_lmer_plot <-
  ggplot(mass_iqr_both_thermo_adj_lmer_tbl, aes(x = term, y = estimate,
                                                color = term)) +
  theme_classic() + 
  geom_hline(yintercept = 0, color = 'red',
             linetype = 2) + # line at null behind coefs
  geom_point(size = 6,
             position=position_dodge(width = 0.0)) +
  geom_errorbar(aes(ymin=(conf.low),
                    ymax=(conf.high)), width=.1,
                position=position_dodge(0.0), size = 1) +
  scale_color_manual(values=c('#440154', '#7ad151','#440154','#7ad151')) +
  facet_grid(cols = vars(model), scales = 'free', space = 'free') +
  theme(strip.background =element_rect(fill= 'white'))+
  theme(strip.text = element_text(colour = 'black')) +
  #coord_flip() + # flip x and y axes
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18)) +
  theme(text = element_text(size=18)) +
  # bold and size title and axes labels
  theme(legend.position = 'none') +
  theme(axis.ticks = element_blank()) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(size = 16, angle=45,
                                   margin = margin(t = 50, r = 0,
                                                   b = 0, l = 0), face = 'bold'), 
        axis.text.y = element_text(size = 16, angle=0,
                                   margin = margin(t = 25, r = 0,
                                                   b = 0, l = 10)),
        axis.title = element_text(size = 18), legend.title=element_blank(),
        legend.text=element_text(size=18),
        legend.position = 'none', #c(0.91, 0.94),
        legend.key = element_blank()) +
  labs(colour = "Parent feeding level") +
  xlab(expression(italic('Standardized IQR of temperature (SD)'))) +
  ylab(expression
       (atop(bold('Beta estimate and 95% CI'),
             paste(italic('Nestling mass (g)')))))
# remove axis ticks

print(mass_iqr_both_thermo_adj_lmer_plot)

## m) Save Plot
# use ggsave to save the plot
ggsave('mass_iqr_both_thermo_adj_lmer_plot.png', plot = mass_iqr_both_thermo_adj_lmer_plot, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 9, 
       height = 6, 
       units = c('in'), dpi = 300, limitsize = TRUE) 




#### Parent feeding level DW plots
## a) Min temp before model: extract estimates and tidy the data frame
mass_min_temp_blups_low_adj_lmer_est <- tidy(mass_min_temp_blups_low_adj_lmer,
                                            conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(nest_min_temp)')

mass_min_temp_blups_med_adj_lmer_est <- tidy(mass_min_temp_blups_med_adj_lmer,
                                             conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(nest_min_temp)')

mass_min_temp_blups_high_adj_lmer_est <- tidy(mass_min_temp_blups_high_adj_lmer,
                                             conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(nest_min_temp)')

mass_min_temp_blups_adj_lmer_est <- tidy(mass_min_temp_blups_adj_lmer,
                                         conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(nest_min_temp):scale(feeding_expontd_blups)')


## b) Rename variables for estimates of interest
mass_min_temp_blups_low_adj_lmer_est['term'][mass_min_temp_blups_low_adj_lmer_est['term'] ==
                                              'scale(nest_min_temp)'] <- 'Low feeding'

mass_min_temp_blups_med_adj_lmer_est['term'][mass_min_temp_blups_med_adj_lmer_est['term'] ==
                                               'scale(nest_min_temp)'] <- 'Avg. feeding'

mass_min_temp_blups_high_adj_lmer_est['term'][mass_min_temp_blups_high_adj_lmer_est['term'] ==
                                               'scale(nest_min_temp)'] <- 'High feeding'

mass_min_temp_blups_adj_lmer_est['term'][mass_min_temp_blups_adj_lmer_est['term'] ==
                                           'scale(nest_min_temp):scale(feeding_expontd_blups)'] <- 'Minimum temp. * feeding level'



## c) Min temp before model: Label the estimates in data frame
mass_min_temp_blups_low_adj_lmer_est$model <- c('Minimum temperature, low feeding model')

mass_min_temp_blups_med_adj_lmer_est$model <- c('Minimum temperature, avg. feeding model')

mass_min_temp_blups_high_adj_lmer_est$model <- c('Minimum temperature, high feeding model')

mass_min_temp_blups_adj_lmer_est$model <- c('Minimum temperature interaction model')



## j) Combine regression estimates into a tidy table
mass_min_feeding_adj_lmer_tbl <- bind_rows(mass_min_temp_blups_low_adj_lmer_est,
                                            mass_min_temp_blups_med_adj_lmer_est,
                                            mass_min_temp_blups_high_adj_lmer_est)

## k) Re-code *nominal* factor (with ordered levels)  
# Set levels (odering) of 'model' variable
mass_min_feeding_adj_lmer_tbl <-
  transform(mass_min_feeding_adj_lmer_tbl,
            model = factor(model,
                           levels = c(
                                      'Low feeding model',
                                      'Avg. feeding model',
                                      'High feeding model')))

mass_min_feeding_adj_lmer_tbl <-
  transform(mass_min_feeding_adj_lmer_tbl,
           term = factor(term,
                         levels = c('Low feeding',
                                    'Avg. feeding',
                                    'High feeding')))

## l) Graph of estimates of associations between maternal prenatal social
# stressors and Pan Tissue EAA in mid-childhood
mass_min_feeding_adj_lmer_tbl_plot <-
  ggplot(mass_min_feeding_adj_lmer_tbl, aes(x = term, y = estimate,
                                                color = term)) +
  theme_classic() + 
  geom_hline(yintercept = 0, color = 'red',
             linetype = 2) + # line at null behind coefs
  geom_point(size = 6,
             position=position_dodge(width = 0.0)) +
  geom_errorbar(aes(ymin=(conf.low),
                    ymax=(conf.high)), width=.1,
                position=position_dodge(0.0), size = 1) +
  scale_color_manual(values=c("#440154", "#31688e", "#35b779")) +
  theme(strip.background =element_rect(fill= 'white'))+
  theme(strip.text = element_text(colour = 'black')) +
  #coord_flip() + # flip x and y axes
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18)) +
  theme(text = element_text(size=18)) +
  # bold and size title and axes labels
  theme(legend.position = 'none') +
  theme(axis.ticks = element_blank()) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(size = 16, angle=45,
                                   margin = margin(t = 50, r = 0,
                                                   b = 0, l = 0), face = 'bold'), 
        axis.text.y = element_text(size = 16, angle=0,
                                   margin = margin(t = 25, r = 0,
                                                   b = 0, l = 10)),
        axis.title = element_text(size = 18), legend.title=element_blank(),
        legend.text=element_text(size=18),
        legend.position = 'none', #c(0.91, 0.94),
        legend.key = element_blank()) +
  labs(colour = "Model term") +
  xlab(expression(italic('Standardized minimum temperature (SD)'))) +
  ylab(expression
       (atop(bold('Beta estimate and 95% CI'),
             paste(italic('Nestling mass (g)')))))
# remove axis ticks

print(mass_min_feeding_adj_lmer_tbl_plot)

## m) Save Plot
# use ggsave to save the plot
ggsave('mass_min_feeding_adj_lmer_tbl_plot.png', plot = mass_min_feeding_adj_lmer_tbl_plot, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 10, 
       height = 10, 
       units = c('in'), dpi = 300, limitsize = TRUE) 



## a) Max temp model: extract estimates and tidy the data frame

mass_max_temp_blups_adj_lmer_est <- tidy(mass_max_temp_blups_adj_lmer,
                                         conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(nest_max_temp)' |
           term == 'scale(nest_max_temp):scale(feeding_expontd_blups)')


## b) Rename variables for estimates of interest
mass_max_temp_blups_adj_lmer_est['term'][mass_max_temp_blups_adj_lmer_est['term'] ==
                                               'scale(nest_max_temp)'] <- 'Maximum temp.'

mass_max_temp_blups_adj_lmer_est['term'][mass_max_temp_blups_adj_lmer_est['term'] ==
                                           'scale(nest_max_temp):scale(feeding_expontd_blups)'] <- 'Maximum temp. * feeding level'





## j) Combine regression estimates into a tidy table
mass_max_feeding_adj_lmer_tbl <- mass_max_temp_blups_adj_lmer_est

## k) Re-code *nominal* factor (with ordered levels)  
# Set levels (odering) of 'model' variable

mass_max_feeding_adj_lmer_tbl <-
  transform(mass_max_feeding_adj_lmer_tbl,
            term = factor(term,
                          levels = c('Maximum temp. * feeding level',
                                     'Maximum temp.')))

## l) Graph of estimates of associations between maternal prenatal social
# stressors and Pan Tissue EAA in mid-childhood
mass_max_feeding_adj_lmer_tbl <-
  ggplot(mass_max_feeding_adj_lmer_tbl, aes(x = term, y = estimate,
                                            color = term)) +
  theme_classic() + 
  geom_hline(yintercept = 0, color = 'red',
             linetype = 2) + # line at null behind coefs
  geom_point(size = 6,
             position=position_dodge(width = 0.0)) +
  geom_errorbar(aes(ymin=(conf.low),
                    ymax=(conf.high)), width=.1,
                position=position_dodge(0.0), size = 1) +
  scale_color_manual(values=c("#440154", "#35b779")) +
  theme(strip.background =element_rect(fill= 'white'))+
  theme(strip.text = element_text(colour = 'black')) +
  #coord_flip() + # flip x and y axes
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18)) +
  theme(text = element_text(size=18)) +
  # bold and size title and axes labels
  theme(legend.position = 'none') +
  theme(axis.ticks = element_blank()) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(size = 16, angle=45,
                                   margin = margin(t = 50, r = 0,
                                                   b = 0, l = 0), face = 'bold'), 
        axis.text.y = element_text(size = 16, angle=0,
                                   margin = margin(t = 25, r = 0,
                                                   b = 0, l = 10)),
        axis.title = element_text(size = 18), legend.title=element_blank(),
        legend.text=element_text(size=18),
        legend.position = 'none', #c(0.91, 0.94),
        legend.key = element_blank()) +
  labs(colour = "Model term") +
  xlab(expression(italic('Model term'))) +
  ylab(expression
       (atop(bold('Beta estimate and 95% CI'),
             paste(italic('Nestling mass (g)')))))
# remove axis ticks

print(mass_max_feeding_adj_lmer_tbl)

## m) Save Plot
# use ggsave to save the plot
ggsave('mass_max_feeding_adj_lmer_tbl.png', plot = mass_max_feeding_adj_lmer_tbl, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 10, 
       height = 10, 
       units = c('in'), dpi = 300, limitsize = TRUE) 






## a) Max temp model: extract estimates and tidy the data frame

mass_iqr_temp_blups_adj_lmer_est <- tidy(mass_iqr_temp_blups_adj_lmer,
                                         conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(nest_iqr_temp)' |
           term == 'scale(nest_iqr_temp):scale(feeding_expontd_blups)')


## b) Rename variables for estimates of interest
mass_iqr_temp_blups_adj_lmer_est['term'][mass_iqr_temp_blups_adj_lmer_est['term'] ==
                                           'scale(nest_iqr_temp)'] <- 'IQR of temp.'

mass_iqr_temp_blups_adj_lmer_est['term'][mass_iqr_temp_blups_adj_lmer_est['term'] ==
                                           'scale(nest_iqr_temp):scale(feeding_expontd_blups)'] <- 'IQR of temp. * feeding level'


## j) Combine regression estimates into a tidy table
mass_iqr_feeding_adj_lmer_tbl <- mass_iqr_temp_blups_adj_lmer_est

## k) Re-code *nominal* factor (with ordered levels)  
# Set levels (odering) of 'model' variable

mass_iqr_feeding_adj_lmer_tbl <-
  transform(mass_iqr_feeding_adj_lmer_tbl,
            term = factor(term,
                          levels = c('IQR of temp. * feeding level',
                                     'IQR of temp.')))

## l) Graph of estimates of associations between maternal prenatal social
# stressors and Pan Tissue EAA in mid-childhood
mass_iqr_feeding_adj_lmer_tbl_plot <-
  ggplot(mass_iqr_feeding_adj_lmer_tbl, aes(x = term, y = estimate,
                                            color = term)) +
  theme_classic() + 
  geom_hline(yintercept = 0, color = 'red',
             linetype = 2) + # line at null behind coefs
  geom_point(size = 6,
             position=position_dodge(width = 0.0)) +
  geom_errorbar(aes(ymin=(conf.low),
                    ymax=(conf.high)), width=.1,
                position=position_dodge(0.0), size = 1) +
  scale_color_manual(values=c("#440154", "#35b779")) +
  theme(strip.background =element_rect(fill= 'white'))+
  theme(strip.text = element_text(colour = 'black')) +
  #coord_flip() + # flip x and y axes
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18)) +
  theme(text = element_text(size=18)) +
  # bold and size title and axes labels
  theme(legend.position = 'none') +
  theme(axis.ticks = element_blank()) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(size = 16, angle=45,
                                   margin = margin(t = 50, r = 0,
                                                   b = 0, l = 0), face = 'bold'), 
        axis.text.y = element_text(size = 16, angle=0,
                                   margin = margin(t = 25, r = 0,
                                                   b = 0, l = 10)),
        axis.title = element_text(size = 18), legend.title=element_blank(),
        legend.text=element_text(size=18),
        legend.position = 'none', #c(0.91, 0.94),
        legend.key = element_blank()) +
  labs(colour = "Model term") +
  xlab(expression(italic('Model term'))) +
  ylab(expression
       (atop(bold('Beta estimate and 95% CI'),
             paste(italic('Nestling mass (g)')))))
# remove axis ticks

print(mass_iqr_feeding_adj_lmer_tbl_plot)

## m) Save Plot
# use ggsave to save the plot
ggsave('mass_iqr_feeding_adj_lmer_tbl_plot.png', plot = mass_iqr_feeding_adj_lmer_tbl_plot, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 10, 
       height = 10, 
       units = c('in'), dpi = 300, limitsize = TRUE) 




## a) Min temp brooding model: extract estimates and tidy the data frame
mass_min_temp_brooding_blups_adj_lmer_est <- tidy(mass_min_temp_brooding_blups_adj_lmer,
                                         conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(nest_min_temp)' |
           term == 'scale(nest_min_temp):scale(brooding_blups)')


## b) Rename variables for estimates of interest
mass_min_temp_brooding_blups_adj_lmer_est['term'][mass_min_temp_brooding_blups_adj_lmer_est['term'] ==
                                           'scale(nest_min_temp)'] <- 'Minimum temp.'

mass_min_temp_brooding_blups_adj_lmer_est['term'][mass_min_temp_brooding_blups_adj_lmer_est['term'] ==
                                           'scale(nest_min_temp):scale(brooding_blups)'] <- 'Minimum temp. * brooding level'


## j) Combine regression estimates into a tidy table
mass_min_temp_brooding_blups_adj_lmer_tbl <- mass_min_temp_brooding_blups_adj_lmer_est

## k) Re-code *nominal* factor (with ordered levels)  
# Set levels (odering) of 'model' variable

mass_min_temp_brooding_blups_adj_lmer_tbl <-
  transform(mass_min_temp_brooding_blups_adj_lmer_tbl,
            term = factor(term,
                          levels = c('Minimum temp. * brooding level',
                                     'Minimum temp.')))

## l) Graph of estimates of associations between maternal prenatal social
# stressors and Pan Tissue EAA in mid-childhood
mass_min_temp_brooding_blups_adj_lmer_plot <-
  ggplot(mass_min_temp_brooding_blups_adj_lmer_tbl, aes(x = term, y = estimate,
                                            color = term)) +
  theme_classic() + 
  geom_hline(yintercept = 0, color = 'red',
             linetype = 2) + # line at null behind coefs
  geom_point(size = 6,
             position=position_dodge(width = 0.0)) +
  geom_errorbar(aes(ymin=(conf.low),
                    ymax=(conf.high)), width=.1,
                position=position_dodge(0.0), size = 1) +
  scale_color_manual(values=c("#440154", "#35b779")) +
  theme(strip.background =element_rect(fill= 'white'))+
  theme(strip.text = element_text(colour = 'black')) +
  #coord_flip() + # flip x and y axes
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18)) +
  theme(text = element_text(size=18)) +
  # bold and size title and axes labels
  theme(legend.position = 'none') +
  theme(axis.ticks = element_blank()) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(size = 16, angle=45,
                                   margin = margin(t = 50, r = 0,
                                                   b = 0, l = 0), face = 'bold'), 
        axis.text.y = element_text(size = 16, angle=0,
                                   margin = margin(t = 25, r = 0,
                                                   b = 0, l = 10)),
        axis.title = element_text(size = 18), legend.title=element_blank(),
        legend.text=element_text(size=18),
        legend.position = 'none', #c(0.91, 0.94),
        legend.key = element_blank()) +
  labs(colour = "Model term") +
  xlab(expression(italic('Model term'))) +
  ylab(expression
       (atop(bold('Beta estimate and 95% CI'),
             paste(italic('Nestling mass (g)')))))
# remove axis ticks

print(mass_min_temp_brooding_blups_adj_lmer_plot)

## m) Save Plot
# use ggsave to save the plot
ggsave('mass_min_temp_brooding_blups_adj_lmer_plot.png', plot = mass_min_temp_brooding_blups_adj_lmer_plot, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 10, 
       height = 10, 
       units = c('in'), dpi = 300, limitsize = TRUE) 



## a) Min temp brooding model: extract estimates and tidy the data frame
mass_max_temp_brooding_blups_adj_lmer_est <- tidy(mass_max_temp_brooding_blups_adj_lmer,
                                                  conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(nest_max_temp)' |
           term == 'scale(nest_max_temp):scale(brooding_blups)')


## b) Rename variables for estimates of interest
mass_max_temp_brooding_blups_adj_lmer_est['term'][mass_max_temp_brooding_blups_adj_lmer_est['term'] ==
                                                    'scale(nest_max_temp)'] <- 'Maximum temp.'

mass_max_temp_brooding_blups_adj_lmer_est['term'][mass_max_temp_brooding_blups_adj_lmer_est['term'] ==
                                                    'scale(nest_max_temp):scale(brooding_blups)'] <- 'Maximum temp. * brooding level'


## j) Combine regression estimates into a tidy table
mass_max_temp_brooding_blups_adj_lmer_tbl <- mass_max_temp_brooding_blups_adj_lmer_est

## k) Re-code *nominal* factor (with ordered levels)  
# Set levels (odering) of 'model' variable

mass_max_temp_brooding_blups_adj_lmer_tbl <-
  transform(mass_max_temp_brooding_blups_adj_lmer_tbl,
            term = factor(term,
                          levels = c('Maximum temp. * brooding level',
                                     'Maximum temp.')))

## l) Graph of estimates of associations between maternal prenatal social
# stressors and Pan Tissue EAA in mid-childhood
mass_max_temp_brooding_blups_adj_lmer_plot <-
  ggplot(mass_max_temp_brooding_blups_adj_lmer_tbl, aes(x = term, y = estimate,
                                                        color = term)) +
  theme_classic() + 
  geom_hline(yintercept = 0, color = 'red',
             linetype = 2) + # line at null behind coefs
  geom_point(size = 6,
             position=position_dodge(width = 0.0)) +
  geom_errorbar(aes(ymin=(conf.low),
                    ymax=(conf.high)), width=.1,
                position=position_dodge(0.0), size = 1) +
  scale_color_manual(values=c("#440154", "#35b779")) +
  theme(strip.background =element_rect(fill= 'white'))+
  theme(strip.text = element_text(colour = 'black')) +
  #coord_flip() + # flip x and y axes
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18)) +
  theme(text = element_text(size=18)) +
  # bold and size title and axes labels
  theme(legend.position = 'none') +
  theme(axis.ticks = element_blank()) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(size = 16, angle=45,
                                   margin = margin(t = 50, r = 0,
                                                   b = 0, l = 0), face = 'bold'), 
        axis.text.y = element_text(size = 16, angle=0,
                                   margin = margin(t = 25, r = 0,
                                                   b = 0, l = 10)),
        axis.title = element_text(size = 18), legend.title=element_blank(),
        legend.text=element_text(size=18),
        legend.position = 'none', #c(0.91, 0.94),
        legend.key = element_blank()) +
  labs(colour = "Model term") +
  xlab(expression(italic('Model term'))) +
  ylab(expression
       (atop(bold('Beta estimate and 95% CI'),
             paste(italic('Nestling mass (g)')))))
# remove axis ticks

print(mass_max_temp_brooding_blups_adj_lmer_plot)

## m) Save Plot
# use ggsave to save the plot
ggsave('mass_max_temp_brooding_blups_adj_lmer_plot.png', plot = mass_max_temp_brooding_blups_adj_lmer_plot, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 10, 
       height = 10, 
       units = c('in'), dpi = 300, limitsize = TRUE) 



## a) Min temp brooding model: extract estimates and tidy the data frame
mass_iqr_temp_brooding_blups_adj_lmer_est <- tidy(mass_iqr_temp_brooding_blups_adj_lmer,
                                                  conf.int = T, conf.level = 0.95) %>%
  filter(term == 'scale(nest_iqr_temp)' |
           term == 'scale(nest_iqr_temp):scale(brooding_blups)')


## b) Rename variables for estimates of interest
mass_iqr_temp_brooding_blups_adj_lmer_est['term'][mass_iqr_temp_brooding_blups_adj_lmer_est['term'] ==
                                                    'scale(nest_iqr_temp)'] <- 'IQR of temp.'

mass_iqr_temp_brooding_blups_adj_lmer_est['term'][mass_iqr_temp_brooding_blups_adj_lmer_est['term'] ==
                                                    'scale(nest_iqr_temp):scale(brooding_blups)'] <- 'IQR of temp. * brooding level'


## j) Combine regression estimates into a tidy table
mass_iqr_temp_brooding_blups_adj_lmer_tbl <- mass_iqr_temp_brooding_blups_adj_lmer_est

## k) Re-code *nominal* factor (with ordered levels)  
# Set levels (odering) of 'model' variable

mass_iqr_temp_brooding_blups_adj_lmer_tbl <-
  transform(mass_iqr_temp_brooding_blups_adj_lmer_tbl,
            term = factor(term,
                          levels = c('IQR of temp. * brooding level',
                                     'IQR of temp.')))

## l) Graph of estimates of associations between maternal prenatal social
# stressors and Pan Tissue EAA in mid-childhood
mass_iqr_temp_brooding_blups_adj_lmer_plot <-
  ggplot(mass_iqr_temp_brooding_blups_adj_lmer_tbl, aes(x = term, y = estimate,
                                                        color = term)) +
  theme_classic() + 
  geom_hline(yintercept = 0, color = 'red',
             linetype = 2) + # line at null behind coefs
  geom_point(size = 6,
             position=position_dodge(width = 0.0)) +
  geom_errorbar(aes(ymin=(conf.low),
                    ymax=(conf.high)), width=.1,
                position=position_dodge(0.0), size = 1) +
  scale_color_manual(values=c("#440154", "#35b779")) +
  theme(strip.background =element_rect(fill= 'white'))+
  theme(strip.text = element_text(colour = 'black')) +
  #coord_flip() + # flip x and y axes
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18)) +
  theme(text = element_text(size=18)) +
  # bold and size title and axes labels
  theme(legend.position = 'none') +
  theme(axis.ticks = element_blank()) +
  theme(axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(size = 16, angle=45,
                                   margin = margin(t = 50, r = 0,
                                                   b = 0, l = 0), face = 'bold'), 
        axis.text.y = element_text(size = 16, angle=0,
                                   margin = margin(t = 25, r = 0,
                                                   b = 0, l = 10)),
        axis.title = element_text(size = 18), legend.title=element_blank(),
        legend.text=element_text(size=18),
        legend.position = 'none', #c(0.91, 0.94),
        legend.key = element_blank()) +
  labs(colour = "Model term") +
  xlab(expression(italic('Model term'))) +
  ylab(expression
       (atop(bold('Beta estimate and 95% CI'),
             paste(italic('Nestling mass (g)')))))
# remove axis ticks

print(mass_iqr_temp_brooding_blups_adj_lmer_plot)

## m) Save Plot
# use ggsave to save the plot
ggsave('mass_iqr_temp_brooding_blups_adj_lmer_plot.png', plot = mass_iqr_temp_brooding_blups_adj_lmer_plot, 
       device = NULL, 
       path = 'Output/blups_nestling_lvl/', scale = 1, width = 10, 
       height = 10, 
       units = c('in'), dpi = 300, limitsize = TRUE) 

