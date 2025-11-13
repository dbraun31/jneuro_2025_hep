rm(list=ls())
library(tidyverse)
library(here)
library(fs)
library(paletteer)
setwd(here())
source('analysis/scripts/helpers/cluster_tools.r')
script_root <- path('writing/figures/revisions/supplement/median_colored/')

# --- IMPORT --- #
m <- readRDS('analysis/data/composite/cluster_result_650ms.rds')
behav <- read.csv('analysis/data/MW_EEG_behavioral.csv')

# --- FORMAT --- #
behav <- behav %>% 
    group_by(subject) %>% 
    mutate(trial = 1:n()) %>% 
    select(subject, trial, arou)

d <- m$eeg
d <- d %>% 
    mutate(trial = as.integer(str_extract(d$probe, 'Probe(\\d+)', group = 1))) %>% 
    relocate(trial, .after = probe) %>% 
    gather(channel, voltage, Fp1:POz) %>% 
    inner_join(behav) %>% 
    filter(!is.na(arou)) %>% 
    group_by(subject) %>% 
    mutate(quarts = case_when(
        arou < quantile(arou, probs = .25) ~ 'Lower quartile',
        arou > quantile(arou, probs = .75) ~ 'Upper quartile',
        TRUE ~ 'Middle'
    )) 

# Visualize split

# Median
scales::show_col(paletteer_d('rcartocolor::Earth'))

pal <- c('#5379A5FF', '#B7AE7AFF')
low <- pal[1]
high <- pal[2]


d <- behav %>% 
    filter(!is.na(arou)) %>% 
    group_by(subject) %>% 
    mutate(condition = ifelse(arou > median(arou), 'Activated', 'Deactivated')) 
    
p <- d %>% 
    ggplot(aes(x = arou)) + 
    geom_histogram(aes(fill = condition)) + 
    facet_wrap(~subject, scales = 'free_y') + 
    scale_fill_manual(values = c(`Deactivated` = low, `Activated` = high)) + 
    labs(
        x = 'Subjective arousal',
        y = 'Frequency',
        fill = ''
    ) + 
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          strip.background = element_rect(fill = NA),
          text = element_text(size = 18),
          legend.title = element_blank(),
          legend.key.size = unit(.4, 'cm'),
          legend.text = element_text(size = 10),
          legend.position = 'bottom',
          axis.text = element_text(size = 8))

p

saveRDS(p, file=path(script_root, '../figure_s1/median.rdata'))
ggsave(path(script_root, 'figures/median_split_subjects.png'), height = 6, width = 11, 
       units = 'in', dpi = 300)
ggsave(path(script_root, 'figures/median_split_subjects.tiff'), height = 6, width = 11, 
       units = 'in', dpi = 300)
    

# --- visualize peak classification --- #

example <- 53

d %>% 
    filter(subject == example) %>% 
    ggplot(aes(x = 1, y = arou)) + 
    geom_boxplot(outliers = FALSE, fill = NA) + 
    geom_jitter(aes(color = condition), width = .08, height = 0) + 
    coord_flip() +
    xlim(.3, 1.8) + 
    labs(
        y = 'Subjective arousal',
        x = '',
        caption = paste0('Example from subject ', example),
        color = ''
    ) + 
    scale_color_manual(values = c(`Deactivated` = low, `Activated` = high)) + 
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          legend.position = c(.9, .9),
          text = element_text(size = 16))

ggsave(path(script_root, 'figures/peak_classification.png'), height = 6, 
       width = 11, units = 'in', dpi = 300)
ggsave(path(script_root, 'figures/peak_classification.tiff'), height = 6, 
       width = 11, units = 'in', dpi = 300)



eeg <- d %>% 
    group_by(subject, trial) %>% 
    summarize(condition_eeg = unique(condition)) 

merge <- behav %>% 
    filter(!subject %in% c(13, 14),
           !is.na(arou)) %>% 
    group_by(subject) %>%
    mutate(condition_b = ifelse(arou > median(arou), 'Activated', 'Deactivated')) %>% 
    inner_join(eeg)
    
    






























