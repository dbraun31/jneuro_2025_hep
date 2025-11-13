rm(list=ls())
library(arrow)
library(tidyverse)
library(data.table)
library(ggpubr)
library(fs)
library(here)
setwd(here())
script_root <- path('writing/figures/revisions/supplement/epoch_segmenting/')
source(path(script_root, 'compute_statistic.r'))
source(path(script_root, 'cluster_tools.r'))

# --- IMPORT AND PROCESS --- #

import_process <- function() {
    # Import segmented data
    d <- read_feather(path(script_root, 'data/epochs_segmented_wide.feather'))
    # Merge in behav data
    behav <- read.csv('analysis/data/MW_EEG_behavioral.csv')
    behav <- behav %>% 
        group_by(subject) %>% 
        mutate(trial = 1:n()) %>% 
        relocate(trial, .after = subject)
    
    # Reshape to long
    d <- d %>% 
        mutate(trial = as.numeric(str_extract(probe, 'Probe(\\d+)', group=1)), 
               subject = as.numeric(subject)) %>% 
        relocate(trial, .after = sequence) %>% 
        inner_join(behav[, c('subject', 'trial', 'arou')]) %>% 
        relocate(arou, .before = Fp1) %>% 
        gather(channel, voltage, Fp1:POz) %>% 
        group_by(subject) %>% 
        mutate(condition = ifelse(arou > median(arou, na.rm=TRUE), 'Activated', 'Deactivated'),
               survey_n = (arou - mean(arou, na.rm=TRUE)) / sd(arou, na.rm=TRUE)) %>% 
        relocate(condition, survey_n, .before = trial) %>% 
        as.data.table()
    
    # rows blew up to 28 mil
    return(d)
}

file <- path(script_root, 'data/epochs_segmented.feather')
if (!file.exists(file)) {
    d <- import_process()
    write_feather(d, file)
} else d <- read_feather(file)


# --- SEGMENT AND MAKE / IMPORT MAPS --- #

# quick changes 7/2
# first <- d[sequence=='first'][, sequence := NULL]
# last <- d[sequence=='last'][, sequence := NULL]

first_map_file <- path(script_root, 'data/first_map.csv')
last_map_file <- path(script_root, 'data/last_map.csv')

maps_exist <- file.exists(first_map_file) & file.exists(last_map_file)

if (! maps_exist) plan(multicore, workers = parallel::detectCores() - 1)

# Compute maps if map files don't exist
if (maps_exist) {
    first_map <- read.csv(first_map_file)
    last_map <- read.csv(last_map_file)
} else {
    first_map <- make_maps(first)
    last_map <- make_maps(last)
}

if (!maps_exist) {
    write.csv(first_map, first_map_file, row.names = FALSE)
    write.csv(last_map, last_map_file, row.names = FALSE)
}
    


# --- RUN CLUSTER TESTS --- #

# med_first <- permutation_cluster_test(first)
# med_last <- permutation_cluster_test(last)
# 
# first$segment <- 'first'
# last$segment <- 'last'
# med <- rbind(first, last)
med <- d
rm(d)
colnames(med)[colnames(med) == 'sequence'] <- 'segment'

num_first <-permutation_cluster_test(first_map, median_s = FALSE)
num_last <-permutation_cluster_test(last_map, median_s = FALSE)

# first_map$segment <- 'first'
# last_map$segment <- 'last'
# num <- rbind(first_map, last_map)


# Only median last has a sig cluster (~0.04)
# maybe double check why the numeric version isn't coming out


# --- VISUALIZE REGRESSION METHOD --- #

num_s <- summary_m(num_last, threshold=.5)
cluster <- num_s[[1]]
channels <- cluster$cluster_info$unique$channels
times <- cluster$cluster_info$unique$times
p <- cluster$p_value
saveRDS(cluster, file=path(script_root, 'data/num_cluster.rds'))

high <- '#FC8D59'
low <- '#91BFDB'

N <- length(unique(num_last$eeg$subject))

pd <- med %>% 
    filter(channel %in% channels) %>% 
    group_by(subject, time, segment, condition) %>% 
    summarize(voltage_ = mean(voltage)) %>% 
    group_by(time, segment, condition) %>% 
    summarize(voltage = mean(voltage_), se = sd(voltage_) / sqrt(N)) %>% 
    filter(!is.na(condition)) %>% 
    mutate(segment = recode(segment, `first` = '-10 to -5 s', `last` = '-5 to 0 s')) 

write.csv(pd, path(script_root, 'data/pd.csv'), row.names=FALSE)

rib <- pd %>% 
    select(-se) %>% 
    filter(time %in% times) %>% 
    spread(condition, voltage) %>% 
    group_by(time, segment) %>% 
    summarize(voltage_min = pmin(Activated, Deactivated), 
              voltage_max = pmax(Activated, Deactivated)) 


pd %>% 
    ggplot(aes(x = time, y = voltage)) + 
    geom_ribbon(data=rib, aes(ymin = voltage_min, ymax = voltage_max), 
                fill = 'darkgreen', alpha = .4, y = 1) + 
    geom_ribbon(aes(fill = condition, ymin = voltage - se, ymax = voltage + se), alpha = .4) + 
    geom_line(aes(color = condition)) + 
    facet_wrap(~segment, nrow = 2) + 
    labs(
        x = 'Time (s) since heartbeat',
        y = latex2exp::TeX('EEG potential ($\\mu~V$)'),
        color = '', 
        fill = '',
        caption = paste0('p = ', p)
    ) + 
    scale_color_manual(values = c(`Activated` = high, `Deactivated` = low)) +
    scale_fill_manual(values = c(`Activated` = high, `Deactivated` = low)) +
    theme_bw() + 
    theme(legend.position = c(.8, .2),
          panel.grid = element_blank(),
          strip.background = element_rect(fill = NA),
          text = element_text(size = 16),
          axis.ticks = element_blank())

ggsave(filename = path(script_root, 'figures/epoch_segmenting_result.png'),
       height = 6, width = 11, units = 'in', dpi = 300)
ggsave(filename = path(script_root, 'figures/epoch_segmenting_result.tiff'),
       height = 6, width = 11, units = 'in', dpi = 300)

# --- VISUALIZE MEDIAN METHOD --- #
N <- length(unique(last_map$subject))

med_ps <- c(min(med_first$result$p_values), min(med_last$result$p_values))
med_ps <- data.frame(time = .55, voltage = -1e-06, segment = c('-10 – -5 seconds preprobe', '-5 – 0 seconds preprobe'),
                     ps = paste0('p = ', med_ps))

channels_all <- med_last$result$channels
times_all <- med_last$result$times
p_value_idx <- which(med_last$result$p_values < .05)

cluster <- med_last$result$clusters[p_value_idx][[1]]
channels <- channels_all[unique(cluster[[2]]) + 1]
times <- times_all[unique(cluster[[1]]) + 1]


low <- '#FC8D59'
high <- '#91BFDB'

pd <- med %>% 
    filter(channel %in% channels, !is.na(condition)) %>% 
    group_by(subject, condition, segment, time) %>% 
    summarize(voltage_ = mean(voltage, na.rm = TRUE)) %>% 
    group_by(condition, segment, time) %>% 
    summarize(voltage = mean(voltage_), se = sd(voltage_) / sqrt(N)) %>% 
    mutate(segment = recode(segment, `first` = '-10 – -5 seconds preprobe',
                            `last` = '-5 – 0 seconds preprobe')) 

rib <- pd %>% 
    filter(time %in% times) %>% 
    select(-se) %>% 
    spread(condition, voltage) %>% 
    mutate(ymin = pmin(Activated, Deactivated),
           ymax = pmax(Activated, Deactivated)) %>% 
    mutate(ymin = ifelse(segment == '-10 – -5 seconds preprobe', NA, ymin),
           ymax = ifelse(segment == '-10 – -5 seconds preprobe', NA, ymax)) 
    

pd %>% 
    ggplot(aes(x = time, y = voltage)) + 
    geom_text(data = med_ps, aes(label = ps)) + 
    geom_ribbon(data = rib, aes(ymin=ymin, ymax=ymax, y=1), fill = 'green', alpha = .4) +
    # geom_vline(xintercept = min(times), linetype = 'dashed', color = 'orange') + 
    # geom_vline(xintercept = max(times), linetype = 'dashed', color = 'orange') + 
    geom_ribbon(aes(ymin = voltage - se, ymax = voltage + se, fill = condition), alpha = .4) +
    geom_line(aes(color = condition)) + 
    facet_wrap(~segment, ncol=1, strip.position='top') + 
    scale_color_manual(values = c(`Deactivated` = low, `Activated` = high)) +
    scale_fill_manual(values = c(`Deactivated` = low, `Activated` = high)) +
    labs(
        x = 'Time (s) since heartbeat',
        y = latex2exp::TeX('EEG potential ($\\mu~V$)'),
        color = '',
        fill = '',
        caption = paste0('Channels: ', paste(channels, collapse = ', '))
    ) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = NA),
          axis.ticks = element_blank(),
          legend.position = c(.35, .63),
          legend.key.size = unit(.5, 'cm'),
          legend.title = element_blank(),
          text = element_text(size = 18))
    
    

ggsave(filename = path(script_root, 'figures/epoch_segmenting_result.png'),
       height = 6, width = 11, units = 'in', dpi = 120)
ggsave(filename = path(script_root, 'figures/epoch_segmenting_result.tiff'),
       height = 6, width = 11, units = 'in', dpi = 120)



