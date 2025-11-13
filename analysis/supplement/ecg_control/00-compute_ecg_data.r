rm(list=ls())
library(tidyverse)
library(reticulate)
library(reshape2)
library(furrr)
library(here)
library(fs)
setwd(here())
use_condaenv('eeg')

# omg when will this be over this is a mess
# make sure the scripts in the writing dir reflect that in the root var

ecg <- read.csv('analysis/scripts/ecg_control/ecg_averaged.csv')
d <- read.csv('analysis/data/MW_EEG_behavioral.csv')
root <- path('writing/figures/revisions/supplement/ecg_control')

# --- DATA FORMATTING --- #
    
d <- d %>% 
    group_by(subject) %>% 
    mutate(trial = 1:n()) %>% 
    relocate(trial, .after = run) %>% 
    select(subject, trial, arou) %>% 
    filter(!(subject %in% c(10, 13, 14)))

ecg <- ecg %>% 
    inner_join(d, by = c('subject', 'trial')) %>% 
    group_by(subject, trial) %>% 
    mutate(time = seq(-.1, .65, 1/250)) %>% 
    relocate(time, .after = sample) %>% 
    group_by(subject) %>% 
    mutate(condition = ifelse(arou > median(arou, na.rm=TRUE), 'Activated', 'Deactivated'),
           arou_n = (arou - mean(arou, na.rm=TRUE)) / sd(arou, na.rm=TRUE)) 


# Write out data for median split cluster test
write.csv(ecg, path(root, 'ecg_median.csv'), row.names=FALSE)


# --- CORRELATION COMPUTATION --- #

flat <- function(data, subject) {
    message(paste0('Subject ', unique(subject)))
    tryCatch({
        v1 <- data$arou_n[!is.na(data$arou_n)]
        v2 <- data$ecg[!is.na(data$arou_n)]
        m <- lm(v2 ~ v1)
        slope <- coef(m)['v1']
        slope
    }, error = function(e) {
        subjects <- unique(data$subject)
        data.frame(subject = subjects, beta = rep(NA_real_, length(subjects)))
        NA
    })
}

plan(multicore, workers = parallel::detectCores() - 1)

flat_maps <- ecg %>% 
    filter(time >= .25 & time <= .45) %>% 
    select(subject, time, ecg, arou_n) %>% 
    group_by(subject, time) %>% 
    nest() %>% 
    mutate(slopes = future_map2(data, subject, flat, .options = furrr_options(seed = TRUE)))

fm <- flat_maps %>% 
    mutate(slope = unlist(slopes)) %>% 
    select(-data, -slopes)

# Fisher transform spearman correlations
fisher_z <- function(v) {
    return((1/2) * log((1 + v) / (1 - v)))
}

# Unbound correlation coefficients 
fm$slope <- fisher_z(fm$slope)
write.csv(fm, path(root, 'fm.csv'), row.names = FALSE)

# Format as array
fm_arr <- acast(fm, subject ~ time, value.var = 'slope')

np <- import('numpy')
np$save(path(root, 'fm_arr.npy'), fm_arr)





