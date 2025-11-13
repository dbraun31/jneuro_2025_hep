library(tidyverse)
library(here)
library(reticulate)
library(fs)
library(reshape2)
library(lme4)
library(furrr)
library(purrr)
setwd(here())



# --- MAIN FUNCTION --- #

make_maps <- function(dat, statistic='beta') {
    
    # Fit model in parallel
    plan(multicore, workers = parallel::detectCores() - 1)
    
    if (! 'survey_n' %in% colnames(dat)) stop('Need normalized survey responses as column survey_n')
    stat_fun <- ifelse(statistic == 'beta', flat, spearman)
    
    maps <- dat %>%
        filter(time >= .25 & time <= .45, !is.na(survey_n)) %>%
        select(subject, channel, time, survey_n, voltage) %>% 
        group_by(subject, channel, time) %>%
        nest() %>% 
        mutate(nested_stats = future_map2(data, subject, stat_fun, .options = furrr_options(seed = TRUE)))
    
    # --- PROCESS RESULT --- #
    
    out <- maps %>% 
        mutate(statistic = unlist(nested_stats)) %>% 
        select(-data, -nested_stats) %>% 
        filter(subject != 10) 
    
    if (statistic != 'beta') {
        out$statistic <- fisher_z(out$statistic)
    }
    
    return(out)
    
}

# --- MODEL DEFINITIONS --- #

# Fisher transform spearman correlations
fisher_z <- function(v) {
    return((1/2) * log((1 + v) / (1 - v)))
}
    
flat <- function(data, subject) {
    message(paste0('Subject ', unique(subject)))
    tryCatch({
        m <- lm(voltage ~ survey_n, data = data)
        beta <- coef(m)['survey_n']
        beta
    }, error = function(e) {
        subjects <- unique(data$subject)
        data.frame(subject = subjects, statistic = rep(NA_real_, length(subjects)))
    })
}

spearman <- function(data, subject) {
    message(paste0('Subject: ', subject))
    tryCatch({
        # Discard NA ES trials
        v1 <- data$survey_n[!is.na(data$survey_n)]
        v2 <- data$voltage[!is.na(data$survey_n)]
        
        # Get correlation
        r = cor(v1, v2, method='spearman')
        r
    }, error = function(e) {
        e
    })
    
}

