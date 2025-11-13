# Import and process data for regression plots
library(arrow)

get_data <- function() {
    behav <- read.csv('analysis/data/MW_EEG_behavioral.csv')
    ecg <- read.csv('analysis/data/MW_ECG_summary_keepIBI.csv')
    alpha <- read_feather('analysis/data/derivatives/spectral/06-psd/occipital_bytrial.feather')
    
    # Format behav
    behav <- behav %>% 
        group_by(subject) %>% 
        mutate(trial = 1:n()) %>% 
        ungroup()
    behav <- get_keep_items(behav)
    
    d <- ecg %>% 
        # Compute trial-level heart metrics
        mutate(trial = as.integer(str_extract(probe, 'Probe(\\d+)', group=1))) %>% 
        select(-probe) %>% 
        mutate(timepoint = sample / 250) %>% 
        group_by(subject, trial) %>% 
        mutate(hr = 60 / (timepoint - lag(timepoint))) %>% 
        ungroup() %>% 
        group_by(subject, trial) %>% 
        summarize(rmssd = sqrt(mean((timepoint - lag(timepoint))^2, na.rm=TRUE)),
                  hr = mean(hr, na.rm=TRUE)) %>% 
        # Merge everything left into ECG
        left_join(behav) %>% 
        left_join(alpha) %>% 
        rename(hrv = rmssd) %>% 
        # Get scaled arou
        group_by(subject) %>% 
        mutate(arou_m = mean(arou, na.rm=TRUE), arou_sd = sd(arou, na.rm=TRUE)) %>% 
        ungroup() %>% 
        mutate(arou_n = (arou - arou_m) / arou_sd) %>% 
        select(-arou_m, -arou_sd) 
    
    # Bads
    bads <- c(10, 13, 14)
    d <- d[!d$subject %in% bads,]
    
    # Remember to drop missing power values for alpha analyses (bad epochs)
    return(d)
    
}


get_keep_items <- function(behav) {
    # Keep the four items with the highest mean individual correlations with arousal
    
    behav <- behav[, !colnames(behav) %in% c('run', 'conf')]
    
    keep <- behav %>% 
        gather(item, response, att:ppl, aff:ling) %>% 
        filter(!is.na(arou), !is.na(response)) %>% 
        group_by(subject, item) %>% 
        summarize(r = cor(response, arou)) %>% 
        group_by(item) %>% 
        summarize(r = mean(r)) %>% 
        arrange(-r) %>% 
        head(n=4) %>% 
        pull(item)
    
    behav <- behav[, c('subject', 'trial', 'arou', keep)]
    return(behav)
    
}