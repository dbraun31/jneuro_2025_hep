library(tidyverse)
library(here)
library(fs)
library(BayesFactor)

setwd(path(here(), 'analysis/scripts/preprocessing/cfa_control'))

hr <- read.csv('heartrate.csv')
amp <- read.csv('ecga.csv')

bads <- c(10, 13, 14)

hr <- hr[!hr$subject %in% bads,]
amp <- amp[!amp$subject %in% bads,]

# Heart rate

hrs <- hr %>% 
    group_by(subject, arou) %>% 
    summarize(hr = mean(heartrate)) %>% 
    spread(arou, hr)

hr_t <- t.test(hrs$activated, hrs$deactivated, paired=TRUE)
bf_h0_hr <- 1 / extractBF(ttestBF(hrs$activated, hrs$deactivated, paired=TRUE))$bf

N <- length(unique(hr$subject))

hr_summary <- hrs %>% 
    gather(arou, hr, activated:deactivated) %>% 
    group_by(arou) %>% 
    summarize(m = mean(hr), se = sd(hr) / sqrt(N))
    
# ECG amplitude

amps <- amp %>% 
    group_by(subject, probe, arou) %>% 
    summarize(peak = mean(peak_ecg)) %>% 
    group_by(subject, arou) %>% 
    summarize(peak_ = mean(peak)) %>% 
    spread(arou, peak_)

amp_t <- t.test(amps$activated, amps$deactivated, paired=TRUE)
bf_h0_amp <- 1/extractBF(ttestBF(amps$activated, amps$deactivated, paired=TRUE))$bf

N <- length(unique(amp$subject))

amp_summary <- amps %>% 
    gather(arou, amp, activated:deactivated) %>% 
    group_by(arou) %>% 
    summarize(m = mean(amp), se = sd(amp) / sqrt(N))


sink(file='summary.txt')
print('HEART RATE')
hr_summary
hr_t
print('BF null')
bf_h0_hr
print('--------------')
print('ECG AMPLITUDE')
amp_summary
amp_t
print('BF null')
bf_h0_amp
sink(NULL)