rm(list=ls())
library(tidyverse)
library(scales)
library(ggsci)
library(ggpubr)
library(fs)
library(here)
setwd(here())
source('analysis/scripts/helpers/cluster_tools.r')

general_text <- 16
script_root <- path('writing/figures/revisions/main_text/figure4_revised')

# ----------------- FIGURE 4 ------------------- #
# Anxiety measures and HEP effect

# ~~ PANEL A: CORRELATION BETWEEN ANX MEASURES AND hep MAGNITUDE ~~ #
s <- read.csv('analysis/data/MW_EEG_survey.csv')
s <- s %>% 
    select(subj_id, GAD7_Score, STAI_Score) %>% 
    rename(subject=subj_id, gad=GAD7_Score, stai=STAI_Score)

# Save out stats from both med and num
med <- readRDS('analysis/data/composite/med_cluster.rds')
num <- readRDS('analysis/data/composite/num_cluster.rds')
methods <- list('m' = summary_m(med)[[1]], 'num' = summary_m(num)[[1]])
labels <- data.frame()

for (method in names(methods)) {

    sm <- methods[[method]]
    channels <- sm$cluster_info$full$channels
    times <- sm$cluster_info$full$times
    code <- paste(times, channels, sep='_')
    
    if (method == 'm') {
        eeg <- med$eeg
        d <- eeg %>% 
            mutate(eeg_code = paste(time, channel, sep = '_')) %>% 
            filter(eeg_code %in% code) %>% 
            select(-eeg_code) %>% 
            group_by(subject, condition) %>% 
            summarize(voltage = mean(voltage)) %>% 
            spread(condition, voltage) %>% 
            mutate(effect = Deactivated - Activated) %>% 
            left_join(s) %>% 
            ungroup()
        d_out <- d
    } else {
        eeg <- num$eeg
        d <- eeg %>% 
            mutate(eeg_code = paste(time, channel, sep = '_')) %>% 
            filter(eeg_code %in% code) %>% 
            select(-eeg_code) %>% 
            group_by(subject) %>% 
            summarize(effect = mean(slope)) %>% 
            left_join(s) %>% 
            ungroup() 
        d_reg <- d
    }
        
    
    # GAD cor
    r <- cor.test(d$gad, d$effect)
    label <- paste0('r(', r$parameter, ') = ', round(r$estimate, 2), ', p = ', round(r$p.value, 3))
    gad_r <- data.frame(x = 5, y = 0.25, label=label, measure = 'GAD')
    
    # STAI cor
    sd <- d[!is.na(d$stai),]
    r <- cor.test(sd$stai, sd$effect)
    label <- paste0('r(', r$parameter, ') = ', round(r$estimate, 2), ', p = ', round(r$p.value, 3))
    stai_r <- data.frame(x = 25, y = .25, label=label, measure = 'STAI-S')
    
    out <- rbind(gad_r, stai_r)
    out$method <- method
    labels <- rbind(labels, out)
}

# Median
#d <- d_out

# Regression
d <- d_reg

select_method <- 'num'

write.csv(labels, path(script_root, 'correlation_results.csv'), row.names=FALSE)

labels <- labels[labels$method==select_method, colnames(labels) != 'method']
labels$x <- c(15, 47)
labels$y <- c(.1, .1)

line_col <- pal_npg()(10)[4]

p1 <- d %>% 
    gather(measure, response, gad, stai) %>% 
    mutate(measure = toupper(measure)) %>% 
    mutate(measure = ifelse(measure == 'STAI', 'STAI-S', measure)) %>% 
    ggplot(aes(x = response, y = effect)) + 
    geom_point() + 
    geom_smooth(method = 'lm', color = line_col) + 
    facet_wrap(~measure, scales='free_x') + 
    geom_label(data = labels, aes(x = x, y = y, label=label)) + 
    labs(
        x = 'Survey Score',
        y = 'Average\nHEP-arousal slope'
    ) + 
    theme_bw() + 
    theme(strip.background = element_rect(fill=NA),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          text = element_text(size = general_text))
        

# ~~ PANEL B: HEP TIMESERIES SPLIT BY GAD ~~ #

eeg <- arrow::read_feather('analysis/data/derivatives/hep/06-evoked-clean/hep_arou.feather')
eeg$subject <- as.numeric(eeg$subject)
eeg <- d %>% 
    filter(!is.na(gad)) %>% 
    select(subject, gad) %>% 
    mutate(gad_d = ifelse(gad < median(gad), 'Low GAD', 'High GAD')) %>% 
    inner_join(eeg)

num_s <- summary_m(num)
channels <- num_s[[1]][['cluster_info']][['unique']][['channels']]
sig_times <- range(num_s[[1]][['cluster_info']][['unique']][['times']])

N <- length(unique(eeg$subject))

pd <- eeg %>% 
    gather(channel, voltage, Fp1:POz) %>% 
    filter(channel %in% channels) %>% 
    group_by(subject, time, condition, gad_d) %>% 
    summarize(voltage_ = mean(voltage)) %>% 
    group_by(time, condition, gad_d) %>% 
    summarize(voltage = mean(voltage_), se = sd(voltage_) / sqrt(N)) 

rib <- pd %>% 
    mutate(gad_d = ifelse(gad_d == 'High GAD', 'High GAD', 'Low GAD')) %>% 
    mutate(gad_d = factor(gad_d, levels = c('Low GAD', 'High GAD'))) %>% 
    filter(time >= sig_times[1], time <= sig_times[2]) %>% 
    group_by(time, gad_d) %>% 
    summarize(ymin = min(voltage), ymax = max(voltage))

high <- '#FC8D59'
low <- '#91BFDB'

p2 <- pd %>% 
    mutate(condition = recode(condition, `Activated` = 'High Arousal',
                              `Deactivated` = 'Low Arousal')) %>% 
    mutate(gad_d = ifelse(gad_d == 'High GAD', 'High GAD', 'Low GAD')) %>% 
    mutate(gad_d = factor(gad_d, levels = c('Low GAD', 'High GAD'))) %>% 
    ungroup() %>% 
    ggplot(aes(x = time, y = voltage)) + 
    geom_ribbon(aes(ymin = voltage - se, ymax = voltage + se, fill = condition), alpha = .3) + 
    geom_line(aes(color = condition)) + 
    geom_ribbon(data = rib, aes(ymin=ymin, ymax=ymax, y = 1), alpha = .4, fill = 'darkgreen') + 
    facet_wrap(~gad_d) +
    labs(
        x = 'Time (s) relative to heartbeat',
        y = latex2exp::TeX('Potential ($\\mu ~V$)')
    ) + 
    scale_color_manual(values = c(`High Arousal` = high, `Low Arousal` = low)) + 
    scale_fill_manual(values = c(`High Arousal` = high, `Low Arousal` = low)) + 
    theme_bw() + 
    theme(
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        strip.background = element_rect(fill = NA),
        text = element_text(size = general_text)
    )



g <- ggarrange(p1, p2, nrow = 2, labels = c('A.', 'B.'))


ggsave(path(script_root, 'figure4_revised_base.png'), height = 6, width = 9, units = 'in', dpi = 300)



p2_zoom <- pd %>% 
    mutate(condition = recode(condition, `Activated` = 'High Arousal',
                              `Deactivated` = 'Low Arousal')) %>% 
    mutate(gad_d = ifelse(gad_d == 'High GAD', 'High GAD', 'Low GAD')) %>% 
    mutate(gad_d = factor(gad_d, levels = c('Low GAD', 'High GAD'))) %>% 
    filter(time > .3, time < .4) %>% 
    ggplot(aes(x = time, y = voltage)) + 
    geom_ribbon(aes(ymin = voltage - se, ymax = voltage + se, fill = condition), alpha = .3) + 
    geom_line(aes(color = condition), size = 2) + 
    #geom_ribbon(data = rib, aes(ymin=ymin, ymax=ymax, y = 1), alpha = .4) + 
    facet_wrap(~gad_d, scales='free') +
    labs(
        x = 'Time (s) relative to heartbeat',
        y = latex2exp::TeX('Potential ($\\mu ~V$)')
    ) + 
    scale_x_continuous(breaks = sig_times, labels = sig_times) + 
    scale_color_manual(values = c(`High Arousal` = high, `Low Arousal` = low)) + 
    scale_fill_manual(values = c(`High Arousal` = high, `Low Arousal` = low)) + 
    theme_bw() + 
    theme(
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 25, angle = 45, hjust=1),
        legend.position = 'bottom',
        legend.title = element_blank(),
        strip.background = element_rect(fill = NA),
        text = element_text(size = general_text)
    )


ggsave(path(script_root, 'figure4_zoom_revised.png'), plot = p2_zoom, 
       height = 6, width = 9, units = 'in', dpi = 300)































