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
    }
        
    
    # GAD cor
    r <- cor.test(d$gad, d$effect)
    label <- paste0('r(', r$parameter, ') = ', round(r$estimate, 2), ', p = ', round(r$p.value, 3))
    gad_r <- data.frame(x = 15, y = 0.9, label=label, measure = 'GAD')
    
    # STAI cor
    sd <- d[!is.na(d$stai),]
    r <- cor.test(sd$stai, sd$effect)
    label <- paste0('r(', r$parameter, ') = ', round(r$estimate, 2), ', p = ', round(r$p.value, 3))
    stai_r <- data.frame(x = 30, y = .9, label=label, measure = 'STAI-S')
    
    out <- rbind(gad_r, stai_r)
    out$method <- method
    labels <- rbind(labels, out)
}

d <- d_out

write.csv(labels, path(script_root, 'correlation_results.csv'), row.names=FALSE)

labels <- labels[labels$method=='m', colnames(labels) != 'method']

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
        y = 'HEP Effect\n(Low - High Arousal)'
    ) + 
    theme_bw() + 
    theme(strip.background = element_rect(fill=NA),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          text = element_text(size = general_text))
        

# ~~ PANEL B: HEP TIMESERIES SPLIT BY STAI ~~ #

eeg <- arrow::read_feather('analysis/data/derivatives/hep/06-evoked-clean/hep_arou.feather')
eeg$subject <- as.numeric(eeg$subject)
eeg <- d %>% 
    filter(!is.na(stai)) %>% 
    select(subject, stai) %>% 
    mutate(stai_d = ifelse(stai < median(stai), 'Low STAI', 'High STAI')) %>% 
    inner_join(eeg)

channels <- c('Fz', 'FC1', 'FC5', 'Fp1', 'F3', 'F7', 'F4', 'Fp2', 'FC2', 'F8')
sig_times <- c(.328, .364)

N <- length(unique(eeg$subject))

pd <- eeg %>% 
    gather(channel, voltage, Fp1:POz) %>% 
    filter(channel %in% channels) %>% 
    group_by(subject, time, condition, stai_d) %>% 
    summarize(voltage_ = mean(voltage)) %>% 
    group_by(time, condition, stai_d) %>% 
    summarize(voltage = mean(voltage_), se = sd(voltage_) / sqrt(N)) 

rib <- pd %>% 
    mutate(stai_d = ifelse(stai_d == 'High STAI', 'High STAI-S', 'Low STAI-S')) %>% 
    mutate(stai_d = factor(stai_d, levels = c('Low STAI-S', 'High STAI-S'))) %>% 
    filter(time >= sig_times[1], time <= sig_times[2]) %>% 
    group_by(time, stai_d) %>% 
    summarize(ymin = min(voltage), ymax = max(voltage))

low <- '#FC8D59'
high <- '#91BFDB'

p2 <- pd %>% 
    mutate(condition = recode(condition, `Activated` = 'High Arousal',
                              `Deactivated` = 'Low Arousal')) %>% 
    mutate(stai_d = ifelse(stai_d == 'High STAI', 'High STAI-S', 'Low STAI-S')) %>% 
    mutate(stai_d = factor(stai_d, levels = c('Low STAI-S', 'High STAI-S'))) %>% 
    ungroup() %>% 
    ggplot(aes(x = time, y = voltage)) + 
    geom_ribbon(aes(ymin = voltage - se, ymax = voltage + se, fill = condition), alpha = .3) + 
    geom_line(aes(color = condition)) + 
    geom_ribbon(data = rib, aes(ymin=ymin, ymax=ymax, y = 1), alpha = .4) + 
    facet_wrap(~stai_d) +
    labs(
        x = 'Time (s) relative to heartbeat',
        y = latex2exp::TeX('Potential ($\\mu ~V$)')
    ) + 
    # To zoom
    #scale_x_continuous(breaks = sig_times, labels = sig_times) + 
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
    mutate(stai_d = ifelse(stai_d == 'High STAI', 'High STAI-S', 'Low STAI-S')) %>% 
    mutate(stai_d = factor(stai_d, levels = c('Low STAI-S', 'High STAI-S'))) %>% 
    filter(time > .3, time < .4) %>% 
    ggplot(aes(x = time, y = voltage)) + 
    geom_ribbon(aes(ymin = voltage - se, ymax = voltage + se, fill = condition), alpha = .3) + 
    geom_line(aes(color = condition), size = 2) + 
    #geom_ribbon(data = rib, aes(ymin=ymin, ymax=ymax, y = 1), alpha = .4) + 
    facet_wrap(~stai_d) +
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
        legend.position = 'bottom',
        legend.title = element_blank(),
        strip.background = element_rect(fill = NA),
        text = element_text(size = general_text)
    )


ggsave(path(script_root, 'figure_4_zoom_revised.png'), plot = p2_zoom, 
       height = 6, width = 9, units = 'in', dpi = 300)









# p1 b

eeg <- num$eeg
d <- eeg %>% 
    mutate(eeg_code = paste(time, channel, sep = '_')) %>% 
    filter(eeg_code %in% code) %>% 
    select(-eeg_code) %>% 
    group_by(subject) %>% 
    summarize(effect = mean(slope)) %>% 
    left_join(s) %>% 
    ungroup() 
    
# GAD cor
r <- cor.test(d$gad, d$effect)
label <- paste0('r(', r$parameter, ') = ', round(r$estimate, 2), ', p = ', round(r$p.value, 3))
gad_r <- data.frame(x = 15, y = 0.25, label=label, measure = 'GAD')

# STAI cor
sd <- d[!is.na(d$stai),]
r <- cor.test(sd$stai, sd$effect)
label <- paste0('r(', r$parameter, ') = ', round(r$estimate, 2), ', p = ', round(r$p.value, 3))
stai_r <- data.frame(x = 30, y = .25, label=label, measure = 'STAI-S')

labels_reg <- rbind(gad_r, stai_r)

p1b <- d %>% 
    gather(measure, response, gad:stai) %>% 
    mutate(measure = recode(measure, `gad` = 'GAD', `stai` = 'STAI-S')) %>% 
    ggplot(aes(x = response, effect)) + 
    geom_point() + 
    geom_smooth(method='lm') + 
    geom_label(data=labels_reg, aes(label=label, x=x, y=y)) + 
    facet_wrap(~measure, scales='free_x') + 
    labs(
        title = 'Regression method',
         x = 'Survey Score',
         y = 'Average slope between\nEEG and arousal'
    ) + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          strip.background = element_rect(fill = NA),
          text = element_text(size = 16))

p1a <- p1 + labs(title = 'Median-split method')

g <- ggarrange(p1a, p1b, nrow = 2)
g

# For screens
ggsave(filename=path(script_root, 'method_compare.png'), plot=g,
       height = 1080, width = 1920, units = 'px', dpi = 120)






















