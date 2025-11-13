rm(list=ls())
library(tidyverse)
library(ggpubr)
library(fs)
library(here)
library(reticulate)
setwd(here())
use_condaenv('eeg')
source('analysis/scripts/helpers/cluster_tools.r')

high <- '#FC8D59'
low <- '#91BFDB'
text_size <- 16
annotation_size <- 4

med <- readRDS('analysis/data/composite/med_cluster.rds')
num <- readRDS('analysis/data/composite/num_cluster.rds')

# PANEL C - MED TIMESERIES
num_cluster <- summary_m(num)[[1]]
times <- num_cluster$cluster_info$full$times
channels <- num_cluster$cluster_info$full$channels
p <- num_cluster$p_value
N <- length(unique(med$eeg$subject))
sig_times <- c(min(num_cluster$cluster_info$unique$times), max(num_cluster$cluster_info$unique$times))


pd <- med$eeg %>% 
    filter(channel %in% channels) %>% 
    group_by(subject, condition, time) %>% 
    summarize(voltage_ = mean(voltage)) %>% 
    group_by(condition, time) %>% 
    summarize(voltage = mean(voltage_), se = sd(voltage_) / sqrt(N)) 


sig_times <- num_cluster$cluster_info$full$times

rib <- pd %>% 
    filter(time %in% sig_times) %>% 
    group_by(time) %>% 
    summarize(ymin = min(voltage), ymax = max(voltage))


panel_c <- pd %>% 
    ggplot(aes(x = time, y = voltage)) + 
    geom_vline(xintercept = 0, color = 'lightgrey') + 
    geom_hline(yintercept = 0, color = 'lightgrey') + 
    geom_ribbon(aes(fill = condition, ymin = voltage - se, ymax = voltage + se), alpha = .4) + 
    geom_line(aes(color = condition)) +
    geom_ribbon(data=rib, aes(x = time, ymax = ymax, ymin = ymin, y = 1), alpha = .4, fill='darkgreen') + 
    annotate('text', x = .42, y= -.2, label = paste0('p = ', round(p, 3)), size = annotation_size) + 
    labs(
        x = 'Time (s) since heartbeat',
        y = latex2exp::TeX('EEG potential $\\mu~V$'),
        color = '',
        fill = ''
    ) +
    scale_color_manual(values = c(`Deactivated` = low, `Activated` = high)) + 
    scale_fill_manual(values = c(`Deactivated` = low, `Activated` = high)) + 
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          #axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 14)),
          legend.position = c(.8, .3),
          text = element_text(size = text_size))
    
    
    
# PANEL D - INDIVIDUALS
sig_code <- paste(num_cluster$cluster_info$full$times, num_cluster$cluster_info$full$channels, sep = '_')

pd <- med$eeg %>% 
    mutate(code = paste(time, channel, sep = '_')) %>% 
    filter(code %in% sig_code) %>% 
    group_by(subject, condition) %>% 
    summarize(voltage = mean(voltage)) 

panel_d <- pd %>% 
    ggplot(aes(x = condition, y = voltage)) + 
    geom_boxplot(aes(fill = condition), outliers = FALSE, alpha = .4) +
    geom_violin(aes(color = condition, fill = condition), alpha = .4) + 
    geom_line(aes(group = subject), linetype = 'dashed', alpha = .6, size = .7) + 
    geom_jitter(width = .01, alpha = .6, size = .7) + 
    labs(
        x = '',
        y = latex2exp::TeX('EEG potential ($\\mu~V$)')
    ) + 
    scale_color_manual(values = c(`Deactivated` = low, `Activated` = high)) + 
    scale_fill_manual(values = c(`Deactivated` = low, `Activated` = high)) + 
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'none',
          text = element_text(size = text_size)
    )
    

# PANEL A - NUM CLUSTER
num_cluster <- summary_m(num)
sig_channels <- num_cluster$`Cluster 5`$cluster_info$unique$channels
sig_times <- num_cluster$`Cluster 5`$cluster_info$unique$times
sig_times <- c(min(sig_times), max(sig_times))
p <- round(num_cluster$`Cluster 5`$p_value, 3)

panel_a <- num$eeg %>% 
    filter(channel %in% sig_channels) %>% 
    group_by(subject, time) %>% 
    summarize(slope_ = mean(slope)) %>% 
    group_by(time) %>% 
    summarize(slope = mean(slope_), se = sd(slope_) / sqrt(N)) %>% 
    mutate(condition = ifelse(slope > 0, 'Activated', 'Deactivated')) %>% 
    spread(condition, slope) %>% 
    gather(condition, slope, Activated:Deactivated) %>% 
    ggplot(aes(x = time, y = slope)) + 
    geom_vline(xintercept = 0, color = 'lightgrey') +
    geom_hline(yintercept = 0, color = 'lightgrey') +
    geom_ribbon(aes(ymin = slope - se, ymax = slope + se, fill = condition), alpha = .4) + 
    geom_line(aes(color = condition)) + 
    geom_vline(xintercept = sig_times[1], linetype = 'dashed', color = 'darkgreen') + 
    geom_vline(xintercept = sig_times[2], linetype = 'dashed', color = 'darkgreen') + 
    scale_color_manual(values = c(`Deactivated` = low, `Activated` = high)) + 
    scale_fill_manual(values = c(`Deactivated` = low, `Activated` = high)) + 
    #annotate('text', x = .43, y = -.1, label = paste0('p = ', p), size = annotation_size) + 
    labs(
        x = 'Time (s) since heartbeat',
        y = 'HEP-arousal slope',
        color = '',
        fill = ''
    ) + 
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'none',
          text = element_text(size = text_size))



# PANEL b - TOPO
panel_b <- plot_topo(num, which(num$result$p_values < .05)[1], n_breaks = 1, nrow=1)




# Composite

lw <- .65
# top <- ggarrange(panel_a, panel_b, labels = c('A.', 'B.'), ncol = 2)
# bottom <- ggarrange(panel_c, panel_d, ncol = 2, widths = c(lw, 1-lw), labels = c('C.', 'D.'))
# g <- ggarrange(top, bottom, nrow=2)

# Swapping order of panels
g <- ggarrange(panel_a, panel_b, panel_c, panel_d,
               labels = c('A.', 'B.', 'C.', 'D.'),
               ncol = 2, nrow=2, widths = c(lw, 1-lw))
g

root <- path('writing/figures/revisions/main_text/figure3_revised/')
ggsave(path(root, 'figure3_revised.png'), plot = g, height = 6, width = 11, units = 'in', dpi = 300)
ggsave(path(root, 'figure3_revised.tiff'), plot = g, height = 6, width = 11, units = 'in', dpi = 300)










































