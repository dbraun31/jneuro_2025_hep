rm(list=ls())
library(tidyverse)
library(ggpubr)
library(here)
library(fs)
setwd(here())
source('analysis/scripts/helpers/cluster_tools.r')

script_root <- path('writing/figures/revisions/supplement/extra_cluster/')

low <- '#FC8D59'
high <- '#91BFDB'
text_size <- 16
label_size <- 4

# Import and test
d <- readRDS('analysis/data/composite/med_cluster.rds')$eeg
m <- permutation_cluster_test(d)


# Parse result
t_obs <- m$result$t_obs
cluster <- summary_m(m)$`Cluster 6`

sig_times <- cluster$cluster_info$unique$times
sig_channels <- cluster$cluster_info$unique$channels



# --- POTENTIAL PLOT --- #


N <- length(unique(d$subject))
pd <- d %>% 
    filter(channel %in% sig_channels) %>% 
    group_by(subject, time, condition) %>% 
    summarize(voltage_ = mean(voltage)) %>% 
    group_by(condition, time) %>% 
    summarize(voltage = mean(voltage_), se = sd(voltage_) / sqrt(N)) 
    
rib <- pd %>% 
    select(-se) %>% 
    filter(time %in% sig_times) %>% 
    spread(condition, voltage) %>% 
    mutate(ymin = pmin(Activated, Deactivated),
           ymax = pmax(Activated, Deactivated)) 
    
p <- cluster$p_value

p1 <- pd %>% 
    ggplot(aes(x = time, y = voltage)) + 
    geom_vline(xintercept = 0, color = 'lightgrey') + 
    geom_hline(yintercept = 0, color = 'lightgrey') + 
    geom_ribbon(aes(fill = condition, ymin = voltage - se, ymax = voltage + se), alpha = .4) + 
    geom_line(aes(color = condition)) + 
    geom_ribbon(data = rib, aes(ymin = ymin, ymax = ymax, y = 1), alpha = .4, fill = 'green') + 
    annotate('text', x = .25, y = -.15, label = paste0('p = ', round(p, 3)), size = label_size) + 
    labs(
        x = 'Time (s) since heartbeat',
        y = latex2exp::TeX('EEG potential ($\\mu~V$)'),
        color = '',
        fill = '',
        caption = paste0('Channels: ', paste(sig_channels, collapse = ', '))
    ) + 
    scale_color_manual(values = c(`Deactivated` = low, `Activated` = high)) + 
    scale_fill_manual(values = c(`Deactivated` = low, `Activated` = high)) + 
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = c(.8, .25),
          text = element_text(size = text_size),
          legend.text = element_text(size = 9),
          legend.key.size = unit(.2, 'cm'))

# --- INDIVIDUAL PLOT --- #

cluster_code <- paste(sig_times, sig_channels, sep = '_')

pd <- d %>% 
    mutate(time_channel = paste(time, channel, sep = '_')) %>% 
    filter(time_channel %in% cluster_code) %>% 
    group_by(subject, condition) %>% 
    summarize(voltage = mean(voltage)) 

p2 <- pd %>% 
    ggplot(aes(x = condition, y = voltage)) + 
    geom_violin(aes(fill = condition), alpha = .6)  +
    geom_boxplot(aes(fill = condition), outliers = FALSE, alpha = .6) + 
    geom_line(aes(group = subject), linetype = 'dashed', alpha = .6) + 
    geom_jitter(width = .01, alpha = .6) +
    labs(
        x = '',
        y = latex2exp::TeX('EEG potential ($\\mu~V$)'),
        fill = ''
    ) + 
    scale_fill_manual(values = c(`Deactivated` = low, `Activated` = high)) +
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          text = element_text(size = text_size),
          legend.position = 'none')




# --- TOPO PLOT --- #

p3 <- plot_topo(m, 6, n_breaks = 4)



# --- COMBINE --- #

g1 <- ggarrange(p1, p2, nrow = 1, labels = c('A.', 'B.'))
g <- ggarrange(g1, p3, nrow = 2, labels = c('', 'C.'))
g

ggsave(filename = path(script_root, 'figures/extra_cluster.png'), plot = g,
      height = 6, width = 11, units = 'in', dpi = 300)

ggsave(filename = path(script_root, 'figures/extra_cluster.tiff'), plot = g,
      height = 6, width = 11, units = 'in', dpi = 300)



































