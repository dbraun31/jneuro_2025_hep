library(tidyverse)
library(lme4)
library(here)
setwd(here())

d <- read.csv('analysis/data/derivatives/spectral/06-psd/alpha-parieto_contrast.csv')

d <- d %>% 
    gather(arousal, power, low, high) %>% 
    mutate(power_log = log(power)) 

d %>% 
    group_by(arousal) %>% 
    summarize(power = mean(power_log), se = sd(power_log) / sqrt(n()), N=n())

t <- t.test(d[d$arousal=='low',]$power_log, d[d$arousal=='high',]$power_log, paired=TRUE)
stat_line <- paste0('t(', t$parameter, ') = ', round(t$statistic, 3), ', p = ', round(t$p.value, 3))
lab <- data.frame(arousal='low', power_log = -20, label = stat_line)

d %>% 
    mutate(arousal = factor(arousal, levels = c('low', 'high'))) %>% 
    ggplot(aes(x = arousal, y = power_log)) + 
    geom_boxplot(fill=NA) + 
    geom_point(aes(group = subject)) + 
    geom_line(aes(group=subject), linetype = 'dotted') + 
    geom_text(data=lab, aes(label=label), hjust = 1.2) + 
    labs(
        x = 'Arousal', 
        y = latex2exp::TeX('Mean $\\alpha$ power (Log transform)'),
        caption = paste0('Averaged over electrodes: ', paste0(c('O1', 'O2', 'Oz', 'P3', 'P4', 'P7', 'P8', 'Pz'), collapse = ', '))
    ) + 
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank())
    
ggsave('analysis/scripts/preprocessing/spectral/08-summarize-psd/alpha_arousal.png',
       height=720, width = 1280, units = 'px', dpi = 120)


m1 <- lmerTest::lmer(power_log ~ arousal + (1 | subject), data = d) 





