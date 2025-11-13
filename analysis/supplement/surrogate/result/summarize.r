rm(list=ls())
library(tidyverse)
library(here)
library(fs)
setwd(here())

null <- scan('result/null_dist_20250627.txt', sep=',')
alt <- as.numeric(readLines('result/alt_sum_t.txt'))
thresholds <- quantile(null, probs = c(.025, .975))

mean(abs(null - mean(null)) >= abs(alt - mean(null)))

p <- (sum(null >= alt) + 1) / (100 + 1)

log_spaced_ints <- function(start, end, n) {
    vals <- exp(seq(log(start), log(end), length.out = n))
    unique(round(vals))
}

breaks <- log_spaced_ints(1, 20, 5)
breaks <- breaks[2:length(breaks)]

null <- data.frame(null = null)

orange <- '#FF7F0E'

p1 <- null %>% 
    ggplot(aes(x = null)) + 
    geom_histogram(fill = 'steelblue', color = 'black', alpha = .7) + 
    geom_vline(xintercept = thresholds[1], linetype = 'dashed', color = 'darkgreen') + 
    geom_vline(xintercept = thresholds[2], linetype = 'dashed', color = 'darkgreen') + 
    geom_point(x = alt, y = 2, shape = 8, size = 7, color = orange) + 
    geom_label(y = 9, x = alt+10, label = paste0('Test statistic: ', round(alt, 3), 
                                              '\np = ', round(p, 3)), hjust=1) + 
    labs(
        x = 'Sum of t values in biggest cluster',
        y = 'Frequency'
    ) + 
    scale_y_continuous(breaks = breaks, labels = breaks) + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank(),
          text = element_text(size = 16), 
          axis.ticks = element_blank())
p1

saveRDS(p1, file='../../../../writing/figures/revisions/supplement/figure_s1/surrogate.rdata')

ggsave('result/surrogate_results.png', width = 11, height = 6, dpi = 120, units = 'in')
ggsave('result/surrogate_results.tiff', width = 11, height = 6, dpi = 300, units = 'in')
