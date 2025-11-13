rm(list=ls())
library(tidyverse)
library(here)
library(ggpubr)
library(fs)
setwd(here())
root <- path('writing/figures/revisions/supplement/figure_s2')

ecg <- readRDS(path(root, 'ecg_control.rdata'))
med <- readRDS(path(root, 'median.rdata'))
surrogate <- readRDS(path(root, 'surrogate.rdata'))

g <- ggarrange(med, surrogate, labels = c('A.', 'B.'),
               heights = c(.75, .25), nrow = 2)

ggsave(filename = path(root, 'figure_s2.png'), plot = g, 
       width = 17.6, height = 25, unit = 'cm', dpi = 300)
ggsave(filename = path(root, 'figure_s2.tiff'), plot = g, 
       width = 17.6, height = 25, unit = 'cm', dpi = 300)
