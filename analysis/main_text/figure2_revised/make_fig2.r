rm(list=ls())
library(tidyverse)
library(here)
library(BayesFactor)
library(scales)
library(arrow)
library(tidytext)
library(ggpubr)
library(fs)
library(ggridges)
setwd(here())
script_root <- path('writing/figures/revisions/main_text/figure2_revised')
#source(path(script_root, 'make_fig2_items.r'))
source(path(script_root, 'import.r'))
source(path(script_root, 'format_posteriors.r'))
source(path(script_root, 'sampling.r'))

legend_title <- 8
legend_text <- 6
x_text <- 8
bf_text <- 3
general_text <- 14

# ---------------------- FIGURE 2 ---------------------------- #

# ~~ PANEL A: BETWEEN SUBJECT VARIABILITY IN AROUSAL ~~ #
d <- read.csv('analysis/data/MW_EEG_behavioral.csv')
bads <- c(10, 13, 14)
d <- d[!d$subject %in% bads,]

N <- length(unique(d$subject))
out <- d %>% 
    group_by(subject) %>% 
    summarize(arou = mean(arou, na.rm=TRUE)) %>% 
    summarize(m = mean(arou), sd = sd(arou))
sink(file = 'writing/figures/stats/fig2_arou_unconditional.txt')
out
sink(NULL)

seed <- 179482
set.seed(seed)
p1 <- d %>% 
    filter(subject %in% sample(subject, 10)) %>% 
    mutate(subject = as.factor(subject)) %>% 
    group_by(subject) %>% 
    mutate(arou_m = mean(arou)) %>% 
    ungroup() %>% 
    ggplot(aes(x = arou, y = reorder(subject, arou_m), fill = ..x..)) + 
    geom_density_ridges_gradient(rel_min_height=.01, alpha = .7) + 
    scale_fill_gradientn(colors = colorRampPalette(brewer_pal(palette='Blues')(9))(100)) + 
    labs(
        x = 'Affective arousal',
        y = 'Participant'
    ) + 
    theme_bw() + 
    theme(
        #panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(size = general_text),
        legend.position = 'none'
    )


# ~~ PANEL B: WITHIN SUBJECT VARIABILITY IN AROUSAL ~~ #

seed <- 22018.79
set.seed(seed)
blue <- '#2171B5'
n = 3
p2 <- d %>% 
    select(subject, arou) %>% 
    group_by(subject) %>% 
    mutate(trial = 1:n()) %>% 
    mutate(trials = max(trial)) %>% 
    ungroup() %>%
    filter(trials == 50) %>% 
    filter(subject %in% sample(subject, n)) %>% 
    ggplot(aes(x = trial, y = arou)) + 
    geom_point(color = blue, size = 1.5) + 
    geom_line(color = blue) + 
    facet_grid(subject~.) + 
    labs(
        x = 'Trial',
        y = 'Affective arousal'
    ) + 
    scale_y_continuous(breaks = c(0, 100), labels = c(0, 100)) + 
    theme_bw() + 
    theme(
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(size = general_text)
    )
    


# ~~ PANEL C: HEART METRICS, TIME ON TASK, AND ALPHA POWER ~~ #
# behav, ecg, and alpha all together
d <- get_data()
posteriors <- get_posteriors(d)
group <- format_group(d, posteriors)
subjects <- format_subjects(d, posteriors)


colors <- brewer_pal(palette = 'RdBu')(10)[c(1, 10, 8, 6)]
colors <- c(colors[1], colors[3])

# Format text
group <- group %>% 
    mutate(p_label = ifelse(p < .001, '< .001', paste0('= ', round(p, 3))),
           label = glue('M = {round(mean,2)} [{round(ci_l,2)}, {round(ci_h,2)}]\np {p_label}'),
           contrast = factor(contrast, levels = unlist(contrast_translate))) %>% 
    arrange(contrast)

bigger <- -30
big <- -20
small <- -.5
smaller <- -.02
group$y <- c(big, big, smaller, small, rep(bigger, 4))

p3 <- subjects %>% 
    filter(subject != 59 | contrast != 'Heart rate') %>% 
    mutate(sig = ifelse(ci_h < 0 | ci_l > 0, 'sig', 'notsig'),
           contrast = factor(contrast, levels = unlist(contrast_translate))) %>% 
    ggplot(aes(x = reorder_within(subject, mean, contrast), y = mean)) + 
    geom_hline(yintercept = 0, linetype = 'dashed') + 
    geom_errorbar(aes(ymin = ci_l, ymax = ci_h, color = sig), width = .2, size = .5) + 
    geom_point(aes(color = sig), size = .5) + 
    geom_text(data = group, aes(x = 45, y = y, label = label), hjust=0, size = 3) + 
    labs(
        x = 'Participant',
        y = 'Slope posterior (outcome ~ affective arousal)'
    ) + 
    scale_x_reordered() + 
    scale_color_manual(values = c(`sig` = colors[2], `notsig` = colors[1])) + 
    coord_flip() + 
    facet_wrap(~contrast, scales = 'free', nrow = 2) + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size=x_text),
          strip.background = element_rect(fill = NA),
          axis.ticks = element_blank(),
          legend.position = 'none',
          text = element_text(size = general_text))
    

top <- ggarrange(p1, p2, nrow=1, labels = c('A.', 'B.'))
g <- ggarrange(top, p3, nrow=2, heights = c(1, 2), labels = c('', 'C.'))

png <- path(script_root, 'figure2_revised.png')
tiff <- path(script_root, 'figure2_revised.tiff')

ggsave(png, plot=g, height = 8, width = 11, units = 'in', dpi = 300)
ggsave(tiff, plot=g, height = 8, width = 11, units = 'in', dpi = 300)

write.csv(group, path(script_root, 'group_summary.csv'), row.names = FALSE)
