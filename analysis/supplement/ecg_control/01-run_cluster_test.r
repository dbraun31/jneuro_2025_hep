rm(list=ls())
library(tidyverse)
library(reticulate)
library(BayesFactor)
library(brms)
library(glue)
library(ggpubr)
library(here)
library(fs)
library(paletteer)

# Config env
use_condaenv('eeg')
setwd(here())
root <- path('writing/figures/revisions/supplement/ecg_control/')
setwd(root)
py_run_string('from cluster_def import one_sample_test, paired_sample_test')
setwd(here())
np <- import('numpy')
one_sample <- py$one_sample_test
paired_sample <- py$paired_sample_test
text <- 18

# Load data 
ecg <- read.csv(path(root, 'ecg_median.csv'))
# This array is already filtered on time >= 0.25
ecg_arr <- np$load(path(root, 'fm_arr.npy'))


# -- RUN CLUSTER TESTS -- #

result_one_sample <- one_sample(ecg_arr)
result_paired_sample <- paired_sample(ecg)



# --- PLOT ECG --- #

N <- length(unique(ecg$subject))

pal <- paletteer_d('nord::lumina')
cs <- c('Activated', 'Deactivated')
lines <- pal[c(2, 4)]
lines <- setNames(lines, cs)
ribbons <- pal[c(1,3)]
ribbons <- setNames(ribbons, cs)

high <- '#FC8D59'
low <- '#91BFDB'


# General timeseries

p0 <- ecg %>% 
    group_by(time, condition) %>% 
    summarize(voltage = mean(ecg), se = sd(ecg) / sqrt(N)) %>% 
    mutate(time = time - 0.1) %>% 
    filter(!is.na(condition)) %>% 
    ggplot(aes(x = time, y = voltage)) + 
    geom_ribbon(aes(ymin = voltage - se, ymax = voltage + se, fill = condition), alpha = .3) + 
    geom_line(aes(color = condition)) + 
    facet_wrap(~condition, nrow=2, strip.position='top') + 
    labs(
        x = 'Time since heartbeat (s)',
        y = 'ECG Voltage',
        color = '',
        fill = ''
    ) + 
    scale_color_manual(values = c(`Deactivated` = low, `Activated` = high)) +
    scale_fill_manual(values = c(`Deactivated` = low, `Activated` = high)) + 
    scale_y_continuous(breaks = c(0, 8e-04)) + 
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'none',
          text = element_text(size = text),
          strip.background = element_rect(fill = NA))


ecg %>% 
    filter(!is.na(arou)) %>% 
    group_by(subject, trial) %>% 
    summarize(ecg = mean(ecg), arou = unique(arou)) %>% 
    head()


# t test
td <- ecg %>% 
    filter(!is.na(arou)) %>% 
    group_by(subject, condition) %>% 
    summarize(ecg = mean(ecg)) %>% 
    spread(condition, ecg)

t_test <- t.test(td$Activated, td$Deactivated, paired=TRUE)
t <- round(t_test$statistic, 2)
df <- t_test$parameter
p <- round(t_test$p.value, 3)
bf_null <- round(1 / extractBF(ttestBF(td$Activated, td$Deactivated, paired=TRUE))$bf, 3)
full_label <- glue('t({df}) = {t}\np = {p}\nBF[01] = {bf_null}')
label <- glue('BF[01] = {bf_null}')
# t(50) = -0.11
# p = 0.916
# BF[01] = 6.526

p1 <- ecg %>% 
    filter(!is.na(arou)) %>% 
    group_by(subject, condition) %>% 
    summarize(ecg = mean(ecg)) %>% 
    ggplot(aes(x = condition, y = ecg)) + 
    geom_boxplot(aes(fill = condition), alpha = .3, outliers = FALSE) + 
    geom_violin(aes(fill = condition), width = .3, trim = FALSE, alpha = .3) + 
    geom_jitter(width = .02, aes(color = condition), size = .5) + 
    scale_fill_manual(values = c(`Activated` = high, `Deactivated` = low)) + 
    scale_color_manual(values = c(`Activated` = high, `Deactivated` = low)) +
    annotate('text', x = 0.5, y = 8e-5, label=label, hjust=0, size = 6) + 
    labs(
        x = 'Affective arousal (median split)',
        y = 'ECG Voltage',
    ) + 
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = 'none',
          text = element_text(size = text))

N <- length(unique(ecg$subject))
ecg %>% 
    filter(!is.na(arou)) %>% 
    group_by(subject, condition) %>% 
    summarize(ecg_ = mean(ecg)) %>% 
    group_by(condition) %>% 
    summarize(ecg = mean(ecg_), se = sd(ecg_) / sqrt(N))

# Means
# condition          ecg         se
# <chr>            <dbl>      <dbl>
#     1 Activated   0.00000389 0.00000306
# 2 Deactivated 0.00000416 0.00000296
    

# LM
rd <- ecg %>% 
    filter(!is.na(arou)) %>% 
    group_by(subject, trial) %>% 
    summarize(ecg = mean(ecg), arou_n = unique(arou_n)) 

if (!file.exists(path(root, 'posterior.rds'))) {
    m0 <- brm(ecg ~ arou_n, data = rd, iter = 10000, chains = 10, cores = 10)
    post <- posterior_samples(m0)
    saveRDS(post, path(root, file='posterior.rds'))
} else {
    post <- readRDS(path(root, 'posterior.rds'))
}

# - GET POSTERIOR PREDICTIVE DISTRIBUTION - #
pred <- data.frame(quants = seq(-3.5, 3.5, .01))

# Mu mean
pred$mu <- sapply(pred$quants, FUN = function(x) mean(post$b_Intercept) + mean(post$b_arou_n) * x)

# Mu uncertainty
mu_dist <- sapply(pred$quants, FUN = function(x) post$b_Intercept + post$b_arou_n * x)
mu_cis <- data.frame(t(apply(mu_dist, 2, FUN = function(x) quantile(x, probs = c(.025, .975)))))
colnames(mu_cis) <- c('mu_l', 'mu_h')
pred <- cbind(pred, mu_cis)

# Response uncertainty
get_ecg_pred <- function(post, arou_n) {
    mus <- post$b_Intercept + post$b_arou_n * arou_n
    ecgs <- rnorm(length(mus), mus, post$sigma)
    return(ecgs)
}

ecg_preds <- sapply(pred$quants, FUN = function(x) get_ecg_pred(post, x))
ecg_cis <- data.frame(t(apply(ecg_preds, 2, FUN = function(x) quantile(x, probs = c(.025, .975)))))
colnames(ecg_cis) <- c('ecg_l', 'ecg_h')
pred <- cbind(pred, ecg_cis)

pred$pred_l <- apply(ecg_preds, 2, FUN = function(x) quantile(x, probs=.025))
pred$pred_h <- apply(ecg_preds, 2, FUN = function(x) quantile(x, probs=.975))

arous <- ecg %>% 
    filter(!is.na(arou)) %>% 
    group_by(subject, trial) %>% 
    summarize(arou = unique(arou))

# Adjust scale
pred$arou <- pred$quants * sd(arous$arou) + mean(arous$arou)
arou_min <- min(arous$arou)
arou_max <- max(arous$arou)
pred <- pred[pred$arou >= arou_min & pred$arou <= arou_max,]

obs <- ecg %>% 
    group_by(subject, trial) %>% 
    summarize(arou = unique(arou), ecg = mean(ecg))
    
# Summarize and plot
b_mean <- round(mean(post$b_arou_n), 8)
b_l <- round(quantile(post$b_arou_n, probs = .025), 8)
b_h <- round(quantile(post$b_arou_n, probs = .975), 8)
m1 <- lm(ecg ~ arou_n, data = rd)
t <- round(summary(m1)$coefficients['arou_n', 't value'], 2)
df <- m1$df.residual
p <- round(summary(m1)$coefficients['arou_n', 'Pr(>|t|)'], 3)
bf_null <- round(1 / extractBF(lmBF(ecg ~ arou_n, data = rd))$bf, 3)


full_label <- glue('Arousal slope = {b_mean}\n95CIs = [{b_l}, {b_h}]\nt({df}) = {t}\np = {p}\nBF[01] = {bf_null}')
label <- glue('BF[01] = {bf_null}')

# Arousal slope = -1.12e-06
# 95CIs = [-3.62e-06, 1.33e-06]
# t(2355) = -0.88
# p = 0.377
# BF[01] = 14.65

p2 <- pred %>% 
    ggplot(aes(x = arou, y = mu)) + 
    geom_ribbon(aes(ymin = mu_l, ymax = mu_h), alpha = .6, fill = 'steelblue') + 
    geom_ribbon(aes(ymin = ecg_l, ymax = ecg_h), alpha = .3, fill = 'steelblue') + 
    geom_line() + 
    annotate('text', label = label, x = 0, y = 1e-3, hjust = 0, size = 6) + 
    geom_point(data=obs, aes(x = arou, y = ecg), alpha = .2, size = .5) + 
    labs(
        x = 'Affective Arousal',
        y = 'ECG Voltage'
    ) + 
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = text))


# Individual slopes (not using)

p3 <- rd %>% 
    group_by(subject) %>% 
    group_modify(~ {
        m2 <- lm(ecg ~ arou_n, data = .x)
        confs <- confint(m2)
        ci_l <- confs['arou_n', 1]
        ci_h <- confs['arou_n', 2]
        data.frame(b = m2$coefficients['arou_n'], ci_l=ci_l, ci_h=ci_h)
    }) %>% 
    ggplot(aes(x = reorder(subject, b), y = b)) + 
    geom_hline(yintercept = 0, linetype = 'dashed') + 
    geom_errorbar(aes(ymin = ci_l, ymax = ci_h), width = 0) + 
    geom_point() +
    coord_flip() + 
    labs(
        x = 'Subject',
        y = 'Affective arousal slope'
    ) + 
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          panel.grid = element_blank(),
          axis.text.y = element_blank(),
          text = element_text(size = 16))



g <- ggarrange(p0, ggarrange(p1, p2, nrow = 1, labels = c('B.', 'C.'), font.label=list(size=18)), 
               nrow = 2, labels = c('A.', ''), font.label=list(size=18))
g


ggsave(filename = path(root, 'figures/ecg_conditions.png'), plot = g, 
       width = 11, height = 10, dpi = 300, units = 'in')
ggsave(filename = path(root, 'figures/ecg_conditions.tiff'), plot = g, 
       width = 11, height = 10, dpi = 300, units = 'in')







# =======================================================================




# --- ORIGINAL OLD PLOTS --- #
dark_blue <- lines[2]
blue <- ribbons[2]

pd <- ecg %>% 
    filter(!is.na(condition)) %>% 
    group_by(time, condition) %>% 
    summarize(voltage = mean(ecg)) %>% 
    spread(condition, voltage) 

p1 <- pd %>% 
    ggplot(aes(x = Activated, y = Deactivated)) + 
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed') + 
    geom_point(color = dark_blue) +
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          text = element_text(size = 16))

p2 <- pd %>% 
    mutate(difference = Deactivated - Activated) %>% 
    ggplot(aes(x = difference)) + 
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    geom_histogram(fill = blue, color = 'black', alpha = .6) +
    labs(
        x = 'Deactivated - Activated',
        y = 'Frequency\nof participants'
    ) + 
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          text = element_text(size = 16))
