rm(list=ls())
library(tidyverse)
library(BayesFactor)
library(here)
library(brms)
library(fs)
library(arrow)
setwd(here())
root <- path('writing/figures/revisions/supplement/ecg_control/')
source('analysis/scripts/helpers/cluster_tools.r')


med <- readRDS('analysis/data/composite/med_cluster.rds')
num <- readRDS('analysis/data/composite/num_cluster.rds')
fm <- read.csv(path(root, 'fm.csv'))

med_times <- range(summary_m(med)[[1]]$cluster_info$full$times)
num_times <- range(summary_m(num)[[1]]$cluster_info$full$times)

# Load data
ecg <- read.csv(path(root, 'ecg_median.csv'))
# Add epochs per trial to ecg
epochs <- read_feather(path(root, 'hep_arou.feather'))
epochs <- epochs %>% 
    rename(trial = probe) %>% 
    mutate(subject = as.integer(subject),
           trial = as.integer(str_extract(trial, 'Probe(\\d+)', group=1))) %>% 
    group_by(subject, trial) %>% 
    summarize(nepochs = unique(epochs_per_probe)) 
ecg <- ecg %>% 
    left_join(epochs)

N <- length(unique(ecg$subject))


# --- EPOCHS PER PROBE --- #

# ~ median ~ #

tdat <- ecg %>% 
    filter(!is.na(arou)) %>% 
    group_by(subject, condition, trial) %>% 
    summarize(nepochs = unique(nepochs)) %>% 
    group_by(subject, condition) %>% 
    summarize(nepochs = mean(nepochs)) %>% 
    spread(condition, nepochs)
t_bf <- ttestBF(tdat$Activated, tdat$Deactivated, paired=TRUE)    
t <- t.test(tdat$Activated, tdat$Deactivated, paired=TRUE)
bfnull <- 1 / extractBF(t_bf)$bf    

t_epochs_med <- paste0('t(', t$parameter, ') = ', round(t$statistic, 2), 
                       ' p = ', round(t$p.value, 3), ' BF01 = ', round(bfnull, 3))
marginal_epochs_med <- tdat %>% 
    gather(condition, nepochs, Activated, Deactivated) %>% 
    group_by(condition) %>% 
    summarize(nepochs_ = mean(nepochs), se = sd(nepochs) / sqrt(N))
    

# ~ regression ~ #
rdat <-  ecg %>% 
    filter(!is.na(arou)) %>% 
    group_by(subject, trial) %>% 
    summarize(nepochs = unique(nepochs), arou = unique(arou)) 
m <- lm(nepochs ~ arou, data = rdat)
cis <- confint(m)[2,]
df <- df.residual(m)
m <- summary(m)$coefficients[2,]
t <- regressionBF(nepochs ~ arou, data=rdat)
bfnull <- 1 / extractBF(t)$bf
r_epochs_num <- paste0('b = ', round(m[1], 3), ' CIs = [', round(cis[1], 3), ', ', round(cis[2], 3),
                       '] t(', df, ') = ', round(m[3], 2), ' p = ', round(m[4], 3),
                       ' BF01 = ', round(bfnull, 3))

# --- HEART RATE --- #

hr <- read.csv(path(root, 'heartrate.csv'))
hr <- hr[hr$subject %in% ecg$subject,]

# Median
tdat <- hr %>% 
    group_by(subject, arou) %>% 
    summarize(heartrate = mean(heartrate)) %>% 
    spread(arou, heartrate)

t_bf <- ttestBF(tdat$activated, tdat$deactivated, paired=TRUE)
bf01 <- round(1 / extractBF(t_bf)$bf, 3)
t <- t.test(tdat$activated, tdat$deactivated, paired=TRUE)
t_val <- round(t$statistic, 3)
df <- t$parameter
cis <- round(t$conf.int, 3)
p <- round(t$p.value, 3)
t_hr_med <- paste0('t(', df, ') = ', t_val, ' CIs = [', cis[1], ', ', cis[2], '] p = ', p, ' BF01 = ', bf01)
marginal_hr_med <- tdat %>% 
    gather(condition, hr, activated, deactivated) %>% 
    group_by(condition) %>% 
    summarize(hr_ = mean(hr), se = sd(hr) / sqrt(N))


# Regression

behav <- ecg %>% 
    filter(!is.na(arou)) %>% 
    group_by(subject, trial) %>% 
    summarize(arou = unique(arou))

hr_num <- hr %>% 
    rename(trial = probe, condition = arou) %>% 
    inner_join(behav) 

quick <- hr_num %>% 
    group_by(trial) %>% 
    summarize(hr = mean(heartrate), arou = mean(arou))


# Frequentist
# Flat
m <- lm(heartrate ~ arou, data = hr_num)
coef <- summary(m)$coefficients[2,]
b <- round(coef[1], 3)
cis <- round(confint(m)[2,], 3)
t <- round(coef[3], 3)
df <- df.residual(m)
p <- ifelse(coef[4] < .001, ' p < .001', paste0(' p = ', round(coef[4], 3)))
r_hr_num_flat <- paste0('b = ', b, ', CIs = [', cis[1], ', ', cis[2], '], t(', df, ') = ', t, p)


# Hierarchical
hr_num$arou_n <- scale(hr_num$arou)[,1]
m1 <- lmerTest::lmer(heartrate ~ arou_n + (1 + arou_n | subject), data = hr_num)
cis <- round(confint(m1)[6,], 3)
coefs <- summary(m1)$coefficients[2,]
df <- round(coefs['df'], 3)
beta <- round(coefs['Estimate'], 3)
t <- round(coefs['t value'], 3)
p <- round(coefs['Pr(>|t|)'], 3)
r_hr_num_hier <- paste0('beta = ', beta, ', CIs = [', cis[1], ', ', cis[2], '], t(', df, ') = ', t, ' p = ', p)

# Bayesian
options(mc.cores = parallel::detectCores() - 1)
m1_b <- brms::brm(heartrate ~ arou_n + (1 + arou_n | subject), 
                  data = hr_num,
                  sample_prior = 'yes',
                  save_pars = save_pars(all=TRUE))
m1_b_null <- brms::brm(heartrate ~ 1 + (1 + arou_n | subject), 
                  data = hr_num,
                  sample_prior = 'yes',
                  save_pars = save_pars(all=TRUE))

posterior <- as_draws_df(m1_b)$b_arou_n
# DOUBLE CHECK THIS
# I think this bf is already 10
bf <- brms::bayes_factor(m1_b, m1_b_null)
cis_b <- round(quantile(posterior, probs = c(.025, .975)), 3)
bf10 <- round(1 / bf$bf, 3)

r_hr_num_hier <- paste0(r_hr_num_hier, ', CIs_b = [', cis_b[1], ', ', cis_b[2], '], BF01 = ', bf10)



# --- MEAN ECG ---#

# ~ median ~ #
tdat <- ecg %>% 
    select(subject, time, condition, ecg) %>% 
    # Time window is 0.328 - 0.364
    filter(time >= med_times[1] & time <= med_times[2]) %>% 
    group_by(subject, condition) %>% 
    summarize(ecg = mean(ecg)) %>% 
    spread(condition, ecg) %>% 
    select(subject:Deactivated) 

t <- t.test(tdat$Activated, tdat$Deactivated, paired = TRUE)
t_val <- round(t$statistic, 3)
df <- t$parameter
cis <- round(t$conf.int, 3)
p <- round(t$p.value, 3)
t_bf <- ttestBF(tdat$Activated, tdat$Deactivated, paired = TRUE)
bf01 <- round(1 / extractBF(t_bf)$bf, 3)
t_mean_med <- paste0('t(', df, ') = ', t_val, ' p = ', p, ' BF01 = ', bf01)
marginal_mean_med <- tdat %>% 
    gather(condition, voltage, Activated, Deactivated) %>% 
    group_by(condition) %>% 
    summarize(voltage_ = mean(voltage), se = sd(voltage) / sqrt(N))

# ~ regression ~ #

rdat <- ecg %>% 
    filter(!is.na(arou),
           time >= num_times[1] & time <= num_times[2]) %>% 
    group_by(subject, trial) %>% 
    summarize(arou = unique(arou), ecg = mean(ecg)) 

m <- lm(ecg ~ arou, data = rdat)
coefs <- summary(m)$coefficients[2, ]
b <- coefs[1]
t <- round(coefs[3], 3)
p <- round(coefs[4], 3)
cis <- confint(m)[2,]
df <- df.residual(m)
reg <- regressionBF(ecg ~ arou, data = rdat)    
bf01 <- round(1 / extractBF(reg)$bf, 3)

r_mean_num <- paste0('b = ', b, ', CIs = [', cis[1], ', ', cis[2], '], t(', df, ') = ', t, ', p = ', p, ', BF01 = ', bf01)


# --- PEAK ECG --- #
# Median
tdat <- ecg %>% 
    select(subject, time, condition, ecg) %>% 
    # Time window is 0.328 - 0.364
    filter(time >= med_times[1] & time <= med_times[2]) %>% 
    group_by(subject, condition) %>% 
    summarize(ecg = max(ecg)) %>% 
    spread(condition, ecg) %>% 
    select(subject:Deactivated) 


t <- t.test(tdat$Activated, tdat$Deactivated, paired = TRUE)
t_val <- round(t$statistic, 3)
df <- t$parameter
cis <- round(t$conf.int, 3)
p <- round(t$p.value, 3)
t_bf <- ttestBF(tdat$Activated, tdat$Deactivated, paired = TRUE)
bf01 <- round(1 / extractBF(t_bf)$bf, 3)

t_peak_med <- paste0('t(', df, ') = ', t_val, ', p = ', p, ', BF01 = ', bf01)

marginal_peak_med <- tdat %>% 
    gather(condition, voltage, Activated, Deactivated) %>% 
    group_by(condition) %>% 
    summarize(voltage_ = mean(voltage), se = sd(voltage) / sqrt(N))

# Regression

rdat <- ecg %>% 
    filter(!is.na(arou),
           time >= num_times[1] & time <= num_times[2]) %>% 
    group_by(subject, trial) %>% 
    summarize(arou = unique(arou), ecg = max(ecg)) 
    
m <- lm(ecg ~ arou, rdat)
coefs <- summary(m)$coefficients[2,]
b <- coefs[1]
cis <- confint(m)[2,]
df <- df.residual(m)
t <- coefs[3]
p <- round(coefs[4], 3)
mb <- regressionBF(ecg ~ arou, rdat)
bf01 <- round(1 / extractBF(mb)$bf, 3)

r_peak_num <- paste0('b = ', b, ', CIs = [', cis[1], ', ', cis[2], '], t(', df, ') = ', t_val, ', p = ', p, ', BF01 = ', bf01)


sink(file=path(root, 'bayes_factor_result.txt'))

print('-------EPOCHS PER PROBE --------')
print(marginal_epochs_med)
print(t_epochs_med)
print(r_epochs_num)

print('')
print('')
print('')
print('')

print('------ HEART RATE ------ ')
print(marginal_hr_med)
print(t_hr_med)
print('FLAT')
print(r_hr_num_flat)
print('HIERARCHICAL')
print(r_hr_num_hier)

print('')
print('')
print('')
print('')

print('-------- ECG MEAN --------')
print(marginal_mean_med)
print(t_mean_med)
print(r_mean_num)
print('')
print('')
print('')
print('')

print('-------- ECG PEAK --------')
print(marginal_peak_med)
print(t_peak_med)
print(r_peak_num)


sink(NULL)


# omg pls be good now
