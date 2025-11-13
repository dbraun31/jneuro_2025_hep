library(tidyverse)
library(here)
# This is what eventually worked to install 
#source("https://raw.githubusercontent.com/apache/arrow/main/r/R/install-arrow.R")
#install_arrow()
library(arrow)
setwd(here())

eeg <- read_feather('analysis/data/derivatives/06-evoked-clean/hep_att.feather')
eeg$subject <- as.numeric(eeg$subject)
d <- read.csv('analysis/data/MW_EEG_behavioral.csv')

N_og <- length(unique(eeg$subject))

d <- d[d$subject %in% eeg$subject,]

d <- d %>% 
    group_by(subject) %>% 
    mutate(probe = 1:n()) %>% 
    ungroup() 


## -- SUBJECT EXCLUSION -- ##
# Up to this point, we're already dropping 6 subjects whose probe count in 
# behavioral data mismatch with the eeg data (for reasons unknown)


# And we're already down two subjects going into 07-averaging.py for reasons
# unknown

# 1. Drop a probe if it has fewer than 5 heartbeats
# 2. Drop a subject if they have fewer than 15 probes per condition

# Faceted histograms by subject
d %>% 
    ggplot(aes(x = att)) +
    geom_histogram(color='black') + 
    facet_wrap(~subject)


## -- FOR PROBE AVERAGED DATA -- ##

# Find all probes for subjects with fewer than five epochs
bad_probes <- eeg %>% 
    group_by(subject, probe) %>% 
    summarize(epochs = unique(epochs_per_probe)) %>% 
    filter(epochs < 5) %>% 
    mutate(probe = as.integer(gsub('Probe', '', probe))) %>% 
    select(subject, probe)

# anti_join is a new one
d <- d %>% 
    anti_join(bad_probes) 

N <- length(unique(d$subject))

bads_count <- d %>% 
    select(subject, run, probe, att) %>% 
    filter(!is.na(att)) %>% 
    group_by(subject) %>% 
    mutate(att_med = median(att, na.rm=TRUE)) %>% 
    ungroup() %>% 
    mutate(condition = ifelse(att < att_med, 'Physical', 'Mental')) %>% 
    group_by(subject, condition) %>% 
    summarize(count = n()) %>% 
    spread(condition, count) %>% 
    mutate_all(~ replace_na(., 0)) %>% 
    mutate(flag = ifelse(Physical < 15 | Mental < 15, as.character(subject), '')) 

# Visualize probes per subject per condition
bads_count %>% 
    gather(condition, count, Mental:Physical) %>% 
    ggplot(aes(x = condition, y = count)) + 
    geom_hline(yintercept = 15, linetype = 'dotted') + 
    geom_boxplot(fill = NA, color = 'steelblue', outlier.shape = NA) + 
    geom_jitter(width = .02, alpha = .7) + 
    ggrepel::geom_label_repel(aes(label = flag)) + 
    geom_line(aes(group = subject), linetype = 'dashed', alpha = .5) + 
    labs(
        x = 'Condition (Median Split)',
        y = 'Total Probes',
        caption = paste0('N = ', N, '\nTrimmed ', nrow(bad_probes), ' probes that had fewer than 5 epochs.')
    ) + 
    theme_bw()

ggsave('analysis/scripts/preprocessing/probes_per_subject.png', height = 1080, width = 1920, units = 'px', dpi=120)

# Bads during preprocessing
bads_count <- as.numeric(bads_count$flag[bads_count$flag != ''])
bads_pre <- c(61, 10, 34, 38, 3, 8, 19, 41, 48)
bads_total <- intersect(bads_count, bads_pre)
print('Total bad subjects')
print(bads_total)

# Dropping 10 (from intersection) and 28 (basically no physical probes) for now
bads_total <- c(bads_total, 28)

write.csv(bads_total, 'analysis/scripts/preprocessing/final_bad_subjects.csv', row.names=FALSE, col.names=FALSE)









## -- WE MIGHT NOT NEED TO FIRST AVERAGE BY PROBES ?? -- ##
N <- length(unique(eeg$subject))

# Visualize epochs per subject per condition
#p <- eeg %>% 
eeg %>% 
    group_by(subject, condition) %>% 
    summarize(n_epochs = unique(n_epochs)) %>% 
    spread(condition, n_epochs) %>% 
    ungroup() %>% 
    mutate(label = ifelse(Mental < mean(Mental) - 2 * sd(Mental) | 
                              Physical < mean(Physical) - 2 * sd(Physical),
                          paste0('sub-', sprintf('%03d', subject)), '')) %>% 
    gather(condition, n_epochs, Mental, Physical) %>% 
    ggplot(aes(x = condition, y = n_epochs)) + 
    geom_boxplot(fill = NA, color = 'steelblue', outlier.shape = NA) + 
    geom_jitter(width = .02, alpha = .7) + 
    ggrepel::geom_label_repel(aes(label=label)) + 
    geom_line(aes(group = subject), linetype = 'dashed', alpha = .5) + 
    labs(
        x = 'Condition (Median split)',
        y = 'Total Epochs',
        caption = paste0('N = ', N)
    ) + 
    ylim(0, 550) +
    theme_bw()

ggsave('analysis/scripts/preprocessing/epochs_per_subject.png', height = 1080, width = 1920, units = 'px', dpi=120)

# Dropping 10 bc bad preprocessing, dropping 28 bc outlier on observations
# also need to drop 3, 13, 48 bc of misalignment

bads <- c(3, 10, 13, 28, 48)

write.csv(bads, 'analysis/scripts/preprocessing/final_bad_subjects.csv', row.names=FALSE, col.names=FALSE)
