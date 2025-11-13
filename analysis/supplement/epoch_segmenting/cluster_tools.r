library(tidyverse)
library(arrow)
if (!'eegUtils' %in% installed.packages()) {
    devtools::install_github('craddm/eegUtils')
}
library(eegUtils)
library(reticulate)
library(RColorBrewer)
use_condaenv('eeg', required=TRUE)
library(here)
library(patchwork)
library(fs)
# Import Python helper functions
setwd(here())
reticulate::py_run_string(
'from writing.figures.revisions.supplement.epoch_segmenting.cluster_tools import (
    get_channel_coordinates,
    permutation_cluster_test_num,
    permutation_cluster_test_med,
    get_channels
)
exit'
)

permutation_cluster_test <- function(eeg,
                               time_window_min=0.25, 
                               bads = c(10, 13, 14), 
                               initial_alpha=.01,
                               median_s = TRUE) {
    # --- PARAMETERS --- #
    # ------------------ #
    # eeg: 
    # time_window_min (double): The minimum time point to use for cluster analysis
    # bads (vector): Vector of subject ids known to have bad data
    # return_results (boolean): 
        # returns as a list the results of the permutation based cluster analysis
    
    # --- RETURNS --- #
    # list of PBCA results
    #--------------------------------------------------------------------------#
    
    item <- 'arou'
    low_anchor <- 'deactivated'
    
    eeg$subject <- as.numeric(eeg$subject)
    eeg <- eeg[!eeg$subject %in% bads,]
    
    
    # --- RUN STATISTICS --- #  
    if (median_s) {
        result <- py$permutation_cluster_test_med(eeg=eeg, time_window_min=time_window_min, initial_alpha=initial_alpha)
    } else {
        result <- py$permutation_cluster_test_num(eeg=eeg, time_window_min=time_window_min, initial_alpha=initial_alpha)
    }
    
    sig_clusters <- length(which(result[['p_values']] < .05))
    
    print('')
    print(paste0('FOUND ', sig_clusters, ' SIGNIFICANT CLUSTER(S)'))
    print('')
    
    return(list(result=result, eeg=eeg, time_window_min=time_window_min))
}


# --- HELPER FUNCTIONS --- #

# PLOTTERS #

plot_topo <- function(m, cluster_id, n_breaks=5, nrow=1) {
    # Function takes in the result of permutation_cluster_test
    # and the idx of *one* cluster
    
    # Will make n_breaks within the time range of the cluster and visualize
      
    result <- m[['result']]
    eeg <- m[['eeg']]
    time_window_min <- m[['time_window_min']]
    
    # Plot topographical map
    channels_total <- get_channels(eeg)
    
    # Get only sig times and channels
    # Get raw values
    times_total <- unique(eeg$time)
    start_time_idx <- which(unique(eeg$time) > time_window_min)[1]
    channels_total <- get_channels(eeg)
    
    # Get indices for t_obs
    time_idxs <- unique(result[['clusters']][[cluster_id]][[1]] + 1)
    channel_idxs <- unique(result[['clusters']][[cluster_id]][[2]] + 1)
    channels <- channels_total[channel_idxs]
    times <- times_total[start_time_idx:(length(times_total))][time_idxs]
    t_obs <- as.data.frame(result[['t_obs']][time_idxs, ])
    colnames(t_obs) <- channels_total
    t_obs <- cbind(data.frame(time=times), t_obs)
    
    # Make time quantiles
    time_points <- quantile(times, probs = seq(0, 1, length.out=n_breaks))
    
    plots <- list()
    count <- 0
    for (time_point in time_points) {
        count <- count + 1
        p <- plot_t_obs(t_obs, time_point = time_point)
        plots[[count]] <- p
    }
    
    out <- purrr::reduce(plots, `+`) + guide_area() + plot_layout(guides='collect', nrow=nrow)
    
    
    return(out)
}

plot_t_obs <- function(t_obs, time_point) {
    # Helper function for plot_topo
    
    # Colors
    colors = rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100))
    
    # Find nearest time point in observations to target time point
    nearest_time_point <- t_obs$time[which.min(abs(t_obs$time - time_point))]
    
    # Get max and min t's (95th percentiles)
    ts <- t_obs %>% 
        gather(channel, t, Fp1:POz) %>% 
        pull(t)
    # Within reasonable limits
    t_min <- max(min(ts), -8)
    t_max <- min(max(ts), 8)
    
    # Objects
    channels_total <- colnames(t_obs)[colnames(t_obs) != 'time']
    ch_pos <- py$get_channel_coordinates(channels_total)
        
    # Plot
    p <- t_obs %>% 
        gather(channel, t, -time) %>% 
        filter(time == nearest_time_point) %>% 
        inner_join(ch_pos) %>% 
        rename(electrode=channel, t_statistic=t) %>% 
        mutate(z=50) %>% 
        ggplot(aes(x = x, y = y, z = z)) +
        geom_topo(chan_markers = 'text', aes(fill = t_statistic, label = electrode)) + 
        labs(
            x = '',
            y = '',
            fill = 't',
            caption = paste0('Time: ', sprintf('%.3f', time_point))
        ) + 
        scale_fill_gradientn(limits = c(t_min, t_max), colors = colors,
                             breaks = c(-3, 4), labels = c('-3 (HA>LA)', '4 (LA>HA)')) + 
        theme_bw() + 
        theme(
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            axis.text = element_blank(),
            text = element_text(size = 14),
            legend.position = 'right'
        )  
    return(p)
}

plot_channels <- function(m, time_range=NA, p_thresh=.05) {
    # Plots all significant clusters
    # m is the output from permutation_cluster_test
    
    result <- m[['result']]
    eeg <- m[['eeg']]
    time_window_min <- m[['time_window_min']]
    
    # For each significant cluster
    sig_clusters <- which(result[['p_values']] < p_thresh)
    
    
    for (cluster_id in sig_clusters) {
    
        # Plot the timeseries of each significant channel
        timeseries <- unique(eeg$time)[unique(eeg$time) >= time_window_min]
        channels_total <- unique(eeg$channel)
        time_idxs <- result[['clusters']][[cluster_id]][[1]]
        channel_idxs <- result[['clusters']][[cluster_id]][[2]]
        channels <- unique(channels_total[channel_idxs+1])
        sig_time_window <- range(timeseries[time_idxs+1])
        
        # Set custom time range for plot
        if (!is.na(time_range)) {
            if (time_range == 'sig') {
                time_min <- sig_time_window[1]
                time_max <- sig_time_window[2]
            } else {
                stopifnot(is.vector(time_range) & length(time_range) == 2, time_range[1] < time_range[2])
                time_min <- time_range[1]
                time_max <- time_range[2]
            }
        } else {
            time_min <- min(eeg$time)
            time_max <- max(eeg$time)
        }
        
        hexs <- scales::brewer_pal(palette = 'Set2')(8)[1:2]
        conditions <- unique(eeg$condition)
        colors <- setNames(hexs, conditions)
        p_value <- result[['p_values']][cluster_id]
        
        p <- eeg %>% 
            filter(channel %in% channels) %>% 
            group_by(subject, condition, channel, time) %>% 
            summarize(voltage = mean(voltage)) %>% 
            group_by(condition, channel, time) %>% 
            summarize(voltage = mean(voltage)) %>% 
            ggplot(aes(x = time, y = voltage)) + 
            geom_hline(yintercept = 0, linetype = 'dashed', alpha = .3) + 
            geom_line(aes(color = condition)) + 
            geom_vline(xintercept = sig_time_window[1], linetype = 'dotted', color = 'steelblue') + 
            geom_vline(xintercept = sig_time_window[2], linetype = 'dotted', color = 'steelblue') + 
            facet_wrap(~channel) + 
            scale_color_manual(values = colors) + 
            labs(
                x = 'Time relative to heartbeat (s)',
                y = latex2exp::TeX('$\\mu v$'),
                color = '',
                caption = paste0('Cluster ID: ', cluster_id, '\np = ', p_value)
            ) + 
            xlim(time_min, time_max) + 
            theme_bw() + 
            theme(
                axis.ticks = element_blank(),
                strip.background = element_rect(fill = NA),
                panel.grid = element_blank(),
                legend.position = 'bottom',
                text = element_text(size = 16)
            )
        print(p)
    }
}

# COMPUTERS #
get_channels <- function(eeg) {
    mask <- sapply(colnames(eeg), grepl, pattern='[A-Z]')
    channels <- colnames(eeg)[mask]
    return(channels)
}

get_significant_cluster_ids <- function(result) {
    # Get a vector of significant cluster ids from result
    return(which(result[['p_values']] < .05))
}

drop_condition_imbalance <- function(eeg, low_anchor, thresh=.9) {
    # Compute which subjects gave very imbalanced responses to a thought probe item
    # !! needs to be updated to work on only one item
    
    # --- PARAMETERS --- #
    # d (df): Data frame of EEG summary data
    # thresh (double): Decimal [0, 1] threshold needed to be exceeded if subject
    #                   is to be labeled bad
    
    # Compute high_anchor
    eeg$condition <- tolower(eeg$condition)
    high_anchor <- unique(eeg$condition)[unique(eeg$condition) != low_anchor]
    
    # Extract imbalanced subjects
    bads <- eeg %>% 
        group_by(subject, condition, probe) %>% 
        summarize(count = n()) %>% 
        group_by(subject, condition) %>% 
        summarize(count = n()) %>% 
        spread(condition, count) %>% 
        mutate({{low_anchor}} := replace_na(.data[[low_anchor]], 0),
               {{high_anchor}} := replace_na(.data[[high_anchor]], 0)) %>% 
        mutate(prop = pmax(.data[[low_anchor]]/(.data[[low_anchor]]+.data[[high_anchor]]),
                           .data[[high_anchor]]/(.data[[low_anchor]]+.data[[high_anchor]]))) %>% 
        filter(prop > thresh) %>% 
        pull(subject)
    
    return(bads)
}



summary_m <- function(m, channels_all=NA, threshold=.05) {
    # Get summary of significant times and channels
    
    sig_clusters <- which(m$result$p_values < threshold)
    eeg <- m$eeg
    
    path <- 'analysis/scripts/hep_controls/surrogate/data/projected'
    if (is.na(channels_all)) channels_all <- py$get_channels(path)
    
    times_all <- unique(eeg$time[eeg$time >= m$time_window_min])
    out <- list()
    
    for (cluster_idx in sig_clusters) {
        cluster <- m$result$clusters[[cluster_idx]]
        
        # Get times
        times_idxs <- cluster[[1]] + 1
        times_idxs_unique <- unique(cluster[[1]] + 1)
        times <- times_all[times_idxs]
        times_unique <- times_all[times_idxs_unique]
        
        # Get channels
        channels_idxs <- cluster[[2]] + 1
        channels_idxs_unique <- unique(cluster[[2]] + 1)
        channels <- channels_all[channels_idxs]
        channels_unique <- channels_all[channels_idxs_unique]
        
        # Sub list
        sublist <- list(full = list(times = times, channels = channels),
                        unique = list(times = times_unique, channels = channels_unique))
        
        # Get p value
        p <- m$result$p_values[cluster_idx]
        
        out[[paste0('Cluster ', cluster_idx)]] <- list(cluster_info = sublist, p_value = p)
        
    }
    
    return(out)
    
}



















