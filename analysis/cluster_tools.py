import mne
from pathlib import Path
import os
import xarray as xr
import pickle
from pyprojroot import here
import numpy as np
import pandas as pd
from mne.channels import find_ch_adjacency
from mne.stats import spatio_temporal_cluster_test
from glob import glob
from mne.stats import spatio_temporal_cluster_1samp_test
from scipy import stats
from collections import OrderedDict

# --- FUNCTIONS THAT GET CALLED FROM R --- #

def get_channel_coordinates(channels):    
    '''
    Gives x, y, z coordinates for EEG channels for plotting
    
    --- PARAMETERS ---
    ------------------
    channels (list of str): 10-20 names of channels to get coordinates for
    '''
    
    # Load montage at assumed location
    file = Path('analysis/scripts/Standard-10-20-Cap81.locs')
    montage = mne.channels.read_custom_montage(here() / file)
    ch_pos = montage.get_positions()['ch_pos']
    ch_pos = OrderedDict((key, value*1000) for key, value in ch_pos.items())
    ch_pos = pd.DataFrame(ch_pos).transpose()
    ch_pos.columns = ['x', 'y', 'z']
    ch_pos.insert(0, 'channel', ch_pos.index)
    ch_pos = ch_pos[ch_pos['channel'].isin(channels)]
    
    return ch_pos

def get_channels(data_root):
    # Data root needs to points to dir with all the subjects
    path = glob(str(Path(data_root) / Path('sub-001/*.fif')))[0]
    raw = mne.io.read_raw_fif(path)
    return raw.info['ch_names']

# STATISTICS #

def permutation_cluster_test_num(eeg, initial_alpha=.01):
    '''
    Conducts a permutation-based clustering analysis across a map of summary
    statistics, performing a one-sample t-test on each element in the input array.

    ** REMOVE BADS BEFORE CALLING THIS FUNCTION **
    
    --- PARAMETERS ---
    ------------------
    eeg: A data frame with subject | time | channel | statistic
        * statistic must be last * !!!
    initial_alpha (float): Alpha for finding initial clusters
    
    --- RETURNS ---
    ---------------
    out (dict) containing results of permutation test:
        t_obs: (N_timepoints x M_channels) matrix with t values as elements
        clusters: list of (array(time_idx, ...), array(channel_idx, ...)) tuples of all found clusters
        p_values: np.array of shape (N_clusters,) where each element is a p value
    '''


    # Make spatial adjacency information
    sample_epochs_p = Path('analysis/data/derivatives/hep/05-epochs-clean/' + 
        'sub-001/sub-001_task-ExperienceSampling_epochs.fif')

    sample_epochs = mne.read_epochs(sample_epochs_p)

    montage = mne.channels.make_standard_montage('standard_1020')
    sample_epochs = sample_epochs.set_montage(montage)
    adjacency, _ = find_ch_adjacency(sample_epochs.info, 'eeg')

    # Format data
    X = _format_num_for_clustering(eeg, sample_epochs)

    n_permutations = 5000
    df = X.shape[0] - 1
    threshold = stats.t.ppf(1 - initial_alpha / 2, df)

    # Run test

    np.random.seed(23423)
    t_obs, clusters, p_values, _ = spatio_temporal_cluster_1samp_test(
            X,
            threshold=threshold,
            n_permutations=n_permutations,
            tail=0,
            adjacency=adjacency
    )

    return {'t_obs': t_obs, 'clusters': clusters, 'p_values': p_values}


def permutation_cluster_test_med(eeg,
                                 initial_alpha=0.01,
                                 path=Path('analysis/data/derivatives/hep/06-evoked-clean')):
                                 
    '''
    Conducts a permutation-based clustering analysis across a *median split* of
    item, analyzing time points in the epoch. 

    ** REMOVE BADS BEFORE CALLING THIS FUNCTION **
    
    --- PARAMETERS ---
    ------------------
    eeg: EEG data as data frame in long
    initial_alpha (float): Alpha for finding initial clusters
    path (pathlib.Path): Path to directory containing data
    
    --- RETURNS ---
    ---------------
    out (dict) containing results of permutation test:
        t_obs: (N_timepoints x M_channels) matrix with t values as elements
        clusters: list of (array(time_idx, ...), array(channel_idx, ...)) tuples of all found clusters
        p_values: np.array of shape (N_clusters,) where each element is a p value
    '''

    item = 'arou'
    low_anchor = 'deactivated'

     
    sample_epochs_p = Path('analysis/data/derivatives/hep/05-epochs-clean/' + 
         'sub-001/sub-001_task-ExperienceSampling_epochs.fif')
    sample_epochs = mne.read_epochs(sample_epochs_p)

    # Get numpy arrays of shape (subjects, time, chans) for each condition
    low, high = _format_med_for_clustering(eeg, sample_epochs)
    
    # Get a sample evoked object for computing distances
    adjacency, _ = find_ch_adjacency(sample_epochs.info, 'eeg')
    
    # Obtain times and channels for output data
    times = sample_epochs.times
    channels = sample_epochs.info['ch_names'][:-1] # Omit ECG
    
    # Format data as list of arrays
    X = [low, high]
        
    # Configure parameters
    df = low.shape[0] - 1
    t_crit = stats.t.ppf(1 - initial_alpha / 2, df)
    tail = 0 
    
    # Run test
    np.random.seed(32423)
    t_obs, clusters, p_values, _ = spatio_temporal_cluster_test(
        X,
        n_permutations=5000,
        threshold=t_crit,
        tail=tail,
        n_jobs=None,
        seed = 1510,
        buffer_size=None,
        adjacency=adjacency,
        stat_fun=_my_t
    )
    
    out = {'t_obs': t_obs, 'clusters': clusters, 'p_values': p_values, 
    'times': times, 'channels': channels}
    
    return out
    
    
# --- GETS CALLED FROM PYTHON ONLY --- #

def _format_num_for_clustering(eeg, sample_epochs):
    '''
    Take in EEG data summarized as dictionary and output in format amenable to 
    permutation-based cluster analysis
    
    --- PARAMETERS ---
    ------------------
    eeg: Data frame with subject | condition | trial | arou | channel | voltage
    sample_epochs: Any epochs object that can give time and channels
    
    --- RETURNS ---
    A numpy array with shape (subjects, timepoints, channels) of summary statistics
    '''

    # Get channels
    channels = sample_epochs.info['ch_names'][:-1] 
    eeg['channel'] = pd.Categorical(eeg['channel'], categories=channels,
                                    ordered=True)
    eeg = eeg.sort_values(['subject', 'time', 'channel'])
    statistic = eeg.columns[-1]

    # Form array
    xarr = (
            eeg
            .set_index(['subject', 'time', 'channel'])
            .to_xarray()
            )[statistic].values

    return xarr


def _format_med_for_clustering(eeg, sample_epochs):
    '''
    Take in EEG data summarized as dictionary and output in format amenable to 
    permutation-based cluster analysis
    
    --- PARAMETERS ---
    ------------------
    eeg: Data frame with subject | condition | trial | survey_n | channel | voltage
    sample_epochs: Any epochs object that can give time and channels
    
    --- RETURNS ---
    Two numpy arrays of shape (subjects, timepoints, channels)
    The first is the array for the low anchor condition
    The second is the array for the high anchor condition
    '''
    
    low_anchor = 'Deactivated'
    high_anchor = 'Activated'

    # Get channel order
    channels = sample_epochs.info['ch_names'][:-1] # Exclude ECG

    # Flag missing values
    if eeg['condition'].isna().any():
        raise ValueError('Clear out NAs from data prior to running cluster test')

    # Average across trials
    eeg = (
            eeg
            .groupby(['subject', 'condition', 'time', 'channel'], as_index=False)
            .agg({'voltage': 'mean'})
          )

    # Preserve channel order
    eeg['channel'] = pd.Categorical(eeg['channel'], categories=channels,
                                    ordered=True)
    eeg = eeg.sort_values(['condition', 'subject', 'time', 'channel'])

    xarr = (
            eeg
            .set_index(['condition', 'subject', 'time', 'channel'])
            .to_xarray()
            )['voltage']

    low = xarr.sel(condition=low_anchor).values
    high = xarr.sel(condition=high_anchor).values

    return low, high
    


def _my_t(a, b):
    # t is positive for first argument being larger
    out = stats.ttest_rel(a, b)
    return out.statistic
