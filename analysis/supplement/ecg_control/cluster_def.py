import mne
from mne.channels import find_ch_adjacency
from mne.stats import permutation_cluster_1samp_test, permutation_cluster_test
from mne.stats import combine_adjacency
from pyprojroot import here
from pathlib import Path
import numpy as np
import os
from scipy import stats
import pickle


def paired_sample_test(ecg, initial_alpha=.01, sig_only=False):
    # Cluster test using median split

    # Format data
    X = _format_input(ecg)

    n_permutations = 5000
    df = X[0].shape[0] - 1
    threshold = stats.t.ppf(1 - initial_alpha / 2, df)

    # Run test
    t_obs, clusters, p_values, _ = permutation_cluster_test(
            X=X,
            threshold=threshold,
            n_permutations=n_permutations,
            tail=0,
            stat_fun=_my_ttest,
            n_jobs=os.cpu_count() - 1,
            seed=32432
    )

    return _format_result(t_obs, clusters, p_values, sig_only)



def one_sample_test(X, initial_alpha=.01, sig_only=False):
    # X needs to be (subjects, time)

    n_permutations = 5000
    df = X.shape[0] - 1
    threshold = stats.t.ppf(1 - initial_alpha / 2, df)

    # Run test

    np.random.seed(23423)
    t_obs, clusters, p_values, _ = permutation_cluster_1samp_test(
            X,
            threshold=threshold,
            n_permutations=n_permutations,
            n_jobs=os.cpu_count() - 1,
            tail=0
    )

    return _format_result(t_obs, clusters, p_values, sig_only)


def _format_result(t_obs, clusters, p_values, sig_only=True):
    # To dict format

    if sum(p_values < .05) == 0 and sig_only:
        return None

    if sig_only:
        clusters = [clusters[i] for i in np.where(p_values<.05)[0]]
    out_clusters = {}

    for number, result in enumerate(zip(clusters, p_values), start=1):
        times = result[0][0]
        channels = result[0][1]
        p_value = result[1]

        out_clusters['Cluster' + str(number)] = {'times_idx': times, 'channels_idx':
                                        channels, 'p_value': p_value}

    out = {'t_obs': t_obs, 'clusters': out_clusters}

    return out


def _format_input(d):
    # Reshape data frame to X
    # X is list of two arrays of (subject, time)

    d = d[d['time'] >= .25]\
            .groupby(['subject', 'time', 'condition'])['ecg']\
            .mean()\
            .reset_index()
       

    # Split by condition and reshape
    high = d[d['condition'] == 'Activated'].drop('condition', axis=1)
    low = d[d['condition'] == 'Deactivated'].drop('condition', axis=1)
    high = high.pivot(index='subject', columns='time', values='ecg').reset_index()
    low = low.pivot(index='subject', columns='time', values='ecg').reset_index()

    # Make sure subject labels are identical
    assert(high['subject'].equals(low['subject']))

    # Convert to array
    high = high.drop('subject', axis=1).to_numpy()
    low = low.drop('subject', axis=1).to_numpy()

    X = [low, high]

    return X


def _my_ttest(a, b):
    out = stats.ttest_rel(a, b)
    return out.statistic

