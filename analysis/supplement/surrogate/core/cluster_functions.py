import mne
from mne.stats import spatio_temporal_cluster_test
from mne.channels import find_ch_adjacency
from pathlib import Path
import numpy as np
from scipy import stats




def cluster_test(d, 
                 data_root, 
                 time_min=.25,
                 time_max=.45,
                 initial_alpha=.01, 
                 nsims=1000):
    
    path = data_root / Path('sub-013/sub-013_task-ExperienceSampling_projected.fif')
    raw = mne.io.read_raw_fif(path)

    # Space
    adjacency, _ = find_ch_adjacency(raw.info, 'eeg')

    # Time
    ds = d[(d['time'] >= time_min) & (d['time'] <= time_max)]

    # Shape data
    low, high = _format_for_clustering(ds, raw)
    X = [low, high]

    # Params
    df = low.shape[0] - 1
    t_crit = stats.t.ppf(1 - initial_alpha / 2, df)
    tail = 0

    t_obs, clusters, p_values, _ = spatio_temporal_cluster_test(
            X,
            n_permutations=nsims,
            threshold=t_crit,
            tail=tail,
            n_jobs=None,
            buffer_size=None,
            adjacency=adjacency,
            stat_fun=_my_t
    )

    return (t_obs, clusters, p_values)


def sum_t(result):
    out = np.array([])

    t_obs = result[0]
    clusters = result[1]

    # Return zero if no clusters are found
    if not clusters:
        return 0

    for cluster in clusters:
        times = cluster[0]
        channels = cluster[1]
        ts = t_obs[times, channels]
        out = np.append(out, sum(ts))

    # Return original element matching max of the abs
    idx = np.where(abs(out) == max(abs(out)))[0][0]
    return out[idx]


def _format_for_clustering(ds, raw):
    
    conditions = ['deactivated', 'activated']
    out = {}

    for condition in conditions:
        dc = ds[ds['condition'].str.lower() == condition]
        dc.drop(['condition', 'time'], axis=1, inplace=True)
        channels = [x for x in dc.columns if x != 'subject']

        grouped = dc.groupby('subject')

        subject_ids = []
        data_list = []

        for subj, group in grouped:
            subject_ids.append(subj)
            data_list.append(group[channels].to_numpy())

        out[condition] = np.stack(data_list)

    return out['deactivated'], out['activated']

def _my_t(a, b):
    # t is positive for first argument being larger
    out = stats.ttest_rel(a, b)
    return out.statistic

