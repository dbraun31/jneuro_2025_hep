import mne
from mne.channels import find_ch_adjacency
from mne.stats import spatio_temporal_cluster_1samp_test
from mne.stats import combine_adjacency
from pyprojroot import here
from pathlib import Path
import numpy as np
import os
from scipy import stats
import pickle

os.chdir(here())

root = Path('analysis/scripts/continuous_arou')

X = np.load(root / Path('fm_arr.npy'))

def run_cluster_test(X, initial_alpha=.01, sig_only=True):

    # Get sample evoked for montage
    evoked_path = Path('analysis/data/derivatives/hep/06-evoked-clean/eeg_dict_arou.pkl')
    with open(evoked_path, 'rb') as file:
        e = pickle.load(file)

    evoked = e['sub-001']['deactivated']['Probe1']
    montage = mne.channels.make_standard_montage('standard_1020')
    evoked = evoked.set_montage(montage)
    adjacency, _ = find_ch_adjacency(evoked.info, 'eeg')

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


