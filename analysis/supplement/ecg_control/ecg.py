import re
import mne
import os
import pickle
import sys
from natsort import natsorted
from pathlib import Path
import pandas as pd
import numpy as np
sys.path.append('../helpers')
from pyprojroot import here
from parsers import get_subject_numbers
os.chdir(here())
import matplotlib.pyplot as plt
import pyarrow.feather as feather

'''
# Exclude subjects who looked bad during preprocessing
# And who have no variability on mental/physical responses
bads = [3, 6, 7, 9, 10, 15, 21, 34, 38, 39, 45, 53, 55, 58, 59, 61, 64]
bads = [str(x).zfill(3) for x in bads]
subjects = [x for x in subjects if x not in bads]
'''
def summarize_eeg_by_probe(in_dir='analysis/data/derivatives/hep/05-epochs-clean',
                           out_dir='analysis/scripts/ecg_control',
                           behav_path='analysis/data/MW_EEG_behavioral.csv'):
    '''
    Summarize the ECG HEP data by trial

    Parameters
    ----------
    in_dir : str
        Local path to directory containing input data.
    out_dir : str
        Local path to save output data.
    behav_path : str
        Local full path to behavioral data.

    Returns
    -------
    Saves a .feather data to file in BIDS naming format.

    '''

    # DIRS
    if not isinstance(in_dir, Path):
        in_dir = Path(in_dir)
    in_file_template = '_task-ExperienceSampling_epochs.fif'
    subjects = get_subject_numbers(str(in_dir))

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # INITS
    bads_log = []
    # Save out one df and one big dict 
    # DF preserves within probe data
    out_df = pd.DataFrame()
    # {subject: {physical: {probe: evoked}, mental: {probe: evoked}}}
    eeg = {}

    for subject in subjects:
        sub_string = 'sub-{}'.format(str(subject).zfill(3))

        file = Path(sub_string + in_file_template)

        # Shape (568, 32, 188)
        # (Heartbeat events, channels, timepoints per epoch)
        epochs = mne.read_epochs(in_dir / Path(sub_string) / file)
        ecg = epochs.get_data()[:, 31, :]

        # Use pandas for grouped average
        trial_dict = {epochs.event_id[k]: k for k in epochs.event_id}
        trials = epochs.events[:, 2]
        d = pd.DataFrame(ecg)
        d['trial'] = [trial_dict[x] for x in trials]
        ecg_avg = d.groupby('trial').mean().reset_index()
        ecg_avg = pd.melt(ecg_avg, id_vars='trial', var_name='sample', value_name='ecg')
        ecg_avg['trial'] = ecg_avg['trial'].apply(lambda x: int(str(x).replace('Probe', '')))
        ecg_avg.insert(0, 'subject', int(subject))
        ecg_avg = ecg_avg.sort_values(by=['trial', 'sample']).reset_index()
        ecg_avg = ecg_avg.drop('index', axis=1)

        out_df = pd.concat([out_df, ecg_avg])

        


    path = Path(out_dir) / Path('ecg_averaged.csv')
    out_df.to_csv(path, index=False)


if __name__ == '__main__':

    summarize_eeg_by_probe()



