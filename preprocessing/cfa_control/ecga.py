import mne
from mne_bids import BIDSPath
from mne.preprocessing import find_ecg_events
import os
import pandas as pd
import numpy as np
from pathlib import Path
from pyprojroot import here
import sys
sys.path.append('../../helpers')
from parsers import get_subject_numbers


def get_peak_amp(epoch, times, win_min=0.328, win_max=0.364):

    mask = np.array([True if e <= win_max and e >= win_min else False for e in times])
    ecg = epoch[0]
    return ecg[mask].max()

def code_condition(d, subject):

    # Extract array
    arou = d[d['subject'] == int(subject)]['arou'].to_numpy()

    high = 'activated'
    low = 'deactivated'
    med = np.nanmedian(arou)

    conditions = [high if x > med else low for x in arou]

    return conditions




if __name__ == '__main__':
    os.chdir(here())

    root_out = Path('analysis/scripts/preprocessing/cfa_control')
    root_in = Path('analysis/data/derivatives/hep/05-epochs-clean')
    subjects = get_subject_numbers(root_in)

    bids_template = BIDSPath(task='ExperienceSampling',
                             suffix='epochs',
                             extension='.fif',
                             root=root_in,
                             check=False)

    d = pd.read_csv('analysis/data/MW_EEG_behavioral.csv')

    result = []

    for subject in subjects:
            
        # Get path
        bids_template.subject = subject
        path = bids_template.fpath

        # Import
        epochs_full = mne.read_epochs(path)
        epochs = epochs_full.pick(['ECG'])
        times = epochs.times
        

        # Get conditions
        conditions = code_condition(d, subject)
        conditions = {i:e for i, e in enumerate(conditions, start=1)}

        probe_ids = np.fromiter(epochs.event_id.values(), dtype='int')
        condition_ids = np.fromiter(conditions.keys(), dtype='int')

        if not np.array_equal(probe_ids, condition_ids):
            conditions = {i: conditions[i] for i in probe_ids if i in condition_ids}

        # Get amplitude
        for i, epoch in enumerate(epochs):
            probe = epochs.events[i, 2]
            condition = conditions.get(probe)
            if not condition:
                continue
            amp = get_peak_amp(epoch, times)
            trial = (int(subject), probe, amp, condition)
            result.append(trial)
            

    result = pd.DataFrame(result, columns=['subject', 'probe', 'peak_ecg', 'arou'])
    result.to_csv(root_out / Path('ecga.csv'), index = False)
