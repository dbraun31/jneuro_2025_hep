import mne
from pyprojroot import here
import os
from pathlib import Path
import numpy as np
import pandas as pd
import sys
sys.path.append('../helpers')
from parsers import get_subject_numbers

os.chdir(here())

'''
Produce HEP evoked data for both early and late in the -10 s trial epoch
window
'''

mne.set_log_level('ERROR')
subjects = get_subject_numbers('analysis/data/derivatives/hep/05-epochs-clean')

script_root = Path('analysis/scripts/epoch_segmenting')
d = pd.DataFrame()


for subject in subjects:

    print(f'\n--PROCESSING SUBJECT {subject}--\n')

    path = Path('analysis/data/derivatives/hep/05-epochs-clean/',
                f'sub-{subject}/sub-{subject}_task-ExperienceSampling_epochs.fif')

    epochs = mne.read_epochs(path)


    trials = np.unique(epochs.events[:,2])
    idxs_first = []
    idxs_last = []

    # Extract indices of early / late epochs
    for trial in trials:
        mask = epochs.events[:, 2] == trial
        trial_idx = np.where(mask)[0]
        half_n = len(trial_idx) // 2
        idxs_first.extend(trial_idx[:half_n])
        idxs_last.extend(trial_idx[half_n:])


    evoked_first = epochs[idxs_first].average(by_event_type=True)
    evoked_last = epochs[idxs_last].average(by_event_type=True)

    channels = evoked_first[0].info['ch_names']

    # Grab those epochs separated by sequence and concat
    for first, last in zip(evoked_first, evoked_last):
        ds = pd.DataFrame(first.data.transpose(), columns=channels)
        ds.insert(0, 'subject', subject)
        ds.insert(1, 'sequence', 'first')
        ds.insert(2, 'probe', first.comment)
        ds.insert(3, 'time', first.times)
        d = pd.concat([d, ds], axis = 0)

        ds = pd.DataFrame(last.data.transpose(), columns=channels)
        ds.insert(0, 'subject', subject)
        ds.insert(1, 'sequence', 'last')
        ds.insert(2, 'probe', last.comment)
        ds.insert(3, 'time', last.times)
        d = pd.concat([d, ds], axis = 0)


d.reset_index().drop('index', axis=1).to_feather(script_root / Path('der_data/epochs_segmented_wide.feather'))
