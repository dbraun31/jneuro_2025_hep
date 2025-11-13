import re
import mne
import os
import pickle
import sys
from pathlib import Path
import pandas as pd
import numpy as np
sys.path.append('../helpers')
from pyprojroot import here
from parsers import get_subject_numbers
os.chdir(here())
import matplotlib.pyplot as plt
import pyarrow.feather as feather


def summarize_eeg_by_probe(probe,
                           low_anchor,
                           high_anchor,
                           in_dir='analysis/data/derivatives/05-epochs-clean',
                           out_dir='analysis/data/derivatives/06-evoked-clean',
                           behav_path='analysis/data/MW_EEG_behavioral.csv'):
    '''
    Summarize the EEG HEP data by responses to one of the thought probes.

    Parameters
    ----------
    probe : str
        The name of the probe (eg, arou)
    low_anchor : str
        Name of the low-anchored side of the scale (eg, deactivated)
    high_anchor : str
        Name of the high-anchored side of the scale (eg, activated)
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

    # INITS
    align_log = []
    # Save out one df and one big dict 
    # DF preserves within probe data
    out_df = pd.DataFrame()
    # {subject: {physical: {probe: evoked}, mental: {probe: evoked}}}
    eeg = {}

    for subject in subjects:
        sub_string = 'sub-{}'.format(str(subject).zfill(3))

        file = Path(sub_string + in_file_template)
        epochs = mne.read_epochs(in_dir / Path(sub_string) / file)

        d_full = pd.read_csv(behav_path)
        d = d_full[d_full['subject']==int(subject)]
        eeg[sub_string] = {low_anchor: {}, high_anchor: {}}

        # Median split based on ES data (omitting NAs)
        probe_list = np.unique(epochs.events[:,2])
        conditions = []
        for response in d[probe]:
            if response < np.nanmedian(d[probe]):
                conditions.append(low_anchor.capitalize())
            else:
                conditions.append(high_anchor.capitalize())

        # Log data alignment and don't continue if misaligned
        entry = {'subject': subject, 'eeg_probes': len(probe_list),
                 'behav_probes': len(conditions)}
        align_log.append(entry)
        if len(probe_list) != len(conditions):
            continue

        condition_dict = {i: conditions[i-1] for i in range(1, len(probe_list)+1)}
        new_descriptions = np.array([condition_dict[x] for x in epochs.events[:,2]])
        new_descriptions = [1 if x == low_anchor.capitalize() else 2 for x in new_descriptions]

        new_events = np.column_stack((epochs.events[:,0], epochs.events[:,1],
                                      new_descriptions))
        new_event_id = {low_anchor.capitalize(): 1, high_anchor.capitalize(): 2}

        epochs.events = new_events
        epochs.event_id = new_event_id

        evokeds = epochs.average(by_event_type=True)

        # Process dict
        low_evokeds = [x for x in evokeds if x.comment == low_anchor.capitalize()][0]
        high_evokeds = [x for x in evokeds if x.comment == high_anchor.capitalize()][0]
        eeg[sub_string][low_anchor] = low_evokeds
        eeg[sub_string][high_anchor] = high_evokeds

        # Process data frame
        dfs = []
        for evoked in evokeds:
            df = evoked.to_data_frame()
            df['subject'] = int(subject)
            df['condition'] = evoked.comment
            df['n_epochs'] = evoked.nave
            dfs.append(df)


        eeg_df = pd.concat(dfs, ignore_index=True)
        head = ['subject', 'condition', 'time']
        chans = [x for x in eeg_df.columns if x not in head]
        eeg_df = eeg_df[head + chans]

        out_df = pd.concat([out_df, eeg_df], ignore_index=True)
        


    if not isinstance(out_dir, Path):
        out_dir = Path(out_dir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # Averaged over heartbeat within subject
    feather.write_feather(out_df, out_dir / Path('hep_{}_by-epoch.feather'.format(probe)))

    # Write eeg dict
    # will parse this before analysis
    with open(out_dir / Path('eeg_dict_{}_by-epoch.pkl'.format(probe)), 'wb') as file:
        pickle.dump(eeg, file)

    file_name = Path('averaging_align_log_{}'.format(probe))
    pd.DataFrame(align_log).to_csv(out_dir / file_name, index=False)

    file.close()



if __name__ == '__main__':
    args = sys.argv[1:]

    if not args:
        raise ValueError('Must supply probe, low_anchor, and high_anchor')

    probe, low_anchor, high_anchor = [x.lower() for x in args]

    summarize_eeg_by_probe(probe, low_anchor, high_anchor)




