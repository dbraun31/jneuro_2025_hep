import re
import mne
import os
import pickle
import sys
from pathlib import Path
import pandas as pd
import numpy as np
sys.path.append('../../helpers')
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
def summarize_eeg_by_probe(item,
                           low_anchor,
                           high_anchor,
                           in_dir='analysis/data/derivatives/hep/05-epochs-clean',
                           out_dir='analysis/data/derivatives/hep/06-evoked-clean',
                           behav_path='analysis/data/MW_EEG_behavioral.csv'):
    '''
    Summarize the EEG HEP data by responses to one of the thought probes.

    Parameters
    ----------
    item : str
        The name of the item (eg, arou)
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

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Patch for loading PCA data
    if 'pc' in item or 'resid' in item:
        split = behav_path.split('.')
        split[0] = split[0] + '_append'
        behav_path = '.'.join(split)

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
        epochs = mne.read_epochs(in_dir / Path(sub_string) / file)

        d_full = pd.read_csv(behav_path)
        d = d_full[d_full['subject']==int(subject)]
        eeg[sub_string] = {low_anchor: {}, high_anchor: {}}

        # Log epochs per probe
        epochs_per_probe = {}
        probes = epochs.event_id.keys()
        for probe in probes:
            epochs_per_probe[probe] = len(epochs[probe])

        # Average over heartbeat events
        evokeds = epochs.average(by_event_type=True)

        # Reorganize the list of evokeds into a dict with
        # {'ProbeID': evoked_object}
        probes = {}
        for evoked in evokeds:
            probes[evoked.comment] = evoked

        # Median split based on ES data (omitting NAs)
        probe_list = list(probes.keys())
        conditions = []
        for response in d[item]:
            if response > np.nanmedian(d[item]):
                conditions.append(high_anchor.capitalize())
            else:
                conditions.append(low_anchor.capitalize())

        # Create dict {probe_number: (evoked_data, condition)}
        if len(probe_list) != len(conditions):
            bads_log.append((sub_string, 'Probes in eeg: {}, probes in behav: {}'\
                    .format(len(probe_list), len(conditions))))
            continue

        for probe, condition, response in zip(probe_list, conditions, d[item]):
            if probe in probes:
                probes[probe] = (probes[probe], condition, response)

        # If at least one low probe observation in conditions
        if any([x == low_anchor.capitalize() for x in conditions]):
            # List of low probe response evokeds
            low_probes = [probes[x] for x in probes if probes[x][1] == low_anchor.capitalize()]
            # Only keep HEP if there are observations
            low_probes = [x for x in low_probes if x[0].nave]
            low_probes_df = pd.DataFrame()
            for low_probe in low_probes:
                evoked = low_probe[0]
                response = low_probe[2]
                dt = evoked.to_data_frame()
                dt['condition'] = low_anchor.capitalize()
                dt['response'] = response
                dt['probe'] = evoked.comment
                low_probes_df = pd.concat([low_probes_df, dt], ignore_index=True)
                eeg[sub_string][low_anchor][evoked.comment] = evoked
        else:
            low_probes_df = pd.DataFrame()

        # If at least one high probe observation in conditions
        if any([x == high_anchor.capitalize() for x in conditions]):
            high_probes = [probes[x] for x in probes if probes[x][1] == high_anchor.capitalize()]
            # Only keep HEP if there are observations
            high_probes = [x for x in high_probes if x[0].nave]
            high_probes_df = pd.DataFrame()
            for high_probe in high_probes:
                evoked = high_probe[0]
                response = high_probe[2]
                dt = evoked.to_data_frame()
                dt['condition'] = high_anchor.capitalize()
                dt['probe'] = evoked.comment
                dt['response'] = response
                high_probes_df = pd.concat([high_probes_df, dt], ignore_index=True)
                eeg[sub_string][high_anchor][evoked.comment] = evoked
        else:
            high_probes_df = pd.DataFrame()

        eeg_df = pd.concat([low_probes_df, high_probes_df], ignore_index=True)
        eeg_df['subject'] = subject
        eeg_df['epochs_per_probe'] = eeg_df['probe'].map(epochs_per_probe)
        head = ['subject', 'condition', 'response', 'probe', 'epochs_per_probe','time']
        chans = [x for x in eeg_df.columns if x not in head]
        eeg_df = eeg_df[head + chans]

        out_df = pd.concat([out_df, eeg_df], ignore_index=True)
        


    # Averaged over heartbeat within probe within subject
    if not isinstance(out_dir, Path):
        out_dir = Path(out_dir)
    feather.write_feather(out_df, out_dir / Path('hep_{}.feather'.format(item)))

    # Write eeg dict
    # will parse this before analysis
    with open(out_dir / Path('eeg_dict_{}.pkl'.format(item)), 'wb') as file:
        pickle.dump(eeg, file)

    log_file = 'analysis/scripts/preprocessing/summaries/averaging_bads_{}.txt'.format(item) 
    with open(log_file, 'w') as file:
        file.write(str(bads_log))

    file.close()


if __name__ == '__main__':
    args = sys.argv[1:]
    if not args:
        raise ValueError('Needs item, low anchor, and high anchor')

    item, low_anchor, high_anchor = [x.lower() for x in args]

    summarize_eeg_by_probe(item, low_anchor, high_anchor)



