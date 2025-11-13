from pyprojroot import here
import os
import pandas as pd
from glob import glob
import numpy as np
import mne
from pathlib import Path
from mne_bids import (
        BIDSPath,
        read_raw_bids,
        write_raw_bids
)
import sys
sys.path.append('../helpers')
from parsers import get_subject_numbers
    

def get_mapping(raw, tentwenty='analysis/scripts/preprocessing/tentwenty.txt'):
    # A function taking in the raw data 
    # and returning a mapping of 
        # channel name (eg, Fp1) to channel type (eg, EEG)
    # relies on there being a txt file with a list of channel names 
    # in 10-20 system 

    tentwenty_names = eval(open(tentwenty, 'r').read())

    mapping = {}
    for name in raw.info['ch_names']:
        if name in tentwenty_names:
            value = 'eeg'
        elif name == 'ECG':
            value = 'ecg'
        else:
            value = 'misc'
        mapping[name] = value
    return mapping


def get_runs(subject):
    # Takes in a subject number as string
    # Ensures we're only operating on ES data
    # Returns all runs for that subject
    files = glob('analysis/data/rawdata/sub-{}/eeg/*.eeg'.format(subject))
    files = [Path(x) for x in files]
    files = [x.name for x in files]
    files = [x.split('_') for x in files]
    runs = []
    for file in files:
        task = [x.split('-')[-1] for x in file if 'task' in x][0]
        if task == 'ExperienceSampling':
            runs.append([x.split('-')[-1] for x in file if 'run' in x][0])

    return sorted(runs)

def get_out_file(out_path, bids_path):
    # Takes in a directory to output files to
    # Builds BIDS-like file name as Path object
    out_path = out_path / Path('sub-{}'.format(subject))
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    new_bids_path = BIDSPath(subject=bids_path.subject,
                             task=bids_path.task,
                             suffix=bids_path.suffix,
                             extension='.fif',
                             root=out_path)
    out_file = out_path / Path(new_bids_path.basename)
    return out_file

def relabel_events(raw, labels_dir=Path('analysis/scripts/preprocessing/labels/')):
    # Automatically relabeling events 
    # Correct csv label file is inferred from number of events
    # If number of events doesn't match one of three numbers, it gets
    # dropped


    target_event_counts = [561, 421, 281]
    suffixes = ['1-14', '15-18', '19on']
    target_event_dict = {k:v for k, v in zip(target_event_counts, suffixes)}

    # If the data event count isn't captured by one of the csv files
    if len(raw.annotations) not in target_event_counts:
        # Handle the deviant runs
        target_event_counts = [562, 253, 225, 279]
        if len(raw.annotations) not in target_event_counts:
            return None

        # Format probe labels
        rep_length = 4 if len(raw.annotations) == 562 else 2
        numbers = np.repeat(np.arange(2, 14), rep_length)
        QS = ['Q{}'.format(x) for x in numbers]
        QS = ['TO']*rep_length + ['Q1_on'] + ['Q1']*(rep_length-1) + QS
        labels_dict = {
                562: ['boundary']*2 + QS * 10,
                253: ['boundary'] + QS * 9,
                225: ['boundary'] + QS * 8, 
                279: (['boundary'] + QS * 10)[:len(raw.annotations)]}

        labels = labels_dict[len(raw.annotations)]

    else:
        # Find the right labels csv for the data
        subject_range = target_event_dict[len(raw.annotations)]
        csv_filename = Path('event_labels_subjects{}.csv'.format(subject_range))
        csv_labels = pd.read_csv(labels_dir / csv_filename, header=None)
        labels = csv_labels.iloc[0, :].values

    # Get old annotation params as lists
    onsets = [x['onset'] for x in raw.annotations]
    durations = [x['duration'] for x in raw.annotations]
    orig_times = [x['orig_time'] for x in raw.annotations]
    descriptions = list(labels) 

    # This only works when you don't specify orig_time
    new_annotations = mne.Annotations(onset=onsets, duration=durations, 
                                      description=descriptions)
    
    raw.set_annotations(new_annotations)

    return raw



if __name__ in '__main__':
    os.chdir(here())

    # Overwrite?
    ans = ''
    while ans not in ['n', 'y']:
        ans = input('\nDo you want to overwrite existing data? (y/n): ')
        ans = ans.strip().lower()
    overwrite = True if ans == 'y' else False

    # Load montage and get paths
    montage = mne.channels.read_custom_montage('analysis/scripts/Standard-10-20-Cap81.locs')
    path = 'analysis/data/rawdata'
    subjects = get_subject_numbers(path)

    # Store bad relabeling data
    bads = []

    # Subject loop
    for subject in subjects:
        runs = get_runs(subject)

        # Build write dir   
        bids_path = BIDSPath(subject=subject,
                             task='ExperienceSampling',
                             suffix='eeg',
                             extension='.vhdr',
                             root='analysis/data/rawdata/')
        out_path = Path('analysis/data/derivatives/00-concatenated')
        out_file = get_out_file(out_path, bids_path)

    
        if overwrite or not os.path.exists(out_file):

            # Build runs 
            raws = []
            for run in runs:
                print('\nProcessing subject {} run {} of {}'.format(subject, run, runs[-1]))
                # Import raw, set channel types and montage
                bids_path.run = run
                raw = read_raw_bids(bids_path, verbose='ERROR')
                raw.set_channel_types(get_mapping(raw))
                raw.set_montage(montage)

                # Downsample
                raw.resample(250)

                # Filter
                eeg_channels = mne.pick_types(raw.info, eeg=True, ecg=False)
                raw.filter(1, 55, picks=eeg_channels)
                raw.notch_filter(freqs = 50, picks=eeg_channels)
                
                # Re reference
                raw = raw.set_eeg_reference('average', projection=True).apply_proj()

                # Drop misc channels
                #raw.load_data()
                eeg_ecg_channels = mne.pick_types(raw.info, eeg=True, ecg=True)
                raw = raw.pick(eeg_ecg_channels)

                # Relabel events
                # Returns None if bad event count
                number_of_events = len(raw.annotations)
                raw = relabel_events(raw)
                if raw:
                    raws.append(raw)
                else:
                    entry = {'subject':subject, 'run':run,
                             'number_of_events': number_of_events}
                    bads.append(entry)

            if len(raws) > 1:
                raw = mne.concatenate_raws(raws)
            elif len(raws) == 1:
                raw = raws[0]
            else:
                raw = None

            if raw and not os.path.exists(out_file.parent):
                os.makedirs(out_file.parent)

            if raw:
                raw.save(out_file, overwrite=overwrite)

    pd.DataFrame(bads).to_csv('analysis/scripts/preprocessing/bad_relabeling.csv', index=False)






