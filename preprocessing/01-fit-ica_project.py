import pandas as pd
import mne
from mne.preprocessing import ICA
from mne_icalabel import label_components
from mne_bids import BIDSPath
from mne_icalabel.iclabel import iclabel_label_components
from concurrent.futures import ThreadPoolExecutor, as_completed
from glob import glob
from pyprojroot import here
from pathlib import Path
import os
import sys
sys.path.append('../helpers')
from tqdm import tqdm
from parsers import get_subject_numbers


def exclude_component(labels, nonbrain_threshold=.6, brain_threshold=.8):
    # Takes in labels as df
    # Thresholds correspond to classification tolerance
    # Returns a vector of length components with 1 as exclude and 0 as keep

    nonbrain = labels[labels.columns[2:]] < nonbrain_threshold
    nonbrain = nonbrain.all(axis=1).astype(int).values
    brain = (labels['brain'] > brain_threshold).astype(int).values

    return 1 - (brain & nonbrain)


def process_subject(subject, path, overwrite=False):
    filename = 'sub-{}_task-ExperienceSampling_eeg.fif'.format(subject) 
    subject_file = '/'.join([path, 'sub-{}'.format(subject), filename])
    raw = mne.io.read_raw_fif(subject_file, preload=True)

    # Make out file names
    out_path = Path(path).parent / Path('01-projected')
    projected_filename = BIDSPath(subject=subject,
                              task='ExperienceSampling',
                              suffix='projected',
                              extension='.fif',
                              root=out_path,
                              check=False)
    ica_filename = projected_filename.copy()
    ica_filename.suffix = 'ica'
    full_out_path = projected_filename.fpath

    if overwrite or not os.path.exists(full_out_path):
        # Fit ICA
        ica = ICA(n_components=20, max_iter='auto', method='infomax',
                  random_state=42)
        ica.fit(raw, picks = ['eeg'])

        # Label components
        columns = ['brain', 'muscle_artifact', 'eye_blink', 'heart_beat',
        'line_noise', 'channel_noise', 'other']
        label_fit = iclabel_label_components(raw, ica)
        label_fit = pd.DataFrame(label_fit, columns = columns)

        # Apply automatic labeling exclusions
        exclude_bin = exclude_component(label_fit)
        ica.exclude = [i for i, e in enumerate(exclude_bin) if e]

        # Apply to data and save
        ica.apply(raw)

        if not os.path.exists(projected_filename.fpath.parent):
            os.makedirs(projected_filename.fpath.parent)
        
        raw.save(projected_filename.fpath, overwrite=overwrite)
        ica.save(ica_filename.fpath, overwrite=overwrite)

def parallel_process(subjects, path, overwrite=False):
    with tqdm(total=len(subjects), desc='Processing') as pbar:
        with ThreadPoolExecutor() as executor:
            futures = {executor.submit(process_subject, subject, path, overwrite): subject for subject in subjects}
            for future in as_completed(futures):
                subject = future.result()
                pbar.update(1)

def serial_process(subjects, path, overwrite):
    for subject in subjects:
        process_subject(subject, path, overwrite)

if __name__ == '__main__':
    # ICA already uses ~ 50% CPU
    parallel = False

    os.chdir(here())
    response = ''
    while response not in ['y','n']:
        response = input('Do you want to overwrite existing data? (y/n) ')
        response = response.strip().lower()
        
    overwrite = False
    if response == 'y':
        overwrite = True

    path = 'analysis/data/derivatives/00-concatenated'
    if not os.path.exists(path):
        print('Data must already be concatenated and be located in\
         analysis/data/derivatives/concatenated/')
        sys.exit(1)

    subjects = get_subject_numbers(path)

    if parallel:
        parallel_process(subjects, path, overwrite=overwrite)
    else:
        serial_process(subjects, path, overwrite)
