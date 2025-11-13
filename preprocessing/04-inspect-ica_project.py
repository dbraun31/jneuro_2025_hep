import pickle
import mne
from mne.preprocessing import ICA
import sys
sys.path.append('../helpers')
from parsers import (get_subject_numbers, get_from_user)
from computers import (exclude_component, get_recent_summary)
import os
from datetime import date
from pathlib import Path
from pyprojroot import here
import pandas as pd
from glob import glob




if __name__ == '__main__':


    # Get overwrite arg from user
    overwrite = get_from_user('Do you want to overwrite existing data?')

    # Whether to manually inspect result
    manual = get_from_user('Do you want to manually inspect the ICA result?')

    # Set up paths
    os.chdir(here())
    data_path = Path('analysis/data/derivatives')
    input_ica_path = Path('02-ICA_result')
    input_raw_path = Path('00-concatenated')
    output_raw_path = Path('03-projected')
    subjects = get_subject_numbers(str(data_path / input_ica_path))
    # Get most recent summary
    bc_dir_files = glob('analysis/scripts/preprocessing/summaries/*_bad_channels.pkl')
    bc_dir = get_recent_summary(bc_dir_files)
    with open(bc_dir, 'rb') as file:
        bad_channels = pickle.load(file)
    ica_drop_log = {}
    ica_drop_log_out = Path(f'analysis/scripts/preprocessing/summaries/ica_drop_log.pkl')
    if os.path.exists(ica_drop_log_out):
        with open(ica_drop_log_out, 'rb') as file:
            ica_drop_log = pickle.load(ica_drop_log_out)

    for subject in subjects:
        print('PROCESSING SUBJECT {}'.format(subject))
        # More paths
        subject_stem = 'sub-{}/sub-{}_task-ExperienceSampling'.format(subject, subject)
        subject_raw = Path(subject_stem + '_eeg.fif')
        subject_ica = Path(subject_stem + '_ica.fif')
        subject_label = Path(subject_stem + '_label.csv')
        subject_projected = Path(subject_stem + '_projected.fif')
        out_path = data_path / output_raw_path / subject_projected

        # Check overwrite
        if not overwrite and os.path.exists(out_path):
            continue

        # Read everything in
        ica = mne.preprocessing.read_ica(data_path / input_ica_path / subject_ica)
        labels = pd.read_csv(data_path / input_ica_path / subject_label)
        raw = mne.io.read_raw_fif(data_path / input_raw_path / subject_raw,
                                  preload = True)

        raw.info['bads'] = bad_channels[f'sub-{subject}']
        if raw.info['bads']:
            raw.interpolate_bads()
        # Apply automatic labeling exclusions
        exclude_bin = exclude_component(labels)
        ica.exclude = [i for i, e in enumerate(exclude_bin) if e]

        if manual:
            # Manually inspect
            ica.plot_sources(raw)
            ica.plot_components()
            response = ''
            while response != 'y':
                print('Press "q" at any time to exit the program.\n')
                print('Components to be rejected:')
                print(ica.exclude)
                response = input('\nAre you ready to move on to the next subject? (y/n) ')
                response = response.lower().strip()
                if response == 'n':
                    response = input('Do you want to start the visualization again? (y/n) ')
                    response = response.lower().strip()
                    if response == 'y':
                        ica.plot_sources(raw)
                        ica.plot_components()
                        response = ''
                            
                if response == 'q':
                    sys.exit(1)

        # Apply ICA to data and save
        ica.apply(raw)
        out_path = data_path / output_raw_path / subject_projected
        if not os.path.exists(out_path):
            os.makedirs(out_path.parent)
        raw.save(out_path, overwrite=overwrite)

        ica_drop_log[subject] = ica.exclude
        with open(ica_drop_log_out, 'wb') as file:
            pickle.dump(ica_drop_log, file)



