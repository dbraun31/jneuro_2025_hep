import mne
from mne.preprocessing import ICA
from mne_icalabel.iclabel import iclabel_label_components
import os
from pyprojroot import here
import sys
sys.path.append('../helpers')
from parsers import (get_subject_numbers, get_overwrite)
from pathlib import Path
import pickle
import pandas as pd
from glob import glob

if __name__ == '__main__':
    os.chdir(here())
    stem = 'analysis/data/derivatives/00-concatenated'
    subjects = get_subject_numbers(stem)
    overwrite = get_overwrite()

    # Load bad channels
    # Set date to match the bad channels file if there's more than one .pkl
    date = ''

    if not date:
        path = 'analysis/scripts/preprocessing/summaries'
        path = glob(path + '/*.pkl')
        if len(path) > 1:
            raise ValueError('Ambiguous .pkl files. Set a date in the script')
        path = path[0]
    else:
        path = 'analysis/scripts/preprocessing/summaries/{}_bad_channels.pkl'.format(date)

    with open(path, 'rb') as file:
        bad_channels = pickle.load(file)

    for subject in subjects:
        print('Subject {} of {}'.format(subject, subjects[-1]))
        read_path = Path(os.path.join(stem, 'sub-{}'.format(subject)))
        read_file = Path('sub-{}_task-ExperienceSampling_eeg.fif'.format(subject))
        out_stem = Path(stem).parent
        out_path = Path(out_stem) / Path('02-ICA_result') / Path('sub-'+subject)
        raw = mne.io.read_raw_fif(read_path / read_file, preload=True)
        out_file = Path('sub-{}_task-ExperienceSampling_ica.fif'.format(subject))
        out_label = Path('sub-{}_task-ExperienceSampling_label.csv'.format(subject))

        if overwrite or not os.path.exists(out_path / out_file):

            to_interpolate = bad_channels['sub-{}'.format(subject)]
            if to_interpolate:
                raw.info['bads'] = to_interpolate
                raw.interpolate_bads()

            ica = ICA(n_components = 20, method='infomax', random_state = 234)

            ica.fit(raw, picks=['eeg'])

            if not os.path.exists(out_path):
                os.makedirs(out_path)

            ica.save(out_path / out_file, overwrite=overwrite)

            # Label components 
            columns = ['brain', 'muscle_artifact', 'eye_blink', 
            'heart_beat', 'line_noise', 'channel_noise', 'other']
            label_fit = iclabel_label_components(raw, ica)
            label_fit = pd.DataFrame(label_fit, columns = columns)
            label_fit.to_csv(out_path / out_label, index=False)
