import mne
import sys
from datetime import datetime
import pickle
import os
from pyprojroot import here
from glob import glob
from pathlib import Path


sys.path.append('../helpers')
from parsers import (get_overwrite, get_subject_numbers)


if __name__ == '__main__':
    os.chdir(here())

    overwrite = get_overwrite()

    data_path = Path('analysis/data/derivatives/01-projected')
    subjects = sorted(get_subject_numbers(str(data_path)))
    date = datetime.today().strftime('%Y-%m-%d')
    out_file = Path('analysis/scripts/preprocessing/summaries/{}_bad_channels.pkl'.format(date))

    # Load output data if present
    if os.path.exists(out_file):
        with open(out_file, 'rb') as file:
            bads_dict = pickle.load(file)
    else:
        bads_dict = {}
    
    # Loop over subjects
    for subject in subjects:
        
        # Skip data if it exists
        if 'sub-{}'.format(subject) in bads_dict and not overwrite:
            continue

        print('Subject: {} of {}\n'.format(subject, subjects[-1]))
        subject_raw = 'sub-{}/sub-{}_task-ExperienceSampling_projected.fif'.format(subject, subject)
        input_path = data_path / Path(subject_raw)

        raw = mne.io.read_raw_fif(input_path)

        raw.plot(n_channels=31, duration=110)
        response = ''
        while response not in ['y', 'n']:
            print('Press "q" at any time to quit the program.\n')
            print('Bad channels:')
            print(raw.info['bads'])
            response = input('\nAre you ready to move on to the next subject?  (y/n) ')
            response = response.strip().lower()

            if response == 'n':
                response = ''
                response = input('\nDo you want to open the visualization again? (y/n) ')
                response = response.strip().lower()
                if response == 'y':
                    raw.plot()
                response = ''
            if response == 'q':
              sys.exit(1)


        bads_dict['sub-{}'.format(subject)] = raw.info['bads']

        if not os.path.exists(out_file.parent):
            os.mkdir(out_file.parent)
            print('der')

        with open(out_file, 'wb') as file:
            pickle.dump(bads_dict, file)


