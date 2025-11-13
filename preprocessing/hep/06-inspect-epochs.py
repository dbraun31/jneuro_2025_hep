import mne
from pathlib import Path
import sys
import os
import pickle
import numpy as np
from autoreject import read_auto_reject
from pyprojroot import here
sys.path.append('../../helpers')
from parsers import (get_overwrite, get_subject_numbers)
from computers import get_hep_epochs

os.chdir(here())

# DIRS
root = Path('analysis/data/derivatives/')
in_root_epochs = Path('hep/04-epochs-raw')
in_root_raw = Path('03-projected')
out_root = Path('hep/05-epochs-clean')
in_file_template_epochs = '_task-ExperienceSampling_epochs.fif'
in_file_template_ar = '_task-ExperienceSampling_autoreject.h5'
in_file_template_raw = '_task-ExperienceSampling_projected.fif'
out_file_template = '_task-ExperienceSampling_epochs.fif'
subjects = get_subject_numbers(str(root / in_root_epochs))
epoch_drop_log = {}
epoch_drop_log_out = Path('analysis/scripts/preprocessing/summaries/epoch_drop_log.pkl')
if os.path.exists(epoch_drop_log_out):
    with open(epoch_drop_log_out, 'rb') as file:
        epoch_drop_log = pickle.load(file)


overwrite = get_overwrite()


for subject in subjects:
        
    print('\n\n')
    print('Processing subject {} of {}'.format(subject, subjects[-1]))
    print('\n\n')

    sub_string = 'sub-{}'.format(subject).zfill(3)

    # epochs dir
    in_file_epochs = Path(sub_string + in_file_template_epochs)
    in_file_epochs = root / in_root_epochs / Path(sub_string) / in_file_epochs
    epochs = mne.read_epochs(in_file_epochs)

    # autoreject dir
    in_file_ar = Path(sub_string + in_file_template_ar)
    in_file_ar = root / in_root_epochs / Path(sub_string) / in_file_ar
    ar = read_auto_reject(in_file_ar)

    # raw dir
    in_file_raw = Path(sub_string + in_file_template_raw)
    in_file_raw = root / in_root_raw / Path(sub_string) / in_file_raw
    raw = mne.io.read_raw_fif(in_file_raw)

    # out dir
    out_dir = root / out_root / Path(sub_string)
    out_file = Path(sub_string + out_file_template)


    if overwrite or not os.path.exists(out_dir / out_file):

        bad_epochs = ar.get_reject_log(epochs).bad_epochs
        black = ['black']*32
        red = ['orange']*32
        colors = [0] * len(bad_epochs)
        ar_bads = []

        for i, e in enumerate(bad_epochs):
            if e:
                ar_bads.append(i)
                colors[i] = red
            else:
                colors[i] = black

        # New file manually inspecting result
        epochs.plot(epoch_colors = colors, n_channels=31, n_epochs=30)

        response = ''
        while response != 'y':
            print('Press "q" at any time to exit the program.\n')
            response = input('\nAre you ready to move on to the next subject? (y/n) ')
            response = response.lower().strip()
            if response == 'n':
                response = input('Do you want to start the visualization again? (y/n) ')
                response = response.lower().strip()
                if response == 'y':
                    epochs.plot(epoch_colors = colors)
                    response = ''
                        
            if response == 'q':
                sys.exit(1)

        # After manual inspection
        manual_bads = [x for x in range(len(bad_epochs)) if x not in epochs.selection]

        # An auto reject bad in manual bads means the user is saying *not* to drop
        # it

        for ar_bad in ar_bads:
            if ar_bad in manual_bads:
                # Remove from ar bads
                ar_bads.remove(ar_bad)
                # Remove from manual bads
                manual_bads.remove(ar_bad)

        # Need to reload to get fresh epochs.selection
        epochs, _ = get_hep_epochs(raw)

        drops = ar_bads + manual_bads
        epochs.drop(drops)
        
        # Save
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        epochs.save(out_dir / out_file, overwrite=overwrite)

        epoch_drop_log[subject] = drops

        
    with open(epoch_drop_log_out, 'wb') as file:
        pickle.dump(epoch_drop_log, file)
