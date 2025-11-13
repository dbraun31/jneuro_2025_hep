import mne
import sys
sys.path.append('../../helpers')
from parsers import get_subject_numbers
from datetime import datetime
from parsers import get_overwrite
from computers import get_hep_epochs
from mne.preprocessing import find_ecg_events
from autoreject import AutoReject
from autoreject import get_rejection_threshold
import os
from pathlib import Path
from pyprojroot import here
import numpy as np
import random

os.chdir(here())

# DIRS
root = Path('analysis/data/derivatives/')
in_root = Path('03-projected')
out_root = Path('hep/04-epochs-raw')
in_file_template = '_task-ExperienceSampling_projected.fif'
out_file_template_ar = '_task-ExperienceSampling_autoreject.h5'
out_file_template_epochs = '_task-ExperienceSampling_epochs.fif'
subjects = get_subject_numbers(str(root / in_root))
hep_drop_summary = {}

overwrite = get_overwrite()

for subject in subjects:

    print('\n\n')
    print('Processing subject {} of {}'.format(subject, sorted(subjects)[-1]))
    print('\n\n')

    sub_string = 'sub-{}'.format(subject.zfill(3))
    in_file = Path(sub_string + in_file_template)
    in_file = root / in_root / Path(sub_string) / in_file
    raw = mne.io.read_raw_fif(in_file)
    out_dir = root / out_root / Path(sub_string)
    out_file_ar = Path(sub_string + out_file_template_ar)
    out_file_epochs = Path(sub_string + out_file_template_epochs)

    if overwrite or not os.path.exists(out_dir / out_file_epochs):

        # Extract epochs around heartbeats nested within probe events
        epochs, prop_dropped = get_hep_epochs(raw)
        # Record number of heartbeats dropped
        hep_drop_summary[subject] = prop_dropped

        # Make it deterministic
        np.random.seed(42)
        random.seed(42)

        ar = AutoReject()
        ar.fit(epochs)
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


        # Save AR and epochs to file
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)


        ar.save(out_dir / out_file_ar, overwrite=overwrite)
        epochs.save(out_dir / out_file_epochs, overwrite=overwrite)


# Write HEP drop summary
date = datetime.today().strftime('%Y-%m-%d')
filename = 'analysis/scripts/preprocessing/summaries/{}_hep_drop_summary.txt'.format(date)
if os.path.exists(filename):
    with open(filename, 'r') as file:
        old_drop_summary = eval(file.read())

    for subject in old_drop_summary:
        if subject not in hep_drop_summary:
            hep_drop_summary[subject] = old_drop_summary[subject]

with open(filename, 'w') as file:
    file.write(str(hep_drop_summary))
file.close()
