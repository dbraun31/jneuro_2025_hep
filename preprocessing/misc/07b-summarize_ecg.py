import mne
import os
from pyprojroot import here
import sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.append('../helpers')
from parsers import get_subject_numbers
os.chdir(here())



def lookup_probe_value(row, event_ids):
    subject = row['subject']
    probe_idx = row['probe_idx']
    return event_ids.get(str(subject).zfill(3), {}).get(probe_idx, 'Unknown')
    



if __name__ == '__main__':


    path = Path('analysis/data/derivatives/hep/05-epochs-clean')
    subjects = get_subject_numbers(str(path))

    event_ids = {}
    d = []

    for subject in subjects:
        file = Path('sub-{}/sub-{}_task-ExperienceSampling_epochs.fif'.format(subject, subject))
        # Inferring R peaks from heartbeat epochs
        # Events are R peaks within probes
        e = mne.read_epochs(path / file)
        event_id_invert = {e.event_id[k]: k for k in e.event_id}

        event_ids[subject] =  event_id_invert

        subject_v = np.full(e.events.shape[0], int(subject))
        new = np.append(subject_v.reshape(-1, 1), e.events[:, [0, 2]], axis=1)

        d.append(new)

    d = pd.DataFrame(np.vstack(d))
    d.columns = ['subject', 'sample', 'probe_idx']

    d['probe'] = d.apply(lookup_probe_value, event_ids=event_ids, axis=1)

    d = d.drop('probe_idx', axis=1)

    d.to_csv('analysis/data/MW_ECG_summary.csv', index=False)

    



