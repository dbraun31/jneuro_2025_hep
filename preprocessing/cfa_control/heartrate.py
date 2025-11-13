import mne
from mne_bids import BIDSPath
from mne.preprocessing import find_ecg_events
import os
import pandas as pd
import numpy as np
from pathlib import Path
from pyprojroot import here
import sys
sys.path.append('../../helpers')
from parsers import get_subject_numbers



def get_heartrate_from_epoch(epoch):
    
    raw = mne.io.RawArray(epoch, epochs.info)
    rpeaks = find_ecg_events(raw)

    times = rpeaks[0][:, 0] / raw.info['sfreq']
    # IBI is how many seconds for one beat
    ibi = np.diff(times)
    # (1 beat / x sec ) * 60
    bpm = 60 / ibi.mean()

    return bpm



def code_condition(d, subject):

    # Extract array
    arou = d[d['subject'] == int(subject)]['arou'].to_numpy()

    high = 'activated'
    low = 'deactivated'
    med = np.nanmedian(arou)

    conditions = [high if x > med else low for x in arou]

    return conditions
    




if __name__ == '__main__':
    os.chdir(here())

    root_out = Path('analysis/scripts/preprocessing/cfa_control')
    root_in = Path('analysis/data/derivatives/spectral/05-epochs-clean')
    subjects = get_subject_numbers(root_in)

    bids_template = BIDSPath(task='ExperienceSampling',
                             suffix='epochs',
                             extension='.fif',
                             root=root_in,
                             check=False)

    d = pd.read_csv('analysis/data/MW_EEG_behavioral.csv')

    result = []

    for subject in subjects:
            
        # Get path
        bids_template.subject = subject
        path = bids_template.fpath

        # Import
        epochs = mne.read_epochs(path)

        hr = [(i, get_heartrate_from_epoch(e)) for i, e in enumerate(epochs, start=1)]
        conditions = code_condition(d, subject)

        out = [(int(subject), hr[0], hr[1], c) for hr, c in zip(hr, conditions)]

        result += out

        
    result = pd.DataFrame(result, columns = ['subject', 'probe', 'heartrate', 'arou'])
    result.to_csv(root_out / Path('heartrate.csv'), index=False)


