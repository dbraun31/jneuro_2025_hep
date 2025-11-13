import mne 
from mne.time_frequency import tfr_morlet
import pickle
from mne_bids import BIDSPath
import sys
sys.path.append('../../../helpers')
from parsers import get_subject_numbers, get_overwrite
from pathlib import Path
from pyprojroot import here
import pandas as pd
import numpy as np
import os

os.chdir(here())


def get_conditions(item, low_anchor, high_anchor, sbehav):
    '''
    Takes in behavioral data from one subject
    Returns array of condition order  

    '''
    codings = np.array(sbehav[item] > np.nanmedian(sbehav[item]))
    conditions = np.where(codings, high_anchor, low_anchor)

    return conditions


def get_psd(low_anchor, high_anchor, epochs, conditions, ROI,
            by_condition=False):
    '''
    Takes in epochs and array of conditions
    Returns a dict {'condition': (3, 7, 2501)}
        array is (trials, channels, freqs (alpha), time points)
    '''

    # Define ROI
    if ROI == 'occipital':
        channels = ['O1', 'O2', 'Oz']
    elif ROI == 'occipital-parietal':
        channels = ['O1', 'O2', 'Oz', 'P3', 'P4', 'P7', 'P8', 'Pz']
    else:
        raise ValueError('Need correctly formatted ROI')

    freqs = np.logspace(np.log10(2), np.log10(30), num=50)
    n_cycles = freqs / 2
    n_jobs = os.cpu_count() - 1
    channels_idx = np.where(np.isin(epochs.info['ch_names'], channels))[0]
    freqs_idx = np.where((freqs >= 8) & (freqs <= 13))[0]

    if not by_condition:
        result = tfr_morlet(epochs, freqs=freqs, n_cycles=n_cycles,
                         average=False, return_itc=False, n_jobs=n_jobs)
        out = result.data[:, channels_idx, :, :]
        out = out[:, :, freqs_idx, :]

    else:
        epochs_low = epochs[np.where(conditions==low_anchor)[0]]
        epochs_high = epochs[np.where(conditions==high_anchor)[0]]
        epochs_d = {'low': epochs_low, 'high': epochs_high}
        out = {}
        for condition in epochs_d:
            result = tfr_morlet(epochs_d[condition], freqs=freqs, n_cycles=n_cycles,
                                 return_itc=False, average=True, n_jobs=n_jobs)
            freqs_idx = np.where((result.freqs >= 8) & (result.freqs <= 13))[0]
            result_data = result.data[:, freqs_idx, :]
            result_data = result_data[channels_idx, :, :]
            out[condition] = result_data


    return out




def main(in_bids_root, script_root, ROI='occipital'):
    '''
    ROI should be 'occipital' 
    '''

    # IO
    subjects = get_subject_numbers(in_bids_root)
    behav_path = Path('analysis/data/MW_EEG_behavioral.csv')
    behav = pd.read_csv(behav_path)
    drop_log_path = Path(script_root.parent) / Path('06-inspect-epochs/drop_log.pkl')
    bad_subjects_path = Path(drop_log_path.parent) / Path('bad_subjects.txt')
    with open(drop_log_path, 'rb') as file:
        drop_log = pickle.load(file)
    '''
    -- Hard coding the bads instead --

    with open(bad_subjects_path, 'r') as file:
        bad_subjects = eval(file.read())
    '''
    bad_subjects = ['010', '013', '014']

    d_list = []

    for subject in subjects:

        print(f'\n\n-------SUBJECT {subject}/{subjects[-1]}--------\n\n')

        if subject in bad_subjects:
            continue
        
        # Keep only relevant subject
        sbehav = behav[behav['subject'] == int(subject)].reset_index()
        sbehav = sbehav.iloc[~np.isin(sbehav.index, drop_log[subject]), :]

        # Load epochs
        in_bids = BIDSPath(subject=subject,
                        task='ExperienceSampling',
                        suffix='epochs',
                        extension='.fif',
                        root=in_bids_root,
                        check=False)
        epochs = mne.read_epochs(in_bids.fpath)

        # Get experience sampling conditions
        conditions = get_conditions('arou', 'deactivated', 'activated', sbehav)

        if len(conditions) != len(epochs):
            raise ValueError('Number of conditions from behavioral data is not equal to number of epochs.')
        if len(np.unique(conditions)) != 2:
            raise ValueError('Subject {} does not have two conditions'.format(subject))

        # Get psd for subject and condition
        # THIS WILL BREAK IF RUNNING BY CONDITION!!
        psd = get_psd('deactivated', 'activated', epochs, conditions, ROI)
        obs = pd.DataFrame({
            'subject': np.full(psd.shape[0], int(subject)),
            'trial': np.arange(1, psd.shape[0] + 1),
            'power': psd.mean(axis=(1,2,3))
        })
        d_list.append(obs)
        
    return d_list




if __name__ == '__main__':

    # THIS WILL BREAK IF RUNNING BY CONDITION !!!

    # IO
    in_root = Path('analysis/data/derivatives/spectral/05-epochs-clean')
    out_root = Path(in_root.parent) / Path('06-psd')
    script_root = Path('analysis/scripts/preprocessing/spectral/07-compute-psd')

    # Run main
    ROI = 'occipital-parietal'
    result = main(in_root, script_root, ROI=ROI)

    # Save result
    d = pd.concat(result, ignore_index=True)
    
    if not os.path.exists(out_root):
        os.makedirs(out_root)
    d.to_feather(out_root / Path('occipital_bytrial.feather'))

