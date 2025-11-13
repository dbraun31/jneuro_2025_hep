import mne
import sys
from pathlib import Path
sys.path.append('../../../helpers/')
from parsers import get_overwrite
from viz import inspect_epochs
from autoreject import read_auto_reject
from mne_bids import BIDSPath
import pickle
import os
from pyprojroot import here
from parsers import get_subject_numbers

os.chdir(here())


def inspect_all_epochs(in_root, out_root, out_log):

    subjects = get_subject_numbers(in_root)

    overwrite = get_overwrite()

    # Import log
    if os.path.exists(out_log):
        with open(out_log, 'rb') as file:
            log = pickle.load(file)
    else:
        log = {}

    for subject in subjects:

        print('\n--------\nPROCESSING SUBJECT {}\n--------\n'.format(subject))
        # IO
        in_bids_epochs = BIDSPath(subject=subject,
                           task='ExperienceSampling',
                           suffix='epochs',
                           extension='.fif',
                           root = in_root,
                           check=False)
        in_bids_raw = in_bids_epochs.copy().update(
                        suffix='projected',
                        root=in_root_raw)

        in_bids_ar = BIDSPath(subject=subject,
                           task='ExperienceSampling',
                           suffix='autoreject',
                           extension='.h5',
                           root=in_root,
                           check=False)

        out_bids = in_bids_epochs.copy().update(root=out_root)

        if overwrite or not os.path.exists(out_bids.fpath):

            epochs = mne.read_epochs(in_bids_epochs.fpath)
            raw = mne.io.read_raw(in_bids_raw.fpath)
            ar = read_auto_reject(in_bids_ar.fpath)

            # If subject is in the drop log, just use it
            if subject in log:
                epochs.drop(log[subject])
            else:
                # Drops are indices of trials dropped
                epochs, drops = inspect_epochs(ar, epochs, raw)
                log[subject] = drops


            if not os.path.exists(out_bids.fpath.parent):
                os.makedirs(out_bids.fpath.parent)

            epochs.save(out_bids.fpath, overwrite=overwrite)



        with open(out_log, 'wb') as file:
            pickle.dump(log, file)


if __name__ == '__main__':

    in_root = Path('analysis/data/derivatives/spectral/04-epochs-raw')
    in_root_raw = Path(in_root.parent.parent) / Path('03-projected')
    out_root = Path(in_root.parent) / Path('05-epochs-clean')
    out_log = Path('analysis/scripts/preprocessing/spectral/06-inspect-epochs/drop_log.pkl')

    inspect_all_epochs(in_root, out_root, out_log)
