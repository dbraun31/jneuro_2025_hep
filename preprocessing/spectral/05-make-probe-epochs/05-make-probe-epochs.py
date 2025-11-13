import mne
from mne_bids import BIDSPath
import os
import sys
sys.path.append('../../../helpers')
from parsers import get_subject_numbers
from parsers import get_overwrite
from pyprojroot import here
from pathlib import Path
from autoreject import AutoReject
from joblib import Parallel, delayed
from tqdm import tqdm
import contextlib
from io import StringIO

os.chdir(here())


def make_probe_epochs(subject, in_root, out_root, overwrite, window_start=-10):
    '''
    in and out_root are pathlib.Paths
    window_start is a float indicating how far back before probe epochs in
    seconds (should be negative)
    should start (they default extend to 0, ie, probe onset) 
    returns subject-wise epoched data as .fifs
    also returns the autoreject object  
    '''
    if window_start > 0:
        raise ValueError('window_start should be <= 0')

        
    # IO
    in_bids = BIDSPath(subject=subject,
                       task='ExperienceSampling',
                       suffix='projected',
                       extension='.fif',
                       check=False)

    # Save epochs and autoreject result
    out_bids = BIDSPath(subject=subject,
                        task='ExperienceSampling',
                        suffix='epochs',
                        extension='.fif',
                        check=False)

    if overwrite or not os.path.exists(out_root / out_bids.fpath):
        # Import and read events
        raw = mne.io.read_raw(in_root / in_bids.fpath)
        events, event_id = mne.events_from_annotations(raw) 


        # Epoch
        tmax = 0
        epochs = mne.Epochs(raw=raw,
                            tmin=window_start,
                            tmax=tmax,
                            event_id=event_id['Q1_on'],
                            events=events,
                            baseline=None,
                            preload=True)

        # Autoreject
        ar = AutoReject()
        # Suppress console logging
        with contextlib.redirect_stdout(StringIO()), contextlib.redirect_stderr(StringIO()):
            ar.fit(epochs)


        if not os.path.exists(out_root / out_bids.fpath.parent):
            os.makedirs(out_root / out_bids.fpath.parent)

        epochs.save(out_root / out_bids.fpath, overwrite=overwrite)
        out_bids.update(suffix = 'autoreject', extension = '.h5')
        ar.save(out_root / out_bids.fpath, overwrite=overwrite)



def run_parallel(in_root, out_root, overwrite, window_start=-10, n_jobs=None):
    # Run in parallel across all subjects

    # Subject numbers as str(\d\d\d)
    subjects = get_subject_numbers(in_root)

    n_jobs = os.cpu_count() - 1

    Parallel(
        n_jobs=n_jobs
    )(
        delayed(make_probe_epochs)(
            subject, 
            in_root, 
            out_root, 
            overwrite,
            window_start
        ) for subject in tqdm(subjects, desc='Processing subjects')
    )



if __name__ == '__main__':

    in_root = Path('analysis/data/derivatives/03-projected')
    out_root = Path(in_root.parent) / Path('spectral/04-epochs-raw')

    overwrite = get_overwrite()

    run_parallel(in_root, out_root, overwrite)

