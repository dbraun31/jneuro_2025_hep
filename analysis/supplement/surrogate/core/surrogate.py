import mne
from pyprojroot import here
import pandas as pd
import pickle
from joblib import Parallel, delayed
from tqdm import tqdm
from glob import glob
import numpy as np
import os
os.chdir(here())
from pathlib import Path
from helpers.parsers import get_subject_numbers
from core.epoch_functions import ParseEpochs
from io import StringIO
import contextlib
from core.cluster_functions import cluster_test, sum_t
from types import SimpleNamespace

class SurrogateTest:

    def __init__(self, subjects, nsims=100, data_root=None, load_existing=False):
        # Assumes script is in root .here dir
        if data_root is None:
            raise ValueError('Point to the data root')

        if not os.path.exists('result'):
            os.mkdir('result')

        if load_existing:
            with open('result/null_dist_cache.pkl', 'rb') as file:
                self.null_dist = pickle.load(file)
            nsims = nsims - len(self.null_dist)
        else:
            self.null_dist = []

        self.nsims = nsims
        self.sim = None 
        self.subjects = subjects
        self.data_root = data_root

    def get_null(self):

        for sim in range(self.nsims):
            print(f'\n\n-------SIMULATION {sim+1}/{self.nsims}-------\n\n')
            self.sim = sim
            self.null_dist.append(self._run_sim())
            with open('result/null_dist_cache.pkl', 'wb') as file:
                pickle.dump(self.null_dist, file)

        return self.null_dist

    def get_alt(self):
        # Return entire result

        d = self._get_cluster_data(pseudo=False)
        result = cluster_test(d, self.data_root)
        return result

    def _run_sim(self):

        # Get data
        d = self._get_cluster_data()
        # Run test
        result = cluster_test(d, self.data_root)
        # Parse result
        return sum_t(result)


    def _get_cluster_data(self, pseudo=True):

        d = pd.DataFrame()

        # Drop subject 10
        subjects = [x for x in self.subjects if x != '010']

        n_jobs = os.cpu_count() - 1

        # Set up tuple of generator functions
        # Each yield of the generator gives _parse_subject with the updated
        # subject number
        jobs = (
                delayed(self._parse_subject)(subject, pseudo=pseudo)
                for subject in tqdm(subjects, desc='Processing subjects')
        )

        results = Parallel(n_jobs=n_jobs)(jobs)

        return pd.concat(results, axis=0)

    def _parse_subject(self, subject, pseudo):
        # Subject as '\d\d\d'
        path = glob(str(self.data_root / Path('sub-' + subject) / Path('*.fif')))[0]
        print(f'\n\n----SIM: {self.sim}/{self.nsims}; SUBJECT: {subject}----\n\n')
        # Suppress console logging
        with contextlib.redirect_stdout(StringIO()), contextlib.redirect_stderr(StringIO()):
            raw = mne.io.read_raw_fif(path)
            pe = ParseEpochs(raw)
            pe.get_labeled_ecg_samples()
            epochs_raw = pe.get_epochs(pseudo=pseudo)
            epochs = pe.clean_epochs(epochs_raw, pseudo=pseudo)
            subject_d = pe.average(epochs)
        return subject_d



def save(obj, name):
    with open(name, 'wb') as file:
        pickle.dump(obj, file)

