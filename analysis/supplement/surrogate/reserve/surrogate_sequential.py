import mne
from pyprojroot import here
import pandas as pd
import pickle
from glob import glob
import numpy as np
import os
os.chdir(here())
from pathlib import Path
from helpers.parsers import get_subject_numbers
from core.epoch_functions import ParseEpochs
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
            print(f'\n\n-------SIMULATION {sim}-------\n\n')
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

        # Make aggregate data
        for subject in subjects:
            print(f'\n\n--------SIM: {self.sim}/{self.nsims}; SUBJECT: {subject}--------\n\n')
            subject_d = self._parse_subject(subject, pseudo=pseudo)
            if subject_d is not None:
                d = pd.concat([d, subject_d], axis=0)

        return d

    def _parse_subject(self, subject, pseudo):
        # Subject as '\d\d\d'
        path = glob(str(self.data_root / Path('sub-' + subject) / Path('*.fif')))[0]
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

def import_d(data_root):
    path = Path(data_root) / Path('sub-001/sub-001_task-ExperienceSampling_projected.fif')
    raw = mne.io.read_raw_fif(path)
    channels = raw.info['ch_names']
    channels = [x for x in channels if x != 'ECG']


    d = pd.read_feather('data/hep_arou.feather')
    d.drop(['probe', 'epochs_per_probe'], axis=1, inplace=True)
    d_long = pd.melt(d,
                     id_vars=['subject', 'condition', 'time'],
                     value_vars=channels,
                     var_name='channel',
                     value_name='voltage')
    d_long = d_long.groupby(['subject', 'condition', 'time',
              'channel']).mean().reset_index()
    d = d_long.pivot(index=['subject', 'condition', 'time'],
                 columns='channel', values='voltage').reset_index()
    d.columns.name = None


    d = d[['subject', 'condition', 'time'] + channels]

    d = d[d['subject'] != '010']

    return d


