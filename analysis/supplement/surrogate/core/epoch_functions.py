import mne
import numpy as np
import pandas as pd
import re
from autoreject import AutoReject
from copy import deepcopy
import pickle
from pyprojroot import here
from mne.preprocessing import find_ecg_events
import os
from io import StringIO
import contextlib
from pathlib import Path
os.chdir(here())
from helpers.parsers import get_subject_numbers
from glob import glob


class ParseEpochs:
    def __init__(self, raw,
                 ecg_dict_path=Path('data/ecg_dict'),
                 behav_path = Path('data/MW_EEG_behavioral.csv')):
        self.ecg_dict_path = ecg_dict_path
        self.subject = re.search(r'(sub-\d+)', str(raw)).group(1)
        self.ecg_samples = None
        self.ecg_d = None
        self.raw = raw
        self.behav = pd.read_csv(behav_path)

    def get_labeled_ecg_samples(self):

        events, event_id = mne.events_from_annotations(self.raw)

        # Get probe onset samples
        q1_samples = events[events[:, 2] == event_id['Q1_on']][:, 0]

        # Get preprobe bins
        sfreq = self.raw.info['sfreq']
        preprobe_bins = np.array([np.array((x - (10 * sfreq), x), dtype=int) for x in q1_samples])
        self.preprobe_bins = preprobe_bins

        # Get true R peaks (read from file first)
        self._get_ecg_events()

        # If subject's ecg samples arent already in the dict
        # Compute and update the written object
        ecg_samples = self.ecg_samples
        if ecg_samples is None:
            ecg_events, _, _ = find_ecg_events(self.raw)
            ecg_samples = ecg_events[:, 0]
            self.ecg_samples = ecg_samples
            self._write_ecg_events()

        # Create probe bins
        bins = np.empty_like(ecg_samples, dtype=object)

        for i, (start, end) in enumerate(preprobe_bins):
            mask = (ecg_samples >= start) & (ecg_samples < end)
            bins[mask] = i

        ecg_mask = [x is not None for x in bins]

        # Obtain set of true R peak samples
        true_rs = ecg_samples[ecg_mask]
        ecg_q1_labels = bins[ecg_mask] + 1

        self.ecg_d = np.column_stack([ecg_q1_labels, true_rs]).astype(int)


    def get_epochs(self, pseudo=True):

        if self.ecg_d is None:
            raise ValueError('Must run get_labeled_ecg_samples first')

        ecg_q1_labels = self.ecg_d[:,0]
        true_rs = self.ecg_d[:, 1]
        # Get number of true heartbeats per probe
        unique, counts = np.unique(ecg_q1_labels, return_counts=True)
        bin_counts = dict(zip(unique, counts))

        # Get pseudo rs k times each

        if pseudo:
            rs = np.array([], dtype=int)
            for i, (start, end) in enumerate(self.preprobe_bins, start=1):
                k = bin_counts[i]

                keep_drawing = True
                while keep_drawing:
                    samples = np.array(range(start, end+1))
                    draws = np.random.choice(samples, size=k, replace=False)

                    if any(np.isin(draws, true_rs)):
                        keep_drawing = False

                rs = np.concatenate([rs, sorted(draws)], axis=0)
        else:
            rs = deepcopy(true_rs)


        values = set(ecg_q1_labels)
        keys = [f'Probe{x}' for x in values]
        event_id = dict(zip(keys, values))
        events = np.column_stack([rs,
                                  np.zeros_like(rs),
                                  ecg_q1_labels]).astype(int)


        tmin = -.1
        tmax = .65
        epochs_raw = mne.Epochs(self.raw,
                            events=events,
                            event_id=event_id,
                            tmin=tmin,
                            tmax=tmax,
                            baseline=None,
                            preload=True)

        return epochs_raw


    def clean_epochs(self, epochs_raw, pseudo=True):
        # If simulation, run autoreject
        if pseudo:
            ar = AutoReject(n_jobs=os.cpu_count()-1)
            epochs_clean = ar.fit_transform(epochs_raw)
            self.epochs = epochs_clean

        # Otherwise, use manually coded data
        else:
            with open('data/epoch_drop_log.pkl', 'rb') as file:
                epoch_drop_log = pickle.load(file)

            subject_num = self.subject.replace('sub-', '')
            bads = epoch_drop_log[subject_num]
            epochs_raw.drop(bads)
            self.epochs = epochs_raw

        return self.epochs

    def average(self, epochs):
        subject_num = int(self.subject.replace('sub-', ''))
        behav = self.behav[self.behav['subject'] == subject_num]
        is_highs = (behav['arou'] > np.nanmedian(behav['arou'])).to_numpy()
        evokeds = epochs.average(by_event_type=True)

        high = []
        low = []
        
        if len(is_highs) != len(evokeds):
            return None

        for is_high, evoked in zip(is_highs, evokeds):
            if np.isnan(evoked.data).all():
                continue
            elif is_high:
                high.append(evoked)
            else:
                low.append(evoked)

        high_evoked = mne.grand_average(high)
        low_evoked = mne.grand_average(low)

        conditions = {'activated': high_evoked, 'deactivated': low_evoked}
        d = pd.DataFrame()
        channels = high_evoked.info['ch_names']

        for condition in conditions:
            evoked = pd.DataFrame(conditions[condition].get_data().transpose(),
                                  columns = channels)
            evoked.insert(0, 'subject', subject_num)
            evoked.insert(1, 'condition', condition)
            evoked.insert(2, 'time', evokeds[0].times)
            d = pd.concat([d, evoked], axis=0)

        return d



    def _get_ecg_events(self):

        read_path = self.ecg_dict_path / Path(self.subject + '_ecg.pkl')
        if os.path.exists(read_path):
            with open(read_path, 'rb') as file:
                self.ecg_samples = pickle.load(file)
        else:
            self.ecg_samples = None



    def _write_ecg_events(self):

        write_path = self.ecg_dict_path / Path(self.subject + '_ecg.pkl')
        with open(write_path, 'wb') as file:
            pickle.dump(self.ecg_samples, file)

