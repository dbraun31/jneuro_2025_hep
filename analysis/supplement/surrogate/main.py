import os
from pyprojroot import here
import mne
os.chdir(here())
from core.surrogate import SurrogateTest, save
import pandas as pd
from helpers.parsers import get_subject_numbers
from core.cluster_functions import cluster_test, sum_t
from pathlib import Path
import sys


if __name__ == '__main__':

    # Arg parse
    args = sys.argv[1:]
    if len(args) > 1 or (args and args[0] != 'full'):
        raise ValueError('Only one arg (full) allowed')

    # Whether to get the null dist
    full = False
    if args:
        full = True

    # Necessary args
    data_root = Path('data/projected')
    subjects = get_subject_numbers(data_root)

    # Init
    surtest = SurrogateTest(subjects, data_root=data_root)

    # Run and save alt
    alt = surtest.get_alt()
    sum_t(alt)
    save(alt, 'result/alt_result.pkl')
    with open('result/alt_sum_t.txt', 'w') as file:
        file.write(f'{sum_t(alt)}\n')
    file.close()

    # Run and save full null simulation
    if full:
        null = surtest.get_null()
        save(null, 'result/null_dist.pkl')

    
