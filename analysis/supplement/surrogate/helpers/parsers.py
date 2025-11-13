from pathlib import Path
from glob import glob

def get_subject_numbers(path):
    # Input is path as pathlib.Path or str to subject directories
    # Output list of string subject numbers

    subject_dirs = [Path(x) for x in glob(str(path) + '/*')]
    subject_dirs = [x for x in subject_dirs if 'sub' in x.name]
    subjects = [x.name.split('-')[1] for x in subject_dirs]
    return sorted(subjects)


def get_overwrite():

    response = ''
    while response not in ['y', 'n']:
        response = input('\nDo you want to overwrite existing data? (y/n) ')
        response = response.strip().lower()
    overwrite = False
    if response == 'y':
        overwrite = True

    return overwrite


def get_from_user(prompt):
    # Returns a boolean from user to a prompt (string)

    response = ''
    while response not in ['y', 'n']:
        response = input('\n'+ prompt + ' (y/n): ')
        response = response.strip().lower()

    if response == 'y':
        return True
    return False


