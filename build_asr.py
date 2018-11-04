import os
import shutil

# config
DIRS = {
    'ASR': 'ASR',
    'MFCC_TRAIN': 'ASR/MFCCs/train',
    'MFCC_TEST': 'ASR/MFCCs/test',
    'LAB_TRAIN': 'ASR/labels/train',
    'LAB_TEST': 'ASR/labels/test',
    'RESULTS': 'ASR/results',
    'LISTS': 'ASR/lists',
    'HMMS': 'ASR/hmms',
    'LIB': 'ASR/lib'
}

TRAIN_DATA = 'data/train/recordings'
TEST_DATA = 'data/test/recordings'

# build new directory tree
if os.path.isdir('ASR'):
    shutil.rmtree('ASR')

for dir in DIRS.values():
    os.makedirs(dir)

# copy training files
for root, dirs, files in os.walk(TRAIN_DATA):
    for file in files:
        fpath = '/' + file
        if file.endswith('.mfcc'):
            shutil.copyfile(root + fpath, DIRS['MFCC_TRAIN'] + fpath)
        if file.endswith('.lab'):
            shutil.copyfile(root + fpath, DIRS['LAB_TRAIN'] + fpath)

# copy testing files
for root, dirs, files in os.walk(TEST_DATA):
    for file in files:
        fpath = '/' + file
        if file.endswith('.mfcc'):
            shutil.copyfile(root + fpath, DIRS['MFCC_TEST'] + fpath)
        if file.endswith('.lab'):
            shutil.copyfile(root + fpath, DIRS['LAB_TEST'] + fpath)