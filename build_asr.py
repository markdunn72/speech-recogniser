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

DICT = '''adrian adrian\nali ali\nandrew andrew\nandy andy\nce ce\nchaorong chaorong\n
          jeremy jeremy\nke ke\nliam liam\nmartinho martinho\nmateusz mateusz\n
          minghong minghong\nnicholas nicholas\nnicole nicole\noliver oliver\n
          sarah sarah\nshaun shaun\ntravis travis\nvincent vincent\n
          vinny vinny\nsil sil'''

GRAM = '''
          '''

TRAIN_DATA = 'data/train/recordings'
TEST_DATA = 'data/test/recordings'


# build new directory tree
if os.path.isdir('ASR'):
    shutil.rmtree('ASR')

for dir in DIRS.values():
    os.makedirs(dir)


# copy training files
trainlist = []

for root, dirs, files in os.walk(TRAIN_DATA):
    for file in files:
        fpath = '/' + file
        if file.endswith('.mfcc'):
            shutil.copyfile(root + fpath, DIRS['MFCC_TRAIN'] + fpath)
            trainlist.append(DIRS['MFCC_TRAIN'][3:] + fpath)
        if file.endswith('.lab'):
            shutil.copyfile(root + fpath, DIRS['LAB_TRAIN'] + fpath)

f = open(DIRS['LISTS'] + '/trainList.txt', 'w+')
[f.write(file + '\n') for file in trainlist]
f.close()


# copy testing files
testlist = []

for root, dirs, files in os.walk(TEST_DATA):
    for file in files:
        fpath = '/' + file
        if file.endswith('.mfcc'):
            shutil.copyfile(root + fpath, DIRS['MFCC_TEST'] + fpath)
            testlist.append(DIRS['MFCC_TEST'][3:] + fpath)
        if file.endswith('.lab'):
            shutil.copyfile(root + fpath, DIRS['LAB_TEST'] + fpath)

f = open(DIRS['LISTS'] + '/testList.txt', 'w+')
[f.write(file + '\n') for file in testlist]
f.close()