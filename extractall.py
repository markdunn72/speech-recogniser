import os
import matlab


def extractall():
    # run mfccextract on all wav files in recordings
    # returns 0 if 1 or more wav files found

    eng = matlab.engine.start_matlab()

    count = 0

    for root, dirs, files in os.walk(os.getcwd() + '/data/recordings'):
        for file in files:
            if file.endswith('.wav'):
                print('extracting features from ' + file + '\n...')
                print(eng.mfccextract(file) + 'created')
                count += 1

    print('finished\nfiles extracted: ' + count)

    return 0 if count else 1
