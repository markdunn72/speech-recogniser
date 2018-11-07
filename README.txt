TO EXTRACT MFCCS
extractall('data/train/recordings');
extractall('data/test/recordings');

TO BUILD ASR LIBRARY
python3 build_asr.py

TO HINIT TRAINING MODELS
cd asr
sh hinit_all.sh
cd ..

TO RUN SPEECH RECOGNITION ON TEST DATA
HVite -T 1 -S lists/testList.txt -d hmms/ -w lib/NET -l results lib/dict lib/words3

TO RUN SPEECH RECOGNITION ON USER DATA
HVite -T 1 -d hmms/ -w lib/NET lib/dict lib/words3 MFCCs/test/demo.mfcc

TO LIST RESULTS
HResults -p -e "???" sil -e "???" sp -L labels/test lib/words3 results/*.rec
