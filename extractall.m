function [extracted, size] = extractall(fDir)
% This function runs mfccextract on all wav files
% in directories within fDir
% returns true if any wav file found else false
sDir = pwd;
cd(fDir);
D = dir;
extracted = false;

for k = 3:length(D)
    currD = D(k).name;
    if ~startsWith(currD, '.')
        cd(fDir);
        cd(currD);
        disp('Checking subdirectory:');
        disp(currD);
        fList = (dir('*.wav'));
        for x = 1:length(fList)
            disp('Extracting features from:');
            f = fList(x).name;
            disp(f);
            size = mfccextract(f);
            extracted = true;
        end
        cd('..');    
    end
cd(sDir);
end
