function runall()

extractall('data/train/recordings');
[~,size] = extractall('data/test/recordings');

cd('data/config');
system(char(strcat('bash hmmprototype.sh', {' '}, num2str(size))));
cd('..');
cd('..');

end