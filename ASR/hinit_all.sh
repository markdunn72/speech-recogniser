for myvar in sil adrian ali andrew andy ce chaorong jeremy ke liam martinho mateusz minghong nicholas nicole oliver sarah shaun travis vincent vinny
do

HInit -S lists/trainList.txt -l ${myvar} -L labels/train -M hmms -o ${myvar} -T 2 lib/hmmprototype.txt

done