for myvar in zero oh one two three four five six seven eight nine
do

HInit -S lists/trainList.txt -l ${myvar} -L labels/train -M hmms -o ${myvar} -T 1 lib/proto4States.txt

done
