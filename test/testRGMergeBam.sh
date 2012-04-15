ERROR=false

../bin/bam rgMergeBam --out results/mergeBam.bam --list testFiles/mergeBam.list
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/mergeBam.bam expected/mergeBam.bam
if [ $? -ne 0 ]
then
    ERROR=true
fi

../bin/bam rgMergeBam -o results/mergeSam.sam -l testFiles/mergeSam.list
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/mergeSam.sam expected/mergeSam.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi

if($ERROR == true)
then
  exit 1
fi
