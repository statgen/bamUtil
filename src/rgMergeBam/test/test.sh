ERROR=false

../../bin/rgMergeBam --out results/merged.bam --list testFiles/merge.list

diff results/merged.bam expected/merged.bam
if [ $? -ne 0 ]
then
    ERROR=true
fi

../../bin/rgMergeBam -o results/merged.sam -l testFiles/mergeSam.list

diff results/merged.sam expected/merged.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi

if($ERROR == true)
then
  exit 1
fi
