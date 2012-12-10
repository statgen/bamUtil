ERROR=false

../bin/bam mergeBam --out results/mergeBam.bam --list testFiles/mergeBam.list
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/mergeBam.bam expected/mergeBam.bam
if [ $? -ne 0 ]
then
    ERROR=true
fi

../bin/bam mergeBam -o results/mergeSam.sam -l testFiles/mergeSam.list
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/mergeSam.sam expected/mergeSam.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Files with same headers, but differ RGs
../bin/bam mergeBam -o results/mergeSam1.sam -i testFiles/sortedSam.sam -i testFiles/sortedSam1.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/mergeSam1.sam expected/mergeSam1.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi


# Files with same headers, but differ RGs, swap the order
../bin/bam mergeBam -o results/mergeSam2.sam -i testFiles/sortedSam1.sam -i testFiles/sortedSam.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/mergeSam2.sam expected/mergeSam2.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi



if($ERROR == true)
then
  exit 1
fi
