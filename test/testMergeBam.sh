ERROR=false

../bin/bam mergeBam --out results/mergeBam.bam --list testFiles/mergeBam.list --noph
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/mergeBam.bam expected/mergeBam.bam
if [ $? -ne 0 ]
then
    ERROR=true
fi

../bin/bam mergeBam -o results/mergeSam.sam -l testFiles/mergeSam.list --noph
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
../bin/bam mergeBam -o results/mergeSam1.sam -i testFiles/sortedSam.sam -i testFiles/sortedSam1.sam --noph
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
../bin/bam mergeBam -o results/mergeSam2.sam -i testFiles/sortedSam1.sam -i testFiles/sortedSam.sam --noph
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/mergeSam2.sam expected/mergeSam2.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi


# Files with same headers, and same RGs (same files)
../bin/bam mergeBam -o results/mergeSam3.sam -i testFiles/sortedSam.sam -i testFiles/sortedSam.sam --noph
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/mergeSam3.sam expected/mergeSam3.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Files with same headers, and same RG ids, but different values.
../bin/bam mergeBam -o results/mergeSam4.sam -i testFiles/sortedSam.sam -i testFiles/sortedSam2.sam --noph 2> results/mergeSam4.log
if [ $? -eq 0 ]
then
    echo "Merge passed when expected to fail."
    ERROR=true
fi

if [ -e results/mergeSam4.sam ]
then
    echo "Unexpected results/mergeSam4.sam file."
    ERROR=true
fi

diff results/mergeSam4.log expected/mergeSam3.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

if($ERROR == true)
then
  echo "Failed testMergeBam.sh"
  exit 1
fi
