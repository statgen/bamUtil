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

diff results/mergeSam4.log expected/mergeSam4.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Files with same headers, and same RG ids, but one has PI
../bin/bam mergeBam -o results/mergeSamPI.sam -i testFiles/sortedSam.sam -i testFiles/sortedSamPI.sam --ignorePI --noph
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/mergeSamPI.sam expected/mergeSam3.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Files with same headers, and same RG ids, but different values and one has PI.
../bin/bam mergeBam -o results/mergeSamPI4.sam -i testFiles/sortedSamPI.sam -i testFiles/sortedSam2.sam -I --noph 2> results/mergeSamPI4.log
if [ $? -eq 0 ]
then
    echo "Merge passed when expected to fail."
    ERROR=true
fi

if [ -e results/mergeSamPI4.sam ]
then
    echo "Unexpected results/mergeSam4.sam file."
    ERROR=true
fi

diff results/mergeSamPI4.log expected/mergeSam4.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Files with same headers, and same RG ids, but different PIs
../bin/bam mergeBam -o results/mergeSamPI1.sam -i testFiles/sortedSamPI1.sam -i testFiles/sortedSamPI.sam -I --noph
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/mergeSamPI1.sam expected/mergeSamPI1.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Files with same headers, and same RG ids, but different values before PI.
../bin/bam mergeBam -o results/mergeSamBeforePIfail.sam -i testFiles/sortedSamPI.sam -i testFiles/sortedSamPI2.sam -I --noph 2> results/mergeSamBeforePIfail.log
if [ $? -eq 0 ]
then
    echo "Merge passed when expected to fail."
    ERROR=true
fi

if [ -e results/mergeSamBeforePIfail.sam ]
then
    echo "Unexpected results/mergeSamBeforePIfail.sam file."
    ERROR=true
fi

diff results/mergeSamBeforePIfail.log expected/mergeSam4.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Files with same headers, and same RG ids, but different values after PI.
../bin/bam mergeBam -o results/mergeSamAfterPIfail.sam -i testFiles/sortedSamPI.sam -i testFiles/sortedSamPI3.sam -I --noph 2> results/mergeSamAfterPIfail.log
if [ $? -eq 0 ]
then
    echo "Merge passed when expected to fail."
    ERROR=true
fi

if [ -e results/mergeSamAfterPIfail.sam ]
then
    echo "Unexpected results/mergeSamAfterPIfail.sam file."
    ERROR=true
fi

diff results/mergeSamAfterPIfail.log expected/mergeSam4.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Files with same headers, and same RG ids, but different values before PI.
../bin/bam mergeBam -o results/mergeSamPIfail.sam -i testFiles/sortedSamPI.sam -i testFiles/sortedSamPI1.sam --noph 2> results/mergeSamPIfail.log
if [ $? -eq 0 ]
then
    echo "Merge passed when expected to fail."
    ERROR=true
fi

if [ -e results/mergeSamPIfail.sam ]
then
    echo "Unexpected results/mergeSamPIfail.sam file."
    ERROR=true
fi

diff results/mergeSamPIfail.log expected/mergeSam4.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

if($ERROR == true)
then
  echo "Failed testMergeBam.sh"
  exit 1
fi
