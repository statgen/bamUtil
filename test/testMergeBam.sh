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

diff results/mergeSamPI4.log expected/mergeSam4PI.log
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

diff results/mergeSamBeforePIfail.log expected/mergeSamPI_2.log
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

diff results/mergeSamAfterPIfail.log expected/mergeSamPI_3.log
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

diff results/mergeSamPIfail.log expected/mergeSamPIfail.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

###############
# TEST REGIONS

### Chr 1
../bin/bam mergeBam -o results/mergeBamReg1List.sam -l testFiles/mergeBam.list --noph -r 1
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/mergeBamReg1List.sam expected/mergeBamReg1.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi

../bin/bam mergeBam -o results/mergeBamReg1File.sam -l testFiles/mergeBam.list --noph -R testFiles/region1.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/mergeBamReg1File.sam expected/mergeBamReg1.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi

### Chr 1,3
../bin/bam mergeBam -o results/mergeBamReg13List.sam -l testFiles/mergeBam.list --noph -r 1,3
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/mergeBamReg13List.sam expected/mergeBamReg13.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi

../bin/bam mergeBam -o results/mergeBamReg13File.sam -l testFiles/mergeBam.list --noph -R testFiles/region13.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/mergeBamReg13File.sam expected/mergeBamReg13.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi

### Chr 1,3 with only parts of chr1
../bin/bam mergeBam -o results/mergeBamReg1P3List.sam -l testFiles/mergeBam.list --noph -r 1:75-1011,1:1750-1750,3
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/mergeBamReg1P3List.sam expected/mergeBamReg1P3.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi

../bin/bam mergeBam -o results/mergeBamReg1P3File.sam -l testFiles/mergeBam.list --noph -R testFiles/region1P3.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/mergeBamReg1P3File.sam expected/mergeBamReg1P3.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi

### Only part of Chr 1
../bin/bam mergeBam -o results/mergeBamReg1PList.sam -l testFiles/mergeBam.list --noph -r 1:1014
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/mergeBamReg1PList.sam expected/mergeBamReg1P.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi

../bin/bam mergeBam -o results/mergeBamReg1PFile.sam -l testFiles/mergeBam.list --noph -R testFiles/region1P.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/mergeBamReg1PFile.sam expected/mergeBamReg1P.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi

###############
# TEST REGIONS ERRORS

../bin/bam mergeBam -o results/mergeErr.sam -l testFiles/mergeBam.list --noph -r 1:1r-5 > results/mergeErr1.err
if [ $? -eq 0 ]
then
    ERROR=true
fi
diff results/mergeErr1.err expected/mergeErr1.err
if [ $? -ne 0 ]
then
    ERROR=true
fi

../bin/bam mergeBam -o results/mergeBamReg1File.sam -l testFiles/mergeBam.list --noph -R testFiles/region1err.txt > results/mergeErr2.err
if [ $? -eq 0 ]
then
    ERROR=true
fi
diff results/mergeErr2.err expected/mergeErr2.err
if [ $? -ne 0 ]
then
    ERROR=true
fi


if($ERROR == true)
then
  echo "Failed testMergeBam.sh"
  exit 1
fi
