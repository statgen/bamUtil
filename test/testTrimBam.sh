ERROR=false

../bin/bam trimBam testFiles/testTrim.sam results/trimSam.sam 2 --noph 2> results/testTrim.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/trimSam.sam expected/trimSam.sam && diff results/testTrim.log expected/testTrim.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

../bin/bam trimBam testFiles/testTrim.sam results/trimSam3.sam 3 --noph 2> results/testTrim3.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/trimSam3.sam expected/trimSam3.sam && diff results/testTrim3.log expected/testTrim3.log
if [ $? -ne 0 ]
then
    ERROR=true
fi


../bin/bam trimBam testFiles/testTrim.sam results/trimSamL1R2.sam -L 1 -R 2 --noph 2> results/testTrimL1R2.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/trimSamL1R2.sam expected/trimSamL1R2.sam && diff results/testTrimL1R2.log expected/testTrimL1R2.log
if [ $? -ne 0 ]
then
    ERROR=true
fi


../bin/bam trimBam testFiles/testTrim.sam results/trimSamL1R2is.sam -L 1 -i -R 2 --noph 2> results/testTrimL1R2is.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/trimSamL1R2is.sam expected/trimSamL1R2is.sam && diff results/testTrimL1R2is.log expected/testTrimL1R2is.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

../bin/bam trimBam testFiles/testTrim.sam results/trimSamL1R2ignoreStrand.sam --ignoreStrand --left 1 --right 2 --noph 2> results/testTrimL1R2ignoreStrand.log
if [ $? -ne 0 ]
then
    ERROright=true
fi

diff results/trimSamL1R2ignoreStrand.sam expected/trimSamL1R2is.sam && diff results/testTrimL1R2ignoreStrand.log expected/testTrimL1R2ignoreStrand.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

##############
# Repeat tests for clipping instead of trimming
../bin/bam trimBam testFiles/testTrim.sam results/clipSam.sam 2 -c --noph 2> results/testClip.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/clipSam.sam expected/clipSam.sam && diff results/testClip.log expected/testClip.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

../bin/bam trimBam testFiles/testTrim.sam results/clipSam3.sam 3 -c --noph 2> results/testClip3.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/clipSam3.sam expected/clipSam3.sam && diff results/testClip3.log expected/testClip3.log
if [ $? -ne 0 ]
then
    ERROR=true
fi


../bin/bam trimBam testFiles/testTrim.sam results/clipSamL1R2.sam -L 1 -R 2 -c --noph 2> results/testClipL1R2.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/clipSamL1R2.sam expected/clipSamL1R2.sam && diff results/testClipL1R2.log expected/testClipL1R2.log
if [ $? -ne 0 ]
then
    ERROR=true
fi


../bin/bam trimBam testFiles/testTrim.sam results/clipSamL1R2is.sam -L 1 -i -R 2 -c --noph 2> results/testClipL1R2is.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/clipSamL1R2is.sam expected/clipSamL1R2is.sam && diff results/testClipL1R2is.log expected/testClipL1R2is.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

../bin/bam trimBam testFiles/testTrim.sam results/clipSamL1R2ignoreStrand.sam --ignoreStrand --left 1 -c --right 2 --noph 2> results/testClipL1R2ignoreStrand.log
if [ $? -ne 0 ]
then
    ERROright=true
fi

diff results/clipSamL1R2ignoreStrand.sam expected/clipSamL1R2is.sam && diff results/testClipL1R2ignoreStrand.log expected/testClipL1R2ignoreStrand.log
if [ $? -ne 0 ]
then
    ERROR=true
fi


if($ERROR == true)
then
  echo "fail testTrimBam.sh"
  exit 1
fi
