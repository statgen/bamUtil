ERROR=false

../bin/bam trimBam testFiles/testTrim.sam results/trimSam.sam 2 2> results/testTrim.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/trimSam.sam expected/trimSam.sam && diff results/testTrim.log expected/testTrim.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

../bin/bam trimBam testFiles/testTrim.sam results/trimSam3.sam 3 2> results/testTrim3.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/trimSam3.sam expected/trimSam3.sam && diff results/testTrim3.log expected/testTrim3.log
if [ $? -ne 0 ]
then
    ERROR=true
fi


../bin/bam trimBam testFiles/testTrim.sam results/trimSamL1R2.sam -L 1 -R 2 2> results/testTrimL1R2.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/trimSamL1R2.sam expected/trimSamL1R2.sam && diff results/testTrimL1R2.log expected/testTrimL1R2.log
if [ $? -ne 0 ]
then
    ERROR=true
fi


../bin/bam trimBam testFiles/testTrim.sam results/trimSamL1R2rev.sam -r -L 1 -R 2 2> results/testTrimL1R2rev.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/trimSamL1R2rev.sam expected/trimSamL1R2rev.sam && diff results/testTrimL1R2rev.log expected/testTrimL1R2rev.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

../bin/bam trimBam testFiles/testTrim.sam results/trimSamL1R2reverse.sam --reverse --left 1 --right 2 2> results/testTrimL1R2reverse.log
if [ $? -ne 0 ]
then
    ERROright=true
fi

diff results/trimSamL1R2reverse.sam expected/trimSamL1R2rev.sam && diff results/testTrimL1R2reverse.log expected/testTrimL1R2reverse.log
if [ $? -ne 0 ]
then
    ERROR=true
fi


if($ERROR == true)
then
  echo "fail testTrimBam.sh"
  exit 1
fi
