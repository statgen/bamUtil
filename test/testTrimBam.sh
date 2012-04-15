ERROR=false

../bin/bam trimBam testFilesLibBam/testSam.sam results/trimSam.sam 2 2> results/testTrim.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/trimSam.sam expected/trimSam.sam && diff results/testTrim.log expected/testTrim.log
if [ $? -ne 0 ]
then
    ERROR=true
fi


if($ERROR == true)
then
  exit 1
fi
