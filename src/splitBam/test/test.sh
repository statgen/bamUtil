ERROR=false

../../bin/splitBam --in testFiles/merged.bam --out results/splitBam 

diff results/splitBam.RG1.bam expected/split.RG1.bam
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/splitBam.RG2.bam expected/split.RG2.bam
if [ $? -ne 0 ]
then
    ERROR=true
fi

../../bin/splitBam -i testFiles/merged.sam -o results/splitSam

diff results/splitSam.RG1.bam expected/split.RG1.bam
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/splitSam.RG2.bam expected/split.RG2.bam
if [ $? -ne 0 ]
then
    ERROR=true
fi

if($ERROR == true)
then
  exit 1
fi
