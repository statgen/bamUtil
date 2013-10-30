ERROR=false

../bin/bam splitBam --in testFiles/splitBam.bam --out results/splitBam --noph

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

../bin/bam splitBam -i testFiles/splitBam.sam -o results/splitSam --noph

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