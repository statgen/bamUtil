ERROR=false

../bin/bam  writeRegion --in testFilesLibBam/sortedBam.bam --out results/regionNeg1.bam --refID -1 2> results/regionNeg1.txt \
&& diff results/regionNeg1.bam expected/regionNeg1.bam && diff results/regionNeg1.txt expected/regionNeg1.txt

if [ $? -ne 0 ]
then
  exit 1
fi

