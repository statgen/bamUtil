ERROR=false

../bin/bam  writeRegion --in testFilesLibBam/sortedBam.bam --out results/regionNeg1.bam --refID -1 2> results/regionNeg1.txt \
&& diff results/regionNeg1.bam expected/regionNeg1.bam && diff results/regionNeg1.txt expected/regionNeg1.txt \
&& \
../bin/bam writeRegion --in testFilesLibBam/sortedBam.bam --out results/regionRead.sam --readName "1:1011:F:255+17M15D20M" 2> results/regionRead.txt \
&& diff results/regionRead.sam expected/regionRead.sam && diff results/regionRead.txt expected/regionRead.txt \
&& \
../bin/bam writeRegion --in testFilesLibBam/sortedBam.bam --out results/regionRead2.sam --readName "1:1011:F:255+17M15D20M" --refName 1 --start 1010 --end 1011 2> results/regionRead2.txt \
&& diff results/regionRead2.sam expected/regionRead2.sam && diff results/regionRead2.txt expected/regionRead2.txt

if [ $? -ne 0 ]
then
  exit 1
fi

