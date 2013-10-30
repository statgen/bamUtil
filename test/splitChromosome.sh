../bin/bam  splitChromosome --in testFilesLibBam/sortedBam.bam --out results/splitSortedBam --noph 2> results/splitChromosome.txt \
&& diff results/splitChromosome.txt expected/splitChromosome.txt \
&& diff results/splitSortedBam1.bam expected/splitSortedBam1.bam \
&& diff results/splitSortedBam2.bam expected/splitSortedBam2.bam \
&& diff results/splitSortedBam3.bam expected/splitSortedBam3.bam \
&& diff results/splitSortedBamUnknownChrom.bam expected/splitSortedBamUnknownChrom.bam \

if [ $? -ne 0 ]
then
  exit 1
fi
