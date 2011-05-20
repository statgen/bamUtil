../bin/bam  splitChromosome --in testFilesLibBam/sortedBam.bam --out results/splitSortedBam 2> results/splitChromosome.txt \
&& diff results/splitChromosome.txt expected/splitChromosome.txt \
&& diff results/splitSortedBam_1.bam expected/splitSortedBam_1.bam \
&& diff results/splitSortedBam_2.bam expected/splitSortedBam_2.bam \
&& diff results/splitSortedBam_3.bam expected/splitSortedBam_3.bam \
&& diff results/splitSortedBam_0.bam expected/splitSortedBam_unknownChrom.bam \

if [ $? -ne 0 ]
then
  exit 1
fi
