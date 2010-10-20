../bam filter --ref testFiles/chr1_partial.fa --in testFiles/testFilter.sam --mis .49 --qu 30  > results/filter.sam ; diff results/filter.sam expected/filter.sam
