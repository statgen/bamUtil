../bin/bam stats --in testFilesLibBam/testSam.sam 2> results/basicStats.txt \
&& diff results/basicStats.txt expected/basicStats.txt \
&& \
../bin/bam stats --in testFilesLibBam/testSam.sam --qual 2> results/qualStats.txt \
&& diff results/qualStats.txt expected/qualStats.txt \
&& \
../bin/bam stats --in testFilesLibBam/testBam.bam --unmappedQual 2> results/unmappedQualStats.txt \
&& diff results/unmappedQualStats.txt expected/unmappedQualStats.txt \
&& \
../bin/bam stats --disableStatistics --in testFilesLibBam/sortedBam.bam --qual 2> results/sortedQualStats.txt \
&& diff results/sortedQualStats.txt expected/sortedQualStats.txt \
&& \
../bin/bam stats --in testFilesLibBam/sortedBam.bam --unmappedQual 2> results/sortedUnmappedQualStats.txt \
&& diff results/sortedUnmappedQualStats.txt expected/sortedUnmappedQualStats.txt \
&& \
../bin/bam stats --in testFilesLibBam/sortedBam.bam --unmappedQual --refFile testFilesLibBam/sortedBam.bam.bai 2> results/sortedUnmappedQualStats1.txt \
&& diff results/sortedUnmappedQualStats1.txt expected/sortedUnmappedQualStats.txt \
&& \
../bin/bam stats --in testFilesLibBam/testSam.sam --qual --maxNumQuals 3 2> results/qualStatsMaxNum.txt \
&& diff results/qualStatsMaxNum.txt expected/qualStatsMaxNum.txt \
&& \
../bin/bam stats --in testFilesLibBam/testBam.bam --unmappedQual --maxNumQuals 3 2> results/unmappedQualStatsMaxNum.txt \
&& diff results/unmappedQualStatsMaxNum.txt expected/unmappedQualStatsMaxNum.txt \
&& \
../bin/bam stats --in testFilesLibBam/sortedBam.bam --unmappedQual --maxNumQuals 3 2> results/sortedUnmappedQualStatsMaxNum.txt \
&& diff results/sortedUnmappedQualStatsMaxNum.txt expected/sortedUnmappedQualStatsMaxNum.txt \

