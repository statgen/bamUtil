../bin/bam stats --in testFilesLibBam/testSam.sam 2> results/stats.txt \
&& diff results/stats.txt expected/stats.txt \
&& \
../bin/bam stats --basic --in testFilesLibBam/testSam.sam 2> results/basicStats.txt \
&& diff results/basicStats.txt expected/basicStats.txt \
&& \
../bin/bam stats --basic --in testFilesLibBam/testSam.sam --qual 2> results/qualStats.txt \
&& diff results/qualStats.txt expected/qualStats.txt \
&& \
../bin/bam stats --basic --in testFilesLibBam/testSam.sam --phred 2> results/phredStats.txt \
&& diff results/phredStats.txt expected/phredStats.txt \
&& \
../bin/bam stats --basic --in testFilesLibBam/testBam.bam --qual 2> results/qualStats2.txt \
&& diff results/qualStats2.txt expected/qualStats.txt \
&& \
../bin/bam stats --basic --in testFilesLibBam/sortedBam.bam 2> results/sortedStats.txt \
&& diff results/sortedStats.txt expected/sortedStats.txt \
&& \
../bin/bam stats --in testFilesLibBam/sortedBam.bam --qual 2> results/sortedQualStats.txt \
&& diff results/sortedQualStats.txt expected/sortedQualStats.txt \
&& \
../bin/bam stats --basic --in testFilesLibBam/sortedBam.bam --unmapped 2> results/unmappedStats.txt \
&& diff results/unmappedStats.txt expected/unmappedStats.txt \
&& \
../bin/bam stats --basic --in testFilesLibBam/sortedBam.bam --unmapped --qual 2> results/unmappedQualStats.txt \
&& diff results/unmappedQualStats.txt expected/unmappedQualStats.txt \
&& \
../bin/bam stats --basic --in testFilesLibBam/sortedBam.bam --unmapped --qual --bamIndex testFilesLibBam/sortedBam.bam.bai 2> results/unmappedQualStats1.txt \
&& diff results/unmappedQualStats1.txt expected/unmappedQualStats.txt \
&& \
../bin/bam stats --basic --in testFilesLibBam/testSam.sam --qual --maxNumReads 3 2> results/qualStatsMaxNum.txt \
&& diff results/qualStatsMaxNum.txt expected/qualStatsMaxNum.txt \
&& \
../bin/bam stats --in testFilesLibBam/sortedBam.bam --unmapped --qual --maxNumReads 1 2> results/unmappedQualStatsMaxNum.txt \
&& diff results/unmappedQualStatsMaxNum.txt expected/unmappedQualStatsMaxNum.txt \
&& \
../bin/bam stats --in testFiles/testStatsBaseQC.sam --cBaseQC results/statsBaseQCsam.txt 2> results/statsBaseQCsam.log \
&& diff results/statsBaseQCsam.txt expected/statsBaseQC.txt && diff results/statsBaseQCsam.log expected/statsBaseQC.log \
&& \
../bin/bam stats --in testFiles/testStatsBaseQC.sam --pBaseQC results/statsBaseQCsamPercent.txt 2> results/statsBaseQCsamPercent.log \
&& diff results/statsBaseQCsamPercent.txt expected/statsBaseQCPercent.txt && diff results/statsBaseQCsamPercent.log expected/statsBaseQC.log \
&& \
../bin/bam stats --in testFiles/testStatsBaseQC.bam --pBaseQC results/statsBaseQCbamPercent.txt 2> results/statsBaseQCbamPercent.log \
&& diff results/statsBaseQCbamPercent.txt expected/statsBaseQCPercent.txt && diff results/statsBaseQCbamPercent.log expected/statsBaseQC.log \
&& \
../bin/bam stats --in testFiles/testStatsBaseQCSorted.bam --cBaseQC results/statsBaseQCreg.txt --regionList testFiles/region.txt 2> results/statsBaseQCreg.log \
&& diff results/statsBaseQCreg.txt expected/statsBaseQCreg.txt && diff results/statsBaseQCreg.log expected/statsBaseQCreg.log \
&& \
../bin/bam stats --in testFiles/testStatsBaseQCSorted.bam --pBaseQC results/statsBaseQCregPercent.txt --regionList testFiles/region.txt 2> results/statsBaseQCregPercent.log \
&& diff results/statsBaseQCregPercent.txt expected/statsBaseQCregPercent.txt && diff results/statsBaseQCregPercent.log expected/statsBaseQCreg.log \
&& \
../bin/bam stats --in testFiles/testStatsBaseQCSorted.bam --cBaseQC results/statsBaseQCregQual2.txt --regionList testFiles/region.txt --minMapQual 2 2> results/statsBaseQCregQual2.log \
&& diff results/statsBaseQCregQual2.txt expected/statsBaseQCregQual2.txt && diff results/statsBaseQCregQual2.log expected/statsBaseQCreg.log \
&& \
../bin/bam stats --in testFiles/testStatsBaseQCSorted.bam --pBaseQC results/statsBaseQCregQual2Percent.txt --regionList testFiles/region.txt --minMapQual 2 2> results/statsBaseQCregQual2Percent.log \
&& diff results/statsBaseQCregQual2Percent.txt expected/statsBaseQCregQual2Percent.txt && diff results/statsBaseQCregQual2Percent.log expected/statsBaseQCreg.log \
&& \
../bin/bam stats --in testFiles/testStatsBaseQCSorted.bam --cBaseQC results/statsBaseQCregQual20.txt --regionList testFiles/region.txt --minMapQual 20 2> results/statsBaseQCregQual20.log \
&& diff results/statsBaseQCregQual20.txt expected/statsBaseQCregQual20.txt && diff results/statsBaseQCregQual20.log expected/statsBaseQCreg.log \
&& \
../bin/bam stats --in testFiles/testStatsBaseQCSorted.bam --pBaseQC results/statsBaseQCregQual20Percent.txt --regionList testFiles/region.txt --minMapQual 20 2> results/statsBaseQCregQual20Percent.log \
&& diff results/statsBaseQCregQual20Percent.txt expected/statsBaseQCregQual20Percent.txt && diff results/statsBaseQCregQual20Percent.log expected/statsBaseQCreg.log \
&& \
../bin/bam stats --in testFiles/testStatsBaseQCSorted.bam --cBaseQC results/statsBaseQCregQual40.txt --regionList testFiles/region.txt --minMapQual 40 2> results/statsBaseQCregQual40.log \
&& diff results/statsBaseQCregQual40.txt expected/statsBaseQCregQual40.txt && diff results/statsBaseQCregQual40.log expected/statsBaseQCreg.log \
&& \
../bin/bam stats --in testFiles/testStatsBaseQCSorted.bam --pBaseQC results/statsBaseQCregQual40Percent.txt --regionList testFiles/region.txt --minMapQual 40 2> results/statsBaseQCregQual40Percent.log \
&& diff results/statsBaseQCregQual40Percent.txt expected/statsBaseQCregQual40Percent.txt && diff results/statsBaseQCregQual40Percent.log expected/statsBaseQCreg.log \
&& \
../bin/bam stats --in testFiles/testStatsBaseQCSorted.bam --cBaseQC results/statsBaseQCregDB.txt --dbsnp testFiles/dbsnp.txt --regionList testFiles/region.txt 2> results/statsBaseQCregDB.log \
&& diff results/statsBaseQCregDB.txt expected/statsBaseQCregDB.txt && diff results/statsBaseQCregDB.log expected/statsBaseQCreg.log \
&& \
../bin/bam stats --in testFiles/testStatsBaseQCSorted.bam --pBaseQC results/statsBaseQCregDBPercent.txt --dbsnp testFiles/dbsnp.txt --regionList testFiles/region.txt 2> results/statsBaseQCregDBPercent.log \
&& diff results/statsBaseQCregDBPercent.txt expected/statsBaseQCregDBPercent.txt && diff results/statsBaseQCregDBPercent.log expected/statsBaseQCreg.log \
&& \
../bin/bam stats --in testFiles/testStatsBaseQC255.sam --cBaseQC results/statsBaseQCsam255.txt 2> results/statsBaseQCsam255.log \
&& diff results/statsBaseQCsam255.txt expected/statsBaseQC255.txt && diff results/statsBaseQCsam255.log expected/statsBaseQC.log \
&& \
../bin/bam stats --in testFiles/testStatsBaseQC255.sam --pBaseQC results/statsBaseQCsam255Percent.txt 2> results/statsBaseQCsam255Percent.log \
&& diff results/statsBaseQCsam255Percent.txt expected/statsBaseQC255Percent.txt && diff results/statsBaseQCsam255Percent.log expected/statsBaseQC.log \
&& \
../bin/bam stats --in testFiles/testStatsBaseQCSorted.bam --pBaseQC results/statsBaseQCregQual2PercentSummary.txt --regionList testFiles/region.txt --minMapQual 2 --baseSum 2> results/statsBaseQCregQual2PercentSummary.log \
&& diff results/statsBaseQCregQual2PercentSummary.txt expected/statsBaseQCregQual2Percent.txt && diff results/statsBaseQCregQual2PercentSummary.log expected/statsBaseQCregQual2PercentSummary.log \
&& \
../bin/bam stats --in testFiles/testStatsBaseQCSorted.bam --cBaseQC results/statsBaseQCregQual20Summary.txt --regionList testFiles/region.txt --minMapQual 20 --baseSum 2> results/statsBaseQCregQual20Summary.log \
&& diff results/statsBaseQCregQual20Summary.txt expected/statsBaseQCregQual20.txt && diff results/statsBaseQCregQual20Summary.log expected/statsBaseQCregQual20Summary.log  \
&& \
../bin/bam stats --in testFiles/testStatsBaseQCSorted.bam --regionList testFiles/region.txt --minMapQual 2 --baseSum 2> results/statsBaseQCregQual2SummaryNoDetail.log \
&& diff results/statsBaseQCregQual2SummaryNoDetail.log expected/statsBaseQCregQual2PercentSummary.log \
&& \
../bin/bam stats --in testFiles/testStatsBaseQCSorted.bam --regionList testFiles/region.txt --minMapQual 20 --baseSum 2> results/statsBaseQCregQual20SummaryNoDetail.log \
&& diff results/statsBaseQCregQual20SummaryNoDetail.log expected/statsBaseQCregQual20Summary.log 

