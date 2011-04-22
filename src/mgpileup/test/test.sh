../../bin/mgpileup -r ../../bam/test/testFiles/chr1_partial.fa -b testFiles/shortSorted.sam > results/mgpileup.vcf 2> results/mgpileup.txt\
 && \
diff results/mgpileup.vcf expected/mgpileup.vcf \
&& \
diff results/mgpileup.txt expected/mgpileup.txt \
&& \

# test pileup with regions.
../../bin/mgpileup -r ../../bam/test/testFiles/chr1_partial.fa -b testFiles/shortSorted.bam -i testFiles/regions.vcf > results/mgpileupRegions.vcf 2> results/mgpileupRegions.txt\
 && \
diff results/mgpileupRegions.vcf expected/mgpileupRegions.vcf \
&& \
diff results/mgpileupRegions.txt expected/mgpileupRegions.txt \
&& \

# test not getting all regions from the file.
../../bin/mgpileup -r ../../bam/test/testFiles/chr1_partial.fa -b testFiles/shortSorted.bam -i testFiles/regionShort.vcf > results/mgpileupRegionShort.vcf 2> results/mgpileupRegionShort.txt\
 && \
diff results/mgpileupRegionShort.vcf expected/mgpileupRegionShort.vcf \
&& \
diff results/mgpileupRegionShort.txt expected/mgpileupRegionShort.txt \
&& \


# test adding/removing with previous vcf.
../../bin/mgpileup -r ../../bam/test/testFiles/chr1_partial.fa -b testFiles/shortSorted.bam -i testFiles/regionMod.vcf -p testFiles/prevPileup.vcf > results/mgpileupRegionMod.vcf 2> results/mgpileupRegionMod.txt\
 && \
diff results/mgpileupRegionMod.vcf expected/mgpileupRegionMod.vcf \
&& \
diff results/mgpileupRegionMod.txt expected/mgpileupRegionMod.txt \
&& \

../../bin/mgpileup -s -r ../../bam/test/testFiles/chr1_partial.fa -b testFiles/shortSorted.sam > results/mgpileupSum.vcf 2> results/mgpileupSum.txt \
&& \
diff results/mgpileupSum.vcf expected/mgpileupSum.vcf \
&& \
diff results/mgpileupSum.txt expected/mgpileupSum.txt 


# test adding/removing with previous vcf and a storage size of 1.
../../bin/mgpileup -r ../../bam/test/testFiles/chr1_partial.fa -b testFiles/shortSorted.bam -i testFiles/regionMod.vcf -p testFiles/prevPileup.vcf -l 1 > results/mgpileupRegionMod.vcf 2> results/mgpileupRegionMod.txt\
 && \
diff results/mgpileupRegionMod.vcf expected/mgpileupRegionMod.vcf \
&& \
diff results/mgpileupRegionMod.txt expected/mgpileupRegionMod.txt \
&& \

../../bin/mgpileup -s -r ../../bam/test/testFiles/chr1_partial.fa -b testFiles/shortSorted.sam > results/mgpileupSum.vcf 2> results/mgpileupSum.txt \
&& \
diff results/mgpileupSum.vcf expected/mgpileupSum.vcf \
&& \
diff results/mgpileupSum.txt expected/mgpileupSum.txt 


if [ $? -ne 0 ]
then
  exit 1
fi

