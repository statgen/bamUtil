../../bin/mgpileup -r ../../bam/test/testFiles/chr1_partial.fa -b testFiles/shortSorted.sam > results/mgpileup.txt \
 && \
diff results/mgpileup.txt expected/mgpileup.txt \
&& \
../../bin/mgpileup -s -r ../../bam/test/testFiles/chr1_partial.fa -b testFiles/shortSorted.sam > results/mgpileupSum.txt \
&& \
diff results/mgpileupSum.txt expected/mgpileupSum.txt 


if [ $? -ne 0 ]
then
  exit 1
fi

