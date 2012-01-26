../bin/bam asp --in testFiles/sortedAsp.sam --out results/asp.txt  --refFile testFilesLibBam/chr1_partial.fa && \
../bin/bam dumpAsp --asp results/asp.txt > results/dumpAsp.txt && \
diff results/dumpAsp.txt expected/dumpAsp.txt && \
../bin/bam dumpAsp --asp results/asp.txt --dataOnly > results/dumpAspDataOnly.txt && \
diff results/dumpAspDataOnly.txt expected/dumpAspDataOnly.txt 