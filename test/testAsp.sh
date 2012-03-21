../bin/bam asp --in testFiles/sortedAsp.sam --out results/asp.asp  --refFile testFilesLibBam/ref_partial.fa 2> results/asp.txt \
diff results/asp.txt expected/asp.txt && \
../bin/bam dumpAsp --asp results/asp.asp > results/dumpAsp.txt && \
diff results/dumpAsp.txt expected/dumpAsp.txt && \
../bin/bam dumpAsp --asp results/asp.asp --dataOnly > results/dumpAspDataOnly.txt && \
diff results/dumpAspDataOnly.txt expected/dumpAspDataOnly.txt 