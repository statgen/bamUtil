TEST_COMMAND = ln -sfT $(LIB_PATH_BAM_UTIL)/bam/test/testFiles testFilesLibBam; \
               mkdir -p results; ./samTests.sh && \
               ./testFilter.sh && ./testSeq.sh && \
               ./testRevert.sh && ./testDiff.sh && \
               ./splitChromosome.sh && ./writeRegion.sh && \
               ./testSqueeze.sh

TEST_CLEAN = rm -f testFilesLibBam

PARENT_MAKE = Makefile.test
include ../Makefile.inc
