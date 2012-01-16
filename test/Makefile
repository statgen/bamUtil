TEST_COMMAND = ln -sfT $(ACTUAL_PATH)/bam/test/testFiles testFilesLibBam; \
               mkdir -p results; ./samTests.sh && \
               ./testFilter.sh && ./testSeq.sh && \
               ./testRevert.sh && ./testDiff.sh && \
               ./splitChromosome.sh && ./writeRegion.sh && \
               ./testSqueeze.sh && ./testCigars.sh && ./testStats.sh && \
               ./testClipOverlap.sh && ./testSplitBam.sh

TEST_CLEAN = rm -f testFilesLibBam

PARENT_MAKE = Makefile.test
include ../Makefile.inc
