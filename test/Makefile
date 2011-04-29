TEST_COMMAND = ln -sf $(TEST_LIB_PATH_BAM_UTIL)/bam/test/testFiles testFilesLibBam; \
               mkdir -p results; ./samTests.sh && \
               ./testFilter.sh && ./testSeq.sh && \
               ./testRevert.sh

TEST_CLEAN = rm -f testFilesLibBam


LIB_PATH_GENERAL ?= ../libStatGen
LIB_PATH_BAM_UTIL ?= $(LIB_PATH_GENERAL)
TEST_LIB_PATH_BAM_UTIL := $(patsubst ../%, ../../%, $(LIB_PATH_BAM_UTIL))

include $(TEST_LIB_PATH_BAM_UTIL)/Makefiles/Makefile.test