SUBDIRS=bam fastqValidator glf glfMerge glfMultiples mgpileup qplot polishBam trimBam rgMergeBam splitBam squeezeBam superDeDuper thunderVCF vcfCooker vcfGenomeStat verifyBamID

# Build in all subdirectories.
#
.PHONY: $(SUBDIRS) all test clean

all: TARGET = all
test: TARGET = test
clean: TARGET = clean
all test clean: $(SUBDIRS)

$(SUBDIRS):
	@$(MAKE) OPTFLAG="$(OPTFLAG)" -C $@ $(TARGET)
