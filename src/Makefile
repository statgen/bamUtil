SUBDIRS=bam fastqValidator glf glfMerge glfMultiples glfSingle karma mgpileup qplot polishBam trimBam rgMergeBam splitBam squeezeBam superDeDuper thunderVCF vcfCooker vcfGenomeStat verifyBamID

# Build in all subdirectories.
#
.PHONY: $(SUBDIRS) all test clean

all: TARGET = all
test: TARGET = test
clean: TARGET = clean
install: TARGET = install
all test clean install: $(SUBDIRS)

help : 
	@echo "Generic Source Distribution"
	@echo " "
	@echo "This Makefile will compile and install" $(TOOL) "on your system"
	@echo " "
	@echo "Type...           To..."
	@echo "make              Compile everything "
	@echo "make help         Display this help screen"
	@echo "make all          Compile everything "
	@echo "make install      Install binaries in $(INSTALLDIR)"
	@echo "make install INSTALLDIR=directory_for_binaries"
	@echo "                  Install binaries in directory_for_binaries"
	@echo "make clean        Delete temporary files"
	@echo "make test         Execute tests (if there are any)"

$(SUBDIRS):
	@$(MAKE) OPTFLAG="$(OPTFLAG)" -C $@ $(TARGET)
