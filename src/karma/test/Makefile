#Include ../Makefile.toolchain

#
# Goncalo's Generic Makefile -- Compiles and installs a Generic Goncalo Tool
# (c) 2000-2007 Goncalo Abecasis
#

# version information
TOOL = karma
# 0.8.7->0.8.8 fix genomeSequence.cpp
# 0.8.8S fix paired end, fix gapped alignment
# 0.9.0 -> new command line. new reference naming.  many map fixes.
#          local (re-)alignments work better.  Tom Mode.  Fix
#          SAM headers.
VERSION = 0.9.0
OPTFLAG?=-O0 -fno-rtti -Werror

# default installation directory
INSTALLDIR=/usr/local/bin

DATE=$(shell date)
NODE=$(shell uname -n)
#
# default compilation flags are 
#
# CFLAGS=-02 -I./libsrc/ -I./$(TOOL)
# 
# The following special options may also be added to the default
# 
#      Option                        Effect
#      -D_FILE_OFFSET_BITS=64        Enables support for swapfiles larger
#                                    than 2GB on supported systems (tested
#                                    on Linux)
#      -D__USE_LONG_INT              Enables support for markers with up
#                                    to 64 alleles (default is 32). Tested
#                                    on systems where gcc supports the long
#                                    long data type and on Windows.
# 
# NB: some of the quality strings in Test.cpp have trigraph character sequences in them,
# for example ??< ({) --- this forces us to require -Wno-trigraphs with g++.
#
# -D__STDC_LIMIT_MACROS is used to include C99 macros for min/max of
#  various integer type sizes (e.g. MAX_UINT32).  These are not
#  part of the C++ standard, but defining this macro includes them.
#
CFLAGS=-pipe -Wall -Wno-trigraphs -ggdb  $(OPTFLAG) -I../libcsg -I../fastQFile -D__ZLIB_AVAILABLE__  -D_FILE_OFFSET_BITS=64 -D__STDC_LIMIT_MACROS

# executable file names and locations
TARGET = $(TOOL)

# Source File Set
# For best results, consider editing this manually ...
TOOLBASE = ReadIndexer WordIndex ReadsProcessor Util MappingStats MapperBase MapperSE MapperPE SamHeader Test MapperSEBaseSpace MapperSEColorSpace MapperPEBaseSpace MapperPEColorSpace UserOptions WordHash DumpInfo
TOOLHDR = $(TOOLBASE:=.h)
TOOLSRC = $(TOOLBASE:=.cpp) Main.cpp
TOOLOBJ = $(TOOLSRC:.cpp=.o)


# Utility Library File Set
LIBRARY_ROOT = ../../../lib
LIBRARY = $(LIBRARY_ROOT)/libStatGen.a $(LIBRARY_ROOT)/samtools/libbam.a 

# private parameters
FETCHDIR=$(HOME)/code
DISTRIBDIR=$(HOME)/code/distrib/$(TOOL)-$(VERSION)

# make everything
all : $(TARGET)

# google test part
GTEST_LIB=../thirdParty/gtest-1.5.0/lib/.libs/libgtest.a
GTEST_INC=../thirdParty/gtest-1.5.0/include

#TOOLTEST=$(TOOLBASE:=_test)
TOOLTEST=gtest Integrity_check MapperBase_test MapperSEBaseSpace_test MapperPEBaseSpace_test MapperSEColorSpace_test MapperPEColorSpace_test 
TOOLTESTHDR=$(TOOLTEST:=.h)
TOOLTESTSRC=$(TOOLTEST:=.cpp)
TOOLTESTOBJ=$(TOOLTESTSRC:.cpp=.o)
TOOLOBJ_NOMAIN=$(TOOLOBJ:Main.o=)
gtest_build: $(TOOLTESTOBJ) $(TOOLOBJ) 
	g++  $(TOOLTESTOBJ) $(TOOLOBJ_NOMAIN) $(GTEST_LIB) $(LIBRARY) -lm -lz -lssl -lpthread -o ./gtest
gtest: gtest_build
	./gtest --gtest_color=yes #| grep -B 1 -A 4 'Failure'
#
# if NCBI36 names an existing reference file, try to run
# karma human test data for it - we get better code
# coverage with human genome sequences than we do with
# phiX.
#
# diff now stops on error, so when diff fails, you must
# either fix the problem or confirm that each diff is
# the result of valid improvements in the code and
# re-commit the *.good file with an explanation.
#
PHIX_REF=phiX.fa
PHIX_REF_CLEAN=phiX-*.um*
NCBI36=/data/local/ref/karma.ref/NCBI36.fa
PHIX_FASTQ=phiX.fastq phiX-pair1.fastq phiX-pair2.fastq phiX-synth-pair1.fastq phiX-synth-pair2.fastq
KARMA_BIN=../../../src/bin/karma

# test: cstest bstest
test: bstest

#bstest: karma $(PHIX_REF) $(PHIX_FASTQ)
bstest: $(PHIX_REF) $(PHIX_FASTQ)
	rm -f phiX-*.um* $(PHIX_REF_CLEAN)
	(ulimit -S -c unlimited ; $(KARMA_BIN) create -w 11 -i $(PHIX_REF))
	(ulimit -S -c unlimited ; $(KARMA_BIN) map -q 15 -w 11 -r $(PHIX_REF) -o phiX.sam -H RG:ID:boingboing -H RG:SM:NA12345 phiX.fastq)
	$(SAMTOOLS_BIN) view -b -S -o phiX.bam phiX.sam
	$(SAMTOOLS_BIN) sort phiX.bam phiX.sorted
	-diff --ignore-matching-lines='^@RG	DT:' --ignore-matching-lines='^@PG	CL:' phiX.sam.good phiX.sam 
	(ulimit -S -c unlimited ; $(KARMA_BIN) map -q 15 -w 11 -r $(PHIX_REF) -o phiX-paired.sam -H RG:ID:boingboing -H RG:SM:NA98765 phiX-pair1.fastq phiX-pair2.fastq)
	$(SAMTOOLS_BIN) view -b -S -o phiX-paired.bam phiX-paired.sam
	$(SAMTOOLS_BIN) sort phiX-paired.bam phiX-paired.sorted
	-diff --ignore-matching-lines='^@RG	DT:' --ignore-matching-lines='^@PG	CL:' phiX-paired.sam.good phiX-paired.sam 
	(ulimit -S -c unlimited ; $(KARMA_BIN) map -a 1000 -w 11 -r $(PHIX_REF) -o phiX-synth-paired.sam -H RG:ID:boingboing -H RG:SM:NA98765 phiX-synth-pair1.fastq phiX-synth-pair2.fastq)
	-diff --ignore-matching-lines='^@RG	DT:' phiX-synth-paired.sam.good phiX-synth-paired.sam 
# useless test - the sole example maps to 4 locations:
#	if [ -f $(NCBI36) ] ; then $(KARMA_BIN) map -a 2000 -o ncbi36.sam -r $(NCBI36) ncbi36_1.fastq ncbi36_2.fastq && diff --ignore-matching-lines='^@RG	DT:' ncbi36.sam ncbi36.sam.good ; else echo "ignored human genome test" ; fi

PHIX_CSFASTQ=$(PHIX_FASTQ:.fastq=.csfastq)
SAMTOOLS_BIN = ../../../lib/samtools/samtools-hybrid
%.csfastq: %.fastq
	../utilities/bs2cs.py $*.fastq
cstest: $(PHIX_REF) $(PHIX_CSFASTQ)
	rm -f phiX-*.um* $(PHIX_REF_CLEAN)
	(ulimit -S -c unlimited ; $(KARMA_BIN) create -w 11 -i $(PHIX_REF))
	(ulimit -S -c unlimited ; $(KARMA_BIN) create -c -w 11 -i $(PHIX_REF))
	(ulimit -S -c unlimited ; $(KARMA_BIN) map -c -q15 -w 11 -r $(PHIX_REF) -o phiX.cs.sam -H RG:ID:boingboing -H RG:SM:NA12345 phiX.csfastq)
	$(SAMTOOLS_BIN) view -b -S -o phiX.cs.bam phiX.cs.sam
	# -diff --ignore-matching-lines='^@RG	DT:' phiX.sam.good phiX.sam 
	$(KARMA_BIN) map -c -q 15 -w 11 -r $(PHIX_REF) -o phiX-paired.cs.sam -H RG:ID:boingboing -H RG:SM:NA98765 phiX-pair1.csfastq phiX-pair2.csfastq
	# -diff --ignore-matching-lines='^@RG	DT:' phiX-paired.sam.good phiX-paired.sam 
	$(KARMA_BIN) map -c -a 1000 -w 11 -r $(PHIX_REF) -o phiX-synth-paired.cs.sam -H RG:ID:boingboing -H RG:SM:NA98765 phiX-synth-pair1.csfastq phiX-synth-pair2.csfastq
	# -diff --ignore-matching-lines='^@RG	DT:' phiX-synth-paired.sam.good phiX-synth-paired.sam 
	# if [ -f $(NCBI36) ] ; then $(KARMA_BIN) map -a 2000 -o ncbi36.sam -r $(NCBI36) ncbi36_1.fastq ncbi36_2.fastq && diff --ignore-matching-lines='^@RG	DT:' ncbi36.sam ncbi36.sam.good ; else echo "ignored human genome test" ; fi

# helpful screen listing available options
help : 
	@echo "Generic Source Distribution"
	@echo " "
	@echo "This Makefile will compile and install" $(TOOL) "on your system"
	@echo " "
	@echo "Type...           To..."
	@echo "make help         Display this help screen"
	@echo "make all          Compile everything "
	@echo "make install      Install binaries in $(INSTALLDIR)"
	@echo "make install INSTALLDIR=directory_for_binaries"
	@echo "                  Install binaries in directory_for_binaries"
	@echo "make clean        Delete temporary files"

# dependencies for executables
$(TARGET) : $(LIBRARY) $(TOOLOBJ) $(BINDIR)
	$(CXX) $(CFLAGS) -o $@ $(TOOLOBJ) $(LIBRARY) -lm -lz -lssl -lstdc++

$(TOOLOBJ): $(TOOLHDR) $(LIBHDR)

clean :
	-rm -f *.o $(TARGET) *~
	-rm -f phiX-bs*.um* phiX.sam phiX.sam.R phiX.sam.stats phiX-paired.sam phiX-paired.sam.R phiX-paired.sam.stats phiX.bam	
	-rm -f phiX-cs*.um* phiX.cs.sam phiX.cs.sam.R phiX.cs.sam.stats phiX-paired.cs.sam phiX-paired.cs.sam.R phiX-paired.cs.sam.stats phiX.cs.bam	
	-rm -f phiX-synth-paired.sam phiX-synth-paired.sam.R phiX-synth-paired.sam.stats
	-rm -f phiX-synth-paired.cs.sam phiX-synth-paired.cs.sam.R phiX-synth-paired.cs.sam.stats
	-rm -f gtest

install : all $(INSTALLDIR)
	@echo " "
	@echo Installing to directory $(INSTALLDIR)
	@echo To select a different directory, run
	@echo " "
	@echo make install INSTALLDIR=your_preferred_dir
	@echo " "
	cp $(EXECUTABLES) $(INSTALLDIR)

$(INSTALLDIR) :
	@echo " "
	@echo Creating directory $(INSTALLDIR)
	@echo " "
	@mkdir -p $(INSTALLDIR)

new-version :
	mkdir -p $(DISTRIBDIR) $(DISTRIBDIR)/$(TOOL)
	mkdir -p $(DISTRIBDIR)/libsrc $(DISTRIBDIR)/examples
	cp ChangeLog LICENSE.twister README $(DISTRIBDIR)
	cp Makefile $(DISTRIBDIR)
	cp -R examples $(DISTRIBDIR)

fetch : 
	cd $(FETCHDIR) ; tcsh -c "cp $(TOOL)/*.{h,cpp} $(DISTRIBDIR)/$(TOOL)"
	cd $(FETCHDIR) ; cp $(LIBSRC) $(LIBHDR) $(DISTRIBDIR)/libsrc
	cd $(DISTRIBDIR); csh ../stamp Goncalo Generic

.c.o :
	$(CXX) $(CFLAGS) -o $@ -c $*.c -DVERSION=\""$(VERSION)\"" -DDATE="\"${DATE}\"" -DNODE="\"${NODE}\"" -DUSER="\"${USER}\"" -I$(GTEST_INC)

.cpp.o : 
	$(CXX) $(CFLAGS) -o $@ -c $*.cpp -DVERSION="\"$(VERSION)\"" -DDATE="\"${DATE}\"" -DNODE="\"${NODE}\"" -DUSER="\"${USER}\"" -I$(GTEST_INC)

archive : clean
	mkdir -p $(TOOL)-$(VERSION)
	cp -R LICENSE.twister README Makefile ChangeLog $(TOOL)-$(VERSION)
	cp -R libsrc pdf examples $(TOOL) $(TOOL)-$(VERSION)
	tar -cvf $(TOOL).tar $(TOOL)-$(VERSION)
	gzip -f --best $(TOOL)-$(VERSION).tar
	rm -rf $(TOOL)-$(VERSION)

distrib : $(EXECUTABLES)
	rm -rf $(TOOL)-$(VERSION)
	mkdir -p $(TOOL)-$(VERSION)
	cp -R LICENSE.twister README ChangeLog $(EXECUTABLES) examples $(TOOL)-$(VERSION)
	tar -cvf `uname`-$(TOOL)-$(VERSION).tar $(TOOL)-$(VERSION)
	gzip -f `uname`-$(TOOL)-$(VERSION).tar
	rm -rf $(TOOL)-$(VERSION)

windowszip : $(EXECUTABLES)
	mkdir -p $(TOOL)-$(VERSION)
	cp -R LICENSE.twister README ChangeLog $(EXECUTABLES) examples  $(TOOL)-$(VERSION)
	zip -r Windows-$(TOOL)-$(VERSION).zip $(TOOL)-$(VERSION)
	rm -rf $(TOOL)-$(VERSION)


.SUFFIXES : .cpp .c .o .X.o $(SUFFIXES)

# for depend target in Makefile.footer:
SOURCES = $(TOOLSRC)

#include ../Makefile.footer

# DO NOT DELETE THIS LINE -- make depend depends on it

ReadIndexer.o: ReadIndexer.h ../libcsg/GenomeSequence.h
ReadIndexer.o: ../libcsg/MemoryMapArray.h ../libcsg/MemoryMap.h
ReadIndexer.o: ../libcsg/StringArray.h ../libcsg/StringBasics.h
ReadIndexer.o: ../libcsg/InputFile.h ../libcsg/FileType.h UserOptions.h
ReadIndexer.o: ../libcsg/SmithWaterman.h ../libcsg/CigarRoller.h
ReadIndexer.o: ../libcsg/Cigar.h ../libcsg/Generic.h WordIndex.h
ReadIndexer.o: ../libcsg/IntHash.h ../libcsg/Sort.h ../libcsg/Constant.h
ReadIndexer.o: WordHash.h ../libcsg/GreedyTupleAligner.h
ReadIndexer.o: ../libcsg/TrimSequence.h
WordIndex.o: ../libcsg/Error.h ../libcsg/InplaceMerge.h Performance.h
WordIndex.o: WordIndex.h ../libcsg/GenomeSequence.h
WordIndex.o: ../libcsg/MemoryMapArray.h ../libcsg/MemoryMap.h
WordIndex.o: ../libcsg/StringArray.h ../libcsg/StringBasics.h
WordIndex.o: ../libcsg/InputFile.h ../libcsg/FileType.h ../libcsg/IntHash.h
WordIndex.o: ../libcsg/Sort.h ../libcsg/Constant.h WordHash.h Util.h
ReadsProcessor.o: ../libcsg/Generic.h MappingStats.h MapperPE.h
ReadsProcessor.o: ../libcsg/GenomeSequence.h ../libcsg/MemoryMapArray.h
ReadsProcessor.o: ../libcsg/MemoryMap.h ../libcsg/StringArray.h
ReadsProcessor.o: ../libcsg/StringBasics.h ../libcsg/InputFile.h
ReadsProcessor.o: ../libcsg/FileType.h MapperBase.h UserOptions.h
ReadsProcessor.o: ../libcsg/Random.h ReadIndexer.h ../libcsg/SmithWaterman.h
ReadsProcessor.o: ../libcsg/CigarRoller.h ../libcsg/Cigar.h WordIndex.h
ReadsProcessor.o: ../libcsg/IntHash.h ../libcsg/Sort.h ../libcsg/Constant.h
ReadsProcessor.o: WordHash.h MapperSE.h Performance.h ReadsProcessor.h
ReadsProcessor.o: SamHeader.h MapperSEBaseSpace.h MapperPEBaseSpace.h
ReadsProcessor.o: ../libcsg/Error.h ../libcsg/MathConstant.h
ReadsProcessor.o: MapperSEColorSpace.h MapperPEColorSpace.h
ReadsProcessor.o: ../libcsg/SimpleStats.h Util.h
Util.o: ../libcsg/Error.h ../libcsg/IntHash.h Util.h ../libcsg/StringBasics.h
Util.o: ../libcsg/InputFile.h ../libcsg/FileType.h
MappingStats.o: MappingStats.h MapperPE.h ../libcsg/GenomeSequence.h
MappingStats.o: ../libcsg/MemoryMapArray.h ../libcsg/MemoryMap.h
MappingStats.o: ../libcsg/StringArray.h ../libcsg/StringBasics.h
MappingStats.o: ../libcsg/InputFile.h ../libcsg/FileType.h MapperBase.h
MappingStats.o: UserOptions.h ../libcsg/Random.h ReadIndexer.h
MappingStats.o: ../libcsg/SmithWaterman.h ../libcsg/CigarRoller.h
MappingStats.o: ../libcsg/Cigar.h ../libcsg/Generic.h WordIndex.h
MappingStats.o: ../libcsg/IntHash.h ../libcsg/Sort.h ../libcsg/Constant.h
MappingStats.o: WordHash.h MapperSE.h Performance.h
MapperBase.o: MappingStats.h MapperPE.h ../libcsg/GenomeSequence.h
MapperBase.o: ../libcsg/MemoryMapArray.h ../libcsg/MemoryMap.h
MapperBase.o: ../libcsg/StringArray.h ../libcsg/StringBasics.h
MapperBase.o: ../libcsg/InputFile.h ../libcsg/FileType.h MapperBase.h
MapperBase.o: UserOptions.h ../libcsg/Random.h ReadIndexer.h
MapperBase.o: ../libcsg/SmithWaterman.h ../libcsg/CigarRoller.h
MapperBase.o: ../libcsg/Cigar.h ../libcsg/Generic.h WordIndex.h
MapperBase.o: ../libcsg/IntHash.h ../libcsg/Sort.h ../libcsg/Constant.h
MapperBase.o: WordHash.h MapperSE.h Performance.h MapperSEColorSpace.h
MapperBase.o: ReadsProcessor.h SamHeader.h MapperSEBaseSpace.h
MapperBase.o: MapperPEBaseSpace.h ../libcsg/Error.h ../libcsg/MathConstant.h
MapperBase.o: Util.h ../bam/SamFlag.h debug.h
MapperSE.o: MapperSE.h ../libcsg/GenomeSequence.h ../libcsg/MemoryMapArray.h
MapperSE.o: ../libcsg/MemoryMap.h ../libcsg/StringArray.h
MapperSE.o: ../libcsg/StringBasics.h ../libcsg/InputFile.h
MapperSE.o: ../libcsg/FileType.h MapperBase.h UserOptions.h
MapperSE.o: ../libcsg/Random.h ReadIndexer.h ../libcsg/SmithWaterman.h
MapperSE.o: ../libcsg/CigarRoller.h ../libcsg/Cigar.h ../libcsg/Generic.h
MapperSE.o: WordIndex.h ../libcsg/IntHash.h ../libcsg/Sort.h
MapperSE.o: ../libcsg/Constant.h WordHash.h MappingStats.h MapperPE.h
MapperSE.o: Performance.h ReadsProcessor.h SamHeader.h MapperSEBaseSpace.h
MapperSE.o: MapperPEBaseSpace.h ../libcsg/Error.h ../libcsg/MathConstant.h
MapperSE.o: Util.h
MapperPE.o: MapperBase.h ../libcsg/GenomeSequence.h
MapperPE.o: ../libcsg/MemoryMapArray.h ../libcsg/MemoryMap.h
MapperPE.o: ../libcsg/StringArray.h ../libcsg/StringBasics.h
MapperPE.o: ../libcsg/InputFile.h ../libcsg/FileType.h UserOptions.h
MapperPE.o: ../libcsg/Random.h ReadIndexer.h ../libcsg/SmithWaterman.h
MapperPE.o: ../libcsg/CigarRoller.h ../libcsg/Cigar.h ../libcsg/Generic.h
MapperPE.o: WordIndex.h ../libcsg/IntHash.h ../libcsg/Sort.h
MapperPE.o: ../libcsg/Constant.h WordHash.h MapperPE.h MapperSE.h
MapperPE.o: MappingStats.h Performance.h ReadsProcessor.h SamHeader.h
MapperPE.o: MapperSEBaseSpace.h MapperPEBaseSpace.h ../libcsg/Error.h
MapperPE.o: ../libcsg/MathConstant.h Util.h
SamHeader.o: SamHeader.h
Test.o: MapperSE.h ../libcsg/GenomeSequence.h ../libcsg/MemoryMapArray.h
Test.o: ../libcsg/MemoryMap.h ../libcsg/StringArray.h
Test.o: ../libcsg/StringBasics.h ../libcsg/InputFile.h ../libcsg/FileType.h
Test.o: MapperBase.h UserOptions.h ../libcsg/Random.h ReadIndexer.h
Test.o: ../libcsg/SmithWaterman.h ../libcsg/CigarRoller.h ../libcsg/Cigar.h
Test.o: ../libcsg/Generic.h WordIndex.h ../libcsg/IntHash.h ../libcsg/Sort.h
Test.o: ../libcsg/Constant.h WordHash.h MapperSEBaseSpace.h
Test.o: MapperSEColorSpace.h MapperPE.h MapperPEBaseSpace.h Test.h
MapperSEBaseSpace.o: ../libcsg/Generic.h MapperSE.h
MapperSEBaseSpace.o: ../libcsg/GenomeSequence.h ../libcsg/MemoryMapArray.h
MapperSEBaseSpace.o: ../libcsg/MemoryMap.h ../libcsg/StringArray.h
MapperSEBaseSpace.o: ../libcsg/StringBasics.h ../libcsg/InputFile.h
MapperSEBaseSpace.o: ../libcsg/FileType.h MapperBase.h UserOptions.h
MapperSEBaseSpace.o: ../libcsg/Random.h ReadIndexer.h
MapperSEBaseSpace.o: ../libcsg/SmithWaterman.h ../libcsg/CigarRoller.h
MapperSEBaseSpace.o: ../libcsg/Cigar.h WordIndex.h ../libcsg/IntHash.h
MapperSEBaseSpace.o: ../libcsg/Sort.h ../libcsg/Constant.h WordHash.h
MapperSEBaseSpace.o: MapperSEBaseSpace.h MappingStats.h MapperPE.h
MapperSEBaseSpace.o: Performance.h ReadsProcessor.h SamHeader.h
MapperSEBaseSpace.o: MapperPEBaseSpace.h ../libcsg/Error.h
MapperSEBaseSpace.o: ../libcsg/MathConstant.h Util.h
MapperSEColorSpace.o: MapperSE.h ../libcsg/GenomeSequence.h
MapperSEColorSpace.o: ../libcsg/MemoryMapArray.h ../libcsg/MemoryMap.h
MapperSEColorSpace.o: ../libcsg/StringArray.h ../libcsg/StringBasics.h
MapperSEColorSpace.o: ../libcsg/InputFile.h ../libcsg/FileType.h MapperBase.h
MapperSEColorSpace.o: UserOptions.h ../libcsg/Random.h ReadIndexer.h
MapperSEColorSpace.o: ../libcsg/SmithWaterman.h ../libcsg/CigarRoller.h
MapperSEColorSpace.o: ../libcsg/Cigar.h ../libcsg/Generic.h WordIndex.h
MapperSEColorSpace.o: ../libcsg/IntHash.h ../libcsg/Sort.h
MapperSEColorSpace.o: ../libcsg/Constant.h WordHash.h MapperSEColorSpace.h
MapperSEColorSpace.o: MappingStats.h MapperPE.h Performance.h
MapperSEColorSpace.o: ReadsProcessor.h SamHeader.h MapperSEBaseSpace.h
MapperSEColorSpace.o: MapperPEBaseSpace.h ../libcsg/Error.h
MapperSEColorSpace.o: ../libcsg/MathConstant.h Util.h
MapperPEBaseSpace.o: MapperPE.h ../libcsg/GenomeSequence.h
MapperPEBaseSpace.o: ../libcsg/MemoryMapArray.h ../libcsg/MemoryMap.h
MapperPEBaseSpace.o: ../libcsg/StringArray.h ../libcsg/StringBasics.h
MapperPEBaseSpace.o: ../libcsg/InputFile.h ../libcsg/FileType.h MapperBase.h
MapperPEBaseSpace.o: UserOptions.h ../libcsg/Random.h ReadIndexer.h
MapperPEBaseSpace.o: ../libcsg/SmithWaterman.h ../libcsg/CigarRoller.h
MapperPEBaseSpace.o: ../libcsg/Cigar.h ../libcsg/Generic.h WordIndex.h
MapperPEBaseSpace.o: ../libcsg/IntHash.h ../libcsg/Sort.h
MapperPEBaseSpace.o: ../libcsg/Constant.h WordHash.h MapperSE.h
MapperPEBaseSpace.o: MappingStats.h Performance.h ReadsProcessor.h
MapperPEBaseSpace.o: SamHeader.h MapperSEBaseSpace.h MapperPEBaseSpace.h
MapperPEBaseSpace.o: ../libcsg/Error.h ../libcsg/MathConstant.h Util.h
MapperPEColorSpace.o: MapperPE.h ../libcsg/GenomeSequence.h
MapperPEColorSpace.o: ../libcsg/MemoryMapArray.h ../libcsg/MemoryMap.h
MapperPEColorSpace.o: ../libcsg/StringArray.h ../libcsg/StringBasics.h
MapperPEColorSpace.o: ../libcsg/InputFile.h ../libcsg/FileType.h MapperBase.h
MapperPEColorSpace.o: UserOptions.h ../libcsg/Random.h ReadIndexer.h
MapperPEColorSpace.o: ../libcsg/SmithWaterman.h ../libcsg/CigarRoller.h
MapperPEColorSpace.o: ../libcsg/Cigar.h ../libcsg/Generic.h WordIndex.h
MapperPEColorSpace.o: ../libcsg/IntHash.h ../libcsg/Sort.h
MapperPEColorSpace.o: ../libcsg/Constant.h WordHash.h MapperSE.h
MapperPEColorSpace.o: MappingStats.h Performance.h ReadsProcessor.h
MapperPEColorSpace.o: SamHeader.h MapperSEBaseSpace.h MapperPEBaseSpace.h
MapperPEColorSpace.o: ../libcsg/Error.h ../libcsg/MathConstant.h Util.h
MapperPEColorSpace.o: MapperSEColorSpace.h MapperPEColorSpace.h
UserOptions.o: UserOptions.h
WordHash.o: Performance.h WordHash.h ../libcsg/GenomeSequence.h
WordHash.o: ../libcsg/MemoryMapArray.h ../libcsg/MemoryMap.h
WordHash.o: ../libcsg/StringArray.h ../libcsg/StringBasics.h
WordHash.o: ../libcsg/InputFile.h ../libcsg/FileType.h Util.h
DumpInfo.o: DumpInfo.h UserOptions.h ../fastQFile/FastQFile.h
DumpInfo.o: ../libcsg/StringBasics.h ../libcsg/InputFile.h
DumpInfo.o: ../libcsg/FileType.h ../fastQFile/BaseComposition.h
DumpInfo.o: ../libcsg/BaseAsciiMap.h ../fastQFile/BaseCount.h
DumpInfo.o: ../fastQFile/FastQStatus.h ../libcsg/GenomeSequence.h
DumpInfo.o: ../libcsg/MemoryMapArray.h ../libcsg/MemoryMap.h
DumpInfo.o: ../libcsg/StringArray.h
Main.o: DumpInfo.h UserOptions.h ../libcsg/Error.h ../libcsg/GenomeSequence.h
Main.o: ../libcsg/MemoryMapArray.h ../libcsg/MemoryMap.h
Main.o: ../libcsg/StringArray.h ../libcsg/StringBasics.h
Main.o: ../libcsg/InputFile.h ../libcsg/FileType.h MapperPE.h MapperBase.h
Main.o: ../libcsg/Random.h ReadIndexer.h ../libcsg/SmithWaterman.h
Main.o: ../libcsg/CigarRoller.h ../libcsg/Cigar.h ../libcsg/Generic.h
Main.o: WordIndex.h ../libcsg/IntHash.h ../libcsg/Sort.h ../libcsg/Constant.h
Main.o: WordHash.h MapperSE.h MapperPEBaseSpace.h ReadsProcessor.h
Main.o: SamHeader.h MapperSEBaseSpace.h MappingStats.h Performance.h Test.h
Main.o: Util.h
