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
#OPTFLAG?=-O0 -fno-rtti -Werror
OPTFLAG?=-O0 -Werror

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
CFLAGS=-pipe -Wall -Wno-trigraphs -ggdb  $(OPTFLAG) -I.. -I../libcsg -I../fastQFile -I../../../lib/general -D__ZLIB_AVAILABLE__  -D_FILE_OFFSET_BITS=64 -D__STDC_LIMIT_MACROS -fopenmp

# executable file names and locations
TARGET = $(TOOL)

# Source File Set
# For best results, consider editing this manually ...
TOOLBASE = ../obj/ColorSpace ../obj/DumpInfo ../obj/FastqReader ../obj/MapperBase ../obj/MapperPEBaseSpace ../obj/MapperPEColorSpace ../obj/MapperPE ../obj/MapperSEBaseSpace ../obj/MapperSEColorSpace ../obj/MapperSE ../obj/MappingStats ../obj/MatchedReadBase ../obj/MatchedReadPE ../obj/MatchedReadSE ../obj/ReadIndexer ../obj/ReadsProcessor ../obj/SamHeader ../obj/Test ../obj/UserOptions ../obj/Util ../obj/WordHash ../obj/WordIndex

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
GTEST_LIB=-lgtest 
#../thirdParty/gtest-1.5.0/lib/.libs/libgtest.a
GTEST_INC=.
#GTEST_INC=../thirdParty/gtest-1.5.0/include

#TOOLTEST=$(TOOLBASE:=_test)
TOOLTEST=gtest Integrity_check MapperBase_test MapperSEBaseSpace_test MapperPEBaseSpace_test MapperSEColorSpace_test MapperPEColorSpace_test MatchedReadBase_test
TOOLTESTHDR=$(TOOLTEST:=.h)
TOOLTESTSRC=$(TOOLTEST:=.cpp)
TOOLTESTOBJ=$(TOOLTESTSRC:.cpp=.o)
TOOLOBJ_NOMAIN=$(TOOLOBJ:Main.o=)
gtest_build: $(TOOLTESTOBJ) #$(TOOLOBJ) 
	(cd .. && make -j )
	g++  $(TOOLTESTOBJ) $(TOOLOBJ_NOMAIN) $(GTEST_LIB) $(LIBRARY) -lm -lz -lssl -lpthread -fopenmp -o ./gtest
gtest: gtest_build
	./gtest --gtest_color=no #| grep -B 1 -A 4 'Failure'
	#./gtest --gtest_color=yes #| grep -B 1 -A 4 'Failure'

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

