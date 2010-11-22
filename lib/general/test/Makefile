PATH_TO_LIB_BASE=../..
PATH_TO_BASE=../$(PATH_TO_LIB_BASE)
include $(PATH_TO_BASE)/Makefile.include

SUBDIRS="inputFileTest cigar"

# try and optimize build on SMP machines
PARALLEL_MAKE+=$(shell if [ X$(OSTYPE) = XFreeBSD ] ; then echo -j `sysctl -n hw.ncpu` ; fi)
PARALLEL_MAKE+=$(shell if [ X$(OSTYPE) = Xlinux ] ; then echo -j `grep -c '^processor' /proc/cpuinfo` ; fi)
PARALLEL_MAKE+=$(shell if [ `uname` = Linux ] ; then echo -j `grep -c '^processor' /proc/cpuinfo` ; fi)

OPTFLAG?=-O0

all test:
	@for i in "$(SUBDIRS)"; do \
		if [ "XXX$$i" = XXX ] ;\
		then \
		    continue; \
		fi;\
		if [ \! -d $$i ] ; \
		then \
		    echo "directory $$i does not exist, skipping." ; \
		    continue ; \
		fi ; \
		(echo "building in directory $$i";cd $$i; $(MAKE) $(PARALLEL_MAKE) OPTFLAG="$(OPTFLAG)" --no-print-directory $@) ; \
		if [ $$? -ne 0 ] ; \
		then \
		    echo "make stopped because of errors." ; \
		    break ; \
		fi \
	done

clean:
	@for i in "$(SUBDIRS)"; do \
		if [ "XXX$$i" = XXX ] ;\
		then \
		    continue; \
		fi;\
		if [ \! -d $$i ] ; \
		then \
		    echo "directory $$i does not exist, skipping." ; \
		    continue ; \
		fi ; \
		(echo "building in directory $$i";cd $$i; $(MAKE) $(PARALLEL_MAKE) --no-print-directory $@) ; \
		if [ $$? -ne 0 ] ; \
		then \
		    echo "make stopped because of errors." ; \
		    break ; \
		fi \
	done
	# clean current directory
	rm -f *.o *~

######################################################################

GTEST_COMPILE_FLAG=-DGTEST_HAS_TR1_TUPLE=0 -D GTEST_HAS_RTTI=0

CFLAGS=-pipe -Wall -Wno-trigraphs -ggdb  $(OPTFLAG) -I. -I.. -I$(INCLUDE_PATH) -D__ZLIB_AVAILABLE__  -D_FILE_OFFSET_BITS=64 -D__STDC_LIMIT_MACROS $(GTEST_COMPILE_FLAG)

LIBRARY=$(REQ_LIBS)

.c.o :
	$(CXX) $(CFLAGS) -o $@ -c $*.c -DVERSION=\""$(VERSION)\"" -DDATE="\"${DATE}\"" -DNODE="\"${NODE}\"" -DUSER="\"${USER}\"" -I$(GTEST_INC)

.cpp.o : 
	$(CXX) $(CFLAGS) -o $@ -c $*.cpp -DVERSION="\"$(VERSION)\"" -DDATE="\"${DATE}\"" -DNODE="\"${NODE}\"" -DUSER="\"${USER}\"" -I$(GTEST_INC)

# google test part
GTEST_LIB=../../thirdParty/gtest-1.5.0/lib/.libs/libgtest.a
GTEST_INC=../../thirdParty/gtest-1.5.0/include

#TOOLTEST=$(TOOLBASE:=_test)
TOOLTEST=gtest TrimSequence_test GreedyTupleAligner_test STLUtilities_test GenomeSequence_test
TOOLTESTHDR=$(TOOLTEST:=.h)
TOOLTESTSRC=$(TOOLTEST:=.cpp)
TOOLTESTOBJ=$(TOOLTESTSRC:.cpp=.o)
TOOLOBJ_NOMAIN=$(TOOLOBJ:Main.o=)

GenomeSequence_test.o: GenomeSequence_test.cpp ../GenomeSequence.h ../GenomeSequence.cpp
	$(CXX) $(CFLAGS) -o $@ -c $< -DVERSION="\"$(VERSION)\"" -DDATE="\"${DATE}\"" -DNODE="\"${NODE}\"" -DUSER="\"${USER}\"" -I$(GTEST_INC)

STLUtilities_test.o: STLUtilities_test.cpp ../STLUtilities.h ../STLUtilities.cpp
	$(CXX) $(CFLAGS) -o $@ -c $< -DVERSION="\"$(VERSION)\"" -DDATE="\"${DATE}\"" -DNODE="\"${NODE}\"" -DUSER="\"${USER}\"" -I$(GTEST_INC)

TrimSequence_test.o: TrimSequence_test.cpp ../TrimSequence.h ../TrimSequence.cpp
	$(CXX) $(CFLAGS) -o $@ -c $< -DVERSION="\"$(VERSION)\"" -DDATE="\"${DATE}\"" -DNODE="\"${NODE}\"" -DUSER="\"${USER}\"" -I$(GTEST_INC)

GreedyTupleAligner_test.o: GreedyTupleAligner_test.cpp ../GreedyTupleAligner.h
	$(CXX) $(CFLAGS) -o $@ -c $< -DVERSION="\"$(VERSION)\"" -DDATE="\"${DATE}\"" -DNODE="\"${NODE}\"" -DUSER="\"${USER}\"" -I$(GTEST_INC)

gtest_build: $(TOOLTESTOBJ) $(TOOLOBJ) 
	g++  $(TOOLTESTOBJ) $(TOOLOBJ_NOMAIN) $(GTEST_LIB) $(LIBRARY) -lm -lz -lssl -lpthread -o ./gtest
gtest: gtest_build
	./gtest --gtest_color=yes

include $(PATH_TO_BASE)/Makefile.footer
