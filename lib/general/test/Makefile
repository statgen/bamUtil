PATH_TO_LIB_BASE=../..
PATH_TO_BASE=../$(PATH_TO_LIB_BASE)
include $(PATH_TO_BASE)/Makefile.include

SUBDIRS=inputFileTest cigar

OPTFLAG?=-O0

.PHONY: $(SUBDIRS) all test clean

all: TARGET = all
test: TARGET = test
clean: TARGET = clean

all test: $(SUBDIRS)

clean: $(SUBDIRS)
	# clean current directory
	rm -f *.o *~

$(SUBDIRS):
	@$(MAKE) OPTFLAG="$(OPTFLAG)" -C $@ $(TARGET)

######################################################################

GTEST_COMPILE_FLAG=-DGTEST_HAS_TR1_TUPLE=0 -D GTEST_HAS_RTTI=0


CFLAGS=-pipe -Wall -Wno-trigraphs -ggdb  $(OPTFLAG) -I. -I.. -I$(INCLUDE_PATH) -D__ZLIB_AVAILABLE__  -D_FILE_OFFSET_BITS=64 -D__STDC_LIMIT_MACROS $(GTEST_COMPILE_FLAG)

LIBRARY=$(REQ_LIBS)

.c.o :
	$(CXX) $(CFLAGS) -o $@ -c $*.c -DVERSION=\""$(VERSION)\"" -DDATE="\"${DATE}\"" -DNODE="\"${NODE}\"" -DUSER="\"${USER}\"" -I$(GTEST_INC)

.cpp.o : 
	$(CXX) $(CFLAGS) -o $@ -c $*.cpp -DVERSION="\"$(VERSION)\"" -DDATE="\"${DATE}\"" -DNODE="\"${NODE}\"" -DUSER="\"${USER}\"" -I$(GTEST_INC)

# Google test part
GTEST_INC = ../../../src/karma/gtest
GTEST_SRC = ../../../src/karma/gtest/gtest/gtest-all.cc
GTEST_OBJ = ../../../src/karma/gtest/gtest/gtest-all.o

$(GTEST_OBJ): $(GTEST_SRC)
	g++ -c -o $@ $< -I$(GTEST_INC)

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

gtest_build: $(TOOLTESTOBJ) $(TOOLOBJ) $(GTEST_OBJ)
	g++  $(TOOLTESTOBJ) $(TOOLOBJ_NOMAIN) $(GTEST_OBJ) $(LIBRARY) -lm -lz -lssl -lpthread -lcrypto -o ./gtest
gtest: gtest_build
	./gtest --gtest_color=yes

include $(PATH_TO_BASE)/Makefile.footer
