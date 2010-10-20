PATH_TO_LIB_BASE=../../..
include $(PATH_TO_LIB_BASE)/Makefile.lib
#
# Goncalo's Generic Makefile -- Compiles and installs a Generic Goncalo Tool
# (c) 2000-2007 Goncalo Abecasis
#

# version information
TOOL = cigarRollerTest
# 0.0.1 - Initial version
VERSION = 0.0.1
SOURCES = CigarRollerTest.cpp
OPTFLAG?=-ggdb -O0 -fno-rtti

# default installation directory
INSTALLDIR=/usr/local/bin

# default C++ compiler
CXX=g++ 

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

CFLAGS= -Wall  $(OPTFLAG) -I../../ -D__ZLIB_AVAILABLE__  -D_FILE_OFFSET_BITS=64

# executable file names and locations
TARGET = $(TOOL)

# Source File Set
# For best results, consider editing this manually ...
TOOLBASE = CigarRollerTest
#TOOLHDR = $(TOOLBASE:=.h)
TOOLSRC = $(TOOLBASE:=.cpp)
TOOLOBJ = $(TOOLSRC:.cpp=.o)


# Utility Library File Set
LIBRARY = $(REQ_LIBS)

# private parameters
#FETCHDIR=$(HOME)/code
#DISTRIBDIR=$(HOME)/code/distrib/$(TOOL)-$(VERSION)

# make everything
all : $(TARGET)

# dependencies for executables
$(TARGET) : $(LIBRARY) $(TOOLOBJ) $(BINDIR)
	$(CXX) $(CFLAGS) -o $@ $(TOOLOBJ) $(LIBRARY) -lm -lz -lssl

$(TOOLOBJ): $(TOOLHDR) $(LIBHDR)

clean :
	-rm -f *.o $(TARGET) *~
	-rm -f results/*

test: all
	mkdir -p results
	./cigarRollerTest

.c.o :
	$(CXX) $(CFLAGS) -o $@ -c $*.c 

.cpp.o : 
	$(CXX) $(CFLAGS) -o $@ -c $*.cpp -DVERSION="\"$(VERSION)\""


DFLAGS=-Y

cleandepend:
	        makedepend -- $(DFLAGS) --

depend:
	        makedepend -- $(DFLAGS) -- $(TOOLSRC) >/dev/null 2>&1

# DO NOT DELETE THIS LINE -- make depend depends on it
