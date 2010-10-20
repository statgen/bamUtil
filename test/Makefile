PATH_TO_SRC_BASE=../..
include ../Makefile.bam
#
# Goncalo's Generic Makefile -- Compiles and installs a Generic Goncalo Tool
# (c) 2000-2007 Goncalo Abecasis
#

# version information
TOOL = 
# 0.0.1 - Initial version
VERSION = 0.0.1
MAIN = 

# default installation directory
#INSTALLDIR=/usr/local/bin


# executable file names and locations
TARGET = $(TOOL)

# Source File Set
# For best results, consider editing this manually ...
TOOLBASE = 
TOOLHDR = $(TOOLBASE:=.h)
TOOLONLYSRC = $(TOOLBASE:=.cpp)
TOOLSRC = $(TOOLBASE:=.cpp) $(MAIN)
TOOLOBJ = $(TOOLSRC:.cpp=.o)
TOOLONLYOBJ = $(TOOLONLYSRC:.cpp=.o)


# Utility Library File Set
LIBRARY = $(REQ_LIBS)
 

# make everything
all : $(TARGET)

# dependencies for executables
$(TARGET) : $(LIBRARY) $(TOOLOBJ) $(BINDIR)
	$(CXX) $(CFLAGS) -o $@ $(TOOLOBJ) $(LIBRARY) -lm -lz -lssl

$(TOOLOBJ): $(TOOLHDR) $(LIBHDR)


test: all
	./samTests.sh
	./testFilter.sh

clean :
	-rm -f *.o $(TARGET) *~

.c.o :
	$(CXX) $(CFLAGS) -o $@ -c $*.c 

.cpp.o : 
	$(CXX) $(CFLAGS) -o $@ -c $*.cpp -DVERSION="\"$(VERSION)\""



.SUFFIXES : .cpp .c .o .X.o $(SUFFIXES)

DFLAGS=-Y

cleandepend:
	        makedepend -- $(DFLAGS) --

depend:
	        makedepend -- $(DFLAGS) -- $(TOOLSRC) >/dev/null 2>&1

# DO NOT DELETE THIS LINE -- make depend depends on it
