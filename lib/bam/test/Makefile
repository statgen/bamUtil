PATH_TO_LIB_BASE=../..
include $(PATH_TO_LIB_BASE)/Makefile.lib
#
# Goncalo's Generic Makefile -- Compiles and installs a Generic Goncalo Tool
# (c) 2000-2007 Goncalo Abecasis
#

# version information
TOOL = samTest
# 0.0.1 - Initial version
VERSION = 0.0.1
MAIN = Main.cpp
OPTFLAG?=-ggdb -O0 -fno-rtti

# default installation directory
#INSTALLDIR=/usr/local/bin


# executable file names and locations
TARGET = $(TOOL)

# Source File Set
# For best results, consider editing this manually ...
TOOLBASE = WriteFiles ValidationTest ReadFiles BamIndexTest ModifyVar Modify SamFileTest Validate
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
	./samTest

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
