PATH_TO_SRC_BASE=..
include Makefile.bam

#
# Goncalo's Generic Makefile -- Compiles and installs a Generic Goncalo Tool
# (c) 2000-2007 Goncalo Abecasis
#

# Source File Set
# For best results, consider editing this manually ...
TOOLBASE = Validate Convert DumpHeader SplitChromosome WriteRegion DumpIndex ReadIndexedBam DumpRefInfo Filter ReadReference
TOOLONLYSRC = Main.cpp

TOOLHDR = $(TOOLBASE:=.h)
TOOLSRC = $(TOOLBASE:=.cpp) $(TOOLONLYSRC)
TOOLOBJ = $(TOOLSRC:.cpp=.o)

EXE=$(BAM_EXE)

# Utility Library File Set
LIBRARY = $(BAM_LIB) $(REQ_LIBS)

OBJECTS=$(patsubst %,$(OBJDIR)/%,$(TOOLOBJ))

# make everything
all : $(EXE)

# dependencies for executables
$(EXE) : $(LIBRARY) $(OBJECTS) $(BINDIR)
	$(CXX) $(CFLAGS) -o $@ $(OBJECTS) $(LIBRARY) -lm -lz -lssl

$(OBJECTS): $(TOOLHDR) $(LIBHDR) $(HDRONLY) | $(OBJDIR)

$(OBJDIR):
	mkdir $(OBJDIR)

clean :
	-rm -f $(OBJDIR)/*.o $(EXE) *~ 
	$(MAKE) $(PARALLEL_MAKE) -C $(BAM_TEST_DIR) OPTFLAG="$(OPTFLAG)" --no-print-directory $@

test: all
	$(MAKE) $(PARALLEL_MAKE) -C $(BAM_TEST_DIR) OPTFLAG="$(OPTFLAG)" --no-print-directory $@

$(OBJDIR)/%.o: %.c
	$(CXX) $(CFLAGS) -o $@ -c $*.c 

$(OBJDIR)/%.o: %.cpp 
	$(CXX) $(CFLAGS) -o $@ -c $*.cpp -DVERSION="\"$(VERSION)\""


.SUFFIXES : .cpp .c .o .X.o $(SUFFIXES)

DFLAGS=-Y

cleandepend:
	        makedepend -- $(DFLAGS) --

depend:
	        makedepend -- $(DFLAGS) -- $(TOOLSRC) >/dev/null 2>&1

# DO NOT DELETE THIS LINE -- make depend depends on it
