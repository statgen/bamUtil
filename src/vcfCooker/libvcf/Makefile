TOOLBASE = VcfFile
SRCONLY = 
HDRONLY = 
NEW_LIB=libvcf.a
VERSION = 0.0.1
PATH_TO_BASE=../../..

include $(PATH_TO_BASE)/Makefile.include

PATH_INCLUDE_TO_LIB_BASE=..
#
# Goncalo's Generic Makefile -- Compiles and installs a Generic Goncalo Tool
# (c) 2000-2007 Goncalo Abecasis
#

# Source File Set
TOOLHDR = $(TOOLBASE:=.h) $(HDRONLY)
TOOLSRC = $(TOOLBASE:=.cpp) $(SRCONLY)
TOOLOBJ = $(TOOLSRC:.cpp=.o)

HEADERS=$(TOOLHDR)
NEWLIBRARY=$(NEW_LIB)
OBJECTS=$(patsubst %,$(OBJDIR)/%,$(TOOLOBJ))

# make everything
all : $(NEWLIBRARY)

$(OBJECTS): $(TOOLHDR) | $(OBJDIR)

$(OBJDIR):
	mkdir $(OBJDIR)

$(NEWLIBRARY) : $(OBJECTS)
	ar -cr $@ $(OBJECTS)
	ranlib $@

clean :
	-rm -f $(OBJDIR)/*.o *~  $(NEW_LIB)

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
