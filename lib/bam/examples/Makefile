PATH_TO_LIB_BASE=../..
include $(PATH_TO_LIB_BASE)/Makefile.lib

TOOL = pileup

HEADERSONLY = 
SOURCES = 
SOURCESONLY = Main.cpp


TARGET = $(TOOL)
TOOLSRC = $(SOURCES:=.cpp) $(SOURCESONLY)
TOOLOBJ = $(TOOLSRC:.cpp=.o)

LIBRARY = $(REQ_LIBS)

all: $(TARGET)

# dependencies for executables
$(TARGET) : $(TOOLOBJ) $(BINDIR)
	$(CXX) $(CFLAGS) -o $@ $(TOOLOBJ) $(LIBRARY) -lm -lz -lssl

prof : $(TOOLOBJ) $(BINDIR)
	$(CXX) $(CFLAGS) -pg -o $(TOOL) $(TOOLOBJ) $(LIBRARY) -lm -lz -lssl

$(TOOLOBJ): $(TOOLHDR)

.cpp.o : 
	$(CXX) $(CFLAGS) -o $@ -c $*.cpp -DVERSION=\"$(VERSION)\"

clean:
	-rm -rf $(TOOL) $(TOOLOBJ)

.SUFFIXES : .cpp .c .o .X.o $(SUFFIXES)

