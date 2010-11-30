DIR_NAME=general
PATH_TO_LIB_BASE=..
PATH_TO_BASE=../$(PATH_TO_LIB_BASE)
include $(PATH_TO_BASE)/Makefile.include

PATH_INCLUDE_TO_LIB_BASE=..

#
# trying to strike the right balance on error/warning
# handling is difficult, but here we make most everthing
# an error, but allow unused results and variables for the
# time being.
#
#CFLAGS=$(OPTFLAG) -Werror -Wno-unused-variable -Wno-unused-result -pipe -Wall -ggdb -D__ZLIB_AVAILABLE__  -D_FILE_OFFSET_BITS=64 -I$(SAMTOOLS_PATH) -D__STDC_LIMIT_MACROS

WARNINGS=-Wno-unused-variable
WARNINGS+=$(shell if [ X$(CCVERSION) \> X4.1.0 ] ; then echo " -Wno-unused-result" ; fi)
CFLAGS=$(OPTFLAG) -Werror $(WARNINGS) -pipe -Wall -ggdb -D__ZLIB_AVAILABLE__  -D_FILE_OFFSET_BITS=64 -I$(SAMTOOLS_PATH) -D__STDC_LIMIT_MACROS

SOURCES_WITH_HEADER=\
	BaseAsciiMap.cpp \
	BasicHash.cpp \
	BaseQualityHelper.cpp \
	BaseUtilities.cpp \
	BgzfFileType.cpp \
	CharBuffer.cpp \
	Cigar.cpp \
	CigarRoller.cpp \
	Error.cpp \
	GenomeSequence.cpp \
	GenotypeLists.cpp \
	GzipFileType.cpp \
	GzipHeader.cpp \
	glfHandler.cpp \
	FastQ.cpp \
	FileType.cpp \
	FortranFormat.cpp \
	Hash.cpp \
	Input.cpp \
	InputFile.cpp \
	IntArray.cpp \
	IntHash.cpp \
	LongHash.cpp \
	LongLongCounter.cpp \
	MapFunction.cpp \
	MathVector.cpp \
	MemoryAllocators.cpp \
	MemoryInfo.cpp \
	MemoryMapArray.cpp \
	MemoryMap.cpp \
	MiniDeflate.cpp \
	Parameters.cpp \
	PedigreeAlleleFreq.cpp \
	Pedigree.cpp \
	PedigreeDescription.cpp \
	PedigreeFamily.cpp \
	PedigreeGlobals.cpp \
	PedigreePerson.cpp \
	QuickIndex.cpp \
	Random.cpp \
	SmithWaterman.cpp \
	Sort.cpp \
	STLUtilities.cpp \
	StringAlias.cpp \
	StringArray.cpp \
	StringBasics.cpp \
	StringHash.cpp \
	StringMap.cpp \
	TrimSequence.cpp \
	UncompressedFileType.cpp \
	WindowsHelper.cpp

SOURCES_NO_HEADER=\
	PedigreeLoader.cpp \
	PedigreeTrim.cpp \
	PedigreeTwin.cpp 

SOURCES=$(SOURCES_WITH_HEADER) $(SOURCES_NO_HEADER) 

HEADERS := $(wildcard *.h)

OBJECT_FILES=$(SOURCES:.cpp=.o)

OBJECTS=$(patsubst %,$(OBJDIR)/%,$(OBJECT_FILES))


all: $(addprefix $(INCLUDE_PATH)/, $(HEADERS)) $(STAT_GEN_LIB)

$(STAT_GEN_LIB) : $(OBJECTS)
	$(AR) cr $@ $(OBJECTS)
	$(RANLIB) $@

$(INCLUDE_PATH)/%.h: %.h
	ln  -s $(PATH_INCLUDE_TO_LIB_BASE)/$(DIR_NAME)/$< $@

$(OBJDIR)/%.o : %.cpp
	$(CXX) $(CFLAGS) -o $@ -c $*.cpp -DVERSION=\"$(VERSION)\"

test : all
	(cd test; $(MAKE) $(PARALLEL_MAKE) OPTFLAG="$(OPTFLAG)" --no-print-directory $@);


CPP_TESTS := $(shell grep -m 1 'if defined(TEST)' *cpp |cut -f1 -d: )
selftest: 
	@for i in $(CPP_TESTS); do \
		if [ "XXX$$i" = XXX ] ;\
		then \
		    continue; \
		fi;\
		(echo "building in self test $$i"; g++ -DTEST -D__STDC_LIMIT_MACROS -fno-rtti -o test_$$i $$i  -L. -lcsg -lm -lz -lssl ; ./test_$$i) ; \
		if [ $$? -ne 0 ] ; \
		then \
		    echo "make stopped because of errors." ; \
		    break ; \
		fi \
	done
clean:
	-rm -rf $(OBJECTS) $(addprefix $(INCLUDE_PATH)/, $(HEADERS))
	$(AR) d $(STAT_GEN_LIB) $(OBJECTS)
	$(RANLIB) $(STAT_GEN_LIB)
	(cd test; $(MAKE) $(PARALLEL_MAKE) OPTFLAG="$(OPTFLAG)" --no-print-directory $@);

.SUFFIXES : .cpp .c .o .X.o $(SUFFIXES)

include $(PATH_TO_BASE)/Makefile.footer

# DO NOT DELETE THIS LINE -- make depend depends on it

$(OBJDIR)/BaseAsciiMap.o: BaseAsciiMap.h StringBasics.h InputFile.h FileType.h
$(OBJDIR)/BasicHash.o: BasicHash.h Error.h
$(OBJDIR)/BaseQualityHelper.o: BaseQualityHelper.h
$(OBJDIR)/BaseUtilities.o: BaseUtilities.h
$(OBJDIR)/BgzfFileType.o: BgzfFileType.h $(SAMTOOLS_PATH)/bgzf.h
$(OBJDIR)/BgzfFileType.o: FileType.h
$(OBJDIR)/CharBuffer.o: CharBuffer.h InputFile.h
$(OBJDIR)/Cigar.o: Cigar.h Generic.h StringBasics.h InputFile.h FileType.h
$(OBJDIR)/Cigar.o: STLUtilities.h
$(OBJDIR)/CigarRoller.o: CigarRoller.h Cigar.h Generic.h StringBasics.h InputFile.h
$(OBJDIR)/CigarRoller.o: FileType.h
$(OBJDIR)/Error.o: Error.h
$(OBJDIR)/GenomeSequence.o: Error.h Generic.h GenomeSequence.h MemoryMapArray.h
$(OBJDIR)/GenomeSequence.o: MemoryMap.h StringArray.h StringBasics.h InputFile.h
$(OBJDIR)/GenomeSequence.o: FileType.h
$(OBJDIR)/GenotypeLists.o: GenotypeLists.h Pedigree.h Constant.h PedigreeAlleles.h
$(OBJDIR)/GenotypeLists.o: LongInt.h PedigreePerson.h PedigreeGlobals.h StringArray.h
$(OBJDIR)/GenotypeLists.o: StringBasics.h InputFile.h FileType.h StringHash.h Hash.h
$(OBJDIR)/GenotypeLists.o: IntArray.h MathVector.h PedigreeFamily.h
$(OBJDIR)/GenotypeLists.o: PedigreeDescription.h PedigreeAlleleFreq.h
$(OBJDIR)/GzipFileType.o: GzipFileType.h
$(OBJDIR)/GzipFileType.o: FileType.h
$(OBJDIR)/GzipHeader.o: GzipHeader.h
$(OBJDIR)/glfHandler.o: glfHandler.h InputFile.h FileType.h StringBasics.h
$(OBJDIR)/glfHandler.o: BaseQualityHelper.h
$(OBJDIR)/FastQ.o: FastQ.h
$(OBJDIR)/FileType.o: FileType.h
$(OBJDIR)/FortranFormat.o: FortranFormat.h StringBasics.h InputFile.h FileType.h
$(OBJDIR)/FortranFormat.o: IntArray.h Error.h
$(OBJDIR)/Hash.o: Hash.h
$(OBJDIR)/Input.o: Input.h Error.h Constant.h
$(OBJDIR)/InputFile.o: InputFile.h FileType.h StringBasics.h GzipHeader.h
$(OBJDIR)/InputFile.o: BgzfFileType.h $(SAMTOOLS_PATH)/bgzf.h
$(OBJDIR)/InputFile.o: GzipFileType.h UncompressedFileType.h
$(OBJDIR)/IntArray.o: IntArray.h Error.h Hash.h Sort.h Constant.h
$(OBJDIR)/IntHash.o: IntHash.h Error.h
$(OBJDIR)/LongLongCounter.o: LongLongCounter.h LongHash.h Error.h
$(OBJDIR)/MapFunction.o: MapFunction.h MathConstant.h
$(OBJDIR)/MathVector.o: MathVector.h StringBasics.h InputFile.h FileType.h MathMatrix.h
$(OBJDIR)/MathVector.o: Error.h MathConstant.h Sort.h Constant.h
$(OBJDIR)/MemoryAllocators.o: MemoryAllocators.h
$(OBJDIR)/MemoryInfo.o: MemoryInfo.h StringBasics.h InputFile.h FileType.h
$(OBJDIR)/MemoryMapArray.o: MemoryMapArray.h MemoryMap.h
$(OBJDIR)/MemoryMap.o: MemoryMap.h
$(OBJDIR)/MiniDeflate.o: MiniDeflate.h
$(OBJDIR)/Parameters.o: Parameters.h StringMap.h StringBasics.h InputFile.h FileType.h
$(OBJDIR)/Parameters.o: Constant.h MathConstant.h Error.h
$(OBJDIR)/PedigreeAlleleFreq.o: PedigreeAlleleFreq.h Pedigree.h Constant.h
$(OBJDIR)/PedigreeAlleleFreq.o: PedigreeAlleles.h LongInt.h PedigreePerson.h
$(OBJDIR)/PedigreeAlleleFreq.o: PedigreeGlobals.h StringArray.h StringBasics.h
$(OBJDIR)/PedigreeAlleleFreq.o: InputFile.h FileType.h StringHash.h Hash.h IntArray.h
$(OBJDIR)/PedigreeAlleleFreq.o: MathVector.h PedigreeFamily.h PedigreeDescription.h
$(OBJDIR)/PedigreeAlleleFreq.o: QuickIndex.h StringMap.h Error.h
$(OBJDIR)/Pedigree.o: Pedigree.h Constant.h PedigreeAlleles.h LongInt.h
$(OBJDIR)/Pedigree.o: PedigreePerson.h PedigreeGlobals.h StringArray.h StringBasics.h
$(OBJDIR)/Pedigree.o: InputFile.h FileType.h StringHash.h Hash.h IntArray.h
$(OBJDIR)/Pedigree.o: MathVector.h PedigreeFamily.h PedigreeDescription.h
$(OBJDIR)/Pedigree.o: PedigreeAlleleFreq.h GenotypeLists.h MemoryInfo.h Error.h Sort.h
$(OBJDIR)/PedigreeDescription.o: PedigreeDescription.h PedigreeGlobals.h Constant.h
$(OBJDIR)/PedigreeDescription.o: StringArray.h StringBasics.h InputFile.h FileType.h
$(OBJDIR)/PedigreeDescription.o: StringHash.h Hash.h IntArray.h MathVector.h
$(OBJDIR)/PedigreeDescription.o: PedigreePerson.h PedigreeAlleles.h LongInt.h
$(OBJDIR)/PedigreeDescription.o: MapFunction.h FortranFormat.h Error.h
$(OBJDIR)/PedigreeFamily.o: Pedigree.h Constant.h PedigreeAlleles.h LongInt.h
$(OBJDIR)/PedigreeFamily.o: PedigreePerson.h PedigreeGlobals.h StringArray.h
$(OBJDIR)/PedigreeFamily.o: StringBasics.h InputFile.h FileType.h StringHash.h Hash.h
$(OBJDIR)/PedigreeFamily.o: IntArray.h MathVector.h PedigreeFamily.h
$(OBJDIR)/PedigreeFamily.o: PedigreeDescription.h PedigreeAlleleFreq.h MathConstant.h
$(OBJDIR)/PedigreeFamily.o: Error.h
$(OBJDIR)/PedigreeGlobals.o: PedigreeGlobals.h Constant.h StringArray.h StringBasics.h
$(OBJDIR)/PedigreeGlobals.o: InputFile.h FileType.h StringHash.h Hash.h IntArray.h
$(OBJDIR)/PedigreeGlobals.o: MathVector.h Sort.h Error.h
$(OBJDIR)/PedigreeLoader.o: Pedigree.h Constant.h PedigreeAlleles.h LongInt.h
$(OBJDIR)/PedigreeLoader.o: PedigreePerson.h PedigreeGlobals.h StringArray.h
$(OBJDIR)/PedigreeLoader.o: StringBasics.h InputFile.h FileType.h StringHash.h Hash.h
$(OBJDIR)/PedigreeLoader.o: IntArray.h MathVector.h PedigreeFamily.h
$(OBJDIR)/PedigreeLoader.o: PedigreeDescription.h PedigreeAlleleFreq.h FortranFormat.h
$(OBJDIR)/PedigreeLoader.o: Error.h
$(OBJDIR)/PedigreePerson.o: PedigreePerson.h Constant.h PedigreeAlleles.h LongInt.h
$(OBJDIR)/PedigreePerson.o: PedigreeGlobals.h StringArray.h StringBasics.h InputFile.h
$(OBJDIR)/PedigreePerson.o: FileType.h StringHash.h Hash.h IntArray.h MathVector.h
$(OBJDIR)/PedigreePerson.o: Error.h
$(OBJDIR)/PedigreeTrim.o: Pedigree.h Constant.h PedigreeAlleles.h LongInt.h
$(OBJDIR)/PedigreeTrim.o: PedigreePerson.h PedigreeGlobals.h StringArray.h
$(OBJDIR)/PedigreeTrim.o: StringBasics.h InputFile.h FileType.h StringHash.h Hash.h
$(OBJDIR)/PedigreeTrim.o: IntArray.h MathVector.h PedigreeFamily.h
$(OBJDIR)/PedigreeTrim.o: PedigreeDescription.h PedigreeAlleleFreq.h
$(OBJDIR)/PedigreeTwin.o: Pedigree.h Constant.h PedigreeAlleles.h LongInt.h
$(OBJDIR)/PedigreeTwin.o: PedigreePerson.h PedigreeGlobals.h StringArray.h
$(OBJDIR)/PedigreeTwin.o: StringBasics.h InputFile.h FileType.h StringHash.h Hash.h
$(OBJDIR)/PedigreeTwin.o: IntArray.h MathVector.h PedigreeFamily.h
$(OBJDIR)/PedigreeTwin.o: PedigreeDescription.h PedigreeAlleleFreq.h Error.h
$(OBJDIR)/$(OBJDIR)/QuickIndex.o: QuickIndex.h MathVector.h StringBasics.h InputFile.h FileType.h
$(OBJDIR)/QuickIndex.o: StringArray.h IntArray.h StringMap.h Error.h
$(OBJDIR)/Random.o: Random.h MathConstant.h Error.h
$(OBJDIR)/SmithWaterman.o: SmithWaterman.h CigarRoller.h Cigar.h Generic.h
$(OBJDIR)/SmithWaterman.o: StringBasics.h InputFile.h FileType.h
$(OBJDIR)/Sort.o: Sort.h Constant.h Error.h
$(OBJDIR)/STLUtilities.o: STLUtilities.h
$(OBJDIR)/StringAlias.o: StringAlias.h
$(OBJDIR)/StringArray.o: StringArray.h StringBasics.h InputFile.h FileType.h Sort.h
$(OBJDIR)/StringArray.o: Constant.h Error.h
$(OBJDIR)/StringBasics.o: StringBasics.h InputFile.h FileType.h Error.h Constant.h
$(OBJDIR)/StringBasics.o: MathConstant.h
$(OBJDIR)/StringHash.o: StringHash.h StringBasics.h InputFile.h FileType.h Constant.h
$(OBJDIR)/StringHash.o: Hash.h Error.h
$(OBJDIR)/StringMap.o: StringMap.h StringBasics.h InputFile.h FileType.h
$(OBJDIR)/TrimSequence.o: TrimSequence.h
$(OBJDIR)/UncompressedFileType.o: UncompressedFileType.h FileType.h
$(OBJDIR)/WindowsHelper.o: WindowsHelper.h
