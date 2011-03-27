EXE=karma
TOOLBASE = DumpInfo MapperBase MapperPEBaseSpace MapperPEColorSpace MapperPE MapperSEBaseSpace MapperSEColorSpace MapperSE MappingStats ReadIndexer ReadsProcessor SamHeader Test UserOptions Util WordHash WordIndex MatchedReadBase ColorSpace MatchedReadSE MatchedReadPE Main FastqReader
SRCONLY = 

DATE=$(shell date)
NODE=$(shell uname -n)
USER=$(shell whoami)
USER_COMPILE_VARS = -DDATE="\"${DATE}\"" -DNODE="\"${NODE}\"" -DUSER="\"${USER}\"" -Wno-trigraphs -fopenmp -lcrypto

# we will export g++ compiling flags to a hidden file
# that helps us identify how KARMA is compiled
COMPILE_FLAGS = .KARMA_COMPILE_FLAGS
%:
	@echo "KARMA is built \non $(DATE)\nby $(USER) \non $(NODE)\n" > $(COMPILE_FLAGS)
	@echo "Compiling flags:\n$(CXXFLAGS) $(CFLAGS)" >> $(COMPILE_FLAGS)
