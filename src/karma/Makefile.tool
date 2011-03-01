EXE=karma
TOOLBASE = DumpInfo MapperBase MapperPEBaseSpace MapperPEColorSpace MapperPE MapperSEBaseSpace MapperSEColorSpace MapperSE MappingStats ReadIndexer ReadsProcessor SamHeader Test UserOptions Util WordHash WordIndex MatchedReadBase ColorSpace MatchedReadSE MatchedReadPE Main
SRCONLY = 

DATE=$(shell date)
NODE=$(shell uname -n)
USER=$(shell whoami)
USER_COMPILE_VARS = -DDATE="\"${DATE}\"" -DNODE="\"${NODE}\"" -DUSER="\"${USER}\"" -Wno-trigraphs -fopenmp

#TESTBASE = MapperBase_test MapperPEBaseSpace_test MapperPEColorSpace_test MapperSEBaseSpace_test MapperSEColorSpace_test 

