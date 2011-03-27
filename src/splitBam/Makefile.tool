EXE = splitBam
SRCONLY = SplitBam.cpp

#OPTFLAG?=-O3
#CFLAGS= -Wall -Wno-write-strings $(OPTFLAG) -I$(SAMTOOLS_PATH) -I$(INCLUDE_PATH) -D__ZLIB_AVAILABLE__  -D_FILE_OFFSET_BITS=64

USER_COMPILE_VARS = -lcrypto -lssl
