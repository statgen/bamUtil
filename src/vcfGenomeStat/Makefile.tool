EXE = vcfGenomeStat
SRCONLY = Main.cpp

#OPTFLAG?=-O3
INCLUDES=-I$(SAMTOOLS_PATH) -I$(TCLAP_PATH)/include/ -I$(INCLUDE_PATH) -I.
CFLAGS=-pipe -Wall -Wno-trigraphs $(OPTFLAG) $(INCLUDES) -D__ZLIB_AVAILABLE__ -D__STDC_LIMIT_MACROS -lcrypto -lssl

