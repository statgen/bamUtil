EXE=mgpileup
TOOLBASE = PileupElementBaseQual
SRCONLY = Main.cpp
HDRONLY = PileupWithGenomeReference.h
USER_LIBS = ../vcfCooker/libvcf/libvcf.a
USER_INCLUDES = -I$(TCLAP_PATH)/include
