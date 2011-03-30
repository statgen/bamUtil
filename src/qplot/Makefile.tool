VPATH=depthGCContent:libsrc:qplot
USER_INCLUDES=-Ilibsrc -IdepthGCContent -Iqplot

EXE=qplot
TOOLBASE=GCContent  GenomeRegionSeqStats QSamFlag Sequence BamQC QCStats
SRCONLY=main.cpp
HDRONLY=FlagDef.h

USER_COMPILE_VARS = -lcrypto
# DO NOT DELETE THIS LINE -- make depend depends on it
