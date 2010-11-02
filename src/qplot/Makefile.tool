VPATH=depthGCContent:libsrc:qplot
CFLAGS+=-Ilibsrc -IdepthGCContent -Iqplot

EXE=qplot
TOOLBASE=GCContent  GenomeRegionSeqStats QSamFlag Sequence BamQC QCStats
SRCONLY=main.cpp
HDRONLY=FlagDef.h

# DO NOT DELETE THIS LINE -- make depend depends on it
