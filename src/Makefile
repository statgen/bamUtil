EXE=bam
TOOLBASE = BamExecutable Validate Convert Diff DumpHeader SplitChromosome WriteRegion DumpIndex ReadIndexedBam DumpRefInfo Filter ReadReference Revert Squeeze FindCigars Stats PileupElementBaseQCStats ClipOverlap MateMapByCoord Asp PileupElementAsp AspReadNameID DumpAsp
SRCONLY = Main.cpp

DATE=$(shell date)
USER=$(shell whoami)
USER_COMPILE_VARS = -DDATE="\"${DATE}\"" -DVERSION="\"${VERSION}\"" -DUSER="\"${USER}\"" -std=c++0x

COMPILE_ANY_CHANGE = BamExecutable

PARENT_MAKE = Makefile.src
include ../Makefile.inc
