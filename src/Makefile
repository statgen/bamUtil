EXE=bam
TOOLBASE = BamExecutable Validate Convert Diff DumpHeader SplitChromosome WriteRegion DumpIndex ReadIndexedBam DumpRefInfo Filter ReadReference Revert Squeeze FindCigars Stats PileupElementBaseQCStats ClipOverlap MateMapByCoord SplitBam TrimBam RGMergeBam PolishBam GapInfo Logger Bam2FastQ Dedup Prediction LogisticRegression MathCholesky HashErrorModel Recab
SRCONLY = Main.cpp
HDRONLY = Covariates.h

DATE=$(shell date)
USER=$(shell whoami)
USER_COMPILE_VARS = -std=c++0x -DDATE="\"${DATE}\"" -DVERSION="\"${VERSION}\"" -DUSER="\"${USER}\""
COMPILE_ANY_CHANGE = BamExecutable

PARENT_MAKE = Makefile.src
include ../Makefile.inc
