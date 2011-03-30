EXE=bam
TOOLBASE = BamExecutable Validate Convert DumpHeader SplitChromosome WriteRegion DumpIndex ReadIndexedBam DumpRefInfo Filter ReadReference Revert
SRCONLY = Main.cpp

VERSION=0.1.1
DATE=$(shell date)
USER=$(shell whoami)
USER_COMPILE_VARS = -DDATE="\"${DATE}\"" -DVERSION="\"${VERSION}\"" -DUSER="\"${USER}\"" -lcrypto

COMPILE_ANY_CHANGE = BamExecutable
