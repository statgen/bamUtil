EXE=bam
TOOLBASE = BamExecutable Validate Convert DumpHeader SplitChromosome WriteRegion DumpIndex ReadIndexedBam DumpRefInfo Filter ReadReference
SRCONLY = Main.cpp

VERSION=0.1.1
DATE=$(shell date)
USER=$(shell whoami)
COMMIT=$(shell git show -s --format="%H by %cn (%ci)")
USER_COMPILE_VARS = -DDATE="\"${DATE}\"" -DVERSION="\"${VERSION}\"" -DUSER="\"${USER}\"" -DCOMMIT="\"${COMMIT}\""
