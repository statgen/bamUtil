EXE=bam
TOOLBASE = BamExecutable Validate Convert Diff DumpHeader SplitChromosome WriteRegion DumpIndex ReadIndexedBam DumpRefInfo Filter ReadReference Revert
SRCONLY = Main.cpp

VERSION=0.1.1
DATE=$(shell date)
USER=$(shell whoami)
USER_COMPILE_VARS = -DDATE="\"${DATE}\"" -DVERSION="\"${VERSION}\"" -DUSER="\"${USER}\"" -lcrypto

COMPILE_ANY_CHANGE = BamExecutable

LIB_PATH_GENERAL ?=../libStatGen
LIB_PATH_BAM_UTIL ?= $(LIB_PATH_GENERAL)
include $(LIB_PATH_BAM_UTIL)/Makefiles/Makefile.src

#download:
#	echo $@
#	if test -d $(LIB_PATH_BAM_UTIL); then echo HI; else echo hi; fi

#$(info $(if $(LIB_PATH_BAM_UTIL), HI, hi))