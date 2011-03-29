EXE=bam
TOOLBASE = BamExecutable Validate Convert DumpHeader SplitChromosome WriteRegion DumpIndex ReadIndexedBam DumpRefInfo Filter ReadReference
SRCONLY = Main.cpp

<<<<<<< HEAD
USER_COMPILE_VARS = -lcrypto -I../../lib/bam

=======
VERSION=0.1.1
DATE=$(shell date)
USER=$(shell whoami)
USER_COMPILE_VARS = -DDATE="\"${DATE}\"" -DVERSION="\"${VERSION}\"" -DUSER="\"${USER}\""

COMPILE_ANY_CHANGE = BamExecutable
>>>>>>> 2745cb85be29b9dceb14ffd00a60f891e19157cf
