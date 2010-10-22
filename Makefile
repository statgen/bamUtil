SUBDIRS="lib src"

# try and optimize build on SMP machines
PARALLEL_MAKE+=$(shell if [ X$(OSTYPE) = XFreeBSD ] ; then echo -j `sysctl -n hw.ncpu` ; fi)
PARALLEL_MAKE+=$(shell if [ X$(OSTYPE) = Xlinux ] ; then echo -j `grep -c '^processor' /proc/cpuinfo` ; fi)
PARALLEL_MAKE+=$(shell if [ `uname` = Linux ] ; then echo -j `grep -c '^processor' /proc/cpuinfo` ; fi)

OPTFLAG?=-O4 -fno-rtti

# redhat gcc version 4.1.2 20070626 doesn't support this flag.
# if important, figure out the right way to detect support for this flag:
# OPTFLAG+=$(shell if [ `uname` = Linux ] ; then echo '-march=native' ; fi)

FASTQ_VERSION=0.0.3
GLF_VERSION=1.0.0
BAM_VERSION=0.0.2

#
# see http://www.gnu.org/software/make/manual/make.html#Phony-Targets
# for a way of improving the following:
#
all test:
	@for i in "$(SUBDIRS)"; do \
		if [ "XXX$$i" = XXX ] ;\
		then \
		    continue; \
		fi;\
		if [ \! -d $$i ] ; \
		then \
		    echo "directory $$i does not exist, skipping." ; \
		    continue ; \
		fi ; \
		(echo "building in directory $$i";cd $$i; $(MAKE) $(PARALLEL_MAKE) OPTFLAG="$(OPTFLAG)" --no-print-directory $@) ; \
		if [ $$? -ne 0 ] ; \
		then \
		    echo "make stopped because of errors." ; \
		    break ; \
		fi \
	done

fastQRelease:
	(cd thirdParty; make clean)
	(cd libcsg; make clean)
	(cd libcsg/test; make clean)
	(cd fastQFile; make clean)
	(cd fastQFile/test; make clean)
	tar cvz --exclude="*~" --exclude-vcs -f fastQFile_$(FASTQ_VERSION).tgz -C .. pipeline/thirdParty pipeline/libcsg pipeline/fastQFile pipeline/Makefile pipeline/Makefile.toolchain


glfRelease:
	(cd thirdParty; make clean)
	(cd libcsg; make clean)
	(cd libcsg/test; make clean)
	(cd glf; make clean)
	tar cvz --exclude="*~" --exclude-vcs -f glfMerge_$(GLF_VERSION).tgz -C .. pipeline/thirdParty pipeline/libcsg pipeline/glf pipeline/Makefile pipeline/Makefile.toolchain


bamRelease:
	(cd thirdParty; make clean)
	(cd libcsg; make clean)
	(cd libcsg/test; make clean)
	(cd bam; make clean)
	tar cvz --exclude="*~" --exclude-vcs -f bam.$(BAM_VERSION).tgz -C .. pipeline/thirdParty pipeline/libcsg pipeline/bam pipeline/Makefile pipeline/Makefile.toolchain


clean:
	@for i in "$(SUBDIRS)"; do \
		if [ "XXX$$i" = XXX ] ;\
		then \
		    continue; \
		fi;\
		if [ \! -d $$i ] ; \
		then \
		    echo "directory $$i does not exist, skipping." ; \
		    continue ; \
		fi ; \
		(echo "building in directory $$i";cd $$i; $(MAKE) $(PARALLEL_MAKE) --no-print-directory $@) ; \
		if [ $$? -ne 0 ] ; \
		then \
		    echo "make stopped because of errors." ; \
		    break ; \
		fi \
	done

release:
	make clean
	find . -name CVS | xargs rm -rf
	find . -name '*~' | xargs rm -rf

