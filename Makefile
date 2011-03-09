SUBDIRS="lib src"

# try and optimize build on SMP machines
PARALLEL_MAKE+=$(shell if [ X$(OSTYPE) = XFreeBSD ] ; then echo -j `sysctl -n hw.ncpu` ; fi)
PARALLEL_MAKE+=$(shell if [ X$(OSTYPE) = Xlinux ] ; then echo -j `grep -c '^processor' /proc/cpuinfo` ; fi)
PARALLEL_MAKE+=$(shell if [ `uname` = Linux ] ; then echo -j `grep -c '^processor' /proc/cpuinfo` ; fi)

OPTFLAG?=-O4

# redhat gcc version 4.1.2 20070626 doesn't support this flag.
# if important, figure out the right way to detect support for this flag:
# OPTFLAG+=$(shell if [ `uname` = Linux ] ; then echo '-march=native' ; fi)

VERSION=0.1.4
RELEASE_FILE=statGen.$(VERSION).tgz

.PHONY: all test clean release help install

#
# see http://www.gnu.org/software/make/manual/make.html#Phony-Targets
# for a way of improving the following:
#
all test install:
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

help : 
	@echo "Generic Source Distribution"
	@echo " "
	@echo "This Makefile will compile and install" $(TOOL) "on your system"
	@echo " "
	@echo "Type...           To..."
	@echo "make              Compile everything "
	@echo "make help         Display this help screen"
	@echo "make all          Compile everything "
	@echo "make install      Install binaries in $(INSTALLDIR)"
	@echo "make install INSTALLDIR=directory_for_binaries"
	@echo "                  Install binaries in directory_for_binaries"
	@echo "make clean        Delete temporary files"
	@echo "make test         Execute tests (if there are any)"

release:
	(make clean)
# the touch gets rid of a tar warning
	touch $(RELEASE_FILE)
	tar cvz --exclude="*~" --exclude=$(RELEASE_FILE) --exclude=mach --exclude-vcs -f $(RELEASE_FILE) ../statgen

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
