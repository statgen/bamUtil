PATH_TO_BASE=..
include $(PATH_TO_BASE)/Makefile.include

TEST_SUBDIRS="general bam fastq glf"
SUBDIRS= "samtools " $(TEST_SUBDIRS)

# Build in all subdirectories.
#
# see http://www.gnu.org/software/make/manual/make.html#Phony-Targets
# for a way of improving the following:
#

all: tclap
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

test: all
	@for i in "$(TEST_SUBDIRS)"; do \
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

#
# from http://tclap.sourceforge.net/
#
tclap: tclap-1.2.0
	ln -s tclap-1.2.0 tclap

#
# tclap is header only - the tests are done using the
# host compiler, but no libraries are used, so no need
# to pass EXPORT_TOOLCHAIN
#
tclap-1.2.0: tclap-1.2.0.tar.gz
	tar xvzf tclap-1.2.0.tar.gz
	(cd tclap-1.2.0; ./configure; make)

clean:
	rm -rf tclap-1.2.0 tclap
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
