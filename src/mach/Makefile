TEST_SUBDIRS=
SUBDIRS="libmach mach1 thunder"

# Build in all subdirectories.
#
# see http://www.gnu.org/software/make/manual/make.html#Phony-Targets
# for a way of improving the following:
#

all: 
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

test: 
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
