SUBDIRS = src
#
#.PHONY: subdirs $(SUBDIRS)
#
#%: $(SUBDIRS)
#
#$(SUBDIRS):
#	$(MAKE) -C $@
#
#all:

PARENT_MAKE := Makefile.tool
include Makefile.inc

#%:
#	@$(MAKE) OPTFLAG="$(OPTFLAG)" -C src $@
#include src/Makefile
