#SUBDIRS = src
#
#.PHONY: subdirs $(SUBDIRS)
#
#%: $(SUBDIRS)
#
#$(SUBDIRS):
#	$(MAKE) -C $@
#
#all:
%:
	@$(MAKE) OPTFLAG="$(OPTFLAG)" -C src $@
#include src/Makefile
