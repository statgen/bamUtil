# Makefile that includes this should specify PARENT_MAKE to be
# the Makefile it should include from libStatGen/Makefiles/$(PARENT_MAKE)

#Update version as necessary
VERSION=1.0.14

CXX11_AVAIL ?= 1

CPP0X=$(shell if [ X$(CCVERSION) \> X4.3.0 ] ; then echo " -std=c++0x" ; fi)

ifeq ($(CXX11_AVAIL), 0)
CPP0X=
endif

CURRENT_PATH := $(dir $(lastword $(MAKEFILE_LIST)))
LIB_PATH_GENERAL ?=../libStatGen
LIB_PATH_BAM_UTIL ?= $(LIB_PATH_GENERAL)

# add any additional ../ as necessary if it is a relative path
ACTUAL_PATH := $(if $(filter-out /%,$(LIB_PATH_BAM_UTIL)), \
        $(addprefix $(CURRENT_PATH), $(LIB_PATH_BAM_UTIL)), \
        $(LIB_PATH_BAM_UTIL))
INCLUDE_MAKE := $(ACTUAL_PATH)/Makefiles/$(PARENT_MAKE)
include $(INCLUDE_MAKE)

.phony: cloneLib

$(INCLUDE_MAKE):
	@echo Unable to locate: $(INCLUDE_MAKE)
	@echo To change the location, set LIB_PATH_GENERAL or LIB_PATH_BAM_UTIL to the appropriate path to libStatGen.  Or specify \"make LIB_PATH_GENERAL=yourPath\" or \"make LIB_PATH_BAM_UTIL=yourPath\"
	@echo Use make cloneLib if you have git and want to clone the current libStatGen at that location.

cloneLib:
	@if test -d $(LIB_PATH_BAM_UTIL); \
	then echo $(LIB_PATH_BAM_UTIL) already exists; \
	else git clone git://github.com/statgen/libStatGen.git $(LIB_PATH_BAM_UTIL); fi
	@echo Call make to compile libStatGen and this tool.

