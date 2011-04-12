EXE=mergeFilterStat
TOOLBASE=FilterStat
SRCONLY = Main.cpp
USER_LIBS = ../vcfCooker/libvcf/libvcf.a /usr/lib/libboost_thread.a -lpthread
USER_INCLUDES = -I../vcfCooker/libvcf
