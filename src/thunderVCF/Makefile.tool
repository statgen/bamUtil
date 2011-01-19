EXE=thunderVCF
TOOLBASE = ShotgunHaplotyper ShotgunManners
SRCONLY = Main.cpp
USER_LIBS = libmach/libmach.a ../vcfCooker/libvcf/libvcf.a
USER_INCLUDES = -Ilibmach -I../vcfCooker/libvcf -I../include
