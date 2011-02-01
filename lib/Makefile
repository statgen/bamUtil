PATH_TO_BASE=..
include $(PATH_TO_BASE)/Makefile.include

SUBDIRS=general bam fastq glf
CLEAN_SUBDIRS= $(patsubst %, %_clean, $(SUBDIRS))

# Build in all subdirectories.
#
# see http://www.gnu.org/software/make/manual/make.html#Phony-Targets
# for a way of improving the following:
#
.PHONY: $(SUBDIRS) all test clean $(CLEAN_SUBDIRS)
all: TARGET = all
test: TARGET = test
clean: TARGET = clean

all test: tclap $(SUBDIRS)

clean: tclap_clean samtools_clean $(CLEAN_SUBDIRS)
	rm -f libStatGen.a

# other subdirectories depend on general
bam fastq glf: general

$(CLEAN_SUBDIRS):  
	@$(MAKE) $(PARALLEL_MAKE) OPTFLAG="$(OPTFLAG)" -C $(patsubst %_clean,%,$@) $(TARGET)

$(SUBDIRS): samtools
	@$(MAKE) $(PARALLEL_MAKE) OPTFLAG="$(OPTFLAG)" -C $@ $(TARGET)


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

tclap_clean:
	rm -rf tclap-1.2.0 tclap

samtools: samtools-0.1.12a
	ln -s samtools-0.1.12a samtools

samtools-0.1.12a: samtools-0.1.12a.tar.bz2
	tar xvf samtools-0.1.12a.tar.bz2 
	@$(MAKE) $(PARALLEL_MAKE) OPTFLAG="$(OPTFLAG)" -C $@

samtools_clean:
	rm -rf samtools-0.1.12a samtools
