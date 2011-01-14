include ../Makefile.src

LIB1_DIR=../vcfCooker/libvcf

LIB2_DIR=libmach

clean: myclean

myclean:
	$(MAKE) $(PARALLEL_MAKE) -C $(LIB1_DIR) OPTFLAG="$(OPTFLAG)" --no-print-directory clean;
	$(MAKE) $(PARALLEL_MAKE) -C $(LIB2_DIR) OPTFLAG="$(OPTFLAG)" --no-print-directory clean;

$(LIB1_DIR)/libvcf.a:
	$(MAKE) $(PARALLEL_MAKE) -C $(LIB1_DIR) OPTFLAG="$(OPTFLAG)" --no-print-directory;

$(LIB2_DIR)/libmach.a:
	$(MAKE) $(PARALLEL_MAKE) -C $(LIB2_DIR) OPTFLAG="$(OPTFLAG)" --no-print-directory;