include ../Makefile.src

LIB_DIR=libvcf

clean: myclean

myclean:
	$(MAKE) $(PARALLEL_MAKE) -C $(LIB_DIR) OPTFLAG="$(OPTFLAG)" --no-print-directory clean;

$(LIB_DIR)/libvcf.a:
	$(MAKE) $(PARALLEL_MAKE) -C $(LIB_DIR) OPTFLAG="$(OPTFLAG)" --no-print-directory;