include ../Makefile.src

LIB_DIR=../vcfCooker/libvcf

clean: myclean

myclean:
	$(MAKE) -C $(LIB_DIR) OPTFLAG="$(OPTFLAG)" --no-print-directory clean;

$(LIB_DIR)/libvcf.a:
	$(MAKE) -C $(LIB_DIR) OPTFLAG="$(OPTFLAG)" --no-print-directory;