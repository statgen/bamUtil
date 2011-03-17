EXE=vcfCooker
SRCONLY = Main.cpp
USER_LIBS = libvcf/libvcf.a
USER_INCLUDES = -Ilibvcf

LIB_DIR=libvcf

maintgt: all

clean: myclean

myclean:
	$(MAKE) -C $(LIB_DIR) OPTFLAG="$(OPTFLAG)" --no-print-directory clean;

$(LIB_DIR)/libvcf.a:
	$(MAKE) -C $(LIB_DIR) OPTFLAG="$(OPTFLAG)" --no-print-directory;