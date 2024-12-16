include $(FEN_DIR)/src/Makefile.defs

override DBG+= -Wall -fcheck=all -g 
ifeq "$(SOURCE)" "methods.f90"
    override CPPDEFS+= -DMPI -DDIM=3
endif
