include $(FEN_DIR)/src/Makefile.defs

DBG += -Wall -g -fcheck=all

code.e: $(SOURCE).o $(OBJ)
	$(FC) $(FLAGS) -o $@ $(SOURCE).o $(OBJ) $(LIBS)

$(SOURCE).o: $(SOURCE).f90 lagrangian_solid.o precision.o constants.o
	$(FC) $(INCLUDE_2DECOMP) $(FLAGS) -c $<
