include $(FEN_DIR)/src/Makefile.defs

override DBG +=-g -Wall -fcheck=bounds -ffpe-trap=invalid,zero,overflow

code.e: $(SOURCE).o $(OBJ)
	$(FC) $(FLAGS) -o $@ $(SOURCE).o $(OBJ) $(LIBS)

$(SOURCE).o: $(SOURCE).f90 precision.o constants.o class_Scalar.o mls.o
	$(FC) $(INCLUDE_2DECOMP) $(FLAGS) -c $<
