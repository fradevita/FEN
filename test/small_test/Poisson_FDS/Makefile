include $(FEN_DIR)/src/Makefile.defs

code.e: $(SOURCE).o $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(SOURCE).o $(OBJ) $(LIBS)

$(SOURCE).o: $(SOURCE).f90 precision.o constants.o class_Grid.o class_Scalar.o Poisson.o
	$(FC) $(INCLUDE_2DECOMP) $(FFLAGS) -c $<
