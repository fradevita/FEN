include $(FEN_DIR)/src/Makefile.defs

override OPT += -O3
override CNNFLAGS += -DIBM

code.e: $(SOURCE).o $(OBJ)
	$(FC) $(FLAGS) -o $@ $(SOURCE).o $(OBJ) $(LIB_2DECOMP) $(LIB_FFTW3)

$(SOURCE).o: $(SOURCE).f90 precision.o constants.o class_Grid.o class_Vector.o navier_stokes.o solver.o ibm.o io.o json.o
	$(FC) $(INCLUDE_2DECOMP) $(FFLAGS) -c $<
