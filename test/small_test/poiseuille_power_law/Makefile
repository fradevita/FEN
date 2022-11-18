include $(FEN_DIR)/src/Makefile.defs

override CNNFLAGS += -DNN

code.e: $(SOURCE).o $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(SOURCE).o $(OBJ) $(LIBS)

$(SOURCE).o: $(SOURCE).f90 precision.o constants.o class_Grid.o fields.o navier_stokes.o solver.o
	$(FC) $(INCLUDE_2DECOMP) $(FFLAGS) -c $<
