include $(FEN_DIR)/src/Makefile.defs

override CNNFLAGS += -DMF
override OPT += -O3

code.e: $(SOURCE).o $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(SOURCE).o $(OBJ) $(LIBS)

$(SOURCE).o: $(SOURCE).f90 precision.o constants.o class_Grid.o volume_of_fluid.o multiphase.o io.o navier_stokes.o solver.o
	$(FC) $(INCLUDE_2DECOMP) $(FFLAGS) -c $<
