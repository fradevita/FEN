include $(FEN_DIR)/src/Makefile.defs

DBG += -Wall -g -fcheck=all
override CNNFLAGS += -DMF

code.e: $(SOURCE).o $(OBJ)
	$(FC) $(FLAGS) -o $@ $(SOURCE).o $(OBJ) $(LIBS)

$(SOURCE).o: $(SOURCE).f90 precision.o constants.o class_Grid.o volume_of_fluid.o multiphase.o navier_stokes.o io.o solver.o json.o
	$(FC) $(INCLUDE_2DECOMP) $(FLAGS) -c $<
