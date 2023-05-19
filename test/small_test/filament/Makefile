include $(FEN_DIR)/src/Makefile.defs

override DBG += -Wall -g -fcheck=all

code.e: $(SOURCE).o $(OBJ)
	$(FC) $(FLAGS) -o $@ $(SOURCE).o $(OBJ) $(LIBS)

$(SOURCE).o: $(SOURCE).f90 precision.o constants.o class_Grid.o lagrangian_solid.o lagrangian_ibm.o navier_stokes.o io.o solver.o class_Vector.o fields.o 
	$(FC) $(INCLUDE_2DECOMP) $(FLAGS) -c $<
