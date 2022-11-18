include $(FEN_DIR)/src/Makefile.defs

override DBG += -Wall -g -fcheck=all

code.e: lid_driven.o $(OBJ)
	$(FC) $(FLAGS) -o $@ lid_driven.o $(OBJ) $(LIBS)

lid_driven.o: lid_driven.f90 precision.o constants.o class_Grid.o class_Scalar.o class_Vector.o fields.o navier_stokes.o solver.o io.o json.o
	$(FC) $(INCLUDE_2DECOMP) $(FLAGS) -c $<
