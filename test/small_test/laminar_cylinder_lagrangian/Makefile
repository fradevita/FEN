include $(FEN_DIR)/src/Makefile.defs

#override DBG += -Wall -g -fcheck=all
override CNNFLAGS += -DIBM
override OPT += -O3

code.e: $(SOURCE).o $(OBJ)
	$(FC) $(FLAGS) -o $@ $(SOURCE).o $(OBJ) $(LIBS)

$(SOURCE).o: $(SOURCE).f90 precision.o constants.o class_Grid.o class_Vector.o eulerian_ibm.o navier_stokes.o solver.o
	$(FC) $(INCLUDE_2DECOMP) $(FLAGS) -c $<
