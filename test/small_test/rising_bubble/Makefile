include $(FEN_DIR)/src/Makefile.defs

override DBG += -g -Wall -fcheck=all -ffpe-trap=invalid,zero,overflow
override CNNFLAGS += -DMF

code.e: $(SOURCE).o $(OBJ)
	$(FC) $(FLAGS) -o $@ $(SOURCE).o $(OBJ) $(LIBS)

$(SOURCE).o: $(SOURCE).f90 precision.o constants.o class_Grid.o fields.o navier_stokes.o solver.o json.o
	$(FC) $(INCLUDE_2DECOMP) $(FLAGS) -c $<
