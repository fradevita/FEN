include $(REO_DIR)/src/Makefile.defs

#override DBG += -Wall -fcheck=all -g
override CNNFLAGS += -DDIM=3
override OPT += -O3

code.e: $(SOURCE).o $(OBJ)
	$(FC) $(FLAGS) -o $@ $(SOURCE).o $(OBJ) $(LIBS)

$(SOURCE).o: $(SOURCE).f90 precision.o constants.o class_Grid.o fields.o navier_stokes.o solver.o json.o
	$(FC) $(INCLUDE_2DECOMP) $(FLAGS) -c $<
