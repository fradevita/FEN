include $(REO_DIR)/src/Makefile.defs

override CNNFLAGS += -DMF

code.e: $(SOURCE).o $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(SOURCE).o $(OBJ) $(LIBS)

$(SOURCE).o: $(SOURCE).f90 precision.o constants.o class_Grid.o volume_of_fluid.o
	$(FC) $(INCLUDE_2DECOMP) $(FFLAGS) -c $<
