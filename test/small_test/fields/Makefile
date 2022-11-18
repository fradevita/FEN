include $(FEN_DIR)/src/Makefile.defs

override DBG += -Wall -g -fcheck=all
override CNNFLAGS += -DDIM=3

code.e: $(SOURCE).o $(OBJ)
	$(FC) $(FLAGS) -o $@ $(SOURCE).o $(OBJ) $(LIBS)

$(SOURCE).o: $(SOURCE).f90 precision.o constants.o io.o class_Grid.o io.o class_Scalar.o class_Vector.o tensors.o
	$(FC) $(INCLUDE_2DECOMP) $(FLAGS) -c $<
