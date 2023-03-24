include $(FEN_DIR)/src/Makefile.defs

override CNNFLAGS += -DIBM
override DBG += -Wall -fcheck=all -g

code.e: $(SOURCE).o $(OBJ)
	$(FC) $(FLAGS) -o $@ $(SOURCE).o $(OBJ) $(LIBS)

$(SOURCE).o: $(SOURCE).f90 precision.o constants.o class_Grid.o polynomial_reconstruction.o utils.o
	$(FC) $(INCLUDE_2DECOMP) $(FLAGS) -c $<
