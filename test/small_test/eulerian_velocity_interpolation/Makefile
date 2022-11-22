include $(FEN_DIR)/src/Makefile.defs

override CNNFLAGS += -DIBM
override DBG += -Wall -fcheck=all -g
override OBJ += class_test_solid.o

code.e: $(SOURCE).o $(OBJ)
	$(FC) $(FLAGS) -o $@ $(SOURCE).o $(OBJ) $(LIBS)

$(SOURCE).o: $(SOURCE).f90 precision.o constants.o class_Grid.o class_Vector.o eulerian_ibm.o class_test_solid.o
	$(FC) $(INCLUDE_2DECOMP) $(FLAGS) -c $<

class_test_solid.o: class_test_solid.f90 class_eulerian_solid.o
	$(FC) $(FLAGS) -c $<
