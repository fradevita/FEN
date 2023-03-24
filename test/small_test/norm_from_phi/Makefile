include $(FEN_DIR)/src/Makefile.defs

override DBG += -Wall -fcheck=all -g

LOBJ = precision.o class_Grid.o class_Scalar.o functions.o io.o halo.o constants.o polynomial_reconstruction.o

code.e: main.o $(LOBJ)
	$(FC) $(FLAGS) -o $@ main.o $(LOBJ) $(LIBS)

main.o: main.f90 $(LOBJ)
	$(FC) $(FLAGS) $(INCLUDE) -c $<
