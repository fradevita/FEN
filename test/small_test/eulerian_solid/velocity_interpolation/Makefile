include $(FEN_DIR)/src/Makefile.defs

override CPPDEFS += -DIBM -DMPI
override DBG += -Wall -fcheck=all -g

test_solid.o: test_solid.f90 precision.o eulerian_solid.o
	$(FC) $(FLAGS) -c $< -o $(ODIR)/$@
