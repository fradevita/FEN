# TEST DESCRIPTION

This directory contains tests for the eulerian solid solver, used for simulation 
of rigid solids defined on a distance field.

## cylinder_Re20

Flow around a cylinder at Reynolds numbe Re = 20. The test compares drag force
with solution available [here](https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark1_re20.html).
Drag is evaluated both with external probes and with eulerian forcing intergral.

## Pan

Migration of cylinder in a periodic Poiseuille flow. The test compares cylinder center of 
mass evolution with solution given in *Pan and Glowinski, Journal of Computational Physics, 181, 2002*.

## tagging

Test to veriy the proper tagging of grid cells.

## velocity_interpolation

The test evaluates convergence order of interpolations subroutine in 2D and 3D.
Expected result is 2nd order convergence.

## Memory

Check that destroying a lagrangian solid correctly free the requested memory.
