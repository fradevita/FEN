# TEST DESCRIPTION

This directory contains tests for the lagrangian solid solver (mass spring network).
The solver can be used both for rigid and deformable solids.

## Cantilever

Clamped filament with a preload. The test compare position in time of the free edge
with solution given in *Heltai et al, Comput. Methods Appl. Mech. Engrg. 316, 2017*.

## cylinder_Re20

Flow around a cylinder at Reynolds numbe Re = 20. The test compares drag force
with solution available [here](https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark1_re20.html).

## Energy

Clamped filament with external load. The test compare the mechanical energy of the solid
with the external work. The objective is to verify energy conservation of the solver.

## Euler-Bernoulli

Comparison with the Euler-Bernoulli solution and convergence order of the method.

## Falling ellipse

An ellipse falling in a closed rectangular box. The test compares center of mass
position, velocity and rotation with solution given in *Xia et al, Fournal of Fluid Mechanics, 625, 2009*.

## Free-filament

Clamped filament oscillating under weigth. The test compare position in time of the
free edge with solution given in *Huang et al, Journal of Computational Physics, 226, 2007*.

## Memory

Check that destroying a lagrangian solid correctly free the requested memory.

## Pan

Migration of cylinder in a periodic Poiseuille flow. The test compares cylinder center of 
mass evolution with solution given in *Pan and Glowinski, Journal of Computational Physics, 181, 2002*.

