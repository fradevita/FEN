# TEST DESCRIPTION

This directory contains tests for the incompressible single-phase Navier-Stokes solver.

## lid_driven

Test for the lid driven cavity. The test compares velocity profiles with solution 
given by *Ghia et al, Journal of computational physics, 48, 1982*.

## poiseuille

Flow between two parallel walls in a periodic domain with applied body force. The test 
compares the solution with the analytical reference solution.

## poiseuille_io

Flow between two parallel walls with inflow velocity profile. The test 
compares the solution with the analytical reference solution.

## taylor_green_vortex

The test evaluates the time decaying solution of a taylor green vortex.
The test compares velocity and pressure solution at a given time with 
analytical reference solution.

