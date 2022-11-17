# FEN

FEN = Fortran Environment for Numerics. It is as a Fortran environment for performing numerical simulations, mainly CFD.

## Overview

- Finite Difference solver for Navier-Stokes equations
- Fast Direct Solver for the solution of the Poisson equation
- Volume of Fluid solver (MTHIN) for multiphase simulations
- Computaional parallelism by MPI

## Installing

### Dependency

FEN uses the [2decomp&FFT](http://www.hector.ac.uk/cse/distributedcse/reports/incompact3d/UserGuide.html) library for the domain decomposition and [FFTW3](http://www.fftw.org/) library for the solution of the Poisson equation.

Add to your `$HOME/.bashrc` the following lines

``` (bash)
export FEN_DIR=...
export FFTW3_DIR=...
export _2DECOMP_DIR=...
```

and replace `...` with proper values for you system. Then type

``` (bash)
source $HOME/.bashrc
cd $FEN_DIR
sh INSTALL.sh
```

After this step the requested libraries should be properly installed.
