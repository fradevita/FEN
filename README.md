# FEN

FEN = Fortran Environment for Numerics. It is a Fortran environment for performing numerical simulations. It is developed for teaching and research activites but with a focus on HPC performances.

## Overview

- Finite difference solver for the incompressible Navier-Stokes equations
- Fast Direct Solver for the solution of the Poisson equation
- Volume of Fluid solver (MTHINC) for multiphase simulations
- Direct forcing Immersed Boundary Method for simulation of flow around solid objects.
- Computaional parallelism by MPI

### Citing

Please kindly cite this publication if the solver helps for yout research
~~~
@article{de2021fully,
  title={A fully Eulerian solver for the simulation of multiphase flows with solid bodies: Application to surface gravity waves},
  author={De Vita, Francesco and De Lillo, Filippo and Verzicco, Roberto and Onorato, Miguel},
  journal={Journal of computational physics},
  volume={438},
  pages={110355},
  year={2021},
  publisher={Elsevier}
}
~~~

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

Possible missing packages: *libblas-dev*, *liblapack-dev*, *gfortran*, *libopenmpi-dev*

## Acknowledgments

FEN is the combinations of several methods / technique that I have learned during past years. The main sources of inspiration are the two open source solvers [Gerris](https://gfs.sourceforge.net/wiki/index.php/Main_Page) and [Basilisk](https://basilisk.fr) and the code used at KTH during my PostDoc, for which many contributions have been given by Pedro Simeon Costa and Marco Edoardo Rosti.
