# Introduction
The program solve the electron density, local density of state, and transmission coefficient of an open system using Non-Equilibrium Green's Function (NEGF) method. In particular, the program computes the quantities based on plane wave basis in transverse direction (xy), and real space grid in transport direction (z).

# Repository Structure
The structure of this repository is
```
├── src
│    ├── main.f90
│    ├── negf.f90
│    ├── grid.f90
│    ├── global.f90
│    ├── tools.f90
│    └── math_kernal.f90
├── build
│    └── makefile
├── doc
└── tools
```
1. src

    All the source code for main code will be contained in this folder.

1. build

    Only `makefile` will be traced by git in this folder, one should always compile the program here by executing
    ```
    $ cd build
    $ make
    ```
1. tools

    It's the folder that contain all the other codes that are not included in main code.
    

Here are some important files in `src/` directory
1. main.f90

    The main program, I call it FIRST layer of the code. Most of the input/output are here. Unit conversion is also done here. I've tried my best to avoid any math formula appear here so that the MPI calculation part is kept as clean and tidy as possible.
1. negf.f90

    This is the SECOND layer of the code. Most of the subroutine are just a "take out" from main program, so that `main.f90` is not going to mess up. One should never try to use any of the subroutine as a tool because they are just too specific.

    Another thing worth notice is that all the physical input parameters (those parameters with unit, except array) are sent in to the subroutine by a derived type variable: `atomic`. This is to avoid hundreds of variable being passed in the subroutine.

1. grid.f90

    This is the THIRD layer of the code. Containing most of the subroutine dealing with basis transformation, construction of Hamiltonian, etc. These subroutines may be use as a tool in other code, provided that you know what they are really doing. Note that all the input parameters in `atomic` are passed from `negf.f90` separately, and they go back to original name (e.g. `atomic%LX` -> `LX`) in all the subroutine here (`grid.f90`). So don't get confused with the name in `main.f90`, where they represent parameters in conventional unit system (Angstrom, eV).

1. math_kernel.f90

    This is the THIRD layer of the code. Containing most of the subroutine that dealing with math formula (e.g. matrix inverse). These subroutines may be use as a tool in other code as well. The logic of the naming is similar to `grid.f90`.

1. global.f90

    This file store the global variable (physical constants) and type definition.

# Installation
Here are some information about environment setting. To successfully compile the program. The system might need to have the following compiler/library installed.

1. gfortran: 

    GNU compiler is suggested for my code. For those who want to use PGI or Intel compiler, some modification of code is required, so I strongly suggest you to use GNU Fortran to compile it.


1. OpenMPI:

    Note that the openmpi should be installed with the default compiler being GNU Fortran mentioned above. One might check this by the following:
    ```
    $ gfortran --version
    GNU Fortran (GCC) x.x.x
    $ mpifort --version
    GNU Fortran (GCC) x.x.x
    ```
    That is, `gfortran --version` and `mpifort --version` should give the same result (this is because mpifort is just a wrapper of gfortran)

1. LAPACK & BLAS & FFTW3

    These libraries should be installed in either `/usr/local/lib` or `/usr/lib`. Also, FFTW3 library requires an include file `fftw3.f03`. Make sure that this file exists in the path specified in `makefile`:
    ```
    # For the file 'fftw3.f03'
    INCLUDE= -I/usr/local/include -I/usr/include
    ```
    If not, you should modify this line in `makefile` to the correct path before compile the code.

# Execution
The program can be executed by `mpirun`. For detailed input parameter setting, please refer to `doc/Manual.pdf`