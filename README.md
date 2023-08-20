# Introduction
This is a program aim for solving a transport system by NEGF method.

# Directory Structure
There are three important directories in this project
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
    
    *  The Fortran codes in this folder are mostly used to test the subroutine in main code. To compile them, one should execute the command specified at the beginning of each code. For example,
        ```
        $ cd tools/
        $ head test_LocalPotential_RtoP.f90 
        program test_LocalPotential_RtoP
            ! To compile, enter the build/ directory, and type:
            ! gfortran ../tools/test_LocalPotential_RtoP.f90 -o test_LocalPotential_RtoP.x grid.o -lfftw3
            use grid
            implicit none
            integer, parameter :: N_x = 160, N_y = 170, Nr = 30
            real*8, parameter :: a_x = 10.D0, a_y = 8.D0, A = 600.D0, R = 2.D0

            real*8 :: sum
            real*8, allocatable :: V_real(:, :)
        $ cd ../build/
        $ gfortran ../tools/test_LocalPotential_RtoP.f90 -o test_LocalPotential_RtoP.x grid.o -lfftw3
        ```
        Then one can execute the program `test_LocalPotential_RtoP.x`
    
    * The Python codes are mostly used to plot graph to visualize the result. To execute the program, go to the proper directory containing the data you want to plot, and execute
        ```
        $ python PythonProgram.py
        ```
# Code structure
There are four important files in `src/` directory
1. main.f90

    The main program, I call it FIRST layer of the code. Most of the input/output are here. Unit conversion is also done here. I've tried my best to avoid any math formula appear here so that the MPI calculation part is kept as clean and tidy as possible.
1. negf.f90

    This is the SECOND layer of the code. Most of the subroutine are just a "take out" from main program, so that `main.f90` is not going to mess up. One should never try to use any of the subroutine as a tool because they are just too specific, if you insist, I promise that it would blow up in the end.

    Another thing worth notice is that all the physical input parameters (those parameters with unit, except array) are sent in to the subroutine by a derived type variable: atomic. This is to avoid hundreds of variable being passed in the subroutine.

1. grid.f90

    This is the THIRD layer of the code. Containing most of the subroutine dealing with basis transformation, construction of Hamiltonian, etc. These subroutines may be use as a tool in other code, provided that you know what they are really doing. Note that all the input parameters in `atomic` are passed from `negf.f90` separately, and they go back to original name (e.g. `atomic%LX` -> `LX`) in all the subroutine here (`grid.f90`). So don't get confused with the name in `main.f90`, where they represent parameters in conventional unit system (Angstrom, eV).

1. math_kernel.f90

    This is the THIRD layer of the code. Containing most of the subroutine that dealing with math formula (e.g. matrix inverse). These subroutines may be use as a tool in other code as well. The logic of the naming is similar to `grid.f90`.

# Compile
Here are some information about environment setting. To successfully compile the program. The system might need to have the following compiler/library installed.

1. gfortran: 

    GNU compiler is suggested for my code. For those who want to use PGI or Intel compiler, some modification is required, so I strongly suggest you to use GNU Fortran to compile it.


1. OpenMPI:

    Note that the openmpi should be installed with the default compiler being GNU Fortran mentioned above. One might check this by the following:
    ```
    $ gfortran --version
    GNU Fortran (GCC) 10.4.0
    Copyright (C) 2020 Free Software Foundation, Inc.
    This is free software; see the source for copying conditions.  There is NO
    warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    $ mpifort --version
    GNU Fortran (GCC) 10.4.0
    Copyright (C) 2020 Free Software Foundation, Inc.
    This is free software; see the source for copying conditions.  There is NO
    warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    ```
    That is, `gfortran --version` and `mpifort --version` should give the same result (this is because mpifort is just a wrapper of gfortran)

1. LAPACK & BLAS & FFTW3

    These libraries should be installed in either `/usr/local/lib` or `/usr/lib`. Also, FFTW3 library requires an include file `fftw3.f03`. Make sure that this file exists in the path specified in `makefile`:
    ```
    # For the file 'fftw3.f03'
    INCLUDE= -I/usr/local/include 
    ```
    If not, you should modify this line in `makefile` to the correct path before compile the code.