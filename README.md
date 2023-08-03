# Introduction
This is a program aim for solving a transport system by NEGF method.

# Structure
There are three important folders in this project
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
# Compile