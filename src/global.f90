module global
    implicit none
    
    ! Derived type definition
    type :: t_parameters
        real*8 :: V_L, V_R, MU, ETA, TEMPERATURE, LX, LY, LZ, ENCUT, GAP, VDS, TRANSMISSION_GRID(2)
    end type
    type :: t_timer
        real :: start, end, sum = 0.0
    end type
    type :: t_kpointmesh
        real*8 :: kx, ky, weight
    end type
    type :: t_transmission
        real*8 :: energy, tau
    end type
    type t_gfunc
        complex*16 :: diagonal, first_column, last_column
    end type


    ! Unit converter, and some physical constant
    !   ( Reference: https://physics.nist.gov/cuu/Constants/ )
    real*8, parameter :: hartree = 27.211386245988D0 ! (in eV)
    !k_B = 8.617333262D-5 is Boltzmann constant (in eV / K)
    real*8, parameter :: k_B = 8.617333262D-5 / hartree ! Boltzmann constant (in atomic unit)
    real*8, parameter :: a_0 = 0.529177210903D0 ! Bohr radius (in Angstrom)
    real*8, parameter :: pi = 3.141592653589793D0

    ! DEBUG switcher
    logical, parameter :: DEBUG = .false.

end module global