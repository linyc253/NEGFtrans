!**********************************************************************
! CALCULATES THE DENSITY AND TRANSMISSION OF A GIVEN POTENTIAL
!   Using NEGF to calculate density and transmission coefficient with 
!   real grid in z-direction, and plane wave basis in xy-direction.
! input files:
!   INPUT
!   POTENTIAL
! output files:
!   OUTPUT
!   DENSITY
!   TRANSMISSION
!   LDOS (optional)
!                                   By YI-CHENG LIN at 2023/6/27
!**********************************************************************
program main
    use anderson_mixing

    implicit none
    real :: V_L, V_R
    namelist /INPUT/ V_L, V_R

    ! Initialize Parameters
    V_L = 1.0
    V_R = 2.0

    ! Read Parameters
    open(unit=15, file="INPUT")
    read(15, INPUT)
    close(15)
    
    ! Print Parameters
    open(unit=16, file="OUTPUT")
    write(16, *) "V_L = ", V_L
    write(16, *) "V_R = ", V_R
    
end program main