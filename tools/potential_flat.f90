program potential_flat
    ! To compile, enter the build/ directory, and type:
    ! gfortran ../tools/potential_flat.f90 -o potential_flat.x
    implicit none
    integer, parameter :: N_x = 40, N_y = 48, N_z = 301

    real*8, allocatable :: V_real(:, :, :)
    integer :: i, j, k

    allocate(V_real(N_x, N_y, N_z))
    V_real(:, :, :) = 0.D0

    open(unit=17, file="POTENTIAL")
    write(17, '(3I5)') N_x, N_y, N_z
    write(17, '(5G17.8)') (((V_real(i, j, k), i=1, N_x), j=1, N_y), k=1, N_z)

end program potential_flat