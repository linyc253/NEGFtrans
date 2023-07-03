program test_LocalPotential
    ! To compile, enter the build/ directory, and type:
    ! gfortran ../tools/test_LocalPotential.f90 -o test_LocalPotential.x grid.o -lfftw3
    use grid
    implicit none
    integer, parameter :: N_x = 160, N_y = 170, N_z = 100
    real*8, parameter :: a_x = 10.D0, a_y = 8.D0, a_z = 6.D0, A = 600.D0, R = 2.D0

    real*8 :: sum
    real*8, allocatable :: V_real(:, :, :)
    complex*16, allocatable :: V_reciprocal(:, :, :)
    integer :: i, j, k

    allocate(V_real(N_x, N_y, N_z))
    allocate(V_reciprocal(N_x, N_y, N_z))

    ! Construct V_real
    V_real(:, :, :) = 0.D0
    sum = 0.D0
    do i=1, N_x
        do j=1, N_y
            do k=1, N_z
                if ((a_x * i / N_x - 5) ** 2 + (a_y * j / N_y - 4) ** 2 + (a_z * k / N_z - 3) ** 2 <= R ** 2) then
                    V_real(i, j, k) = A / R ** 2 * (&
                    (a_x * i / N_x - 5) ** 2 + (a_y * j / N_y - 4) ** 2 + (a_z * k / N_z - 3) ** 2) - A
                end if
                sum = sum + V_real(i, j, k)
            end do
        end do
    end do

    !print *, sum / (N_x * N_y * N_z)

    call LocalPotential_RealToPlanewave(sum / (N_x * N_y * N_z))
    
end program test_LocalPotential