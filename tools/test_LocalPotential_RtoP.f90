program test_LocalPotential_RtoP
    ! To compile, enter the build/ directory, and type:
    ! gfortran ../tools/test_LocalPotential_RtoP.f90 -o test_LocalPotential_RtoP.x grid.o -lfftw3
    use grid
    implicit none
    integer, parameter :: N_x = 160, N_y = 170, Nr = 30
    real*8, parameter :: a_x = 10.D0, a_y = 8.D0, A = 600.D0, R = 2.D0

    real*8 :: sum
    real*8, allocatable :: V_real(:, :)
    complex*16, allocatable :: V_reciprocal(:, :)
    integer :: i, j

    allocate(V_real(N_x, N_y))
    allocate(V_reciprocal(-Nr/2: Nr/2, -Nr/2: Nr/2))

    ! Construct V_real
    V_real(:, :) = 0.D0
    sum = 0.D0
    do i=1, N_x
        do j=1, N_y
            if ((a_x * i / N_x - 5) ** 2 + (a_y * j / N_y - 4) ** 2 <= R ** 2) then
                V_real(i, j) = A / R ** 2 * (&
                (a_x * i / N_x - 5) ** 2 + (a_y * j / N_y - 4) ** 2) - A
            end if
            sum = sum + V_real(i, j)
        end do
    end do

    print *, sum / (N_x * N_y)
    V_reciprocal(:, :) = 1.D0
    call LocalPotential_RealToPlanewave(V_real, V_reciprocal)
    print *, V_reciprocal(0, 0)
    print *, V_reciprocal(1, 1)
    print *, V_reciprocal(-1, -1)

end program test_LocalPotential_RtoP