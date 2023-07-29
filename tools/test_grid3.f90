program test_grid3
    ! To compile, enter the build/ directory, and type:
    ! gfortran ../tools/test_grid3.f90 -o test_grid3.x grid.o -lfftw3
    use grid
    implicit none
    integer, parameter :: N_x = 50, N_y = 60, Nr_x = 24, Nr_y = 30
    real*8, parameter :: Lx = 10.D0, Ly = 8.D0, A = 5.D0, R = 2.D0   

    real*8 :: sum
    real*8, allocatable :: V_loc_real(:, :)
    complex*16, allocatable :: V_real(:, :, :, :), V_reciprocal(:, :, :, :), V_loc_reciprocal(:, :)
    integer :: i, j
    integer, allocatable :: nx_grid(:), ny_grid(:)

    include 'constant.f90'

    allocate(V_real(N_x, N_y, N_x, N_y), V_loc_real(N_x, N_y))
    allocate(V_reciprocal(-Nr_x/2: Nr_x/2, -Nr_y/2: Nr_y/2, -Nr_x/2: Nr_x/2, -Nr_y/2: Nr_y/2))
    allocate(V_loc_reciprocal(-Nr_x/2: Nr_x/2, -Nr_y/2: Nr_y/2))

    ! Construct V_real, V_loc_real
    V_real(:, :, :, :) = 0.D0
    do j=1, N_y
        do i=1, N_x
            if ((Lx * i / N_x - 5) ** 2 + (Ly * j / N_y - 4) ** 2 <= R ** 2) then
                V_loc_real(i, j) = A / R ** 2 * (&
                (Lx * i / N_x - 5) ** 2 + (Ly * j / N_y - 4) ** 2) - A
                V_real(i, j, i, j) = V_loc_real(i, j)
            end if
        end do
    end do

    ! Perform RealToPlanewave transformation with two different method
    call LocalPotential_RealToPlanewave(V_loc_real, V_loc_reciprocal)
    call NonLocalPotential_RealToPlanewave(V_real, V_reciprocal)
    
    sum = 0.D0
    do j=-Nr_y / 2, Nr_y / 2
        do i=-Nr_x / 2, Nr_x / 2
            V_loc_reciprocal(i, j) = V_loc_reciprocal(i, j) / (N_x * N_y)
            sum = sum + abs(V_reciprocal(i, j, 0, 0) - V_loc_reciprocal(i, j)) ** 2
        end do
    end do

    print *, "Error =", sqrt(sum)
    open(unit=16, file="grid3_slice.dat")
    write(16, *) "    X               V_reciprocal      V_loc_reciprocal"

    do j=-Nr_y / 2, Nr_y / 2
        write(16, '(3(G12.5,7X))') j * Ly / N_y, real(V_reciprocal(0, j, 0, 0)), real(V_loc_reciprocal(0, j))
    end do

end program test_grid3