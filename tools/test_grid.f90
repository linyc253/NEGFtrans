program test_grid
    ! To compile, enter the build/ directory, and type:
    ! gfortran ../tools/test_grid.f90 -o test_grid.x grid.o -lfftw3
    use grid
    implicit none
    integer, parameter :: N_x = 50, N_y = 60, Nr = 50
    real*8, parameter :: Lx = 10.D0, Ly = 8.D0, A = 5.D0, R = 2.D0, ENCUT = 180.D0

    real*8 :: sum, trans_factor
    real*8, allocatable :: V_real(:, :), V_new_real(:, :)
    complex*16, allocatable :: V_reciprocal(:, :), V_new(:, :), subPotential(:, :)
    integer :: i, j, N
    integer, allocatable :: nx_grid(:), ny_grid(:)

    allocate(V_real(N_x, N_y), V_new_real(N_x, N_y))
    allocate(V_reciprocal(-Nr/2: Nr/2, -Nr/2: Nr/2), V_new(N_x, N_y))

    ! Construct V_real
    V_real(:, :) = 0.D0
    do i=1, N_x
        do j=1, N_y
            if ((Lx * i / N_x - 5) ** 2 + (Ly * j / N_y - 4) ** 2 <= R ** 2) then
                V_real(i, j) = A / R ** 2 * (&
                (Lx * i / N_x - 5) ** 2 + (Ly * j / N_y - 4) ** 2) - A
            end if
        end do
    end do

    ! Grid layout
    call LocalPotential_RealToPlanewave(V_real, V_reciprocal)
    N = PlaneWaveBasis_construction_findsize(ENCUT, Lx, Ly)
    allocate(nx_grid(N), ny_grid(N), subPotential(N, N))
    print *, N
    call PlaneWaveBasis_construction(ENCUT, Lx, Ly, nx_grid, ny_grid)

    do j=1, N
        do i=1, N
            subPotential(i, j) =&
             V_reciprocal(nx_grid(i) - nx_grid(j), ny_grid(i) - ny_grid(j))
        end do
    end do
    print *, "subPotential"
    call print_c_matrix(subPotential)

    call Greensfunction_PlanewaveToReal(subPotential, nx_grid, ny_grid, V_new)

    sum = 0.D0
    do j=1, N_y
        do i=1, N_x
            V_new_real(i, j) = real(V_new(i, j))
        end do
    end do

    trans_factor = maxval(V_new_real)
    do j=1, N_y
        do i=1, N_x
            V_new_real(i, j) = V_new_real(i, j) - trans_factor
            sum = sum + (V_real(i, j) - V_new_real(i, j)) ** 2
        end do
    end do
    print *, "Error =", sqrt(sum)
    print *, "Min(V_real) =", minval(V_real)
    print *, "Min(V_new_real) =", minval(V_new_real)
    open(unit=16, file="slice.dat")
    write(16, *) "   X                 V_real             V_new_real"
    do i=1, N_x
        write(16, '(3(G12.5,7X))') i * Lx / N_x, V_real(i, N_y / 2), V_new_real(i, N_y / 2)
    end do

end program test_grid