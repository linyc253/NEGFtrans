program test_grid2
    ! To compile, enter the build/ directory, and type:
    ! gfortran ../tools/test_grid2.f90 -o test_grid2.x grid.o -lfftw3
    use grid
    implicit none
    integer, parameter :: N_x = 50, N_y = 60, Nr = 50
    real*8, parameter :: Lx = 10.D0, Ly = 8.D0, A = 5.D0, R = 2.D0, ENCUT = 180.D0

    real*8 :: sum
    complex*16, allocatable :: V_real(:, :, :, :), V_reciprocal(:, :, :, :), V_new_diag(:, :) &
    ,subPotential(:, :)
    integer :: i, j, ii, jj, N
    integer, allocatable :: nx_grid(:), ny_grid(:)

    include 'constant.f90'

    allocate(V_real(N_x, N_y, N_x, N_y), V_new_diag(N_x, N_y))
    allocate(V_reciprocal(-Nr/2: Nr/2, -Nr/2: Nr/2, -Nr/2: Nr/2, -Nr/2: Nr/2))

    ! Construct V_real
    do jj=1, N_y
        do ii=1, N_x
            do j=1, N_y
                do i=1, N_x
                    V_real(i, j, ii, jj) = &
                    exp(cos(2 * pi * (real(i) / N_x - N_x / 2)) * cos(2 * pi * (real(i) / N_x - real(ii) / N_x))) &
                    + exp(cos(2 * pi * (real(j) / N_y - N_y / 2)) * cos(2 * pi * (real(j) / N_y - real(jj) / N_y)))
                end do
            end do
        end do
    end do

    ! Grid layout
    call NonLocalPotential_RealToPlanewave(V_real, V_reciprocal)
    N = PlaneWaveBasis_construction_findsize(ENCUT, Lx, Ly)
    allocate(nx_grid(N), ny_grid(N), subPotential(N, N))
    print *, "N = ", N
    call PlaneWaveBasis_construction(ENCUT, Lx, Ly, nx_grid, ny_grid)

    do j=1, N
        do i=1, N
            subPotential(i, j) =&
             V_reciprocal(nx_grid(i), ny_grid(i), nx_grid(j), ny_grid(j))
        end do
    end do
    print *, "subPotential"
    call print_c_matrix(subPotential)

    call Greensfunction_PlanewaveToReal(subPotential, nx_grid, ny_grid, V_new_diag)

    sum = 0.D0
    do j=1, N_y
        do i=1, N_x
            sum = sum + abs(V_real(i, j, i, j) - V_new_diag(i, j)) ** 2
        end do
    end do

    print *, "Error =", sqrt(sum)
    open(unit=16, file="grid2_slice.dat")
    write(16, *) "   X                 V_real             V_new_diag"
    ! do i=1, N_x
    !     write(16, '(3(G12.5,7X))') i * Lx / N_x, real(V_real(i, N_y / 2, i, N_y / 2)), real(V_new_diag(i, N_y / 2))
    ! end do

    do j=1, N_y
        write(16, '(3(G12.5,7X))') j * Ly / N_y, real(V_real(N_x / 2, j, N_x / 2, j)), real(V_new_diag(N_x / 2, j))
    end do

end program test_grid2