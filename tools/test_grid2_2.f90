program test_grid2
    ! To compile, enter the build/ directory, and type:
    ! gfortran ../tools/test_grid2_2.f90 -o test_grid2_2.x grid.o -lfftw3
    use grid
    implicit none
    integer, parameter :: N_x = 50, N_y = 60, Nr_x = 24, Nr_y = 30, N_trial = 20
    real*8, parameter :: Lx = 10.D0, Ly = 8.D0, A = 5.D0, R = 2.D0    

    real*8 :: sum, ENCUT
    complex*16, allocatable :: V_real(:, :, :, :), V_reciprocal(:, :, :, :), V_new_diag(:, :) &
    ,subPotential(:, :)
    integer :: i, j, ii, jj, N, trial
    integer, allocatable :: nx_grid(:), ny_grid(:)

    include 'constant.f90'

    allocate(V_real(N_x, N_y, N_x, N_y), V_new_diag(N_x, N_y))
    allocate(V_reciprocal(-Nr_x/2: Nr_x/2, -Nr_y/2: Nr_y/2, -Nr_x/2: Nr_x/2, -Nr_y/2: Nr_y/2))

    ! Construct V_real
    do jj=1, N_y
        do ii=1, N_x
            do j=1, N_y
                do i=1, N_x
                    V_real(i, j, ii, jj) = &
                    complex(exp(cos(2 * pi * (real(i) / N_x)) * cos(2 * pi * (real(i) / N_x - real(ii) / N_x))) &
                    , exp(cos(2 * pi * (real(j) / N_y) - 0.3) * cos(2 * pi * (real(j) / N_y - real(jj) / N_y))))
                end do
            end do
        end do
    end do

    ! Perform RealToPlanewave transformation
    call NonLocalPotential_RealToPlanewave(V_real, V_reciprocal)

    open(unit=16, file="encut_test.dat")
    write(16, *) "   ENCUT                 Error"
    do trial=0, N_trial
        ENCUT = 20.D0 + real(trial) / N_trial * (200.D0 - 20.D0)

        ! Grid layout
        N = PlaneWaveBasis_construction_findsize(ENCUT / hartree, Lx, Ly)
        allocate(nx_grid(N), ny_grid(N), subPotential(N, N))
        ! print *, "N = ", N
        call PlaneWaveBasis_construction(ENCUT / hartree, Lx, Ly, nx_grid, ny_grid)

        do j=1, N
            do i=1, N
                subPotential(i, j) =&
                V_reciprocal(nx_grid(i), ny_grid(i), nx_grid(j), ny_grid(j))
            end do
        end do
        ! print *, "subPotential"
        ! call print_c_matrix(subPotential)

        ! Perform PlanewaveToReal transformation
        call Greensfunction_PlanewaveToReal(subPotential, nx_grid, ny_grid, V_new_diag)

        sum = 0.D0
        do j=1, N_y
            do i=1, N_x
                sum = sum + abs(V_real(i, j, i, j) - V_new_diag(i, j)) ** 2
            end do
        end do
        write(16, *) ENCUT, sqrt(sum)

        deallocate(nx_grid, ny_grid, subPotential)
    end do

end program test_grid2