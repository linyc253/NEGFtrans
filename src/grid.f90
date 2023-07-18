module grid
    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'
    private
    public LocalPotential_RealToPlanewave
contains
    subroutine LocalPotential_RealToPlanewave(V_real, V_reciprocal)
        ! Transform the potential in real space to the matrix element in plane wave basis
        real*8, intent(in) :: V_real(:, :)
        complex*16, allocatable, intent(inout) :: V_reciprocal(:, :)
        !                ^                 ^
        ! Becuase we want to pass the boundary of V_reciprocal to this subroutine

        real*8, allocatable :: work(:, :)
        complex*16, allocatable :: work_out(:, :)
        integer :: N_x, N_y, i, j, N_prod, Nr_x, Nr_y, STATUS
        type(C_PTR) :: plan

        N_x = size(V_real, 1)
        N_y = size(V_real, 2)
        N_prod = N_x * N_y
        allocate(work(N_x, N_y))
        allocate(work_out(N_x / 2 + 1, N_y))

        ! Check the form of V_reciprocal
        Nr_x = size(V_reciprocal, 1) - 1
        Nr_y = size(V_reciprocal, 2) - 1
        if((Nr_x > N_x) .or. (Nr_y > N_y)) then
            print *, "ERROR: 'The size of V_reciprocal minus one' should not be"
            print *, "       larger than 'the size of V_real'"
            call exit(STATUS)
        end if
        if((-lbound(V_reciprocal, 1) /= ubound(V_reciprocal, 1)) .or. &
        (-lbound(V_reciprocal, 2) /= ubound(V_reciprocal, 2))) then
            print *, "ERROR: The index boundary of V_reciprocal must be symmetric"
            print *, "       about zero"
            call exit(STATUS)
        end if

        ! Times phase factor
        do i=1, N_x
            do j=1, N_y         
                work(i, j) = V_real(i, j) * (-1) ** (j - 1)
            end do
        end do


        ! Perform FFT
        plan = fftw_plan_dft_r2c_2d(N_y, N_x, work, work_out, FFTW_ESTIMATE)
        call fftw_execute_dft_r2c(plan, work, work_out)

        ! Construct V_reciprocal
        do i=0, Nr_x / 2
            do j=-Nr_y / 2, Nr_y / 2
                V_reciprocal(i, j) = work_out(i + 1, j + 1 + N_y / 2) / N_prod
            end do
        end do

        do i=-Nr_x / 2, -1
            do j=-Nr_y / 2, Nr_y / 2
                V_reciprocal(i, j) = conjg(V_reciprocal(-i, -j))
            end do
        end do
    end subroutine LocalPotential_RealToPlanewave

    subroutine sub_Hamiltonian(nx_grid, ny_grid, Lx, Ly, kx, ky, V_reciprocal, Hamiltonian)
        ! Construct the Hamiltonian in transverse direction (xy-direction) for each z
        integer, intent(in) :: nx_grid(:), ny_grid(:)
        real*8, allocatable, intent(in) :: V_reciprocal(:, :)
        real*8, intent(in) :: Lx, Ly, kx, ky
        real*8, intent(out) :: Hamiltonian(:, :)

        integer :: i, j, N

        include 'constant.f90'

        N = size(nx_grid)

        do j=1, N
            do i=1, N
                Hamiltonian(i, j) =&
                 V_reciprocal(nx_grid(i) - nx_grid(j), ny_grid(i) - ny_grid(j))
            end do
        end do

        do i=1, N
            Hamiltonian(i, i) = Hamiltonian(i, i) + &
             ((2 * nx_grid(i) * pi / Lx + kx) ** 2 + (2 * ny_grid(i) * pi / Ly + ky) ** 2) / 2
        end do

    end subroutine sub_Hamiltonian

    function PlaneWaveBasis_construction_findsize(ENCUT, Lx, Ly) result(N)
        ! Calculate the size (N) of the planewave grid by ENCUT
        ! Use this before you call the subroutine "PlaneWaveBasis_construction"
        ! Then allocate the array by: allocate(nx_grid(N), ny_grid(N))
        real*8, intent(in) :: ENCUT, Lx, Ly
        integer :: N

        integer :: i, j, nx_max, ny_max

        include 'constant.f90'

        N = 0
        nx_max = ceiling(Lx / pi * sqrt(ENCUT / 2))
        ny_max = ceiling(Ly / pi * sqrt(ENCUT / 2))
        do j=-ny_max, ny_max
            do i= -nx_max, nx_max
                if (((2 * i * pi / Lx) ** 2 + (2 * j * pi / Ly) ** 2) / 2 &
                 <= ENCUT) then
                    N = N + 1
                end if
            end do
        end do

    end function PlaneWaveBasis_construction_findsize

    subroutine PlaneWaveBasis_construction(ENCUT, Lx, Ly, nx_grid, ny_grid)
        ! Construction of the planewave basis, limited by ENCUT
        real*8, intent(in) :: ENCUT, Lx, Ly
        integer, intent(out) :: nx_grid(:), ny_grid(:)

        integer :: i, j, nx_max, ny_max, n

        include 'constant.f90'

        n = 1
        nx_max = ceiling(Lx / pi * sqrt(ENCUT / 2))
        ny_max = ceiling(Ly / pi * sqrt(ENCUT / 2))
        do j=-ny_max, ny_max
            do i= -nx_max, nx_max
                if (((2 * i * pi / Lx) ** 2 + (2 * j * pi / Ly) ** 2) / 2 &
                 <= ENCUT) then
                    nx_grid(n) = i
                    ny_grid(n) = j
                    n = n + 1
                end if
            end do
        end do

    end subroutine PlaneWaveBasis_construction

end module grid