module grid
    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'
    private
    public LocalPotential_RealToPlanewave, PlaneWaveBasis_construction_findsize, &
    Greensfunction_PlanewaveToReal, PlaneWaveBasis_construction, sub_Hamiltonian, &
    print_c_matrix, NonLocalPotential_RealToPlanewave
contains
    subroutine LocalPotential_RealToPlanewave(V_real, V_reciprocal)
        ! Transform a local potential in real space into the matrix element in plane wave basis
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

    subroutine Greensfunction_PlanewaveToReal(subGreenFunc, nx_grid, ny_grid, subGreenDiag)
        ! Transform a Green's function in plane wave basis into real space, output diagonal elements only
        integer, intent(in) :: nx_grid(:), ny_grid(:)
        complex*16, intent(in) :: subGreenFunc(:, :)
        complex*16, intent(out) :: subGreenDiag(:, :)

        integer :: i, j, ii, jj, N, N_x, N_y
        complex*16, allocatable :: GreenTensor(:, :, :, :)
        complex*16, allocatable :: work(:, :), work_out(:, :)
        type(C_PTR) :: plan, iplan

        include 'constant.f90'

        N = size(nx_grid)
        N_x = size(subGreenDiag, 1)
        N_y = size(subGreenDiag, 2)
        allocate(GreenTensor(N_x, N_y, N_x, N_y))
        allocate(work(N_x, N_y), work_out(N_x, N_y))
        
        ! Construct GreenTensor
        GreenTensor(:, :, :, :) = complex(0.D0, 0.D0)
        do j=1, N
            do i=1, N
                GreenTensor(nx_grid(i) + N_x / 2, ny_grid(i) + N_y / 2, &
                 nx_grid(j) + N_x / 2, ny_grid(j) + N_y / 2)&
                 = subGreenFunc(i, j)
            end do
        end do

        plan = fftw_plan_dft_2d(N_y, N_x, work, work_out, FFTW_FORWARD, FFTW_ESTIMATE)
        iplan = fftw_plan_dft_2d(N_y, N_x, work, work_out, FFTW_BACKWARD, FFTW_ESTIMATE)

        ! Perform FFT on second index
        do j=1, N_y
            do i=1, N_x

                do jj=1, N_y
                    do ii=1, N_x
                        work(ii, jj) = GreenTensor(i, j, ii, jj)
                    end do
                end do
                call fftw_execute_dft(plan, work, work_out)
                do jj=1, N_y
                    do ii=1, N_x
                        GreenTensor(i, j, ii, jj) = work_out(ii, jj)
                    end do
                end do
                
            end do
        end do

        ! Perform IFFT on first index
        do jj=1, N_y
            do ii=1, N_x

                do j=1, N_y
                    do i=1, N_x
                        work(i, j) = GreenTensor(i, j, ii, jj)
                    end do
                end do
                call fftw_execute_dft(iplan, work, work_out)
                subGreenDiag(ii, jj) = work_out(ii, jj)

            end do
        end do

    end subroutine Greensfunction_PlanewaveToReal

    subroutine NonLocalPotential_RealToPlanewave(V_real, V_reciprocal)
        ! Transform a non-local potential in real space into the matrix element in plane wave basis
        complex*16, intent(in) :: V_real(:, :, :, :)
        complex*16, allocatable, intent(inout) :: V_reciprocal(:, :, :, :)
        !                ^                 ^
        ! Becuase we want to pass the boundary of V_reciprocal to this subroutine

        integer :: i, j, ii, jj, N_x, N_y, Nr_x, Nr_y, STATUS
        real*8 :: N_prod
        complex*16, allocatable :: V_work(:, :, :, :)
        complex*16, allocatable :: work(:, :), work_out(:, :)
        type(C_PTR) :: plan, iplan

        include 'constant.f90'

        N_x = size(V_real, 1)
        N_y = size(V_real, 2)
        allocate(V_work(N_x, N_y, N_x, N_y))
        allocate(work(N_x, N_y), work_out(N_x, N_y))

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
        do jj=1, N_y
            do ii=1, N_x
                do j=1, N_y
                    do i=1, N_x
                        V_work(i, j, ii, jj) = V_real(i, j, ii, jj) * (-1) ** (i + j - ii - jj)
                    end do
                end do
            end do
        end do

        plan = fftw_plan_dft_2d(N_y, N_x, work, work_out, FFTW_FORWARD, FFTW_ESTIMATE)
        iplan = fftw_plan_dft_2d(N_y, N_x, work, work_out, FFTW_BACKWARD, FFTW_ESTIMATE)

        ! Perform IFFT on second index
        do j=1, N_y
            do i=1, N_x

                do jj=1, N_y
                    do ii=1, N_x
                        work(ii, jj) = V_work(i, j, ii, jj)
                    end do
                end do
                call fftw_execute_dft(iplan, work, work_out)
                do jj=1, N_y
                    do ii=1, N_x
                        V_work(i, j, ii, jj) = work_out(ii, jj)
                    end do
                end do
                
            end do
        end do

        ! Perform FFT on first index
        do jj=1, N_y
            do ii=1, N_x

                do j=1, N_y
                    do i=1, N_x
                        work(i, j) = V_work(i, j, ii, jj)
                    end do
                end do
                call fftw_execute_dft(plan, work, work_out)
                do j=1, N_y
                    do i=1, N_x
                        V_work(i, j, ii, jj) = work_out(i, j)
                    end do
                end do    

            end do
        end do

        ! Devide by ((N_x * N_y) ** 2), and rearrange index
        N_prod = real((N_x * N_y) ** 2)
        do jj=-Nr_y / 2, Nr_y / 2
            do ii=-Nr_x / 2, Nr_x / 2
                do j=-Nr_y / 2, Nr_y / 2
                    do i=-Nr_x / 2, Nr_x / 2
                        V_reciprocal(i, j, ii, jj) = &
                        V_work(i + 1 + N_x / 2, j + 1 + N_y / 2, ii + 1 + N_x / 2, jj + 1 + N_y / 2) / N_prod
                    end do
                end do
            end do
        end do
    end subroutine NonLocalPotential_RealToPlanewave

    subroutine sub_Hamiltonian(nx_grid, ny_grid, Lx, Ly, kx, ky, V_reciprocal, Hamiltonian)
        ! Construct the Hamiltonian in transverse direction (xy-direction) for each z
        integer, intent(in) :: nx_grid(:), ny_grid(:)
        complex*16, allocatable, intent(in) :: V_reciprocal(:, :)
        real*8, intent(in) :: Lx, Ly, kx, ky
        complex*16, intent(out) :: Hamiltonian(:, :)

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

        integer :: i, j, nx_max, ny_max, nn

        include 'constant.f90'

        nn = 1
        nx_max = ceiling(Lx / pi * sqrt(ENCUT / 2))
        ny_max = ceiling(Ly / pi * sqrt(ENCUT / 2))
        do j=-ny_max, ny_max
            do i= -nx_max, nx_max
                if (((2 * i * pi / Lx) ** 2 + (2 * j * pi / Ly) ** 2) / 2 &
                 <= ENCUT) then
                    nx_grid(nn) = i
                    ny_grid(nn) = j
                    nn = nn + 1
                end if
            end do
        end do

    end subroutine PlaneWaveBasis_construction

    subroutine print_c_matrix(Matrix)
        complex*16, intent(in) :: Matrix(:, :)

        real*8, allocatable :: Matrix_abs(:, :)
        integer :: i, j, N, M
        N = size(Matrix, 1)
        M = size(Matrix, 2)

        allocate(Matrix_abs(N, M))
        do j=1, M
            do i=1, N
                Matrix_abs(i, j) = abs(Matrix(i, j))
            end do
        end do
        print *, "(", N, ",", M, ")"
        print *, maxval(Matrix_abs)
        print *, "=============================="
    end subroutine print_c_matrix

    subroutine print_c_tensor(Tensor)
        complex*16, intent(in) :: Tensor(:, :, :, :)

        real*8, allocatable :: Tensor_abs(:, :, :, :)
        integer :: i, j, ii, jj, N, M
        N = size(Tensor, 1)
        M = size(Tensor, 3)

        allocate(Tensor_abs(N, N, M, M))
        do j=1, M
            do i=1, N
                do jj=1, M
                    do ii=1, N
                        Tensor_abs(i, ii, j, jj) = abs(Tensor(i, ii, j, jj))
                    end do
                end do
            end do
        end do
        print *, "(", N, ",", N, ",", M, ",", M, ")"
        print *, maxval(Tensor_abs)
        print *, "=============================="
    end subroutine print_c_tensor
end module grid