module grid
    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'
    private
    public LocalPotential_RealToPlanewave
contains
    subroutine LocalPotential_RealToPlanewave(V_real, V_reciprocal)
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

end module grid