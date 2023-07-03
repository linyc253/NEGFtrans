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
        integer :: N_x, N_y, i, j, N_prod
        type(C_PTR) :: plan

        N_x = size(V_real, 1)
        N_y = size(V_real, 2)
        N_prod = N_x * N_y
        allocate(work(N_x, N_y))
        allocate(work_out(N_x / 2 + 1, N_y))

        ! Times phase factor
        do i=1, N_x
            do j=1, N_y         
                work(i, j) = V_real(i, j) * (-1) ** (j - 1)
            end do
        end do


        ! Perform FFT
        plan = fftw_plan_dft_r2c_2d(N_y, N_x, work, work_out, FFTW_ESTIMATE)
        call fftw_execute_dft_r2c(plan, work, work_out)

        ! Post Process
        do i=0, N_x / 2
            do j=-N_y / 2, N_y / 2
                V_reciprocal(i, j) = work_out(i + 1, j + 1 + N_y / 2) / N_prod
            end do
        end do

        do i=-N_x / 2, -1
            do j=-N_y / 2, N_y / 2
                V_reciprocal(i, j) = conjg(V_reciprocal(-i, -j))
            end do
        end do
    end subroutine LocalPotential_RealToPlanewave

end module grid