module tools
    implicit none
    private
    public LDOS_size, LDOS_grid_converter
contains
    function LDOS_size(LDOS_GRID_value, N_i) result(NN)
        integer, intent(in) :: LDOS_GRID_value, N_i
        integer :: NN

        if(LDOS_GRID_value == 0) then
            NN = N_i
        else
            NN = 1
        end if
        
    end function LDOS_size

    subroutine LDOS_grid_converter(i, LDOS_GRID_value, N_i, ii, rescale_factor)
        integer, intent(in) :: i, LDOS_GRID_value, N_i
        integer, intent(out) :: ii
        real*8, intent(inout) :: rescale_factor

        if(LDOS_GRID_value == 0) then
            ii = i
        else
            ii = 1
        end if

        if(LDOS_GRID_value == 0) then
            rescale_factor = rescale_factor
        else if(LDOS_GRID_value == -1) then
            rescale_factor = rescale_factor / N_i
        else if(LDOS_GRID_value == i) then
            rescale_factor = rescale_factor
        else
            rescale_factor = 0.D0
        end if

    end subroutine LDOS_grid_converter

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

end module tools