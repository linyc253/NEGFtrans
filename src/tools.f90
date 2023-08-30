module tools
    implicit none
    private
    public LDOS_size, LDOS_grid_converter, mpi_distribution, print_introduction
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

    subroutine LDOS_grid_converter(i, LDOS_GRID_value, delta_i, ii, rescale_factor)
        integer, intent(in) :: i, LDOS_GRID_value
        real*8, intent(in) :: delta_i
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
            rescale_factor = rescale_factor * delta_i
        else if(LDOS_GRID_value == i) then
            rescale_factor = rescale_factor
        else
            rescale_factor = 0.D0
        end if

    end subroutine LDOS_grid_converter

    function mpi_distribution(mpi_size, N_E, Nk) result(BEST_size_per_energy)
        integer, intent(in) :: mpi_size, N_E, Nk
        integer :: size_per_energy, BEST_size_per_energy, energy_per_grid, waste, min_waste, k_loop

        if(mpi_size > N_E * Nk) then
            BEST_size_per_energy = mpi_size / N_E 
        else
            min_waste = 1000000000
            do size_per_energy=1, min(Nk, mpi_size)
                energy_per_grid = min(mpi_size / size_per_energy, 1)
                k_loop = ceiling(real(Nk) / size_per_energy)
                waste = max(mpi_size - size_per_energy * energy_per_grid, mpi_size - N_E) * k_loop +&
                (k_loop * size_per_energy - Nk) * size_per_energy
                waste = waste * ceiling(real(N_E) / energy_per_grid)
                print *, size_per_energy, waste
                if(waste < min_waste) then
                    BEST_size_per_energy = size_per_energy
                    min_waste = waste
                end if
            end do
        end if

    end function

    subroutine print_introduction(unit)
        integer, intent(in) :: unit
        write(unit, *) "**********************************************************************"
        write(unit, *) " CALCULATES THE DENSITY AND TRANSMISSION OF A GIVEN POTENTIAL"
        write(unit, *) "   Using NEGF to calculate density and transmission coefficient with "
        write(unit, *) "   real grid in z-direction, and plane wave basis in xy-direction."
        write(unit, *) " input files:"
        write(unit, *) "   INPUT"
        write(unit, *) "   POTENTIAL"
        write(unit, *) " output files:"
        write(unit, *) "   OUTPUT"
        write(unit, *) "   DENSITY"
        write(unit, *) "   TRANSMISSION"
        write(unit, *) "   LDOS"
        write(unit, *) "                                   By YI-CHENG LIN at 2023/6/27"
        write(unit, *) "**********************************************************************"
        write(unit, *) ""
    end subroutine print_introduction
    
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