module math_kernel
    use global
    implicit none
    private
    public GreensFunction_tri_solver, coefficient_simpson, value_simpson, fermi_func
    logical :: reset = .true.
contains
    subroutine inverse(Matrix, inv_Matrix)
        complex*16, intent(in) :: Matrix(:, :)
        complex*16, intent(out) :: inv_Matrix(:, :)
        complex*16, allocatable, save :: work(:)
        integer, allocatable, save :: IPIV(:)
        integer, save :: N, NB
        integer :: i, j, STATUS, ilaenv
        if(reset) then
            N = size(Matrix, 1)
            NB = ilaenv(1, "zgetri", "", N, N, -1, -1)
            if(NB < 1) NB = N
            if(allocated(IPIV)) deallocate(IPIV)
            if(allocated(work)) deallocate(work)
            allocate(IPIV(N), work(NB * N))
            reset = .false.
        end if

        do j=1, N
            do i=1, N
                inv_Matrix(i, j) = Matrix(i, j)
            end do
        end do
        
        call zgetrf(N, N, inv_Matrix, N, IPIV, STATUS)
        call zgetri(N, inv_Matrix, N, IPIV, work, NB * N, STATUS)
        
    end subroutine

    subroutine matrix_product(A, B, C)
        complex*16, intent(in) :: A(:, :), B(:, :)
        complex*16, intent(out) :: C(:, :)
        integer :: N
        N = size(A, 1)
        call zgemm('N', 'N', N, N, N, complex(1.D0, 0.D0), A, N, B, N, complex(0.D0, 0.D0), C, N)
    end subroutine

    subroutine GreensFunction_tri_solver(G_inv_D, G_Function, DFL)
        !********************************************************
        ! This is a triangular matrix inverse solver, see
        !       "Journal of Applied Physics 81, 7845 (1997)"
        ! for details (Recursive Green function algorithm)
        !
        ! Assume the off-diagonal blocks are identities
        !********************************************************
    
            implicit none
            complex*16, intent(in) :: G_inv_D(:, :, :)
            character(*), intent(in) :: DFL
            type(t_gfunc), intent(out) :: G_Function(:, :, :)
            ! Structure of G_Function:
            !   G_Function%diagonal     : diagonal element of the inverse matrix, calculated if "D" contained in DFL
            !   G_Function%first_column : first column element the inverse matrix, calculated if "F" contained in DFL
            !   G_Function%last_column  : last column element the inverse matrix, calculated if "L" contained in DFL

            complex*16, allocatable :: generator_l(:, :, :), generator_r(:, :, :), Matrix(:, :)
            logical :: if_D, if_F, if_L
            
            integer :: N_z, N, i, j, k

            N_z = size(G_inv_D, 3)
            N = size(G_inv_D, 1)
            allocate(generator_l(N, N, N_z - 1))
            allocate(generator_r(N, N, N_z - 1))
            allocate(Matrix(N, N))
            if_D = scan(DFL, "Dd") /= 0
            if_F = scan(DFL, "Ff") /= 0
            if_L = scan(DFL, "Ll") /= 0
            G_Function(:, :, :)%diagonal = 0.D0
            G_Function(:, :, :)%first_column = 0.D0
            G_Function(:, :, :)%last_column = 0.D0
    
    
            ! Calculate generator_r()
            if(if_D .or. if_F) then
                call inverse(G_inv_D(:, :, N_z), generator_r(:, :, N_z - 1))
                do k=2, N_z - 1
                    do j=1, N
                        do i=1, N
                            Matrix(i, j) = G_inv_D(i, j, N_z - k + 1) - generator_r(i, j, N_z - k + 1)
                        end do
                    end do
                    call inverse(Matrix, generator_r(:, :, N_z - k))
                end do
            end if
    
            ! Get G_Function%diagonal (i.e. diagonal element)
            if(if_D) then
                do j=1, N
                    do i=1, N
                        Matrix(i, j) = G_inv_D(i, j, 1) - generator_r(i, j, 1)
                    end do
                end do
                call inverse(Matrix, G_Function(:, :, 1)%diagonal)
            
                do k=2, N_z
                    call matrix_product(G_Function(:, :, k - 1)%diagonal, generator_r(:, :, k - 1), Matrix)
                    do i=1, N
                        Matrix(i, i) = Matrix(i, i) + 1.D0
                    end do
                    call matrix_product(generator_r(:, :, k - 1), Matrix, G_Function(:, :, k)%diagonal)
                end do
            end if
    
            ! Get G_Function%first_column (i.e. first column element)
            if(if_F) then
                do j=1, N
                    do i=1, N
                        Matrix(i, j) = G_inv_D(i, j, 1) - generator_r(i, j, 1)
                    end do
                end do
                call inverse(Matrix, G_Function(:, :, 1)%first_column)
                
                do k=2, N_z
                    call matrix_product(generator_r(:, :, k - 1), G_Function(:, :, k - 1)%first_column, &
                    G_Function(:, :, k)%first_column)
                    do j=1, N
                        do i=1, N
                            G_Function(i, j, k)%first_column = - G_Function(i, j, k)%first_column
                        end do
                    end do
                end do
            end if

            ! Calculate generator_l()
            if(if_L) then
                call inverse(G_inv_D(:, :, 1), generator_l(:, :, 1))
                do k=2, N_z - 1
                    do j=1, N
                        do i=1, N
                            Matrix(i, j) = G_inv_D(i, j, k) - generator_l(i, j, k - 1)
                        end do
                    end do
                    call inverse(Matrix, generator_l(:, :, k))
                end do
            end if
    
            ! Get G_Function%last_column (i.e. last column element)
            if(if_L) then
                do j=1, N
                    do i=1, N
                        Matrix(i, j) = G_inv_D(i, j, N_z) - generator_l(i, j, N_z - 1)
                    end do
                end do
                call inverse(Matrix, G_Function(:, :, N_z)%last_column)

                do k=1, N_z - 1
                    call matrix_product(generator_l(:, :, N_z - k), G_Function(:, :, N_z - k + 1)%last_column, &
                    G_Function(:, :, N_z - k)%last_column)
                    do j=1, N
                        do i=1, N
                            G_Function(i, j, N_z - k)%last_column = - G_Function(i, j, N_z - k)%last_column
                        end do
                    end do
                end do
            end if
            
            reset = .true.
        end subroutine GreensFunction_tri_solver

        function coefficient_simpson(i, N, a, b) result(coefficient)
            ! Provide the scalar factor (c_i) of simpson integral method
            ! \[\int_{a}^{b} f(x) dx = \sum_i c_i f(x_i)\]
            ! The index "i" start from 1 ~ N (Fortran convention)
            integer, intent(in) :: i, N
            real*8, intent(in) :: a, b
            real*8 :: coefficient

            if(mod(N, 2) == 1) then
                ! Simpson 1/3 rule apply to all
                if((i == 1) .or. (i == N)) then
                    coefficient = 1.D0 / 3
                else
                    if(mod(i, 2) == 0) coefficient = 4.D0 / 3
                    if(mod(i, 2) == 1) coefficient = 2.D0 / 3
                end if

            else
                ! Simpson 1/3 rule apply to 1 ~ (N-3)
                ! Simpson 3/8 rule apply to (N-3) ~ N
                if(i == 1) then
                    coefficient = 1.D0 / 3
                else if(i == N - 3) then
                    coefficient = 1.D0 / 3 + 1.D0 * 3 / 8
                else if((i == N - 2) .or. (i == N - 1)) then
                    coefficient = 3.D0 * 3 / 8
                else if(i == N) then
                    coefficient = 1.D0 * 3 / 8
                else
                    if(mod(i, 2) == 0) coefficient = 4.D0 / 3
                    if(mod(i, 2) == 1) coefficient = 2.D0 / 3
                end if

            end if

            coefficient = coefficient * (b - a) / (N - 1)

        end function coefficient_simpson

        function value_simpson(i, N, a, b) result(val)
            ! Provide the grid value (x_i) of simpson integral method
            ! \[\int_{a}^{b} f(x) dx = \sum_i c_i f(x_i)\]
            ! The index "i" start from 1 ~ N (Fortran convention)
            integer, intent(in) :: i, N
            real*8, intent(in) :: a, b
            real*8 :: val

            val = (i - 1) * (b - a) / (N - 1) + a

        end function value_simpson

        function fermi_func(energy, mu, temperature) result(f_E)
            implicit none
            real*8, intent(in) :: energy, mu, temperature
            real*8 :: f_E
            f_E = 1.D0 / (exp((energy - mu) / (k_B * temperature)) + 1)
        end function fermi_func

        
end module math_kernel