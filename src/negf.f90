module negf
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

    subroutine GreensFunction_tri_solver(G_inv_D, G_Function)
        !********************************************************
        ! This is a triangular matrix inverse solver, see
        !       "Journal of Applied Physics 81, 7845 (1997)"
        ! for details (Recursive Green function algorithm)
        !
        ! Assume G_inv_DU, G_inv_DL are identities
        !********************************************************
    
            implicit none
            complex*16, intent(in) :: G_inv_D(:, :, :)
            complex*16, intent(out) :: G_Function(:, :, :, :)
            ! Structure of G_Function is (N, N, N_z, 3), with fourth dimension store:
            !   1: diagonal element
            !   2: first column element
            !   3: last column element
            complex*16, allocatable :: generator_l(:, :, :), generator_r(:, :, :), Matrix(:, :)
            
            integer :: N_z, N, i, j, k
    
            N_z = size(G_inv_D, 3)
            N = size(G_inv_D, 1)
            allocate(generator_l(N, N, N_z - 1))
            allocate(generator_r(N, N, N_z - 1))
            allocate(Matrix(N, N))
    
            ! Calculate generator_l()
            ! Check the order if you want to generalize the code into block matrix
            call inverse(G_inv_D(:, :, 1), generator_l(:, :, 1))
            do k=2, N_z - 1
                do j=1, N
                    do i=1, N
                        Matrix(i, j) = G_inv_D(i, j, k) - generator_l(i, j, k - 1)
                    end do
                end do
                call inverse(Matrix, generator_l(:, :, k))
            end do
    
            ! Calculate generator_r()
            call inverse(G_inv_D(:, :, N_z), generator_r(:, :, N_z - 1))
            do k=2, N_z - 1
                do j=1, N
                    do i=1, N
                        Matrix(i, j) = G_inv_D(i, j, N_z - k + 1) - generator_r(i, j, N_z - k + 1)
                    end do
                end do
                call inverse(Matrix, generator_r(:, :, N_z - k))
            end do
    
            ! Get G_Function(:, 1) (i.e. diagonal element)
            do j=1, N
                do i=1, N
                    Matrix(i, j) = G_inv_D(i, j, 1) - generator_r(i, j, 1)
                end do
            end do
            call inverse(Matrix, G_Function(:, :, 1, 1))

            do k=2, N_z
                call matrix_product(G_Function(:, :, k - 1, 1), generator_r(:, :, k - 1), Matrix)
                do i=1, N
                    Matrix(i, i) = Matrix(i, i) + 1.D0
                end do
                call matrix_product(generator_r(:, :, k - 1), Matrix, G_Function(:, :, k, 1))
            end do
    
            ! Get G_Function(:, 2) (i.e. first column element)
            do j=1, N
                do i=1, N
                    G_Function(i, j, 1, 2) = G_Function(i, j, 1, 1)
                end do
            end do
            
            do k=2, N_z
                call matrix_product(generator_r(:, :, k - 1), G_Function(:, :, k - 1, 2), G_Function(:, :, k, 2))
                do j=1, N
                    do i=1, N
                        G_Function(i, j, k, 2) = - G_Function(i, j, k, 2)
                    end do
                end do
            end do
    
            ! Get G_Function(:, 3) (i.e. last column element) 
            ! Check the order if you want to generalize the code into block matrix
            do j=1, N
                do i=1, N
                    G_Function(i, j, N_z, 3) = G_Function(i, j, N_z, 1)
                end do
            end do
            do k=1, N_z - 1
                call matrix_product(generator_l(:, :, N_z - k), G_Function(:, :, N_z - k + 1, 3), G_Function(:, :, N_z - k, 3))
                do j=1, N
                    do i=1, N
                        G_Function(i, j, N_z - k, 3) = - G_Function(i, j, N_z - k, 3)
                    end do
                end do
            end do
            
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
            include 'constant.f90'
            f_E = 1D0 / (exp((energy - mu) / (k_B * temperature)) + 1)
        end function fermi_func

        
end module