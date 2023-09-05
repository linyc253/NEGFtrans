program test_inverse
    ! To compile, enter the build/ directory, and type:
    ! gfortran ../tools/test_inverse.f90 -o test_inverse.x math_kernel.o global.o -llapack -lblas
    use math_kernel
    use global
    implicit none
    integer, parameter :: N = 30, N_z = 50
    real*8 :: rand_real, rand_imag, start_time, end_time, sum
    complex*16, allocatable :: E_minus_H(:, :), Blocks(:, :, :), Greensfunc(:, :), work(:)
    type(t_gfunc), allocatable :: G_function(:, :, :)
    integer, allocatable :: IPIV(:)
    integer :: i, j, z, N_tot, STATUS, NB, ilaenv

    N_tot = N * N_z
    allocate(E_minus_H(N_tot, N_tot), Blocks(N, N , N_z), Greensfunc(N_tot, N_tot))
    allocate(G_function(N, N, N_z))

    do z=1, N_z
        do j=1, N
            do i=1, N
                call random_number(rand_real)
                call random_number(rand_imag)
                Blocks(i, j, z) = complex(rand_real, rand_imag)
                E_minus_H(i + N * (z - 1), j + N * (z - 1)) = Blocks(i, j, z)
            end do
        end do
    end do

    do z=1, N_z - 1
        do i=1, N
                E_minus_H(i + N * (z - 1), i + N * (z - 1 + 1)) = 1.D0
                E_minus_H(i + N * (z - 1 + 1), i + N * (z - 1)) = 1.D0
        end do
    end do
    ! do i=1, N_tot
    !     print '(12f8.2)', (real(E_minus_H(i, j)), j=1, N_tot)
    ! end do

    do j=1, N_tot
        do i=1, N_tot
            Greensfunc(i, j) = E_minus_H(i, j)
        end do
    end do

    ! Inverse the matrix, Method 1
    call cpu_time(start_time)
    NB = ilaenv(1, "zgetri", "", N_tot, N_tot, -1, -1)
    if(NB < 1) NB = N_tot
    allocate(IPIV(N_tot), work(NB * N_tot))
    call zgetrf(N_tot, N_tot, Greensfunc, N_tot, IPIV, STATUS)
    call zgetri(N_tot, Greensfunc, N_tot, IPIV, work, NB * N_tot, STATUS)
    call cpu_time(end_time)
    print *, "Method 1", end_time - start_time, "s"

    ! do j=1, N
    !     print '(10f8.2)', (real(Greensfunc(i, j)), i=1, N)
    ! end do

    ! Inverse the matrix, Method 2
    call cpu_time(start_time)
    call GreensFunction_tri_solver(Blocks, G_function, "FLD")
    call cpu_time(end_time)
    print *, "Method 2", end_time - start_time, "s"

    ! do j=1, N
    !     print '(10f8.2)', (real(G_function(i, j, 1, 1)), i=1, N)
    ! end do

    ! Check
    sum = 0.D0
    do z=1, N_z
        do j=1, N
            do i=1, N
                sum = sum + &
                abs(Greensfunc(i + N * (z - 1), j + N * (z - 1)) - G_function(i, j, z)%diagonal) ** 2
            end do
        end do
    end do
    print *, "Error1 =", sqrt(sum)

    sum = 0.D0
    do z=1, N_z
        do j=1, N
            do i=1, N
                sum = sum + &
                abs(Greensfunc(i + N * (z - 1), j) - G_function(i, j, z)%first_column) ** 2
            end do
        end do
    end do
    print *, "Error2 =", sqrt(sum)

    sum = 0.D0
    do z=1, N_z
        do j=1, N
            do i=1, N
                sum = sum + &
                abs(Greensfunc(i + N * (z - 1), j + N * (N_z - 1)) - G_function(i, j, z)%last_column) ** 2
            end do
        end do
    end do
    print *, "Error3 =", sqrt(sum)

    
end program test_inverse