program test_simpson
    ! To compile, enter the build/ directory, and type:
    ! gfortran ../tools/test_simpson.f90 -o test_simpson.x math_kernel.o -llapack -lblas
    use math_kernel
    implicit none
    integer, parameter :: N_integral = 20
    real*8, parameter :: l_bound = 0.D0, u_bound = 1.D0

    integer :: i
    real*8 :: sum

    sum = 0.D0
    do i=1, N_integral
        sum = sum + coefficient_simpson(i, N_integral, l_bound, u_bound) * &
            exp(- value_simpson(i, N_integral, l_bound, u_bound))
    end do

    print '(f35.30)', 1.D0 - exp(-1.D0)
    print '(f35.30)', sum
    print *, "ERROR =", 1.D0 - exp(-1.D0) - sum

    
end program test_simpson