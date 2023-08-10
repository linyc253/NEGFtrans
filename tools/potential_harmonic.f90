program potential_harmonic
    ! To compile, enter the build/ directory, and type:
    ! gfortran ../tools/potential_harmonic.f90 -o potential_harmonic.x
    implicit none
    integer, parameter :: N_x = 50, N_y = 56, N_z = 100

    ! Unit: Angstrom, eV
    real*8, parameter :: Lx = 10.D0, Ly = 8.D0, Lz = 15.D0, A = 97.9D0, R = 3.D0
    real*8 :: omega
    real*8, allocatable :: V_real(:, :, :)
    integer :: i, j, k
    include 'constant.f90'

    allocate(V_real(N_x, N_y, N_z))
    V_real(:, :, :) = 0.D0
    do k=1, N_z
        do j=1, N_y
            do i=1, N_x
                if ((Lx * i / N_x - Lx / 2) ** 2 + (Ly * j / N_y - Ly / 2) ** 2 +&
                (Lz * k / N_z - Lz / 2) ** 2 <= R ** 2) then
                    V_real(i, j, k) = A / R ** 2 * (&
                    (Lx * i / N_x - Lx / 2) ** 2 + (Ly * j / N_y - Ly / 2) ** 2 +&
                    (Lz * k / N_z - Lz / 2) ** 2) - A
                end if
            end do
        end do
    end do

    omega = sqrt(2 * (A / hartree) / (R / a_0) ** 2) ! atomic unit
    write(*, *) "omega =", omega, "(atomic unit)"
    write(*, *) "(E_0 + E_1) / 2 =", 2 * omega * hartree - A, "eV"
    write(*, *) "(E_1 + E_2) / 2 =", 3 * omega * hartree - A, "eV"
    write(*, *) "(E_2 + E_3) / 2 =", 4 * omega * hartree - A, "eV"

    open(unit=17, file="POTENTIAL")
    write(17, '(3I5)') N_x, N_y, N_z
    write(17, '(5G17.8)') (((V_real(i, j, k), i=1, N_x), j=1, N_y), k=1, N_z)

end program potential_harmonic