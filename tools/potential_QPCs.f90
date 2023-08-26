program potential_QPCs
    ! To compile, enter the build/ directory, and type:
    ! gfortran ../tools/potential_QPCs.f90 -o potential_QPCs.x global.o
    use global
    implicit none
    integer, parameter :: N_x = 50, N_y = 90, N_z = 300

    ! Unit: Angstrom, eV
    real*8, parameter :: Lx = 8.D0, Ly = 18.D0, Lz = 40.D0, Az = 10.D0, Height = 21.D0, sigma = 1.5D0
    real*8, parameter :: gap_x = 12.D0, gap_y = 2.D0 ! eV
    real*8, allocatable :: V_real(:, :, :)
    integer :: i, j, k
    real*8 :: x, y, z, f, omega_x, omega_y
    omega_x = gap_x / hartree ! by \Delta E = \hbar\omega  => \Delta E = \omega in atomic unit
    omega_y = gap_y / hartree

    allocate(V_real(N_x, N_y, N_z))
    do k=1, N_z
        do j=1, N_y
            do i=1, N_x
                x = i * Lx / N_x
                y = j * Ly / N_y
                z = k * Lz / N_z
                f = (tanh((z - (Lz / 2 - Az / 2)) / sigma) - tanh((z - (Lz / 2 + Az / 2)) / sigma)) / 2 ! no unit
                V_real(i, j, k) = omega_x ** 2 / 2 * ((x - Lx / 2) / a_0) ** 2 + f * omega_y ** 2 / 2 * ((y - Ly / 2) / a_0) ** 2 
                V_real(i, j, k) = V_real(i, j, k) * hartree ! atomic unit (hartree) -> eV
                if(V_real(i, j, k) > Height) V_real(i, j, k) = Height
            end do
        end do
    end do

    write(*, *) "Max(POTENTIAL) =", maxval(V_real)
    write(*, *) "n_x = 1"
    do i=0, 6
        write(*, *) "  n_y =", i, "           E =", (gap_x * 0.5D0) + (gap_y * (i + 0.5D0)), "eV"
    end do
    write(*, *) "n_x = 2 comes in at", (gap_x * 1.5D0) + (gap_y * 0.5D0), "eV"

    open(unit=17, file="POTENTIAL")
    write(17, '(3I5)') N_x, N_y, N_z
    write(17, '(5G17.8)') (((V_real(i, j, k), i=1, N_x), j=1, N_y), k=1, N_z)

end program potential_QPCs