program potential_QPC
    ! Can be compiled by
    ! gfortran potential_QPC.f90 -o potential_QPC.x
    implicit none
    integer, parameter :: N_x = 50, N_y = 90, N_z = 200

    ! Unit: Angstrom, eV
    real*8, parameter :: Lx = 8.D0, Ly = 15.D0, Lz = 30.D0, Ax = 2.D0, Ay = 8.D0, Az = 10.D0, Height = 27.D0
    real*8, parameter :: a_0 = 0.529177210903D0 ! Bohr radius (in Angstrom)
    real*8, parameter :: pi = 3.141592653589793D0
    real*8, parameter :: hartree = 27.211386245988D0 ! (in eV)
    real*8, allocatable :: V_real(:, :, :)
    integer :: i, j, k

    allocate(V_real(N_x, N_y, N_z))
    V_real(:, :, :) = 0.D0
    do k=1, N_z
        do j=1, N_y
            do i=1, N_x
                if((i * Lx / N_x > Ax) .or. ((j * Ly / N_y > Ay) .and. &
                 (k * Lz / N_z > Lz / 2 - Az / 2) .and. (k * Lz / N_z <= Lz / 2 + Az / 2))) V_real(i, j, k) = Height
            end do
        end do
    end do

    do i=1, 6
        write(*, *) "n =", i, "           E =", &
        (pi ** 2 / (2 * (Ax / a_0) ** 2) + (i * pi) ** 2 / (2 * (Ay / a_0) ** 2)) * hartree, "eV"
    end do

    open(unit=17, file="POTENTIAL")
    write(17, '(3I5)') N_x, N_y, N_z
    write(17, '(5G17.8)') (((V_real(i, j, k), i=1, N_x), j=1, N_y), k=1, N_z)

end program potential_QPC