!**********************************************************************
! CALCULATES THE DENSITY AND TRANSMISSION OF A GIVEN POTENTIAL
!   Using NEGF to calculate density and transmission coefficient with 
!   real grid in z-direction, and plane wave basis in xy-direction.
! input files:
!   INPUT
!   POTENTIAL
! output files:
!   OUTPUT
!   DENSITY
!   TRANSMISSION
!   LDOS (optional)
!                                   By YI-CHENG LIN at 2023/6/27
!**********************************************************************
program main
    use negf
    use grid

    implicit none
    real*8 :: V_L, V_R, MU, ETA, TEMPERATURE, LX, LY, LZ, ENCUT
    integer :: NGX, NGY, N_circle, N_line_per_eV
    namelist /INPUT/ V_L, V_R, MU, ETA, TEMPERATURE, LX, LY, LZ, NGX, NGY, ENCUT, N_circle, N_line_per_eV

    integer :: N_x, N_y, N_z, N, i, j, k, ii, STATUS, N_line, i_job
    real :: start_time, end_time
    real*8 :: kx, ky, circle_L, circle_R, line_L, line_R, E_c, E_R, theta, sum_density, rescale_factor
    complex*16 :: energy
    real*8, allocatable :: V_real(:, :, :), Density(:, :, :)
    complex*16, allocatable :: V_reciprocal(:, :), Hamiltonian(:, :, :), E_minus_H(:, :, :), G_Function(:, :, :, :)&
    , GreenDiag(:, :, :)
    integer, allocatable :: nx_grid(:), ny_grid(:)

    include 'constant.f90'

    open(unit=16, file="OUTPUT")
    call cpu_time(start_time)
    write(16, *) "**********************************************************************"
    write(16, *) " CALCULATES THE DENSITY AND TRANSMISSION OF A GIVEN POTENTIAL"
    write(16, *) "   Using NEGF to calculate density and transmission coefficient with "
    write(16, *) "   real grid in z-direction, and plane wave basis in xy-direction."
    write(16, *) " input files:"
    write(16, *) "   INPUT"
    write(16, *) "   POTENTIAL"
    write(16, *) " output files:"
    write(16, *) "   OUTPUT"
    write(16, *) "   DENSITY"
    write(16, *) "   TRANSMISSION"
    write(16, *) "   LDOS (optional)"
    write(16, *) "                                   By YI-CHENG LIN at 2023/6/27"
    write(16, *) "**********************************************************************"
    write(16, *) ""

    !===Initialize Parameters===
    V_L = 0.D0 ! eV
    V_R = 0.D0 ! eV
    MU = 5.568D0 ! eV (correspond to RS=3.0, i.e. Au)
    ETA = 1.D-5 ! Hartree
    TEMPERATURE = 300.D0 ! K
    LX = 0.D0 ! Angstrom
    LY = 0.D0 ! Angstrom
    LZ = 0.D0 ! Angstrom
    ENCUT = 100.D0 ! eV
    
    NGX = 0
    NGY = 0
    N_circle = 40 ! (points/pi)
    N_line_per_eV = 500 ! (points/eV)

    !===Read Files===
    write(16, *) "Read in INPUT..."
    open(unit=15, file="INPUT")
    read(15, INPUT)
    close(15)
    if((LX == 0.D0) .or. (LY == 0.D0) .or. (LZ == 0.D0)) then
        write(16, *) "ERROR: LX, LY, LZ must be specified in INPUT"
        call exit(STATUS)
    end if
    write(16, *) "-success-"

    write(16, *) "----------"
    write(16, *) "Read in POTENTIAL..."
    open(unit=17, file="POTENTIAL")
    read(17, *) N_x, N_y, N_z
    allocate(V_real(N_x, N_y, N_z))
    read(17, *) (((V_real(i, j, k), i=1, N_x), j=1, N_y), k=1, N_z)
    close(17)
    if((mod(N_x, 2) == 1) .or. (mod(N_y, 2) == 1)) then
        write(16, *) "ERROR: The number of grid point in x & y direction must be even number"
        call exit(STATUS)
    end if
    write(16, *) "-success-"
    write(16, *) "The size of POTENTIAL is:"
    write(16, '(5X, 3I5)') N_x, N_y, N_z

    if(NGX == 0) NGX = (N_x / 3) * 2
    if(NGY == 0) NGY = (N_y / 3) * 2
    ! UNDONE: Check NGX, NGY

    write(16, *) "The parameters are:"
    write(16, INPUT)

    !===Convert unit===
    V_L = V_L / hartree ! eV -> hartree
    V_R = V_R / hartree ! eV -> hartree
    MU = MU / hartree ! eV -> hartree
    ETA = ETA ! Hartree
    TEMPERATURE = TEMPERATURE ! K
    LX = LX / a_0 ! Angstrom -> bohr
    LY = LY / a_0 ! Angstrom -> bohr
    LZ = LZ / a_0 ! Angstrom -> bohr
    ENCUT = ENCUT / hartree ! eV -> hartree
    do k=1, N_z
        do j=1, N_y
            do i=1, N_x
                V_real(i, j, k) = V_real(i, j, k) / hartree ! eV -> hartree
            end do
        end do
    end do



    !===Grid layout===
    write(16, *) "----------"
    write(16, *) "Layout the grid..."

    allocate(V_reciprocal(-NGX/2: NGX/2, -NGY/2: NGY/2))
    N = PlaneWaveBasis_construction_findsize(ENCUT, Lx, Ly)
    allocate(nx_grid(N), ny_grid(N), Hamiltonian(N, N, N_z), E_minus_H(N, N, N_z))
    allocate(G_function(N, N, N_z, 3), GreenDiag(N_x, N_y, N_z), Density(N_x, N_y, N_z))
    call PlaneWaveBasis_construction(ENCUT, Lx, Ly, nx_grid, ny_grid)

    write(16, *) "Grid size in x, y direction is: N =", N
    write(16, *) "Grid size in z direction is: N_z =", N_z
    write(16, *) "The calculation time is proportional to (N^3 N_z)"
    write(16, *) "----------"
    write(16, *) "NEGF calculation start..."

    !============================================================================================
    !==========================================START=============================================
    !===Construct the Hamiltonian===
    kx = 0.D0 ! Temporary
    ky = 0.D0 ! Temporary
    do k=1, N_z
        call LocalPotential_RealToPlanewave(V_real(:, :, k), V_reciprocal)
        call sub_Hamiltonian(nx_grid, ny_grid, Lx, Ly, kx, ky, V_reciprocal, Hamiltonian(:, :, k))
    end do
    
    !===NEGF integration===
    
    circle_L = minval(V_real) - 1.D0 / hartree
                               !  ^ This improve precision significantly 
    circle_R = MU - 6 * k_B * temperature
    line_L = circle_R
    line_R = MU + 6 * k_B * temperature

    E_c = (circle_R + circle_L) / 2
    E_R = (circle_R - circle_L) / 2
    N_line = IDNINT(N_line_per_eV * ((line_R - line_L) * hartree))
    write(16, *) "Total energy points for integration:", N_circle + N_line
    write(16, *) "circle_L : circle_R = ", circle_L, ":", circle_R!!
    write(16, *) "line_L : line_R = ", line_L, ":", line_R!!
    close(16)
    Density(:, :, :) = 0.D0
    sum_density = 0.D0
    do i_job=1, N_circle + N_line
        if(i_job <= N_circle) then
            ! CASE1
            theta = value_simpson(i_job, N_circle, 0.D0, pi)
            energy = E_c + E_R * exp(complex(0.D0, theta))
        else
            ! CASE2
            energy = value_simpson(i_job - N_circle, N_line, line_L, line_R)
        end if

        ! Main calculation
        call E_minus_H_construction(Hamiltonian, energy + complex(0.D0, ETA), nx_grid, ny_grid, Lx, Ly, Lz, &
         kx, ky, V_L, V_R, E_minus_H)
        call GreensFunction_tri_solver(E_minus_H, G_Function)
        do ii=1, 3
            do k=1, N_z
                do j=1, N
                    do i=1, N
                        ! Rescale to compensate the extra factor in E_minus_H
                        G_Function(i, j, k, ii) = G_Function(i, j, k, ii) * 2.D0 * (LZ / N_z) ** 2
                    end do
                end do
            end do
        end do
        do k=1, N_z
            call Greensfunction_PlanewaveToReal(G_Function(:, :, k, 1), nx_grid, ny_grid, GreenDiag(:, :, k))
        end do

        
        if(i_job <= N_circle) then
            ! CASE1
            do k=1, N_z
                do j=1, N_y
                    do i=1, N_x
                        Density(i, j, k) = Density(i, j, k) + &
                         2.D0 / pi * coefficient_simpson(i_job, N_circle, 0.D0, pi) * &
                         aimag(complex(0.D0, 1.D0) * E_R * exp(complex(0.D0, theta)) * GreenDiag(i, j, k))
                    end do
                end do
            end do
        else
            ! CASE2
            do k=1, N_z
                do j=1, N_y
                    do i=1, N_x
                        Density(i, j, k) = Density(i, j, k) - &
                         2.D0 / pi * coefficient_simpson(i_job - N_circle, N_line, line_L, line_R) * &
                         aimag(fermi_func(dble(energy), MU, TEMPERATURE) * GreenDiag(i, j, k))
                    end do
                end do
            end do

        end if
        
        open(unit=16, file="OUTPUT", status="old", position="append")
        write(16, '(A2, I4, A1, I4)') "=>", i_job, "/", N_circle + N_line
        close(16)
    end do

    !===========================================END==============================================
    !============================================================================================

    ! Rescale　density
    rescale_factor = dble(N_z) / (LX * LY * LZ * a_0 ** 3)
    do k=1, N_z
        do j=1, N_y
            do i=1, N_x
                Density(i, j, k) = Density(i, j, k) * rescale_factor ! unit: 1/Angstrom^3
            end do
        end do
    end do

    !　Write out data
    open(unit=16, file="OUTPUT", status="old", position="append")
    write(16, *) "NEGF calculation DONE"
    write(16, *) "Total charge is:", sum(Density) * (LX * LY * LZ * a_0 ** 3) / (N_x * N_y * N_z)
    write(16, *) "Writing DENSITY..."

    open(unit=18, file="DENSITY")
    write(18, '(3I5)') N_x, N_y, N_z
    write(18, '(5G17.8)') (((Density(i, j, k), i=1, N_x), j=1, N_y), k=1, N_z)
    close(18)

    call cpu_time(end_time)
    write(16,*) "CPU time: ", end_time - start_time, "s"
    write(16,*) "********************"
    write(16,*) "*                  *"
    write(16,*) "*      DONE        *"
    write(16,*) "*                  *"
    write(16,*) "********************"
    close(16)
end program main