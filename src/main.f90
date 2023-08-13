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
    use math_kernel
    use grid
    use negf

    implicit none
    real*8 :: V_L, V_R, MU, ETA, TEMPERATURE, LX, LY, LZ, ENCUT, GAP
    integer :: NGX, NGY, N_circle, N_line_per_eV
    namelist /INPUT/ V_L, V_R, MU, ETA, TEMPERATURE, LX, LY, LZ, NGX, NGY, ENCUT, N_circle, N_line_per_eV, GAP
    
    ! All the input parameters will be converted to atomic unit and stored in "atomic"
    type(t_parameters) :: atomic
    integer :: N_x, N_y, N_z, N, i, j, k, STATUS, N_line, i_job
    type(t_timer) :: total_time, inverse_time, RtoP_time, PtoR_time
    real*8 :: rescale_factor
    real*8, allocatable :: V_real(:, :, :), Density(:, :, :)
    complex*16, allocatable :: V_reciprocal(:, :), V_reciprocal_all(:, :, :)
    integer, allocatable :: nx_grid(:), ny_grid(:)

    include 'constant.f90'

    open(unit=16, file="OUTPUT")
    call cpu_time(total_time%start)
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
    GAP = 6.D0 ! k_B * TEMPERATURE
    
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
    !XXXXXXXXXXXXXXXXXXXXXX UNDONE: Check NGX, NGY

    write(16, *) "The parameters are:"
    write(16, INPUT)

    !===Convert unit===
    ! Note that the name is changed from "parameter" to "atomic%parameter"
    atomic%V_L = V_L / hartree ! eV -> hartree
    atomic%V_R = V_R / hartree ! eV -> hartree
    atomic%MU = MU / hartree ! eV -> hartree
    atomic%ETA = ETA ! Hartree
    atomic%TEMPERATURE = TEMPERATURE ! K
    atomic%LX = LX / a_0 ! Angstrom -> bohr
    atomic%LY = LY / a_0 ! Angstrom -> bohr
    atomic%LZ = LZ / a_0 ! Angstrom -> bohr
    atomic%ENCUT = ENCUT / hartree ! eV -> hartree
    atomic%GAP = GAP * k_B * atomic%TEMPERATURE ! (k_B T) -> hartree
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

    allocate(Density(N_x, N_y, N_z))
    allocate(V_reciprocal(-NGX/2: NGX/2, -NGY/2: NGY/2), V_reciprocal_all(-NGX/2: NGX/2, -NGY/2: NGY/2, 1: N_z))

    N = PlaneWaveBasis_construction_findsize(atomic%ENCUT, atomic%LX, atomic%LY)
    allocate(nx_grid(N), ny_grid(N))
    call PlaneWaveBasis_construction(atomic%ENCUT, atomic%LX, atomic%LY, nx_grid, ny_grid)

    write(16, *) "Grid size in x, y direction is: N =", N
    write(16, *) "Grid size in z direction is: N_z =", N_z
    write(16, *) "The 'Inverse Matrix' time is proportional to (N^3 N_z)"
    write(16, *) "----------"
    write(16, *) "NEGF calculation start..."
    close(16)

    !============================================================================================
    !==========================================START=============================================

    ! STEP 1. Construct 'V_reciprocal_all', single-thread
    call cpu_time(RtoP_time%start)
    do k=1, N_z
        call LocalPotential_RealToPlanewave(V_real(:, :, k), V_reciprocal)
        do j= -NGY/2, NGY/2
            do i= -NGX/2, NGX/2
                V_reciprocal_all(i, j, k) = V_reciprocal(i, j)
            end do
        end do
    end do
    call cpu_time(RtoP_time%end)
    RtoP_time%sum = RtoP_time%sum + RtoP_time%end - RtoP_time%start

    ! STEP 2. Use NEGF to find 'Density', multi-thread
    Density(:, :, :) = 0.D0
    N_line = IDNINT(N_line_per_eV * ((2 * atomic%GAP) * hartree))
    do i_job=1, N_circle + N_line
        ! Density = Density + Density_contribution_of_that_energy_point
        call Equilibrium_Density(i_job, V_reciprocal_all, nx_grid, ny_grid, N_circle, N_line, minval(V_real), atomic, Density,&
        inverse_time, PtoR_time)

        open(unit=16, file="OUTPUT", status="old", position="append")
        write(16, '(A2, I4, A1, I4)') "=>", i_job, " /", N_circle + N_line
        close(16)
    end do
    !===========================================END==============================================
    !============================================================================================
    

    ! Rescale　density. (Unit: 1/Angstrom^3)
    rescale_factor = dble(N_z) / (LX * LY * LZ)
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
    write(16, *) "Total charge is:", sum(Density) * (LX * LY * LZ) / (N_x * N_y * N_z)
    write(16, *) "Writing DENSITY..."

    open(unit=18, file="DENSITY")
    write(18, '(3I5)') N_x, N_y, N_z
    write(18, '(5G17.8)') (((Density(i, j, k), i=1, N_x), j=1, N_y), k=1, N_z)
    close(18)

    call cpu_time(total_time%end)
    total_time%sum = total_time%sum + total_time%end - total_time%start
    write(16,*) "CPU time (Total): ", total_time%sum, "s"
    write(16,*) "  CPU time (Real_space to Plane_wave): ", RtoP_time%sum, "s"
    write(16,*) "  CPU time (Inverse Matrix): ", inverse_time%sum, "s"
    write(16,*) "  CPU time (Plane_wave to Real_space): ", PtoR_time%sum, "s"
    write(16,*) "  CPU time (Others): ", total_time%sum - RtoP_time%sum - inverse_time%sum - PtoR_time%sum, "s"
    write(16,*) "********************"
    write(16,*) "*                  *"
    write(16,*) "*      DONE        *"
    write(16,*) "*                  *"
    write(16,*) "********************"
    close(16)
end program main