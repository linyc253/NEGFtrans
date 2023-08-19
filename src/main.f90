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
    use mpi !%%
    use global

    implicit none
    real*8 :: V_L, V_R, MU, ETA, TEMPERATURE, LX, LY, LZ, ENCUT, GAP
    integer :: NGX, NGY, N_circle, N_line_per_eV, NKX, NKY
    namelist /INPUT/ V_L, V_R, MU, ETA, TEMPERATURE, LX, LY, LZ, NGX, NGY, ENCUT, &
     N_circle, N_line_per_eV, GAP, NKX, NKY
    
    ! All the input parameters will be converted to atomic unit and stored in "atomic"
    type(t_parameters) :: atomic
    integer :: N_x, N_y, N_z, N, i, j, k, STATUS, N_line, i_job, N_job, rank = 0, mpi_size = 1
    type(t_timer) :: total_time, inverse_time, RtoP_time, PtoR_time
    type(t_kpointmesh), allocatable :: kpoint(:)
    real*8 :: rescale_factor
    real*8, allocatable :: V_real(:, :, :), Density(:, :, :)
    complex*16, allocatable :: V_reciprocal(:, :), V_reciprocal_all(:, :, :)
    integer, allocatable :: nx_grid(:), ny_grid(:)

    call MPI_INIT(STATUS) !%%
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, STATUS) !%%
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, STATUS) !%%

    if(rank == 0) then
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
    end if

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
    NKX = 0
    NKY = 0
    N_circle = 40 ! (points/pi)
    N_line_per_eV = 500 ! (points/eV)

    !===Read Files===
    ! Read INPUT
    if(rank == 0) write(16, *) "Read in INPUT..."
    open(unit=15, file="INPUT")
    read(15, INPUT)
    close(15)

    ! Modify NGX, NGY
    if(NGX <= ceiling((LY / a_0) / pi * sqrt(ENCUT / hartree / 2) * 2)) then
        NGX = ceiling((LX / a_0) / pi * sqrt(ENCUT / hartree / 2) * 2) + 1
    end if
    if(NGY <= ceiling((LY / a_0) / pi * sqrt(ENCUT / hartree / 2) * 2)) then
        NGY = ceiling((LY / a_0) / pi * sqrt(ENCUT / hartree / 2) * 2) + 1
    end if

    if(rank == 0) then
        ! Check LX, LY, LZ
        if((LX == 0.D0) .or. (LY == 0.D0) .or. (LZ == 0.D0)) then
            write(16, *) "ERROR: LX, LY, LZ must be specified in INPUT"
            call exit(STATUS)
        end if
        if((NKX == 0) .or. (NKY == 0)) then
            write(16, *) "ERROR: NKX, NKY must be specified in INPUT"
            call exit(STATUS)
        end if
        write(16, *) "-success-"
        write(16, *) "The parameters are:"
        write(16, INPUT)
    end if

    ! Read POTENTIAL
    if(rank == 0) write(16, *) "----------"
    if(rank == 0) write(16, *) "Read in POTENTIAL..."
    open(unit=17, file="POTENTIAL")
    read(17, *) N_x, N_y, N_z
    allocate(V_real(N_x, N_y, N_z))
    read(17, *) (((V_real(i, j, k), i=1, N_x), j=1, N_y), k=1, N_z)
    close(17)

    if(rank == 0) then
        ! Check N_x, N_y
        if((mod(N_x, 2) == 1) .or. (mod(N_y, 2) == 1)) then
            write(16, *) "ERROR: The number of grid point in x & y direction must be even number"
            call exit(STATUS)
        end if
        write(16, *) "-success-"
        write(16, *) "The size of POTENTIAL is:", N_x, N_y, N_z
        write(16, *) "----------"
    end if

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

    !===K-point layout===
    allocate(kpoint(NKX * NKY))
    call Kpoint_mesh_construction(NKX, NKY, atomic%LX, atomic%LY, kpoint)
    if(rank == 0) then
        write(16, *) "Number of K-point grid is:", size(kpoint)
        write(16, *) "   kx          ky           weight"
        write(16, '(3G12.3)') (kpoint(i)%kx, kpoint(i)%ky, kpoint(i)%weight, i=1, size(kpoint))
        write(16, *) "----------"
    end if

    !===Grid layout===
    allocate(Density(N_x, N_y, N_z))
    allocate(V_reciprocal(-NGX: NGX, -NGY: NGY), V_reciprocal_all(-NGX: NGX, -NGY: NGY, 1: N_z))

    N = PlaneWaveBasis_construction_findsize(atomic%ENCUT, atomic%LX, atomic%LY)
    allocate(nx_grid(N), ny_grid(N))
    call PlaneWaveBasis_construction(atomic%ENCUT, atomic%LX, atomic%LY, nx_grid, ny_grid)

    if(rank == 0) then
        write(16, *) "Layout the Grid..."
        write(16, *) "Grid size in x, y direction is: N =", N
        write(16, *) "Grid size in z direction is: N_z =", N_z
        write(16, *) "The 'Inverse Matrix' time is proportional to (N^3 N_z)"
        write(16, *) "----------"
        write(16, *) "Transform Potential to plane wave basis..."
        write(16, *) "----------"
        close(16)
    end if

    !============================================================================================
    !==========================================START=============================================

    ! STEP 1. Construct 'V_reciprocal_all'
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

    ! STEP 2. Use NEGF to find 'Density', parallelized
    if(rank == 0) then
        open(unit=16, file="OUTPUT", status="old", position="append")
        write(16, *) "Density calculation start..."
        close(16)
    end if
    Density(:, :, :) = 0.D0
    N_line = IDNINT(N_line_per_eV * ((2 * atomic%GAP) * hartree))

    N_job = (N_circle + N_line) * size(kpoint)
    i_job = 1 + rank
    do while(i_job <= N_job)
        ! Density = Density + Density_contribution_of_that_energy_point
        call Equilibrium_Density(i_job, V_reciprocal_all, nx_grid, ny_grid, N_circle, N_line, minval(V_real), atomic, &
        kpoint, Density, inverse_time, PtoR_time)

        i_job = i_job + mpi_size
        
        if(rank == 0) then
            open(unit=16, file="OUTPUT", status="old", position="append")
            write(16, '(A2, I4, A2, I4)') "=>", min(i_job - 1, N_job), " /", N_job
            close(16)
        end if
    end do

    ! Collect data from different rank
    if(rank == 0) then !%%
        call MPI_Reduce(MPI_IN_PLACE, Density, size(Density), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, STATUS) !%%
    else !%%
        call MPI_Reduce(Density, Density, size(Density), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, STATUS) !%%
    end if !%%

    !===========================================END==============================================
    !============================================================================================
    if(rank == 0) then
        ! Rescale density. (Unit: 1/Angstrom^3)
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
        write(16, *) "----------"
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
    end if

    call MPI_FINALIZE(STATUS) !%%
end program main