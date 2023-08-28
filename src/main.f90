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
!   LDOS
!                                   By YI-CHENG LIN at 2023/6/27
!**********************************************************************
program main
    use math_kernel
    use grid
    use negf
    use mpi !%%
    use global
    use tools

    implicit none
    real*8 :: V_L, V_R, MU, ETA, TEMPERATURE, LX, LY, LZ, ENCUT, GAP, VDS, TRANSMISSION_GRID(3), &
     LDOS_ENERGY_GRID(3)
    integer :: NGX, NGY, N_circle, N_line_per_eV, NKX, NKY, LDOS_GRID(3)
    logical :: if_density, if_transmission, if_LDOS
    namelist /INPUT/ V_L, V_R, MU, ETA, TEMPERATURE, LX, LY, LZ, NGX, NGY, ENCUT, &
     N_circle, N_line_per_eV, GAP, NKX, NKY, VDS, TRANSMISSION_GRID, if_density, if_transmission, if_LDOS, &
     LDOS_ENERGY_GRID, LDOS_GRID
    
    ! All the input parameters will be converted to atomic unit and stored in "atomic"
    type(t_parameters) :: atomic
    integer :: N_x, N_y, N_z, N, i, j, k, i_E, subsize, STATUS, N_line, i_job, N_job, rank = 0, mpi_size = 1, N_E
    integer :: ii, jj, kk
    type(t_timer) :: total_time, inverse_time, RtoP_time, PtoR_time, trans_time, LDOS_inverse_time, LDOS_PtoR_time
    type(t_kpointmesh), allocatable :: kpoint(:)
    type(t_transmission), allocatable :: Transmission(:)
    real*8 :: rescale_factor
    real*8, allocatable :: V_real(:, :, :), Density(:, :, :), LDOS(:, :, :, :)
    complex*16, allocatable :: V_reciprocal(:, :), V_reciprocal_all(:, :, :), Density_Matrix(:, :, :), &
     sendbuf(:, :, :), recvbuf(:, :, :), LDOS_Matrix(:, :, :, :)
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
        write(16, *) "   LDOS"
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
    VDS = 0.D0 ! eV
    TRANSMISSION_GRID(:) = 0.D0
    LDOS_ENERGY_GRID(:) = 0.D0
    LDOS_GRID(:) = -10
    
    NGX = 0
    NGY = 0
    NKX = 0
    NKY = 0
    N_circle = 40 ! (points/pi)
    N_line_per_eV = 500 ! (points/eV)

    if_density = .false.
    if_transmission = .false.
    if_LDOS = .false.

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
    
    ! Check INPUT parameters
    if(rank == 0) then
        write(16, *) "Checking Parameters..."
        if((LX == 0.D0) .or. (LY == 0.D0) .or. (LZ == 0.D0)) then
            write(16, *) "ERROR: LX, LY, LZ must be specified in INPUT"
            call exit(STATUS)
        end if
        if((NKX == 0) .or. (NKY == 0)) then
            write(16, *) "ERROR: NKX, NKY must be specified in INPUT"
            call exit(STATUS)
        end if
        if(.not. (if_density .or. if_transmission .or. if_LDOS)) then
            write(16, *) "ERROR: At least one of (if_xxxxx) be .true."
            write(16, *) "Nothing to do..."
            call exit(STATUS)
        end if
        if(if_transmission .and. (TRANSMISSION_GRID(3) == 0.D0)) then
            write(16, *) "ERROR: TRANSMISSION_GRID must be specified for (if_transmission = .true.)"
            call exit(STATUS)
        end if
        if(if_LDOS .and. (LDOS_ENERGY_GRID(3) == 0.D0)) then
            write(16, *) "ERROR: LDOS_ENERGY_GRID must be specified for (if_LDOS = .true.)"
            call exit(STATUS)
        end if
        if(if_LDOS .and. ((LDOS_GRID(1) < -1) .or. (LDOS_GRID(2) < -1) .or. (LDOS_GRID(3) < -1))) then
            write(16, *) "ERROR: LDOS_GRID must be specified for (if_LDOS = .true.)"
            call exit(STATUS)
        end if
        write(16, *) "-success-"
        write(16, *) "The Parameters are:"
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
    atomic%VDS = VDS / hartree ! eV -> hartree
    atomic%TRANSMISSION_GRID(1) = TRANSMISSION_GRID(1) / hartree ! eV -> hartree
    atomic%TRANSMISSION_GRID(2) = TRANSMISSION_GRID(2) / hartree ! eV -> hartree
    atomic%LDOS_ENERGY_GRID(1) = LDOS_ENERGY_GRID(1) / hartree ! eV -> hartree
    atomic%LDOS_ENERGY_GRID(2) = LDOS_ENERGY_GRID(2) / hartree ! eV -> hartree
    
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
    allocate(nx_grid(N), ny_grid(N), Density_Matrix(N, N, N_z))
    call PlaneWaveBasis_construction(atomic%ENCUT, atomic%LX, atomic%LY, nx_grid, ny_grid)

    if(rank == 0) then
        write(16, *) "Layout the Grid..."
        write(16, *) "Grid size in x, y direction is: N =", N
        write(16, *) "Grid size in z direction is: N_z =", N_z
        write(16, *) "The 'Inverse Matrix' time is proportional to (N^3 N_z)"
        write(16, *) "----------"
        write(16, *) "Transform Potential to plane wave basis..."
        close(16)
    end if

    !============================================================================================
    !==========================================START=============================================

    ! ===STEP 1. Construct 'V_reciprocal_all'===
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

    if(rank == 0) then
        open(unit=16, file="OUTPUT", status="old", position="append")
        write(16, *) "DONE"
        write(16, *) "----------"
        close(16)
    end if

    if(if_density) then
        ! ===STEP 2. Use NEGF to find 'Density_Matrix', parallelized===
        N_line = IDNINT(N_line_per_eV * ((2 * atomic%GAP) * hartree))
        if(rank == 0) then
            open(unit=16, file="OUTPUT", status="old", position="append")
            write(16, *) "Density calculation start..."
            write(16, *) "Number of energy points for integration:", N_circle + N_line
            write(16, *) "  ├── hemi-circle integral:", N_circle
            write(16, *) "  └── line integral:", N_line
            write(16, *) "Get Green's Function (Inverse Matrix):"
            close(16)
        end if

        Density_Matrix(:, :, :) = 0.D0
        N_job = (N_circle + N_line) * size(kpoint)
        i_job = 1 + rank
        do while(i_job <= N_job)
            ! Density_Matrix = Density_Matrix + [Density_Matrix contribution of that energy point]
            call Equilibrium_Density(i_job, V_reciprocal_all, nx_grid, ny_grid, N_circle, N_line, minval(V_real), atomic, &
            kpoint, Density_Matrix, inverse_time)

            i_job = i_job + mpi_size
            
            if(rank == 0) then
                open(unit=16, file="OUTPUT", status="old", position="append")
                write(16, '(A2, I6, A2, I6)') "=>", min(i_job - 1, N_job), " /", N_job
                close(16)
            end if
        end do

        ! Reduce from different rank
        if(rank == 0) then !%%
            call MPI_Reduce(MPI_IN_PLACE, Density_Matrix, size(Density_Matrix), &
            MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, STATUS) !%%
        else !%%
            call MPI_Reduce(Density_Matrix, Density_Matrix, size(Density_Matrix), &
            MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, STATUS) !%%
        end if !%%

        ! ===STEP 3. Transform Density_Matrix from plane wave basis to real space, parallelized===
        if(rank == 0) then
            open(unit=16, file="OUTPUT", status="old", position="append")
            write(16, *) "Transform Green's Function from plane wave basis to real space:"
            close(16)
        end if

        ! Scatter to different rank
        call cpu_time(PtoR_time%start)
        subsize = ((size(Density_Matrix, 3) - 1) / mpi_size + 1)
        allocate(sendbuf(N, N, subsize * mpi_size), recvbuf(N, N, subsize))
        if(rank == 0) then
            sendbuf(:, :, :) = 0.D0
            do k=1, N_z
                do j=1, N
                    do i=1, N
                        sendbuf(i, j, k) = Density_Matrix(i, j, k)
                    end do
                end do
            end do
        end if
        call MPI_Scatter(sendbuf, size(recvbuf), MPI_DOUBLE_COMPLEX, &
                        recvbuf, size(recvbuf), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, STATUS) !%%
        deallocate(sendbuf)
        allocate(sendbuf(N_x, N_y, subsize))

        ! Calculation
        do k=1, subsize
            call Greensfunction_PlanewaveToReal(recvbuf(:, :, k), nx_grid, ny_grid, sendbuf(:, :, k))

            if(rank == 0) then
                open(unit=16, file="OUTPUT", status="old", position="append")
                write(16, '(A2, I6, A2, I6)') "=>", k * mpi_size, " /", subsize * mpi_size
                close(16)
            end if
        end do

        ! Gather from different rank
        deallocate(recvbuf)
        allocate(recvbuf(N_x, N_y, subsize * mpi_size))
        call MPI_Gather(sendbuf, size(sendbuf), MPI_DOUBLE_COMPLEX, &
                        recvbuf, size(sendbuf), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, STATUS) !%%
        if(rank == 0) then
            do k=1, N_z
                do j=1, N_y
                    do i=1, N_x
                        Density(i, j, k) = aimag(recvbuf(i, j, k))
                    end do
                end do
            end do
        end if

        ! Done
        call cpu_time(PtoR_time%end)
        PtoR_time%sum = PtoR_time%sum + PtoR_time%end - PtoR_time%start

        if(rank == 0) then
            open(unit=16, file="OUTPUT", status="old", position="append")
            write(16, *) "Density calculation DONE"
            write(16, *) "----------"
            close(16)
        end if
    end if

    ! ===STEP 4. Use NEGF to find 'Transmission coefficient', parallelized===
    if(if_transmission) then
        ! Transmission-energy-grid layout
        N_E = nint(TRANSMISSION_GRID(3))
        allocate(Transmission(N_E))
        do i=1, N_E
            Transmission(i)%energy = atomic%TRANSMISSION_GRID(1) + &
            (atomic%TRANSMISSION_GRID(2) - atomic%TRANSMISSION_GRID(1)) * (i - 1) / N_E
            !                                                               ^ Last point is not included
        end do

        if(rank == 0) then
            open(unit=16, file="OUTPUT", status="old", position="append")
            write(16, *) "Transmission calculation start..."
            write(16, *) "Number of energy grid points:", N_E
            close(16)
        end if

        Transmission(:)%tau = 0.D0
        N_job = N_E * size(kpoint)
        i_job = 1 + rank
        do while(i_job <= N_job)
            call Transmission_Coefficient(i_job, V_reciprocal_all, nx_grid, ny_grid, atomic,&
            kpoint, Transmission, trans_time)

            i_job = i_job + mpi_size

            if(rank == 0) then
                open(unit=16, file="OUTPUT", status="old", position="append")
                write(16, '(A2, I6, A2, I6)') "=>", min(i_job - 1, N_job), " /", N_job
                close(16)
            end if
        end do

        ! Reduce from different rank
        if(rank == 0) then !%%
            call MPI_Reduce(MPI_IN_PLACE, Transmission%tau, size(Transmission%tau), &
            MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, STATUS) !%%
        else !%%
            call MPI_Reduce(Transmission%tau, Transmission%tau, size(Transmission%tau), &
            MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, STATUS) !%%
        end if !%%

        if(rank == 0) then
            open(unit=16, file="OUTPUT", status="old", position="append")
            write(16, *) "Transmission calculation DONE"
            write(16, *) "----------"
            close(16)
        end if
    end if

    
    if(if_LDOS) then
        ! ===STEP 5. Use NEGF to find 'Local Density Of State', parallelized===
        ! Layout LDOS grid
        N_E = nint(LDOS_ENERGY_GRID(3))
        allocate(LDOS_Matrix(N, N, N_z, N_E))
        allocate(LDOS(LDOS_size(LDOS_GRID(1), N_x), LDOS_size(LDOS_GRID(2), N_y), &
         LDOS_size(LDOS_GRID(3), N_z), N_E))

        if(rank == 0) then
            open(unit=16, file="OUTPUT", status="old", position="append")
            write(16, *) "LDOS calculation start..."
            write(16, *) "Number of energy grid points:", N_E
            close(16)
        end if
    
        ! MPI calculation
        LDOS_Matrix(:, :, :, :) = 0.D0
        N_job = N_E * size(kpoint)
        i_job = 1 + rank
        do while(i_job <= N_job)
            call Local_Density_Of_State(i_job, V_reciprocal_all, nx_grid, ny_grid, atomic,&
            kpoint, LDOS_Matrix, LDOS_inverse_time)

            i_job = i_job + mpi_size
            
            if(rank == 0) then
                open(unit=16, file="OUTPUT", status="old", position="append")
                write(16, '(A2, I6, A2, I6)') "=>", min(i_job - 1, N_job), " /", N_job
                close(16)
            end if
        end do

        ! Reduce from different rank
        if(rank == 0) then !%%
            call MPI_Reduce(MPI_IN_PLACE, LDOS_Matrix, size(LDOS_Matrix), &
            MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, STATUS) !%%
        else !%%
            call MPI_Reduce(LDOS_Matrix, LDOS_Matrix, size(LDOS_Matrix), &
            MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, STATUS) !%%
        end if !%%

        ! ===STEP 6. Transform LDOS_Matrix from plane wave basis to real space, parallelized===
        if(rank == 0) then
            open(unit=16, file="OUTPUT", status="old", position="append")
            write(16, *) "Transform LDOS from plane wave basis to real space:"
            close(16)
        end if

        ! Scatter to different rank
        call cpu_time(LDOS_PtoR_time%start)
        subsize = ((size(LDOS_Matrix, 3) + size(LDOS_Matrix, 4) - 1) / mpi_size + 1)
        allocate(sendbuf(N, N, subsize * mpi_size), recvbuf(N, N, subsize))
        if(rank == 0) then
            sendbuf(:, :, :) = 0.D0
            do i_E=1, N_E
                do k=1, N_z
                    do j=1, N
                        do i=1, N
                            sendbuf(i, j, k + (i_E - 1) * N_z) = LDOS_Matrix(i, j, k, i_E)
                        end do
                    end do
                end do
            end do
        end if
        call MPI_Scatter(sendbuf, size(recvbuf), MPI_DOUBLE_COMPLEX, &
                        recvbuf, size(recvbuf), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, STATUS) !%%
        deallocate(sendbuf)
        allocate(sendbuf(N_x, N_y, subsize))

        ! Calculation
        do k=1, subsize
            call Greensfunction_PlanewaveToReal(recvbuf(:, :, k), nx_grid, ny_grid, sendbuf(:, :, k))

            if(rank == 0) then
                open(unit=16, file="OUTPUT", status="old", position="append")
                write(16, '(A2, I6, A2, I6)') "=>", k * mpi_size, " /", subsize * mpi_size
                close(16)
            end if
        end do

        ! Gather from different rank
        deallocate(recvbuf)
        allocate(recvbuf(N_x, N_y, subsize * mpi_size))
        call MPI_Gather(sendbuf, size(sendbuf), MPI_DOUBLE_COMPLEX, &
                        recvbuf, size(sendbuf), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, STATUS) !%%
        
        LDOS(:, :, :, :) = 0.D0
        if(rank == 0) then
            do i_E=1, N_E
                do k=1, N_z
                    do j=1, N_y
                        do i=1, N_x
                            rescale_factor = 1.D0
                            call LDOS_grid_converter(i, LDOS_GRID(1), N_x, ii, rescale_factor)
                            call LDOS_grid_converter(j, LDOS_GRID(2), N_y, jj, rescale_factor)
                            call LDOS_grid_converter(k, LDOS_GRID(3), N_z, kk, rescale_factor)

                            LDOS(ii, jj, kk, i_E) = LDOS(ii, jj, kk, i_E) + &
                             aimag(recvbuf(i, j, k + (i_E - 1) * N_z)) * rescale_factor
                        end do
                    end do
                end do
            end do
        end if

        ! Done
        call cpu_time(LDOS_PtoR_time%end)
        LDOS_PtoR_time%sum = LDOS_PtoR_time%sum + LDOS_PtoR_time%end - LDOS_PtoR_time%start

        if(rank == 0) then
            open(unit=16, file="OUTPUT", status="old", position="append")
            write(16, *) "LDOS calculation DONE"
            write(16, *) "----------"
            close(16)
        end if

    end if

    !===========================================END==============================================
    !============================================================================================
    
    if(rank == 0) then

        !　Write out data
        open(unit=16, file="OUTPUT", status="old", position="append")
        write(16, *) ""
        write(16, *) "===== ALL NEGF calculation DONE ====="
        write(16, *) ""
        
        if(if_density) then
            ! Rescale Density (Unit: 1/Angstrom^3)
            rescale_factor = dble(N_z) / (LX * LY * LZ)
            do k=1, N_z
                do j=1, N_y
                    do i=1, N_x
                        Density(i, j, k) = Density(i, j, k) * rescale_factor ! unit: 1/Angstrom^3
                    end do
                end do
            end do
            ! Write
            write(16, *) "Total charge is:", sum(Density) * (LX * LY * LZ) / (N_x * N_y * N_z)
            write(16, *) "Writing DENSITY..."
            open(unit=18, file="DENSITY")
            write(18, '(3I5)') N_x, N_y, N_z
            write(18, '(5G17.8)') (((Density(i, j, k), i=1, N_x), j=1, N_y), k=1, N_z)
            close(18)
        end if

        if(if_transmission) then
            ! Convert unit
            do i= 1, size(Transmission)
                Transmission(i)%energy = Transmission(i)%energy * hartree ! unit: eV
            end do
            ! Write
            write(16, *) "Writing TRANSMISSION..."
            open(unit=19, file="TRANSMISSION")
            write(19, *) "   Energy (eV)    Transmission_coefficient"
            write(19, '(2G17.8)') (Transmission(i)%energy, Transmission(i)%tau, i=1, size(Transmission))
            close(19)
        end if

        if(if_LDOS) then
            ! Rescale LDOS (Unit: 1/Angstrom^3)
            rescale_factor = dble(N_z) / (LX * LY * LZ)
            do i_E=1, size(LDOS, 4)
                do k=1, size(LDOS, 3)
                    do j=1, size(LDOS, 2)
                        do i=1, size(LDOS, 1)
                            LDOS(i, j, k, i_E) = LDOS(i, j, k, i_E) * rescale_factor ! unit: 1/Angstrom^3
                        end do
                    end do
                end do
            end do
            ! Write
            write(16, *) "Writing LDOS..."
            open(unit=20, file="LDOS")
            write(20, '(4I5)') size(LDOS, 1), size(LDOS, 2), size(LDOS, 3), size(LDOS, 4)
            write(20, '(5G17.8)') ((((LDOS(i, j, k, i_E), i=1, size(LDOS, 1)), j=1, size(LDOS, 2)),&
             k=1, size(LDOS, 3)), i_E=1, size(LDOS, 4))
            close(20)
        end if

        call cpu_time(total_time%end)
        total_time%sum = total_time%sum + total_time%end - total_time%start
        
        write(16,*) "CPU time (Total): ", total_time%sum, "s"
        write(16,*) " ├── CPU time (Real space to Plane wave): ", RtoP_time%sum, "s"
        write(16,*) " ├── CPU time (Density): ", inverse_time%sum + PtoR_time%sum, "s"
        write(16,*) " │   ├── CPU time (Inverse Matrix): ", inverse_time%sum, "s"
        write(16,*) " │   └── CPU time (Plane wave to Real space): ", PtoR_time%sum, "s"
        write(16,*) " ├── CPU time (Transmission (Inverse Matrix)): ", trans_time%sum, "s"
        write(16,*) " ├── CPU time (LDOS): ", LDOS_inverse_time%sum + LDOS_PtoR_time%sum, "s"
        write(16,*) " │   ├── CPU time (Inverse Matrix): ", LDOS_inverse_time%sum, "s"
        write(16,*) " │   └── CPU time (Plane wave to Real space): ", LDOS_PtoR_time%sum, "s"
        write(16,*) " └── CPU time (Others): ", &
         total_time%sum - RtoP_time%sum - inverse_time%sum - PtoR_time%sum - &
         trans_time%sum - LDOS_inverse_time%sum - LDOS_PtoR_time%sum, "s"
        write(16,*) "********************"
        write(16,*) "*                  *"
        write(16,*) "*      DONE        *"
        write(16,*) "*                  *"
        write(16,*) "********************"
        close(16)
    end if

    call MPI_FINALIZE(STATUS) !%%
end program main