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
    real*8 :: V_L, V_R, EPSILON, TEMPERATURE, LX, LY, LZ, ENCUT
    integer :: NGX, NGY
    namelist /INPUT/ V_L, V_R, EPSILON, TEMPERATURE, LX, LY, LZ, NGX, NGY, ENCUT

    integer :: N_x, N_y, N_z, N, i, j, k, STATUS
    real*8 :: kx, ky
    real*8, allocatable :: V_real(:, :, :)
    complex*16, allocatable :: V_reciprocal(:, :), Hamiltonian(:, :, :)
    integer, allocatable :: nx_grid(:), ny_grid(:)

    include 'constant.f90'

    open(unit=16, file="OUTPUT")
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
    EPSILON = 1.D-5 ! Hartree
    TEMPERATURE = 300.D0 ! K
    LX = 0.D0 ! Angstrom
    LY = 0.D0 ! Angstrom
    LZ = 0.D0 ! Angstrom
    
    NGX = 0
    NGY = 0

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
    EPSILON = EPSILON ! Hartree
    TEMPERATURE = TEMPERATURE ! K
    LX = LX / a_0 ! Angstrom -> bohr
    LY = LY / a_0 ! Angstrom -> bohr
    LZ = LZ / a_0 ! Angstrom -> bohr



    !===Grid layout===
    write(16, *) "----------"
    write(16, *) "Layout the grid..."

    allocate(V_reciprocal(-NGX/2: NGX/2, -NGY/2: NGY/2))
    N = PlaneWaveBasis_construction_findsize(ENCUT, Lx, Ly)
    allocate(nx_grid(N), ny_grid(N), Hamiltonian(N, N, N_z))
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
    
    !===Main NEGF calculation===

    


    !===========================================END==============================================
    !============================================================================================

    close(16)
end program main