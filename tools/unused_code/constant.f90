    ! Unit converter, and some physical constant
    !   ( Reference: https://physics.nist.gov/cuu/Constants/ )
    real*8 :: hartree, k_B, a_0, pi
    hartree = 27.211386245988D0 ! in eV
    k_B = 8.617333262D-5 ! in eV
    k_B = k_B / hartree ! convert to atomic unit
    a_0 = 0.529177210903D0 ! Bohr radius (in Angstrom)
    pi = 3.141592653589793D0


