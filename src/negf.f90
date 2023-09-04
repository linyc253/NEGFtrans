module negf
    use math_kernel
    use grid
    use global
    implicit none
    private
    public Equilibrium_Density, Transmission_Coefficient, Local_Density_Of_State
contains
    subroutine Equilibrium_Density(i_job, V_reciprocal_all, nx_grid, ny_grid, N_circle, N_line, min_V, atomic,&
        kpoint, Density_Matrix, inverse_time)
        real*8, intent(in) :: min_V
        complex*16, intent(in), allocatable :: V_reciprocal_all(:, :, :)
        integer, intent(in) :: i_job, N_circle, N_line, nx_grid(:), ny_grid(:)
        type(t_parameters), intent(in) :: atomic
        type(t_timer), intent(inout) :: inverse_time
        type(t_kpointmesh) :: kpoint(:)
        complex*16, intent(inout) :: Density_Matrix(:, :, :)

        integer :: N_z, N, i, j, k, i_integral, i_kpoint, N_integral
        complex*16, allocatable :: Hamiltonian(:, :, :), E_minus_H(:, :, :)
        type(t_gfunc), allocatable :: G_Function(:, :, :)
        real*8 :: circle_L, circle_R, line_L, line_R, E_c, E_R, theta, kx, ky
        complex*16 :: energy

        ! Distribute i_job
        N_integral = N_circle + N_line
        i_integral = 1 + mod(i_job - 1, N_integral)
        i_kpoint = 1 + (i_job - 1) / N_integral ! fractional part (remainder) is discarded

        ! Allocate arrays
        N_z = size(Density_Matrix, 3)
        N = size(nx_grid)
        allocate(Hamiltonian(N, N, N_z), E_minus_H(N, N, N_z))
        allocate(G_function(N, N, N_z))
        
        ! Determine integral parameters
        circle_L = min_V - 1.D0 / hartree
                        !  ^ This improve precision significantly 
        circle_R = atomic%MU - atomic%GAP
        line_L = atomic%MU - atomic%GAP
        line_R = atomic%MU + atomic%GAP

        E_c = (circle_R + circle_L) / 2
        E_R = (circle_R - circle_L) / 2
        
        if((i_job == 1) .and. DEBUG) then
            open(unit=16, file="OUTPUT", status="old", position="append")
            write(16, *) "DEBUG: circle_L : circle_R = ", circle_L, ":", circle_R
            write(16, *) "DEBUG: line_L : line_R = ", line_L, ":", line_R
            close(16)
        end if

        ! Construct Hamiltonian
        kx = kpoint(i_kpoint)%kx
        ky = kpoint(i_kpoint)%ky
        call Hamiltonian_construction(nx_grid, ny_grid, atomic%LX, atomic%LY, kx, ky, V_reciprocal_all, Hamiltonian)

        if(i_integral <= N_circle) then
            ! CASE1
            theta = value_simpson(i_integral, N_circle, 0.D0, pi)
            energy = E_c + E_R * exp(complex(0.D0, theta))
        else
            ! CASE2
            energy = value_simpson(i_integral - N_circle, N_line, line_L, line_R)
        end if

        ! Construct (EI - Hamiltonian)
        call E_minus_H_construction(Hamiltonian, energy + complex(0.D0, atomic%ETA), nx_grid, ny_grid, atomic%LX, &
         atomic%LY, atomic%LZ, kx, ky, atomic%V_L, atomic%V_R, E_minus_H)

        ! Solve for G_Function = (EI - Hamiltonian)^{-1}
        call cpu_time(inverse_time%start)
        call GreensFunction_tri_solver(E_minus_H, G_Function)
        call cpu_time(inverse_time%end)
        inverse_time%sum = inverse_time%sum + inverse_time%end - inverse_time%start

        ! Rescale to compensate the extra factor in E_minus_H
        do k=1, N_z
            do j=1, N
                do i=1, N
                    G_Function(i, j, k)%diagonal = G_Function(i, j, k)%diagonal         * 2.D0 * (atomic%LZ / N_z) ** 2
                    G_Function(i, j, k)%first_column = G_Function(i, j, k)%first_column * 2.D0 * (atomic%LZ / N_z) ** 2
                    G_Function(i, j, k)%last_column = G_Function(i, j, k)%last_column   * 2.D0 * (atomic%LZ / N_z) ** 2
                end do
            end do
        end do
        
        if(i_integral <= N_circle) then
            ! CASE1
            do k=1, N_z
                do j=1, N
                    do i=1, N
                        Density_Matrix(i, j, k) = Density_Matrix(i, j, k) + &
                         kpoint(i_kpoint)%weight * &
                         2.D0 / pi * coefficient_simpson(i_integral, N_circle, 0.D0, pi) * &
                         complex(0.D0, 1.D0) * E_R * exp(complex(0.D0, theta)) * G_Function(i, j, k)%diagonal
                    end do
                end do
            end do
        else
            ! CASE2
            do k=1, N_z
                do j=1, N
                    do i=1, N
                        Density_Matrix(i, j, k) = Density_Matrix(i, j, k) - &
                         kpoint(i_kpoint)%weight * &
                         2.D0 / pi * coefficient_simpson(i_integral - N_circle, N_line, line_L, line_R) * &
                         fermi_func(dble(energy), atomic%MU, atomic%TEMPERATURE) * G_Function(i, j, k)%diagonal
                    end do
                end do
            end do

        end if

    end subroutine Equilibrium_Density
    
    subroutine Transmission_Coefficient(i_job, V_reciprocal_all, nx_grid, ny_grid, atomic,&
        kpoint, Transmission, trans_time)
        complex*16, intent(in), allocatable :: V_reciprocal_all(:, :, :)
        integer, intent(in) :: i_job, nx_grid(:), ny_grid(:)
        type(t_parameters), intent(in) :: atomic
        type(t_timer), intent(inout) :: trans_time
        type(t_kpointmesh) :: kpoint(:)
        type(t_transmission), intent(inout) :: Transmission(:)

        integer :: i, j, k, N, N_z, N_E, i_energy, i_kpoint
        complex*16, allocatable :: Hamiltonian(:, :, :), E_minus_H(:, :, :)
        type(t_gfunc), allocatable :: G_Function(:, :, :)
        real*8 :: kx, ky, tau
        complex*16 :: energy

        ! Distribute i_job
        N_E = size(Transmission)
        i_energy = 1 + mod(i_job - 1, N_E)
        i_kpoint = 1 + (i_job - 1) / N_E ! fractional part (remainder) is discarded

        ! Allocate arrays
        N = size(nx_grid)
        N_z = size(V_reciprocal_all, 3)
        energy = Transmission(i_energy)%energy
        allocate(Hamiltonian(N, N, N_z), E_minus_H(N, N, N_z))
        allocate(G_function(N, N, N_z))

        ! Construct Hamiltonian
        kx = kpoint(i_kpoint)%kx
        ky = kpoint(i_kpoint)%ky
        call Hamiltonian_construction(nx_grid, ny_grid, atomic%LX, atomic%LY, kx, ky, V_reciprocal_all, Hamiltonian)


        ! Construct (EI - Hamiltonian)
        call E_minus_H_construction(Hamiltonian, energy + complex(0.D0, atomic%ETA), nx_grid, ny_grid, atomic%LX, &
         atomic%LY, atomic%LZ, kx, ky, atomic%V_L, atomic%V_R, E_minus_H)

        ! Solve for G_Function = (EI - Hamiltonian)^{-1}
        call cpu_time(trans_time%start)
        call GreensFunction_tri_solver(E_minus_H, G_Function)
        call cpu_time(trans_time%end)
        trans_time%sum = trans_time%sum + trans_time%end - trans_time%start

        ! Rescale to compensate the extra factor in E_minus_H
        do k=1, N_z
            do j=1, N
                do i=1, N
                    G_Function(i, j, k)%diagonal = G_Function(i, j, k)%diagonal         * 2.D0 * (atomic%LZ / N_z) ** 2
                    G_Function(i, j, k)%first_column = G_Function(i, j, k)%first_column * 2.D0 * (atomic%LZ / N_z) ** 2
                    G_Function(i, j, k)%last_column = G_Function(i, j, k)%last_column   * 2.D0 * (atomic%LZ / N_z) ** 2
                end do
            end do
        end do


        ! Calculate Transmission
        call Transmission_solver(G_Function, Transmission(i_energy)%energy, nx_grid, ny_grid, atomic%LX, &
        atomic%LY, atomic%LZ, kx, ky, atomic%V_L, atomic%V_R, tau)

        Transmission(i_energy)%tau = Transmission(i_energy)%tau + tau * kpoint(i_kpoint)%weight
        
    end subroutine Transmission_Coefficient

    subroutine Local_Density_Of_State(i_kpoint, i_energy, N_E, V_reciprocal_all, nx_grid, ny_grid, atomic,&
        kpoint, LDOS_Matrix, LDOS_inverse_time)
        complex*16, intent(in), allocatable :: V_reciprocal_all(:, :, :)
        integer, intent(in) :: i_kpoint, i_energy, N_E, nx_grid(:), ny_grid(:)
        type(t_parameters), intent(in) :: atomic
        type(t_timer), intent(inout) :: LDOS_inverse_time
        type(t_kpointmesh) :: kpoint(:)
        complex*16, intent(inout) :: LDOS_Matrix(:, :, :)

        integer :: N_z, N, i, j, k
        complex*16, allocatable :: Hamiltonian(:, :, :), E_minus_H(:, :, :)
        type(t_gfunc), allocatable :: G_Function(:, :, :)
        real*8 :: kx, ky
        complex*16 :: energy

        ! Allocate arrays
        N_z = size(V_reciprocal_all, 3)
        N = size(nx_grid)
        allocate(Hamiltonian(N, N, N_z), E_minus_H(N, N, N_z))
        allocate(G_function(N, N, N_z))

        ! Construct Hamiltonian
        kx = kpoint(i_kpoint)%kx
        ky = kpoint(i_kpoint)%ky
        call Hamiltonian_construction(nx_grid, ny_grid, atomic%LX, atomic%LY, kx, ky, V_reciprocal_all, Hamiltonian)

        ! Get energy value
        energy = atomic%LDOS_ENERGY_GRID(1) + &
        (atomic%LDOS_ENERGY_GRID(2) - atomic%LDOS_ENERGY_GRID(1)) * (i_energy - 1) / N_E
        !                                                               ^ Last point is not included

        ! Construct (EI - Hamiltonian)
        call E_minus_H_construction(Hamiltonian, energy + complex(0.D0, atomic%ETA), nx_grid, ny_grid, atomic%LX, &
         atomic%LY, atomic%LZ, kx, ky, atomic%V_L, atomic%V_R, E_minus_H)

        ! Solve for G_Function = (EI - Hamiltonian)^{-1}
        call cpu_time(LDOS_inverse_time%start)
        call GreensFunction_tri_solver(E_minus_H, G_Function)
        call cpu_time(LDOS_inverse_time%end)
        LDOS_inverse_time%sum = LDOS_inverse_time%sum + LDOS_inverse_time%end - LDOS_inverse_time%start

        ! Rescale to compensate the extra factor in E_minus_H
        do k=1, N_z
            do j=1, N
                do i=1, N
                    G_Function(i, j, k)%diagonal = G_Function(i, j, k)%diagonal         * 2.D0 * (atomic%LZ / N_z) ** 2
                    G_Function(i, j, k)%first_column = G_Function(i, j, k)%first_column * 2.D0 * (atomic%LZ / N_z) ** 2
                    G_Function(i, j, k)%last_column = G_Function(i, j, k)%last_column   * 2.D0 * (atomic%LZ / N_z) ** 2
                end do
            end do
        end do

        ! Add contribution to LDOS_Matrix
        do k=1, N_z
            do j=1, N
                do i=1, N
                    LDOS_Matrix(i, j, k) = LDOS_Matrix(i, j, k) - &
                     kpoint(i_kpoint)%weight * 1.D0 / pi  * G_Function(i, j, k)%diagonal
                end do
            end do
        end do
        
    end subroutine Local_Density_Of_State

end module negf