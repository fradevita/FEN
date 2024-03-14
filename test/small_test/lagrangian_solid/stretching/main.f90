program main

    ! The script reproduces figures 4.11 and 4.12 of 
    ! Tanaka, Wada and Nakamura, Computational Biomechanics (2012)

    use precision_mod
    use lagrangian_solid_mod
    use lagrangian_solid_2D_mod
    use euclidean_mod

    implicit none

    ! Parameters
    real(dp), parameter :: L0 = 20.0e-3_dp         ! Lenght of the specimen [m]
    real(dp), parameter :: D0 = 5.0e-3_dp          ! Width of the specimen  [m]
    real(dp), parameter :: eps_x = 0.02_dp         ! stretching in x dir     []
    real(dp), parameter :: L = L0*(1.0_dp + eps_x) ! stretched length       [m]
    real(dp), parameter :: dt = 1.0e-2_dp          ! timestep               [s]
    real(dp), parameter :: tol = 1.0e-12_dp        ! tolerance for steady state

    ! Variables
    integer                   :: n, step, nbcl, nbcr, nbct, nbcb, e_id, f_id, d_id
    integer     , allocatable :: ibcl(:), ibcr(:), ibct(:), ibcb(:)
    real(dp)                  :: ks, ka, W_o, W, ymax, ymin, D
    type(lagrangian_solid_2D) :: S
    character(len=1 )         :: arg1, arg2

    ! Load the solid from the file
    call S%loadFromFile('mesh.txt')
    S%is_open = .true.
    S%is_deformable = .true.

    ! Get command line arguments (ka and ks)
    call get_command_argument(1, arg1)
    call get_command_argument(2, arg2)
 
    select case(arg1)
    case('0')
        ks = 1.0_dp
    case('1')
        ks = 1.1_dp
    case('2')
        ks = 1.2_dp
    case('3')
        ks = 1.3_dp
    case('4')
        ks = 1.4_dp
    case('5')
        ks = 1.5_dp
    end select

    select case(arg2)
    case('0')
        ka = 0.0_dp
    case('1')
        ka = 0.2_dp
    case('2')
        ka = 0.4_dp
    end select

    ! Set material properties
    call S%setMass(1.0_dp*S%number_of_mass_points)
    if (ks > 0._dp) call S%setStretchingConstant(ks)
    if (ka > 0._dp) call S%setAreaConstant(ka)

    ! Set damping using damping ratio equal to 1 (critically damped system)
    S%gamma = 2.0_dp*sqrt(S%mass_points(1)%m*ks)

    ! This requires the explicit implicit solver for the deformation
    S%advance_deformation => expImp

    ! Find boundary points
    nbcl = 0
    nbcr = 0
    nbct = 0
    nbcb = 0
    do n = 1,S%number_of_mass_points
        if (S%mass_points(n)%X(1) == 0.0_dp) nbcl = nbcl + 1
        if (S%mass_points(n)%X(1) == L0) nbcr = nbcr + 1
        if (S%mass_points(n)%X(2) == 0.0_dp) nbcb = nbcb + 1
        if (S%mass_points(n)%X(2) == D0) nbct = nbct + 1
    end do
    allocate(ibcl(nbcl))
    allocate(ibcr(nbcr))
    allocate(ibcb(nbcb))
    allocate(ibct(nbct))

    nbcl = 0
    nbcr = 0
    nbct = 0
    nbcb = 0
    do n = 1,S%number_of_mass_points
        if (S%mass_points(n)%X(1) == 0.0_dp) then
            nbcl = nbcl + 1
            ibcl(nbcl) = n
        end if
        if (S%mass_points(n)%X(1) == L0) then
            nbcr = nbcr + 1
            ibcr(nbcr) = n
        endif
        if (S%mass_points(n)%X(2) == 0.0_dp) then
            nbcb = nbcb + 1
            ibcb(nbcb) = n
        end if
        if (S%mass_points(n)%X(2) == D0) then
            nbct = nbct + 1
            ibct(nbct) = n
        endif
    end do
    S%apply_Constraints => test_constraints

    ! Apply the stretching
    do n = 1,S%number_of_mass_points
        S%mass_points(n)%X(1) = S%mass_points(n)%X(1)*(1.0_dp + eps_x)
    end do
    call S%apply_Constraints()
    call S%update()

    ! Open output file
    open(newunit = e_id, file = 'data/energy.dat')
    open(newunit = f_id, file = 'data/forces.dat')

    ! Print initial configuration
    call S%writeSTL('data/ini.stl')
    call S%writeGNU('data/ini.gnu')

    ! Initial potential energy
    W_o = S%getPotentialEnergy()

    ! Start time loop
    time_loop: do step = 1,2000000
        
        ! Update solution
        call S%advance(dt, [0.0_dp, 0.0_dp, 0.0_dp])
        
        ! Output energy and forces
        W = S%getPotentialEnergy()
        if (mod(step,10) == 0) then
            write(e_id,*) step*dt, S%getKineticEnergy(), W
            write(f_id,*) step*dt, S%internalForcesIntegral()
            flush(e_id)
            flush(f_id)
        endif

        ! Check steady state on potential energy value
        if (abs(W_o - W)/W_o < tol) exit time_loop
        print *, step, abs(W - W_o)/W_o
        W_o = W

    end do time_loop

    ! Compute Poisson ration for a y-section at L/2
    ymin = 1.0_dp
    ymax = -1.0_dp
    do n = 1,S%number_of_mass_points
        if (abs(S%mass_points(n)%X(1) - L/2.0_dp) < 1.0e-6_dp) then

            if (S%mass_points(n)%X(2) > ymax) ymax = S%mass_points(n)%X(2)
            if (S%mass_points(n)%X(2) < ymin) ymin = S%mass_points(n)%X(2)
        end if
    end do
    D = ymax - ymin

    ! Compute Poisson ratio using the mean value of ymax and ymin
    ymin = 0.0_dp
    ymax = 0.0_dp
    do n = 1,nbcb
        ymin = ymin + S%mass_points(ibcb(n))%X(2)
    end do
    ymin = ymin/nbcb
    do n = 1,nbct
        ymax = ymax + S%mass_points(ibct(n))%X(2)
    end do
    ymax = ymax/nbct

    open(newunit = d_id, file = 'data/D.dat')
    write(d_id,*) ymax, ymin, D, (D0 - D)/D0, (D0 - ymax+ymin)/D0

    ! Print final configuration
    call S%writeSTL('data/end.stl')
    call S%writeGNU('data/end.gnu')

    ! Free memory
    call S%destroy()

contains

    !===========================================================================
    subroutine test_constraints(self)

        class(lagrangian_solid), intent(inout), target :: self

        ! x = 0
        do n = 1,nbcl
            self%mass_points(ibcl(n))%X(1) = 0.0_dp
        end do

        ! x = L
        do n = 1,nbcr
            self%mass_points(ibcr(n))%X(1) = L
        end do

    end subroutine
    !===========================================================================

end program main
