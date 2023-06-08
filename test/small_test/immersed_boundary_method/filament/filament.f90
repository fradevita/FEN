program spring_mass

    ! Test for the filament in a uniform current

    use mpi
    use precision_mod        , only : dp
    use global_mod           , only : ierror, myrank, pi
    use grid_mod
    use vector_mod
    use ibm_mod              , only : lagrangian_solid_list
    use navier_stokes_mod    , only : viscosity, g, v, set_timestep, dt_o
    use IO_mod               , only : stdout
    use solver_mod
    use fields_mod           , only : curl
    use lagrangian_solid_mod
    use lagrangian_ibm_mod   , only : Nfstep
    
    implicit none

    ! Parameters
    real(dp), parameter :: U = 1.0_dp          ! Uniform velocity
    real(dp), parameter :: L = 1.0_dp          ! Filament length
    real(dp), parameter :: rhof = 1.0_dp       ! Fluid density
    real(dp), parameter :: Re = 200.0_dp       ! Reynolds number = UL/nu
    real(dp), parameter :: Fr = 1.414_dp       ! Froud nmber = U/sqrt(gL)
    real(dp), parameter :: gamma = 150.0_dp    ! density ratio (solid/fluid)
    real(dp), parameter :: eps = 25000000.0_dp ! Non-dimensional elastic modulus
    real(dp), parameter :: beta = 0.0015_dp    ! Non-dimensional bending stiffness
    real(dp), parameter :: h = 0.01_dp         ! Filament thickness
    real(dp), parameter :: w = 1.0_dp          ! Unit width (2D simulation)

    ! Get derived parameters
    real(dp), parameter :: nu = U*L/Re             ! Fluid kinematic viscosity
    real(dp), parameter :: gravity = (U/Fr)**2/L   ! Gravitational acceleration
    real(dp), parameter :: rhos = gamma*rhof       ! Solid density
    real(dp), parameter :: E = eps*rhof*U**2       ! Elastic modulus
    real(dp), parameter :: B = beta*rhof*U**2*L**3 ! Bending modulus

    ! Variables
    integer             :: Nx, Ny, Nz, step, out_id, fid, lid
    real(dp)            :: Lx, Ly, Lz, Ep0, L0, dt, time, Lt
    character(len=7)    :: sn
    character(len=16)   :: filename
    type(grid)          :: comp_grid
    type(solid), target :: S
    type(bc_type)       :: bc(4)
    type(vector)        :: omega

    ! Initialize MPI
    call mpi_init(ierror)
    call mpi_comm_rank(mpi_comm_world, myrank, ierror)

    ! The domain is a square box of size 8
    Lx = 8.0_dp
    Ly = 8.0_dp

    ! Set the resolution
    Nx = 1024
    Ny = 1024

    ! Since 2D
    Nz = 1
    Lz = Lx*float(Nz)/float(Nx)

    bc(1)%s = 'Wall'
    bc(2)%s = 'Wall'
    bc(3)%s = 'Inflow'
    bc(4)%s = 'Outflow'

    ! Create the grid
    call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, [0.0_dp, 0.0_dp, 0.0_dp], 8, 1, bc)

    ! Create the filament
    ! Set the mass of the filament
    ! rhos = gamma*rhof
    ! volume = h*L*w
    S%Vol = L*h*w
    S%rho = rhos
    S%M = rhos*L*h*w

    ! Turn on the deformable flag
    S%is_deformable = .true.

    ! The filament is an open geometry
    S%is_open = .true.

    ! Create the lagrangian solid from the mesh file
    call S%create('mesh.txt')

    ! Half the last mass since it is connected to one single edge
    S%mass_points(S%number_of_mass_points)%M = S%mass_points(S%number_of_mass_points)%M*0.5_dp

    ! Set the elastic in-plane constant of the springs
    ! eq. 25 of de Tullio and Pascazio JCP 2016.
    ! Eh = eps*rhof*U**2*L
    ! Ke = Eh*SumAi/l**2
    ! A = l*1
    ! l = L/number_of_mass_points
    S%ke = E*h*real(S%number_of_edges, dp)/L
    S%edges%ke = E*h*real(S%number_of_edges, dp)/L

    ! Set the bending constant
    ! eq. 30 of de Tullio and Pascazio JCP 2016.
    ! kb = B*2/sqrt(3)
    ! B = beta*rhof*U**2*L**3
    S%kb = B*real(S%number_of_mass_points, dp)

    ! Compute the initial length of the filament
    L0 = S%get_total_length()

    ! Compute initial potential energy
    Ep0 = S%get_potential_energy()

    ! Set the constraints procedure
    S%apply_constraints => test_constraints

    ! Allocate the array of solid
    allocate(lagrangian_solid_list(1))
    lagrangian_solid_list(1)%pS => S

    ! Set the viscosity
    viscosity = rhof*nu

    ! Set the body force
    g(2) = gravity

    ! Initialize the solver
    call init_solver(comp_grid)
    step = 0
    time = 0.0_dp

    ! Explicitly set the timestep
    dt = 1.0e-3_dp
    dt_o = dt
    
    ! Set boundary conditions
    v%y%bc%bottom = U
    v%y%bc%left = U
    v%y%bc%right = U
    
    ! Print the poitns
    call S%print_configuration(step)

    ! Open output file
    open(newunit = out_id, file = 'out.txt')
    open(newunit = fid   , file = 'forces.txt')
    open(newunit = lid   , file = 'length.txt')

    ! Allocate the vorticity vector
    call omega%allocate(comp_grid, 1)
    
    Nfstep = 3

    !==== Start Time loop =========================================================================
    time_loop: do while(time < 30)

        step = step + 1
        time = time + dt

        ! Advance in time the solution
        call advance_solution(comp_grid, step, dt)

        ! Advance solver status to log file
        call print_solver_status(stdout, step, time, dt)

        ! Output
        if (myrank == 0) write(out_id,*) time, S%mass_points(S%number_of_mass_points)%X(1:2)
        if (mod(step,100) == 0) then
            call S%integrate_hydrodynamic_forces
            if (myrank == 0) write(fid,*) time, S%center_of_mass%Fh(:)
            if (myrank == 0) call S%print_configuration(step)
            Lt = S%get_total_length()
            if (myrank == 0) write(lid,*) time, Lt, L0
            write(sn,'(I0.7)') step
            call curl(v, omega)
            filename = 'data/vrt_'//sn
            call omega%x%write(filename)
        end if

    end do time_loop

    ! Free memory
    call destroy_solver

    ! Finalize the simulation
    call MPI_FINALIZE(ierror)

contains
    
    !========================================================================================
    subroutine test_constraints(self)
        
        ! In/Out variables
        class(solid), intent(inout) :: self

        ! Clamped condition at X = 0
        self%mass_points(1)%X(1) = Lx/2.0_dp
        self%mass_points(1)%X(2) = 2.0_dp
        self%mass_points(1)%V = 0.0_dp
        self%mass_points(1)%A = 0.0_dp

    end subroutine test_constraints
    !========================================================================================

end program spring_mass
