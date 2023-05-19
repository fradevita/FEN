program spring_mass

    ! Test for the filament in a uniform current

    use precision        , only : dp
    use constants        , only : pi
    use class_Grid       , only : base_grid, bc_type
    use lagrangian_ibm   , only : number_of_solid_bodies, solid_array
    use navier_stokes    , only : viscosity, g, v, set_timestep, dt_o
    use io               , only : stdout
    use solver
    use class_Vector
    use fields           , only : curl
    use lagrangian_solid , only : apply_constraints
    
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
    integer  :: ierror, Nx, Ny, Nz, step, out_id, fid, lid
    real(dp) :: Lx, Ly, Lz, Ep0, L0, dt, time, Lt
    character(len=7) :: sn
    character(len=16) :: filename
    type(bc_type) :: bc(4)
    type(vector) :: omega

    ! Initialize MPI
    call mpi_init(ierror)

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
    call base_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, [0.0_dp, 0.0_dp, 0.0_dp], 8, 1, bc)

    ! Set one solid body
    number_of_solid_bodies = 1

    ! Allocate the array of solid
    allocate(solid_array(number_of_solid_bodies))

    ! Set the mass of the filament
    ! rhos = gamma*rhof
    ! volume = h*L*w
    solid_array(1)%Vol = L*h*w
    solid_array(1)%rho = rhos
    solid_array(1)%M = rhos*L*h*w

    ! Turn on the deformable flag
    solid_array(1)%is_deformable = .true.

    ! The filament is an open geometry
    solid_array(1)%is_open = .true.

    ! Create the lagrangian solid from the mesh file
    call solid_array(1)%create('mesh.txt')

    solid_array(1)%mass_points(solid_array(1)%number_of_mass_points)%M = &
        solid_array(1)%mass_points(solid_array(1)%number_of_mass_points)%M*0.5_dp

    ! Set the elastic in-plane constant of the springs
    ! eq. 25 of de Tullio and Pascazio JCP 2016.
    ! Eh = eps*rhof*U**2*L
    ! Ke = Eh*SumAi/l**2
    ! A = l*1
    ! l = L/number_of_mass_points
    !solid_array(1)%edges%ke = eps*rhof*U**2*L*2.0_dp*real(solid_array(1)%number_of_mass_points, dp)
    solid_array(1)%edges%ke = E*h*real(solid_array(1)%number_of_edges, dp)/L
    solid_array(1)%ke = E*h*real(solid_array(1)%number_of_edges, dp)/L


    ! Set the bending constant
    ! eq. 30 of de Tullio and Pascazio JCP 2016.
    ! kb = B*2/sqrt(3)
    ! B = beta*rhof*U**2*L**3
    !solid_array(1)%kb = beta*rhof*U**2*L**3*2.0_dp/sqrt(3.0_dp)
    solid_array(1)%kb = B*real(solid_array(1)%number_of_mass_points, dp)

    ! Compute the initial length of the filament
    L0 = solid_array(1)%get_total_length()

    ! Compute initial potential energy
    Ep0 = solid_array(1)%get_potential_energy()

    ! Set the constraints procedure
    apply_constraints => test_constraints

    ! Set the viscosity
    viscosity = rhof*nu

    ! Set the body force
    g(2) = gravity

    ! Initialize the solver
    call init_solver
    step = 0
    time = 0.0_dp

    ! Explicitly set the timestep
    dt = 1.0e-3_dp
    dt_o = dt
    
    ! Set boundary conditions
    v%y%bc%b = U
    v%y%bc%l = U
    v%y%bc%r = U
    
    ! Print the poitns
    call solid_array(1)%print_configuration(step)

    ! Open output file
    open(newunit = out_id, file = 'out.txt')
    open(newunit = fid   , file = 'forces.txt')
    open(newunit = lid   , file = 'length.txt')

    ! Allocate the vorticity vector
    call omega%allocate()
    
    !==== Start Time loop ===================================================================
    time_loop: do while(time < 30)

        step = step + 1
        time = time + dt

        ! Advance in time the solution
        call advance_solution(step, dt)

        ! Advance solver status to log file
        call print_solver_status(stdout, step, time, dt)

        ! Output
        if (base_grid%rank == 0) write(out_id,*) time, &
            solid_array(1)%mass_points(solid_array(1)%number_of_mass_points)%X(1:2)
        if (mod(step,100) == 0) then
            call solid_array(1)%integrate_hydrodynamic_forces
            if (base_grid%rank == 0) write(fid,*) time, solid_array(1)%center_of_mass%Fh(:)
            if (base_grid%rank == 0) call solid_array(1)%print_configuration(step)
            Lt = solid_array(1)%get_total_length()
            if (base_grid%rank == 0) write(lid,*) time, Lt, L0
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
    subroutine test_constraints(obj)

        use lagrangian_solid, only : solid
        
        type(solid), intent(inout) :: obj

        obj%mass_points(1)%X(1) = Lx/2.0_dp
        obj%mass_points(1)%X(2) = 2.0_dp
        obj%mass_points(1)%V = 0.0_dp
        obj%mass_points(1)%A = 0.0_dp
        
    end subroutine test_constraints
    !========================================================================================
  
end program spring_mass
