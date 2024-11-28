program Pan

    ! Particle migration in a channel flow at Re 12

    use mpi
    use precision_mod       , only : dp
    use global_mod          , only : ierror, myrank
    use IO_mod              , only : stdout
    use grid_mod
    use eulerian_circle_mod , only : Circle
    use ibm_mod             , only : solid_list
    use eulerian_ibm_mod
    use navier_stokes_mod   , only : set_timestep, g, viscosity
    use solver_mod

    implicit none

    ! Parameters
    real(dp), parameter :: radius = 0.125_dp    !< cylinder radius

    ! Variables
    integer              :: Nx, Ny, Nz, step, out_id
    real(dp)             :: Lx, Ly, Lz, time, dt, origin(3), deltap, mu
    type(grid), target   :: comp_grid
    type(bc_type)        :: bc(4)
    type(Circle), target :: C
    character(len=1)     :: case

    ! Initialize MPI
    call mpi_init(ierror)
    call mpi_comm_rank(mpi_comm_world, myrank, ierror)

    ! Select case setup
    call get_command_argument(1, case)
    select case(case)
    case('1')
        Deltap = 1.763e-3_dp
        mu = 3.2498036e-3_dp
    case('2')
        Deltap = 8.167e-4_dp
        mu = 1.5e-3_dp
    case('3')
        Deltap = 2.337e-4_dp
        mu = 4.2834760e-4_dp
    end select
    
    ! The domain is a unit square
    Lx = 1.0_dp
    Ly = 1.0_dp
    
    ! Set the resolution
    Nx = 96
    Ny = 96

    ! Since 2D
    Nz = 1
    Lz = Lx*float(Nz)/float(Nx)
    origin = [0.0_dp, 0.0_dp, 0.0_dp]

    ! Set bc
    bc(1)%s = 'Periodic'
    bc(2)%s = 'Periodic'
    bc(3)%s = 'Wall'
    bc(4)%s = 'Wall'

    ! Create the grid
    call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 4, 1, bc)

    ! Set body force
    g(1) = Deltap

    ! Create the solid body
    C%name = 'Circle'
    C%R = radius
    C%G => comp_grid
    C%density = 1.0_dp
    call C%setup()
    C%use_probes = .true.
    call C%load_surface_points('mesh.txt')
    
    ! Set initial position
    C%center_of_mass%X(1:2) = [0.5_dp, 0.4_dp]
    
    ! Set body force on the solid
    C%center_of_mass%Fe(1) = g(1)*C%mass

    ! Setup the solid list 
    allocate(solid_list(1))
    solid_list(1)%pS => C

    ! Set the viscoisty
    viscosity = mu

    ! Initialize the solver
    call init_solver(comp_grid)
    step = 0
    time = 0.0_dp
    call set_timestep(comp_grid, dt, 1.0_dp)
    dt = dt/2.0_dp

    open(newunit = out_id, file = 'case'//case//'.csv')
    if (myrank == 0) write(out_id, '(A29)') 't,x,y,omega,Fx,Fy,Fz,Mx,My,Mz'

    !==== Start Time loop ===================================================================
    time_loop: do while(time < 1.0e+3_dp)

        step = step + 1
        time = time + dt

        ! Advance in time the solution
        call advance_solution(comp_grid, step, dt)

        ! Advance solver status to log file
        call print_solver_status(stdout, step, time, dt)

        ! Print position
        if ((mod(step,10) == 0) .and. myrank == 0) then
            write(out_id,'(*(E16.8,:,","))') time, C%center_of_mass%X(1:2), C%center_of_mass%V(6), &
                                             C%center_of_mass%Fh
            flush(out_id)
        endif

    end do time_loop

    ! free memory
    call destroy_solver
    call comp_grid%destroy()

    ! Finalize the simulation
    call MPI_FINALIZE(ierror)

end program Pan
