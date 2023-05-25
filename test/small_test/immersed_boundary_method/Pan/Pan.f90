program Pan

    ! Particle migration in a channel flow at Re 12

    use mpi
    use precision_mod       , only : dp
    use global_mod          , only : ierror, myrank, pi
    use IO_mod              , only : stdout
    use grid_mod
    use eulerian_circle_mod , only : Circle
    use ibm_mod             , only : Eulerian_Solid_list
    use eulerian_ibm_mod
    use navier_stokes_mod   , only : set_timestep, g, viscosity
    use solver_mod

    implicit none

    ! Parameters
    real(dp), parameter :: radius = 0.125_dp        !< cylinder radius
    real(dp), parameter :: Deltap = 1.763e-3_dp     !< external body force
    real(dp), parameter :: mu = 3.2498036e-3_dp     !< fluid viscosity

    ! Variables
    integer              :: Nx, Ny, Nz, step
    real(dp)             :: Lx, Ly, Lz, time, dt, origin(3)
    type(grid), target   :: comp_grid
    type(bc_type)        :: bc(4)
    type(Circle), target :: C

    ! Initialize MPI
    call mpi_init(ierror)
    call mpi_comm_rank(mpi_comm_world, myrank, ierror)

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
    call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 1, 1, bc)

    ! Set body force
    g(1) = Deltap

    ! Create the solid body
    C = circle(X = [0.5_dp, 0.4_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], R = radius, name = 'Circle')
    C%G => comp_grid
    allocate(Eulerian_Solid_list(1))
    Eulerian_Solid_list(1)%pS => C
    call C%load_surface_points('mesh.txt')
    
    ! Set body force on the solid
    C%eF(1:3) = [g(1)*C%mass, g(2)*C%mass, 0.0_dp]

    ! Set the viscoisty
    viscosity = mu

    ! Initialize the solver
    call init_solver(comp_grid)
    step = 0
    time = 0.0_dp
    call set_timestep(comp_grid, dt, 1.0_dp)

    !==== Start Time loop ===================================================================
    time_loop: do while(time < 1.0e+3_dp)

        step = step + 1
        time = time + dt

        ! Advance in time the solution
        call advance_solution(comp_grid, step, dt)

        ! Advance solver status to log file
        call print_solver_status(stdout, step, time, dt)

        ! Print position
        if (myrank == 0) then
            call C%print_csv(time)
        endif

    end do time_loop

    ! free memory
    call destroy_solver
    call comp_grid%destroy()

    ! Finalize the simulation
    call MPI_FINALIZE(ierror)

end program Pan
