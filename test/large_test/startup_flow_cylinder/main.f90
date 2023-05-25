program main

    ! Startup flow around a fixed cylinder at Re = 1000

    use mpi
    use precision_mod        , only : dp
    use global_mod           , only : ierror, myrank, pi
    use grid_mod
    use vector_mod           , only : vector
    use eulerian_circle_mod  , only : circle
    use navier_stokes_mod    , only : v, set_timestep, viscosity, p, mu, rho, g
    use fields_mod           , only : curl
    use ibm_mod
    use eulerian_ibm_mod     , only : compute_hydrodynamic_loads
    use solver_mod
    use IO_mod               , only : stdout
    
    implicit none

    ! Parameters
    real(dp), parameter :: D = 1.0_dp        !< Reference lenght is the cylinder diameter
    real(dp), parameter :: radius = D/2.0_dp !< cylinder radius
    real(dp), parameter :: U = 1.0_dp        !< reference velocity
    real(dp), parameter :: Re = 1.0e+3_dp    !< Reynolds number

    ! Variables
    integer              :: Nx, Ny, Nz, step, out_id
    real(dp)             :: Lx, Ly, Lz, time, dt, origin(3), Cd
    type(grid), target   :: comp_grid          
    type(vector)         :: omega
    type(bc_type)        :: bc(4)
    type(circle), target :: C
    character(len=7 )    :: ss
    character(len=20)    :: filename

    ! Initialize MPI
    call mpi_init(ierror)
    call mpi_comm_rank(mpi_comm_world, myrank, ierror)

    ! The domain has size [18D, 12D]
    Lx = 12.0_dp*D
    Ly = 18.0_dp*D

    ! Set the resolution to 512 points per diameter
    Nx = 512*12/2
    Ny = 512*18/2

    ! Since 2D
    Nz = 1
    Lz = Lx*float(Nz)/float(Nx)

    ! Create the grid
    bc(1)%s = 'Wall'
    bc(2)%s = 'Wall'
    bc(3)%s = 'Inflow'
    bc(4)%s = 'Outflow'
    origin = [0.0_dp, 0.0_dp, 0.0_dp]
    call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 8, 1, bc)

    ! Create the circle
    C = circle(X = [6.0_dp, 6.0_dp, 6.0_dp], R = radius, name = 'C')
    C%G => comp_grid
    allocate(Eulerian_Solid_list(1))
    Eulerian_Solid_list(1)%pS => C
    ! Need to create the surface mesh to use probes for hydrodynamic forces component 
    ! evaluation
    call C%load_surface_points('mesh.txt')

    ! Set the viscoisty
    viscosity = 1.0_dp/Re
    
    ! Initialize the solver
    call init_solver(comp_grid)
    step = 0
    time = 0.0_dp
    call set_timestep(comp_grid, dt, 3*U)

    ! Set the inflow boundary condition
    v%y%bc%left = U
    v%y%bc%right = U
    v%y%bc%bottom = U

    ! Open output file
    open(newunit = out_id, file = 'out.txt')

    call save_fields(0)

    ! Allocate memory for vorticity field (for output only)
    call omega%allocate(comp_grid)

    !==== Start Time loop ===================================================================
    time_loop: do while (time < 3.0_dp)

        step = step + 1
        time = time + dt

        ! Advance in time the solution
        call advance_solution(comp_grid, step, dt)

        ! Advance solver status to log file
        call print_solver_status(stdout, step, time, dt)

        ! Output fields
        if (mod(step,131) == 0) then
            call save_fields(step)
            ! Compute vorticity
            call curl(v, omega)
            write(ss,'(I0.7)') step
            filename = 'data/vr_'//ss//'.raw'
            call omega%x%write(filename)     
        endif

        ! Evaluate hydrodynamic forces
        if (mod(step,100) == 0) then
            ! Drag coefficient can be computed integrating the eulerain Forcing Fe
            Cd = 2.0_dp*Fe%y%integral()/U**2/D/Lz
            ! Must use the probes to evaluate viscous contribution and 
            call compute_hydrodynamic_loads(C, v, p, mu, rho, g)
            ! Write output
            if (myrank == 0) write(out_id, *) time, Cd, 2.0_dp*C%hFv(2)/U**2/D, 2.0_dp*C%hFp(2)/U**2/D
        endif

    end do time_loop

    ! free memory
    call destroy_solver

    ! Finalize the simulation
    call MPI_FINALIZE(ierror)

end program main
