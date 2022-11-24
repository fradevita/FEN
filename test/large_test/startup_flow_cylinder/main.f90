program main

    ! Startup flow around a fixed cylinder at Re = 1000

    use mpi
    use precision            , only : dp
    use constants            , only : pi
    use class_Grid           , only : base_grid, bc_type
    use class_Vector         , only : vector
    use class_eulerian_circle, only : circle
    use navier_stokes        , only : v, set_timestep, viscosity
    use fields               , only : curl
    use ibm
    use solver
    use io                   , only : stdout
    use json
    
    implicit none

    ! Parameters
    real(dp), parameter :: D = 1.0_dp        !< Reference lenght is the cylinder diameter
    real(dp), parameter :: radius = D/2.0_dp !< cylinder radius
    real(dp), parameter :: U = 1.0_dp        !< reference velocity
    real(dp), parameter :: Re = 1.0e+3_dp    !< Reynolds number

    ! Variables
    integer              :: ierror, Nx, Ny, Nz, step, out_id
    real(dp)             :: Lx, Ly, Lz, time, dt, origin(3), Cd, Cl
    type(vector)         :: omega
    type(bc_type)        :: bc(4)
    type(circle), target :: C
    character(len=7 )    :: ss
    character(len=20)    :: filename

    ! Initialize MPI
    call mpi_init(ierror)

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
    call base_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 8, 1, bc)

    ! Create the circle
    C = circle(X = [6.0_dp, 6.0_dp, 6.0_dp], R = radius)
    allocate(Eulerian_Solid_list(1))
    Eulerian_Solid_list(1)%pS => C

    ! Set the viscoisty
    viscosity = 1.0_dp/Re
    
    ! Initialize the solver
    call init_solver
    step = 0
    time = 0.0_dp
    call set_timestep(dt, 3*U)

    ! Set the inflow boundary condition
    v%y%bc%l = U
    v%y%bc%r = U
    v%y%bc%b = U

    ! Open output file
    open(newunit = out_id, file = 'out.txt')

    call print_setup_json(dt)
    call save_fields(0)

    ! Allocate memory for vorticity field (for output only)
    call omega%allocate()

    !==== Start Time loop ===================================================================
    time_loop: do while (time < 3.0_dp)

        step = step + 1
        time = time + dt

        ! Advance in time the solution
        call advance_solution(step, dt)

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
            Cd = 2.0_dp*Fe%y%integral()/U**2/D/Lz
            Cl = 2.0_dp*Fe%x%integral()/U**2/D/Lz
            if (base_grid%rank == 0) write(out_id, *) time, Cd, Cl
        endif

    end do time_loop

    ! free memory
    call destroy_solver

    ! Finalize the simulation
    call MPI_FINALIZE(ierror)

end program main
