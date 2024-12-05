program test_io_fsi_e

    use mpi
    use precision_mod
    use global_mod          , only : ierror, myrank, pi
    use grid_mod
    use scalar_mod
    use eulerian_circle_mod , only : Circle
    use ibm_mod             , only : solid_list
    use eulerian_ibm_mod
    use navier_stokes_mod   , only : set_timestep, g, viscosity
    use solver_mod
    use IO_mod

    implicit none

    ! Parameters
    real(dp), parameter :: radius = 0.125_dp    !< cylinder radius
    real(dp), parameter :: deltap = 2.337e-4_dp !< pressure gradient (body force)
    real(dp), parameter :: mu = 4.2834760e-4_dp !< fluid viscosity

    ! Variables
    integer              :: Nx, Ny, Nz, step, i, j, k
    real(dp)             :: Lx, Ly, Lz, time, dt, origin(3)
    type(grid), target   :: comp_grid
    type(bc_type)        :: bc(4)
    type(Circle), target :: C
    character(len=1)     :: case

    ! Init MPI
    call mpi_init(ierror)
    call mpi_comm_rank(mpi_comm_world, myrank, ierror)

    ! Create the numerical grid
    Nx = 96
    Ny = 96
    Nz = 1
    Lx = 1.0_dp
    Ly = Lx
    Lz = Lx*float(Nz)/float(Nx)

    ! Set bc
    bc(1)%s = 'Periodic'
    bc(2)%s = 'Periodic'
    bc(3)%s = 'Wall'
    bc(4)%s = 'Wall'
    call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, [0.0_dp, 0.0_dp, 0.0_dp], 1, 1)
    
    ! Set body force
    g(1) = Deltap

    ! Create the solid body
    C%name = 'Circle'
    C%R = radius
    C%G => comp_grid
    C%density = 1.0_dp
    call C%setup()
    C%use_probes = .true.

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

    ! Select run case
    call get_command_argument(1, case)
    select case(case)
    case('1')
        ! Set initial position
        C%center_of_mass%X(1:2) = [0.5_dp, 0.4_dp]

        ! Load surface points
        call C%load_surface_points('mesh.txt')
    
        ! Perform 10 timestep
        do i = 1,10
            step = step + 1
            time = time + dt
            call advance_solution(comp_grid, step, dt)

            ! Save simulation state
            if (i == 5) then
                call save_state(step)
            endif 
        end do
        call save_state(step)
    case('2')
        ! Now load state 5 and perform again 4 timestep
        step = 5
        time = step*dt
        call load_state(step)
        
        do i = 6, 10
            step = step + 1
            time = time + dt
            call advance_solution(comp_grid, step, dt)
        end do
        call save_state(step*2)
    end select

    ! Clean memory
    call destroy_solver()
    call comp_grid%destroy()
    call mpi_finalize(ierror)

end program test_io_fsi_e
