program main

    ! Flow around a fixed cylinder at Re = 20

    use mpi
    use precision       , only : dp
    use constants       , only : pi
    use class_Grid      , only : base_grid, bc_type
    use ibm             , only : Lagrangian_Solid_list
    use lagrangian_solid, only : solid
    use solver          , only : init_solver, destroy_solver, advance_solution, print_solver_status
    use navier_stokes   , only : v, p, set_timestep, g, viscosity, mu, rho
    use lagrangian_ibm  , only : compute_hydrodynamic_loads
    use io              , only : stdout
    
    implicit none

    ! Parameters
    real(dp), parameter :: radius = 0.05_dp !< Cylinder radius
    real(dp), parameter :: U = 0.3_dp       !< Maximum inflow velocity

    ! Variables
    integer             :: ierror, Nx, Ny, Nz, step, i
    real(dp)            :: Lx, Ly, Lz, time, dt, origin(3)
    character(len=3)    :: arg
    type(bc_type)       :: bc(4)
    type(solid), target :: C 

    ! Initialize MPI
    call mpi_init(ierror)

    ! The program takes as input argument the number of point in the x direction
    call get_command_argument(1, arg)

    ! The domain
    Lx = 0.41_dp
    Ly = Lx*5.0_dp

    ! Set the resolution
    read(arg,'(I03)') Nx
    Ny = Nx*5

    ! Since 2D
    Nz = 1
    Lz = Lx*float(Nz)/float(Nx)

    origin = [0.0_dp, 0.0_dp, 0.0_dp]
    bc(1)%s = 'Wall'
    bc(2)%s = 'Wall'
    bc(3)%s = 'Inflow'
    bc(4)%s = 'Outflow'
    ! Create the grid
    call base_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 8, 1, bc)

    ! Create the cylinder
    call C%create('mesh.txt', name = trim('C_'//arg))

    ! Allocate the array of solid
    allocate(lagrangian_solid_list(1))
    Lagrangian_Solid_list(1)%pS => C
    
    ! Set the viscoisty
    viscosity = 0.001_dp
    
    ! Initialize the solver
    call init_solver
    step = 0
    time = 0.0_dp
    call set_timestep(dt, U)
    
    ! Set the inflow boundary condition
    do i = base_grid%lo(1),base_grid%hi(1)
        v%y%bc%b(i,:) = 4.0_dp*U*base_grid%x(i)*(Lx - base_grid%x(i))/Lx**2
    end do

    ! Output
    if (base_grid%rank == 0) call C%write_csv(time)

    !==== Start Time loop ===================================================================
    time_loop: do while (time < 3.0_dp)

        step = step + 1
        time = time + dt

        ! Advance in time the solution
        call advance_solution(step, dt)

        ! Advance solver status to log file
        call print_solver_status(stdout, step, time, dt)

        ! Compute the forces on the cylinder
        call compute_hydrodynamic_loads(lagrangian_solid_list(1)%pS, v, p, mu, rho, g)
        call lagrangian_solid_list(1)%pS%integrate_hydrodynamic_forces()

        ! Output forces
        if (base_grid%rank == 0) call C%write_csv(time)
    end do time_loop

    ! free memory
    call destroy_solver

    ! Finalize the simulation
    call MPI_FINALIZE(ierror)

end program main
