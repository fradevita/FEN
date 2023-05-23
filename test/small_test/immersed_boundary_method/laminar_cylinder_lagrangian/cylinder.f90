program main

    ! Flow around a fixed cylinder at Re = 20

    use mpi
    use precision_mod        , only : dp
    use global_mod           , only : ierror, pi
    use grid_mod
    use ibm_mod              , only : Lagrangian_Solid_list
    use lagrangian_solid_mod , only : solid
    use solver_mod        
    use navier_stokes_mod    , only : v, p, set_timestep, g, viscosity, mu, rho
    use lagrangian_ibm_mod   , only : compute_hydrodynamic_loads
    use IO_mod               , only : stdout
    
    implicit none

    ! Parameters
    real(dp), parameter :: radius = 0.05_dp !< Cylinder radius
    real(dp), parameter :: U = 0.3_dp       !< Maximum inflow velocity

    ! Variables
    integer             :: Nx, Ny, Nz, step, i
    real(dp)            :: Lx, Ly, Lz, time, dt, origin(3)
    character(len=3)    :: arg
    type(grid)          :: comp_grid
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
    call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 4, 1, bc)

    ! Create the cylinder
    call C%create('mesh.txt', name = trim('C_'//arg))

    ! Allocate the array of solid
    allocate(lagrangian_solid_list(1))
    Lagrangian_Solid_list(1)%pS => C
    
    ! Set the viscoisty
    viscosity = 0.001_dp
    
    ! Initialize the solver
    call init_solver(comp_grid)
    step = 0
    time = 0.0_dp
    call set_timestep(comp_grid, dt, U)
    
    ! Set the inflow boundary condition
    do i = comp_grid%lo(1),comp_grid%hi(1)
        v%y%bc%bottom(i,:) = 4.0_dp*U*comp_grid%x(i)*(Lx - comp_grid%x(i))/Lx**2
    end do

    ! Output
    if (comp_grid%rank == 0) call C%write_csv(time)

    !==== Start Time loop ===================================================================
    time_loop: do while (time < 3.0_dp)

        step = step + 1
        time = time + dt

        ! Advance in time the solution
        call advance_solution(comp_grid, step, dt)

        ! Advance solver status to log file
        call print_solver_status(stdout, step, time, dt)

        ! Compute the forces on the cylinder
        call compute_hydrodynamic_loads(lagrangian_solid_list(1)%pS, v, p, mu, rho, g)
        call lagrangian_solid_list(1)%pS%integrate_hydrodynamic_forces()

        ! Output forces
        if (comp_grid%rank == 0) call C%write_csv(time)
    end do time_loop

    ! free memory
    call destroy_solver

    ! Finalize the simulation
    call MPI_FINALIZE(ierror)

end program main
