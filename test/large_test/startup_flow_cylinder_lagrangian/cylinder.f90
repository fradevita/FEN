program main

    ! Flow around a fixed cylinder at Re = 20

    use mpi
    use precision_mod        , only : dp
    use global_mod           , only : ierror, myrank, pi
    use grid_mod
    use ibm_mod              , only : Lagrangian_Solid_list
    use lagrangian_solid_mod , only : solid
    use solver_mod        
    use navier_stokes_mod    , only : v, p, set_timestep, g, viscosity, mu, rho, dt_o
    use lagrangian_ibm_mod   , only : compute_hydrodynamic_loads
    use IO_mod               , only : stdout
    
    implicit none

    ! Parameters
    real(dp), parameter :: D = 1.0_dp        !< Reference lenght is the cylinder diameter
    real(dp), parameter :: radius = D/2.0_dp !< cylinder radius
    real(dp), parameter :: U = 1.0_dp        !< reference velocity
    real(dp), parameter :: Re = 1.0e+3_dp    !< Reynolds number

    ! Variables
    integer             :: Nx, Ny, Nz, step
    real(dp)            :: Lx, Ly, Lz, time, dt, origin(3)
    type(grid)          :: comp_grid
    type(bc_type)       :: bc(4)
    type(solid), target :: C 

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

    origin = [0.0_dp, 0.0_dp, 0.0_dp]
    bc(1)%s = 'Wall'
    bc(2)%s = 'Wall'
    bc(3)%s = 'Inflow'
    bc(4)%s = 'Outflow'
    ! Create the grid
    comp_grid%name = 'grid'
    call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 8, 1, bc)

    ! Create the cylinder
    call C%create('mesh.txt', name = 'C')

    ! Allocate the array of solid
    allocate(lagrangian_solid_list(1))
    Lagrangian_Solid_list(1)%pS => C
    
    ! Set the viscoisty
    viscosity = 1.0_dp/Re
    
    ! Initialize the solver
    call init_solver(comp_grid)
    step = 0
    time = 0.0_dp
    call set_timestep(comp_grid, dt, 2*U)
    dt = dt/2.0_dp
    dt_o = dt
 
    ! Set the inflow boundary condition
    v%y%f = U
    v%y%bc%left = U
    v%y%bc%right = U
    v%y%bc%bottom = U

    ! Output
    if (comp_grid%rank == 0) call C%write_csv(time)
    call save_fields(step)
    
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
        if (comp_grid%rank == 0) then
            call C%write_csv(time)
            call C%print_configuration(step)
        endif
        if (mod(step,10) == 0)  call save_fields(step)
    end do time_loop

    ! free memory
    call destroy_solver

    ! Finalize the simulation
    call MPI_FINALIZE(ierror)

end program main
