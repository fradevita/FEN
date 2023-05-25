program main

    ! Flow around a fixed cylinder at Re = 20

    use mpi
    use precision_mod        , only : dp
    use global_mod           , only : ierror, pi
    use grid_mod
    use ibm_mod              , only : Lagrangian_Solid_list
    use lagrangian_solid_mod , only : solid
    use solver_mod        
    use navier_stokes_mod    , only : set_timestep, g, viscosity
    use lagrangian_ibm_mod   , only : compute_hydrodynamic_loads
    use IO_mod               , only : stdout
    
    implicit none

    ! Parameters
    real(dp), parameter :: radius = 0.125_dp        !< cylinder radius
    real(dp), parameter :: Deltap = 1.763e-3_dp     !< external body force
    real(dp), parameter :: mu = 3.2498036e-3_dp     !< fluid viscosityy

    ! Variables
    integer             :: Nx, Ny, Nz, step
    real(dp)            :: Lx, Ly, Lz, time, dt, origin(3)
    type(grid)          :: comp_grid
    type(bc_type)       :: bc(4)
    type(solid), target :: C 

    ! Initialize MPI
    call mpi_init(ierror)

    ! The domain
    Lx = 1.0_dp
    Ly = 1.0_dp

    ! Set the resolution
    Nx = 96
    Ny = 96

    ! Since 2D
    Nz = 1
    Lz = Lx*float(Nz)/float(Nx)

    origin = [0.0_dp, 0.0_dp, 0.0_dp]
    bc(1)%s = 'Periodic'
    bc(2)%s = 'Periodic'
    bc(3)%s = 'Wall'
    bc(4)%s = 'Wall'
    ! Create the grid
    call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 1, 1, bc)

    ! Create the cylinder
    ! Set the mass of the solid body
    C%M(1:2) = pi*0.125_dp**2
    C%M(3) = 0.5_dp*pi*0.125_dp**4
    ! Create the lagrangian solid
    call C%create('mesh.txt', name = 'C')
    ! Allocate the array of solid
    allocate(lagrangian_solid_list(1))
    Lagrangian_Solid_list(1)%pS => C

    ! Set body force
    g(1) = Deltap

    ! Set the viscoisty
    viscosity = mu
    
    ! Initialize the solver
    call init_solver(comp_grid)
    step = 0
    time = 0.0_dp
    call set_timestep(comp_grid, dt, 1.0_dp)
    
    ! Output
    if (comp_grid%rank == 0) call C%write_csv(time)
    call save_fields(0)
    
    !==== Start Time loop ===================================================================
    time_loop: do while (time < 400.0_dp)

        step = step + 1
        time = time + dt

        ! Advance in time the solution
        call advance_solution(comp_grid, step, dt)

        ! Advance solver status to log file
        call print_solver_status(stdout, step, time, dt)

        ! Output forces
        if (comp_grid%rank == 0) call C%write_csv(time)

        if (mod(step,100) == 0) then
            call save_fields(step)
            if (comp_grid%rank == 0) call C%print_configuration(step)
        endif
    end do time_loop

    ! free memory
    call destroy_solver

    ! Finalize the simulation
    call MPI_FINALIZE(ierror)

end program main
