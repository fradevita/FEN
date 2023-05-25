program main

    use mpi
    use precision_mod       , only : dp
    use global_mod          , only : ierror, myrank, pi
    use grid_mod
    use eulerian_sphere_mod
    use navier_stokes_mod   , only : g, set_timestep, dt_o
    use ibm_mod             , only : Eulerian_Solid_list
    use solver_mod
    use IO_mod

    implicit none

    ! Parameters
    real(dp), parameter :: D = 1.0_dp   ! Sphere diameter
    integer , parameter :: N = 16       ! Grid point in the diameter
    real(dp), parameter :: L = 2.0_dp*D ! Domain size

    ! Variables
    integer              :: step
    real(dp)             :: time, dt
    type(grid), target   :: comp_grid
    type(sphere), target :: S

    call mpi_init(ierror)
    call mpi_comm_rank(mpi_comm_world, myrank, ierror)

    ! Create the grid
    call comp_grid%setup(2*N, 2*N, 2*N, L, L, L, [-L/2.0_dp,-L/2.0_dp,-L/2.0_dp], 2, 4)

    ! Create the sphere
    S = sphere(R = D/2.0_dp, name = 'S')
    S%G => comp_grid
    call S%setup()
    allocate(Eulerian_Solid_list(1))
    Eulerian_Solid_list(1)%pS => S

    ! Set body force
    g(1) = 0.1_dp

    ! Setup the solver
    call init_solver(comp_grid)
    step = 0
    time = 0.0_dp
    call set_timestep(comp_grid, dt, 1.0_dp)

    call save_fields(step)

    !==== Start Time loop ===================================================================
    time_loop: do while(step < 100)

        step = step + 1
        time = time + dt

        ! Advance in time the solution
        call advance_solution(comp_grid, step, dt)

        ! Advance solver status to log file
        call print_solver_status(stdout, step, time, dt)

        if (mod(step,10) == 0) call save_fields(step)

    end do time_loop

    ! free memory
    call destroy_solver
    call comp_grid%destroy()

    call mpi_finalize(ierror)

end program