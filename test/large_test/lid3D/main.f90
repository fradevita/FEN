program main

    use mpi
    use precision_mod 
    use global_mod              , only : ierror, myrank
    use grid_mod
    use vector_mod
    use navier_stokes_mod       , only : set_timestep, viscosity, v
    use solver_mod
    use IO_mod                  , only : stdout

    implicit none

    ! Physical Parameters
    real(dp), parameter :: L = 1.0_dp     ! Cube edges
    real(dp), parameter :: U = 1.0_dp     ! Wall velocity
    real(dp), parameter :: Re = 1000.0_dp ! Reynolds number = U L / nu
    real(dp), parameter :: Lx = L
    real(dp), parameter :: Ly = L
    real(dp), parameter :: Lz = L

    ! Computational parameters
    integer , parameter :: Nx = 64
    integer , parameter :: Ny = 64
    integer , parameter :: Nz = 64
    real(dp), parameter :: tol = 1.0e-8_dp
    
    ! Variables
    integer       :: i, step
    real(dp)      :: dt, time, diff_x, diff_y, diff_z, maxdiff
    type(grid)    :: comp_grid
    type(bc_type) :: bc(6)
    type(vector)  :: v_old

    ! Setup MPI
    call mpi_init(ierror)
    call mpi_comm_rank(mpi_comm_world, myrank, ierror)

    !**** Define the computational grid ************************************************************
    bc(1)%s = 'Wall'
    bc(2)%s = 'Wall'
    bc(3)%s = 'Wall'
    bc(4)%s = 'Wall'
    bc(5)%s = 'Wall'
    bc(6)%s = 'Wall'
    call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, [0.0_dp, 0.0_dp, 0.0_dp], 2, 4, bc)
    call comp_grid%print_json()
    !***********************************************************************************************

    !**** Initilize the solver *********************************************************************
    viscosity = U*L/Re
    call init_solver(comp_grid)
    v%x%bc%top = U
    step = 0
    time = 0.0_dp
    call set_timestep(comp_grid, dt, U)
    dt = dt/2.0_dp
    !***********************************************************************************************

    !**** Start Time loop **************************************************************************
    v_old%x%f = v%x%f
    v_old%y%f = v%y%f
    v_old%z%f = v%z%f
    time_loop: do while(time < 2000.0_dp)

        step = step + 1
        time = time + dt

        ! Advance in time the solution
        call advance_solution(comp_grid, step, dt)

        ! Advance solver status to log file
        call print_solver_status(stdout, step, time, dt)

        ! Check for steady state
        diff_x = maxval(abs(v%x%f - v_old%x%f))
        diff_y = maxval(abs(v%y%f - v_old%y%f))
        diff_z = maxval(abs(v%z%f - v_old%z%f))
        call mpi_allreduce(mpi_in_place, diff_x, 1, mpi_real8, mpi_max, mpi_comm_world, ierror)
        call mpi_allreduce(mpi_in_place, diff_y, 1, mpi_real8, mpi_max, mpi_comm_world, ierror)
        call mpi_allreduce(mpi_in_place, diff_z, 1, mpi_real8, mpi_max, mpi_comm_world, ierror)
        maxdiff = max(diff_x, diff_y, diff_z)
        if (myrank == 0) print *, time, maxdiff
        v_old%x%f = v%x%f
        v_old%y%f = v%y%f
        v_old%z%f = v%z%f
        if (maxdiff < tol) then  
            call save_fields(0)
            call v%x%write('u.raw')
            call v%y%write('v.raw')
            exit time_loop
        endif

    end do time_loop
    !***********************************************************************************************
    
    ! free memory
    call destroy_solver
    call comp_grid%destroy()

    call mpi_finalize(ierror)

end program main
