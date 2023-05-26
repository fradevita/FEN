program main

    use mpi
    use precision_mod       , only : dp
    use global_mod          , only : ierror, myrank, pi
    use grid_mod
    use scalar_mod
    use eulerian_sphere_mod
    use navier_stokes_mod   , only : g, set_timestep, v, viscosity
    use ibm_mod             , only : Eulerian_Solid_list
    use solver_mod
    use IO_mod
    use fields_mod

    implicit none

    ! Parameters
    real(dp), parameter :: D = 1.0_dp   ! Sphere diameter
    integer , parameter :: N = 16       ! Grid point in the diameter
    real(dp), parameter :: L = 2.0_dp*D ! Domain size

    ! Variables
    integer              :: step
    real(dp)             :: time, dt
    type(grid), target   :: comp_grid
    type(scalar)         :: div
    type(sphere), target :: S
    type(bc_type)        :: bc(6)   

    call mpi_init(ierror)
    call mpi_comm_rank(mpi_comm_world, myrank, ierror)

    ! Create the grid
    bc(1)%s = 'Periodic'
    bc(2)%s = 'Periodic'
    bc(3)%s = 'Periodic'
    bc(4)%s = 'Periodic'
    bc(5)%s = 'Wall'
    bc(6)%s = 'Wall'
    call comp_grid%setup(2*N, 2*N, 2*N, L, L, L, [0.0_dp, 0.0_dp, 0.0_dp], 2, 2, bc)

    ! Create the sphere
    S = sphere(X = [1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], R = D/2.0_dp, name = 'S')
    S%G => comp_grid
    call S%setup()
    allocate(Eulerian_Solid_list(1))
    Eulerian_Solid_list(1)%pS => S

    viscosity = 0.1_dp

    ! Setup the solver
    call init_solver(comp_grid)
    step = 0
    time = 0.0_dp
    call set_timestep(comp_grid, dt, 1.0_dp)

    g(1) = 1.0_dp

    call save_fields(step)
 
    !v%x%bc%back = 1.0_dp
    call div%allocate(comp_grid)
    !==== Start Time loop ===================================================================
    time_loop: do while(step < 1000)

        step = step + 1
        time = time + dt

        ! Advance in time the solution
        call advance_solution(comp_grid, step, dt)

        ! Advance solver status to log file
        call print_solver_status(stdout, step, time, dt)


        call divergence(v, div)
        block
            integer :: i ,j ,k, imax, jmax, kmax
            real(dp) :: maxdiv
            maxdiv = 0.0_dp
            do k = 1,N
                do j = 1,N
                    do i = 1,N
                        if (abs(div%f(i,j,k)) > maxdiv) then
                            maxdiv = abs(div%f(i,j,k))
                            imax = i
                            jmax = j
                            kmax = k
                        endif
                    end do
                end do
            end do
            write(99,*) maxdiv, imax, jmax, kmax
        end block

        if (mod(step,1) == 0) call save_fields(step)

    end do time_loop

    ! free memory
    call destroy_solver
    call comp_grid%destroy()

    call mpi_finalize(ierror)

end program
