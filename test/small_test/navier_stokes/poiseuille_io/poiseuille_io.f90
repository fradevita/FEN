program test_Poiseuille_io

    use mpi
    use IO_mod
    use precision_mod     , only : dp
    use grid_mod     
    use scalar_mod
    use solver_mod        , only : init_solver, advance_solution, print_solver_status, destroy_solver
    use navier_stokes_mod , only : v, set_timestep

    implicit none

    integer          :: ierror, Nx, Ny, Nz, n, step, i
    real(dp)         :: Lx, Ly, Lz , dt, time, origin(3), max_r
    type(grid)       :: comp_grid
    type(bc_type)    :: bc(4)
    type(scalar)     :: uo
    character(len=3) :: sn 
    character(len=5) :: filename

    ! Set the tolerance
    real(dp) :: tol = 1.0d-8

    ! Initialize MPI
    call mpi_init(ierror)

    ! The physical size in y is 1 and in x is 4
    Lx = 1.0_dp
    Ly = 4.0_dp
    origin = [0.0_dp, 0.0_dp, 0.0_dp]
    bc(1)%s = 'Wall'
    bc(2)%s = 'Wall'
    bc(3)%s = 'Inflow'
    bc(4)%s = 'Outflow'

    ! Solve the problem for different levels or refinement
    refinement_loop: do n = 1,4

        Nx = 8*2**(n-1)
        Ny = 4*Nx
        Nz = 1
        Lx = Ly*float(Nx)/float(Ny)
        Lz = Ly*float(Nz)/float(Ny)

        ! Create the grid
        call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 4, 1, bc)

        ! Old u to check steady state
        call uo%allocate(comp_grid, 1)

        ! Initialize the solver
        call init_solver(comp_grid)
        step = 0
        time = 0.0_dp
        call set_timestep(comp_grid, dt, 1.0_dp)

        ! Set Inflow bc
        do i = comp_grid%lo(1),comp_grid%hi(1)
            v%y%bc%bottom(i,:) = -0.5_dp*(comp_grid%x(i)**2 - comp_grid%x(i))
        end do     

        !==== Start Time loop ===================================================================
        time_loop: do

            step = step + 1
            time = time + dt

            ! Advance in time the solution
            call advance_solution(comp_grid, step, dt)

            ! Print solver status
            call print_solver_status(stdout, step, time, dt)

            ! Check for steady state
            max_r = maxval(v%y%f - uo%f)
            call mpi_allreduce(mpi_in_place,max_r,1,mpi_real8,mpi_max,mpi_comm_world,ierror)
            if (max_r < tol .and. step > 2) then
                exit time_loop
            endif
            uo = v%y
        end do time_loop

        ! Write u
        write(sn,'(I0.3)') nx
        filename = 'v_'//sn
        call v%y%write(filename)

        ! Free the memory allocated by the solver
        call destroy_solver()
        call uo%destroy()
        call comp_grid%destroy

    end do refinement_loop

    ! Finalize the simulation
    call MPI_FINALIZE(ierror)

end program
