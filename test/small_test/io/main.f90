program test_io_ns

    use mpi
    use precision_mod
    use global_mod          , only : ierror, myrank, pi
    use grid_mod
    use scalar_mod
    use navier_stokes_mod
    use solver_mod
    use IO_mod

    implicit none

    integer    :: Nx, Ny, Nz, step, i, j, k
    real(dp)   :: Lx, Ly, Lz, time, dt
    type(grid) :: comp_grid
    type(scalar) :: usave, vsave, wsave, psave
    
    ! init MPI
    call mpi_init(ierror)
    call mpi_comm_rank(mpi_comm_world, myrank, ierror)

    ! Create the numerical grid
    Nx = 4
    Ny = 4
    Lx = 2.0_dp*pi
    Ly = Lx
#if DIM==3
    Nz = 4
    Lz = Lx
#else
    Nz = 1
    Lz = Lx*float(Nz)/float(Nx)
#endif
    call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, [0.0_dp, 0.0_dp, 0.0_dp], 1, 1)

    ! Initialize the solver
    call init_solver(comp_grid)
    step = 0
    time = 0.0_dp

    ! Set initial fields
    call init_fields

    call set_timestep(comp_grid, dt, 2.0_dp)

    ! Perform 2 timestep
    do i = 1,2
        call advance_solution(comp_grid, step, dt)

        ! Print solver status
        step = step + 1
        time = time + dt
        !call print_solver_status(stdout, step, time, dt)

        ! Save simulation state
        call save_state(step)

        if (i == 1) then
            usave%f = v%x%f
            vsave%f = v%y%f
            wsave%f = v%z%f
            psave%f = p%f
        endif
    end do

    ! Now load state 1 and perform again one timestep
    call load_state(1)
    call advance_solution(comp_grid, step, dt)
    !call print_solver_status(stdout, step, time, dt)
    call save_state(3)

    ! Clean memory
    call destroy_solver()
    call comp_grid%destroy()
    call mpi_finalize(ierror)

contains

    !===============================================================================================
    subroutine init_fields

        integer :: i, j, k
        real(dp) :: x, y, z

        do k = comp_grid%lo(3),comp_grid%hi(3)
            do j = comp_grid%lo(2),comp_grid%hi(2)
                do i = comp_grid%lo(1),comp_grid%hi(1)
                    x = float(i)*comp_grid%delta
                    y = (float(j) - 0.5_dp)*comp_grid%delta
                    v%x%f(i,j,k) = -cos(x)*sin(y)
                    x = (float(i) - 0.5_dp)*comp_grid%delta
                    y = float(j)*comp_grid%delta
                    v%y%f(i,j,k) = sin(x)*cos(y)
#if DIM==3
                    x = (float(i) - 0.5_dp)*comp_grid%delta
                    z = float(k)*comp_grid%delta
                    v%z%f(i,j,k) = sin(x)*cos(z)
#endif
                    x = (float(i) - 0.5_dp)*comp_grid%delta
                    y = (float(j) - 0.5_dp)*comp_grid%delta
                    p%f(i,j,k) = -0.25_dp*(cos(2.0_dp*x) + cos(2.0_dp*y))
                end do
            end do
        end do

        ! Update halos and bc
        call p%update_ghost_nodes
        call v%update_ghost_nodes

    end subroutine init_fields
    !===============================================================================================

end program test_io_ns