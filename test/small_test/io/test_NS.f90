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

    integer          :: Nx, Ny, Nz, step, i, j, k
    real(dp)         :: Lx, Ly, Lz, time, dt
    type(grid)       :: comp_grid
    character(len=1) :: case
    
    ! init MPI
    call mpi_init(ierror)
    call mpi_comm_rank(mpi_comm_world, myrank, ierror)

    ! Create the numerical grid
    Nx = 16
    Ny = 16
    Nz = 16
    Lx = 2.0_dp*pi
    Ly = Lx
    Lz = Lx
    call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, [0.0_dp, 0.0_dp, 0.0_dp], 1, 1)

    ! Initialize the solver
    call init_solver(comp_grid)
    step = 0
 
    ! Set initial fields
    call init_fields
    call set_timestep(comp_grid, dt, 2.0_dp)

    ! Read command line argument
    call get_command_argument(1, case)

    select case(case)
    case('1')        
        ! Perform 10 timestep
        do i = 1,10
            step = step + 1
            call advance_solution(comp_grid, step, dt)
            
            ! Save simulation state after 5 steps
            if (i == 5) call save_state(step)
        end do
        call save_state(step)
    case('2')
        ! Load state 5 and perform again 5 timestep
        step = 5
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