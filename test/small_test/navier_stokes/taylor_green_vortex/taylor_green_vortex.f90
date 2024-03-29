program test_TGV

    use mpi
    use IO_mod
    use precision_mod     , only : dp
    use global_mod        , only : pi
    use grid_mod
    use solver_mod        , only : init_solver, advance_solution, destroy_solver, print_solver_status
    use navier_stokes_mod , only : v, set_timestep, p

    implicit none

    integer          :: ierror, nx, ny, nz, step, n
    real(dp)         :: Lx, Ly, Lz, dt, time, origin(3)
    type(grid)       :: comp_grid
    type(bc_type)    :: bc(4)
    character(len=3) :: sn 
    character(len=5) :: filename

    ! The simulation is performed up to t = 0.3
    real(dp), parameter :: Tmax = 0.3_dp

    ! Initialize MPI
    call mpi_init(ierror)

    ! The domain is a periodic squared box of size 2 pi
    Lx = 2.0_dp*pi
    Ly = Lx
    origin = [0.0_dp, 0.0_dp, 0.0_dp]
    bc(1)%s = 'Periodic'
    bc(2)%s = 'Periodic'
    bc(3)%s = 'Periodic'
    bc(4)%s = 'Periodic'

    ! Solve the problem for several different levels resolution
    refinement_loop: do n = 1,5

        ! Setup the grid
        Nx = 16*2**(n-1)
        Ny = Nx
        Nz = 1
        Lz = Lx*float(Nz)/float(Nx)
        call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 4, 1, bc)
    
        ! Initialize the solver
        call init_solver(comp_grid)
        step = 0
        time = 0.0_dp

        ! Set initial fields
        call init_fields

        call set_timestep(comp_grid, dt, 2.0_dp)

        !==== Start Time loop ===================================================================
        time_loop: do while(time < Tmax)

            step = step + 1
            time = time + dt

            ! Advance in time the solution
            call advance_solution(comp_grid, step, dt)

            ! Print solver status
            call print_solver_status(stdout, step, time, dt)

        end do time_loop

        ! Write u
        write(sn,'(I0.3)') ny
        filename = 'u_'//sn
        call v%x%write(filename)

        ! Free the memory allocated by the solver
        call destroy_solver()
        call comp_grid%destroy
  
    end do refinement_loop

    ! Finalize the simulation
    call MPI_FINALIZE(ierror)

contains

    !=======================================================================================
    subroutine init_fields

        integer :: i, j, k
        real(dp) :: x, y

        do k = comp_grid%lo(3),comp_grid%hi(3)
        do j = comp_grid%lo(2),comp_grid%hi(2)
            do i = comp_grid%lo(1),comp_grid%hi(1)
            x = float(i)*comp_grid%delta
            y = (float(j) - 0.5_dp)*comp_grid%delta
            v%x%f(i,j,k) = -cos(x)*sin(y)
            x = (float(i) - 0.5_dp)*comp_grid%delta
            y = float(j)*comp_grid%delta
            v%y%f(i,j,k) = sin(x)*cos(y)
            x = (float(i) - 0.5_dp)*comp_grid%delta
            y = (float(j) - 0.5_dp)*comp_grid%delta
            p%f(i,j,k) = -0.25_dp*(cos(2.0_dp*x) + cos(2.0_dp*y))
            end do
        end do
        end do

        ! Update halos and bc
        call p%update_ghost_nodes
        call v%update_ghost_nodes

        !call compute_explicit_terms(dv_o)

    end subroutine init_fields
    !=======================================================================================

end program test_TGV
