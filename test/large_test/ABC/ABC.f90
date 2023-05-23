program ABC

    ! Test case for the Arnold-Beltrami-Childress flow.

    use mpi
    use precision_mod , only : dp
    use global_mod    , only : pi, ierror
    use grid_mod 
    use solver_mod
    use navier_stokes_mod , only : viscosity, set_timestep, v, v, p
    use IO_mod

    implicit none

    ! Variables
    integer           :: Nx, Ny, Nz, step, e_id, r
    real(dp)          :: Lx, Ly, Lz, time, dt
    type(grid)        :: comp_grid 
    type(bc_type)     :: bc(6)
  
    ! Initialize MPI
    call mpi_init(ierror)

    ! The domain is a squared box of size 2 pi
    Lx = 2.0_dp*pi
    Ly = 2.0_dp*pi
    Lz = 2.0_dp*pi

    open(newunit = e_id, file = 'error')

    resolution_cycle: do r = 1,3
    
        ! Set the resolution
        Nx = 32*2**(r-1)
        Ny = 32*2**(r-1)
        Nz = 32*2**(r-1)

        ! Set boundary conditions
        bc(1)%s = 'Periodic'
        bc(2)%s = 'Periodic'
        bc(3)%s = 'Periodic'
        bc(4)%s = 'Periodic'
        bc(5)%s = 'Periodic'
        bc(6)%s = 'Periodic'
     
        ! Create the grid
        call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, [-pi, -pi, -pi], 8, 1, bc)
        
        ! Initialize the solver
        call init_solver(comp_grid)
        step = 0
        time = 0.0_dp
        
        ! Set the viscosity
        viscosity = 1.0e-1_dp

        ! Compute the timestep
        call set_timestep(comp_grid, dt, 2.0_dp)
        
        ! Initial condition
        call init_fields

        !==== Start Time loop ===================================================================
        time_loop: do while(time <= 0.1_dp)

            step = step + 1
            time = time + dt

            ! Advance in time the solution
            call advance_solution(comp_grid, step, dt)

            ! Advance solver status to log file
            call print_solver_status(stdout, step, time, dt)

        end do time_loop

        call compute_error()

        call destroy_solver
        call comp_grid%destroy()

    end do resolution_cycle
    
  ! Finalize the simulation
  call MPI_FINALIZE(ierror)

contains

    !=======================================================================================
    subroutine init_fields

        ! Initialize the ABC flow field.

        ! Local variables
        integer :: i, j, k

        ! Set the ABC flow
        do k = comp_grid%lo(3),comp_grid%hi(3)
        do j = comp_grid%lo(2),comp_grid%hi(2)
            do i = comp_grid%lo(1),comp_grid%hi(1)
                v%x%f(i,j,k) = cos(comp_grid%y(j)) + sin(comp_grid%z(k))
                v%y%f(i,j,k) = sin(comp_grid%x(i)) + cos(comp_grid%z(k))
                v%z%f(i,j,k) = cos(comp_grid%x(i)) + sin(comp_grid%y(j))
                p%f(i,j,k) = -(cos(comp_grid%x(i))*sin(comp_grid%y(j)) + &
                                sin(comp_grid%x(i))*cos(comp_grid%z(k)) + &
                                cos(comp_grid%y(j))*sin(comp_grid%z(k)))
            end do
        end do
        end do

        ! Update halos and bc
        call v%update_ghost_nodes()
        call p%update_ghost_nodes()

    end subroutine init_fields
    !=======================================================================================

    !=======================================================================================
    subroutine compute_error

        ! Evaluate error wrt the analytical solution

        ! Local variables
        integer  :: i, j, k
        real(dp) :: sol, e, L2(3), Linf(3)

        Linf = 0.0_dp
        L2 = 0.0_dp

        do k = comp_grid%lo(3),comp_grid%hi(3)
        do j = comp_grid%lo(2),comp_grid%hi(2)
            do i = comp_grid%lo(1),comp_grid%hi(1)
                sol = (cos(comp_grid%y(j)) + sin(comp_grid%z(k)))*exp(-viscosity*0.1_dp)
                e = abs(v%x%f(i,j,k) - sol)
                if (e > Linf(1)) Linf(1) = e
                L2(1) = L2(1) + e**2
                
                
                sol = (sin(comp_grid%x(i)) + cos(comp_grid%z(k)))*exp(-viscosity*0.1_dp)
                e = abs(v%y%f(i,j,k) - sol)
                if (e > Linf(2)) Linf(2) = e
                L2(2) = L2(2) + e**2
                
                sol = (cos(comp_grid%x(i)) + sin(comp_grid%y(j)))*exp(-viscosity*0.1_dp)
                e = abs(v%z%f(i,j,k) - sol)
                if (e > Linf(3)) Linf(3) = e
                L2(3) = L2(3) + e**2
            end do
        end do
        end do

        call mpi_allreduce(mpi_in_place,L2,3,mpi_real8,mpi_sum,mpi_comm_world,ierror)
        call mpi_allreduce(mpi_in_place,Linf,3,mpi_real8,mpi_max,mpi_comm_world,ierror)
        
        if (comp_grid%rank == 0) then
            write(e_id,*) Nx, Linf, sqrt(L2)
            flush(e_id)
        endif

    end subroutine compute_error
    !=======================================================================================
  
end program ABC
