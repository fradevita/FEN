program ABC

  ! Test case for the Arnold-Beltrami-Childress flow.

  use mpi
  use precision    , only : dp
  use constants    , only : pi
  use class_Grid   , only : base_grid, bc_type
  use solver
  use navier_stokes, only : viscosity, set_timestep, v, constant_CFL, CFL
  use io
  use json

  implicit none

  ! Variables
  integer  :: ierror, Nx, Ny, Nz, step, e_id, r, error 
  real(dp) :: Lx, Ly, Lz, time, dt
  type(bc_type) :: bc(6)
  character(len=7) :: sn
  character(len=14) :: filename
  
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
     call base_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, [-pi, -pi, -pi], 8, 1, bc)
     
     ! Initialize the solver
     call init_solver
     step = 0
     time = 0.0_dp
     
     ! Set the viscosity
     viscosity = 1.0_dp

     ! Compute the timestep
     call set_timestep(dt, 2.0_dp)
     
     ! Initial condition
     call init_fields
     
     call print_setup_json(dt)
     call save_fields(step)

     !==== Start Time loop ===================================================================
     time_loop: do while(time <= 0.1_dp)

        step = step + 1
        time = time + dt

        ! Advance in time the solution
        call advance_solution(step, dt)

        ! Advance solver status to log file
        call print_solver_status(stdout, step, time, dt)

     end do time_loop

     call save_fields(step)
     call compute_error()

     call destroy_solver

  end do resolution_cycle
  
  ! Finalize the simulation
  call MPI_FINALIZE(ierror)

contains

  !=======================================================================================
  subroutine init_fields

    ! Initialize the ABC flow field.
    
    use navier_stokes, only : v, p

    ! Local variables
    integer :: i, j, k

    ! Set the ABC flow
    do k = base_grid%lo(3),base_grid%hi(3)
       do j = base_grid%lo(2),base_grid%hi(2)
          do i = base_grid%lo(1),base_grid%hi(1)
             v%x%f(i,j,k) = cos(base_grid%y(j)) + sin(base_grid%z(k))
             v%y%f(i,j,k) = sin(base_grid%x(i)) + cos(base_grid%z(k))
             v%z%f(i,j,k) = cos(base_grid%x(i)) + sin(base_grid%y(j))
               p%f(i,j,k) = -(cos(base_grid%x(i))*sin(base_grid%y(j)) + &
                              sin(base_grid%x(i))*cos(base_grid%z(k)) + &
                              cos(base_grid%y(j))*sin(base_grid%z(k)))
          end do
       end do
    end do

    ! Update halos and bc
    call v%apply_bc()
    call p%apply_bc()

  end subroutine init_fields
  !=======================================================================================

  !=======================================================================================
  subroutine compute_error

    ! Evaluate error wrt the analytical solution

    use navier_stokes, only : v, p

    ! Local variables
    integer  :: i, j, k, imax, jmax, kmax
    real(dp) :: sol, e, eu_max, ev_max, ew_max, L2(3), Linf(3)

    Linf = 0.0_dp
    L2 = 0.0_dp

    do k = base_grid%lo(3),base_grid%hi(3)
       do j = base_grid%lo(2),base_grid%hi(2)
          do i = base_grid%lo(1),base_grid%hi(1)
             sol = (cos(base_grid%y(j)) + sin(base_grid%z(k)))*exp(-viscosity*0.1_dp)
             e = abs(v%x%f(i,j,k) - sol)
             if (e > Linf(1)) Linf(1) = e
             L2(1) = L2(1) + e**2
             
             
             sol = (sin(base_grid%x(i)) + cos(base_grid%z(k)))*exp(-viscosity*0.1_dp)
             e = abs(v%y%f(i,j,k) - sol)
             if (e > Linf(2)) Linf(2) = e
             L2(2) = L2(2) + e**2
             
             sol = (cos(base_grid%x(i)) + sin(base_grid%y(j)))*exp(-viscosity*0.1_dp)
             e = abs(v%z%f(i,j,k) - sol)
             if (e > Linf(3)) Linf(3) = e
             L2(3) = L2(3) + e**2
          end do
       end do
    end do

    call mpi_allreduce(mpi_in_place,L2,3,mpi_real8,mpi_sum,mpi_comm_world,error)
    call mpi_allreduce(mpi_in_place,Linf,3,mpi_real8,mpi_max,mpi_comm_world,error)
    
    if (base_grid%rank == 0) write(e_id,*) Nx, Linf, sqrt(L2)

  end subroutine compute_error
  !=======================================================================================
  
end program ABC
