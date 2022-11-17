program isotropic

  ! Test case for the forces isotropic turbulance, comparison with data from basilisk
  ! and hit3d.

  use mpi
  use precision    , only : dp
  use constants    , only : pi
  use class_Grid   , only : base_grid, bc_type
  use solver
  use navier_stokes, only : viscosity, set_timestep, v, constant_CFL, CFL, dt_o
  use io
  use json

  implicit none

  ! Variables
  integer  :: ierror, Nx, Ny, Nz, step, oid
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

  ! Set the resolution
  Nx = 128
  Ny = 128
  Nz = 128

  ! Set boundary conditions
  bc(1)%s = 'Periodic'
  bc(2)%s = 'Periodic'
  bc(3)%s = 'Periodic'
  bc(4)%s = 'Periodic'
  bc(5)%s = 'Periodic'
  bc(6)%s = 'Periodic'
  
  ! Create the grid
  call base_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, [0.0_dp, 0.0_dp, 0.0_dp], 2, 4, bc)
  
  ! Set the viscosity
  viscosity = 1.0e-2_dp

  ! Initialize the solver
  call init_solver
  step = 0
  time = 0.0_dp

  ! Compute the timestep
  call set_timestep(dt, 1.0_dp)
  constant_CFL = .true.
  CFL = 0.3_dp
  
  ! Initial condition
  call init_fields

  call print_setup_json(dt)
  call save_fields(step)

  ! Open output file
  open(newunit = oid, file = 'output.txt')

  !==== Start Time loop ===================================================================
  time_loop: do

     step = step + 1
     time = time + dt

     ! Apply the linear forcing
     call forcing()

     ! Advance in time the solution
     call advance_solution(step, dt)

     ! Advance solver status to log file
     call print_solver_status(stdout, step, time, dt)

     ! Compute output
     call output
     if (mod(step,1000)==0) call save_fields(step)

  end do time_loop

contains

  !=======================================================================================
  subroutine init_fields

    use navier_stokes, only : v

    integer :: i, j, k
    real(dp) :: x, y, z

    ! Set the ABC flow
    do k = base_grid%lo(3),base_grid%hi(3)
       do j = base_grid%lo(2),base_grid%hi(2)
          do i = base_grid%lo(1),base_grid%hi(1)
             v%x%f(i,j,k) = cos(base_grid%y(j)) + sin(base_grid%z(k)) + rand()*1.0d-2
             v%y%f(i,j,k) = sin(base_grid%x(i)) + cos(base_grid%z(k)) + rand()*1.0d-2
             v%z%f(i,j,k) = cos(base_grid%x(i)) + sin(base_grid%y(j)) + rand()*1.0d-2
          end do
       end do
    end do

    ! Update halos and bc
    call v%apply_bc()

  end subroutine init_fields
  !=======================================================================================

  !=======================================================================================
  subroutine forcing

    use navier_stokes, only : v, S

    integer  :: i, j, k
    real(dp) :: ubar, vbar, wbar 

    ! Compute average velocity in each direction
    ubar = 0.0_dp
    vbar = 0.0_dp
    wbar = 0.0_dp
    do k = base_grid%lo(3),base_grid%hi(3)
       do j = base_grid%lo(2),base_grid%hi(2)
          do i = base_grid%lo(1),base_grid%hi(1)
             ubar = ubar + v%x%f(i,j,k)
             vbar = vbar + v%y%f(i,j,k)
             wbar = wbar + v%z%f(i,j,k)
          end do
       end do
    end do

    call mpi_allreduce(mpi_in_place,ubar,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
    call mpi_allreduce(mpi_in_place,vbar,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
    call mpi_allreduce(mpi_in_place,wbar,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
    
    ubar = ubar/(base_grid%nx*base_grid%ny*base_grid%nz)
    vbar = vbar/(base_grid%nx*base_grid%ny*base_grid%nz)
    wbar = wbar/(base_grid%nx*base_grid%ny*base_grid%nz)
   
    ! Force the velocity field
    do k = base_grid%lo(3),base_grid%hi(3)
       do j = base_grid%lo(2),base_grid%hi(2)
          do i = base_grid%lo(1),base_grid%hi(1)
             S%x%f(i,j,k) = 0.1_dp*(v%x%f(i,j,k) - ubar)
             S%y%f(i,j,k) = 0.1_dp*(v%y%f(i,j,k) - vbar)
             S%z%f(i,j,k) = 0.1_dp*(v%z%f(i,j,k) - wbar)
          end do
       end do
    end do

  end subroutine forcing
  !=======================================================================================
    
  !=======================================================================================
  subroutine output

    use navier_stokes, only : v

    ! Local variables
    integer  :: i, j, k
    real(dp) :: ubar, vbar, wbar, ke
    
    ! Compute average velocity in each direction
    ubar = 0.0_dp
    vbar = 0.0_dp
    wbar = 0.0_dp
    do k = base_grid%lo(3),base_grid%hi(3)
       do j = base_grid%lo(2),base_grid%hi(2)
          do i = base_grid%lo(1),base_grid%hi(1)
             ubar = ubar + v%x%f(i,j,k)
             vbar = vbar + v%y%f(i,j,k)
             wbar = wbar + v%z%f(i,j,k)
          end do
       end do
    end do

    call mpi_allreduce(mpi_in_place,ubar,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
    call mpi_allreduce(mpi_in_place,vbar,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
    call mpi_allreduce(mpi_in_place,wbar,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
    
    ubar = ubar/(base_grid%nx*base_grid%ny*base_grid%nz)
    vbar = vbar/(base_grid%nx*base_grid%ny*base_grid%nz)
    wbar = wbar/(base_grid%nx*base_grid%ny*base_grid%nz)

    ! Compute kinetic energy and viscous dissipation
    ke = 0.0_dp
    do k = base_grid%lo(3),base_grid%hi(3)
       do j = base_grid%lo(2),base_grid%hi(2)
          do i = base_grid%lo(1),base_grid%hi(1)
             ke = ke + (v%x%f(i,j,k) - ubar)**2 + (v%y%f(i,j,k) - vbar)**2 + &
                       (v%z%f(i,j,k) - wbar)**2
          end do
       end do
    end do

    call mpi_allreduce(mpi_in_place,ke,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)

    if (base_grid%rank == 0) then
       write(oid,*) time, 0.5_dp*ke/float(nx*ny*nz) 
       flush(oid)
    endif
    
  end subroutine output
  !=======================================================================================
  
end program isotropic
