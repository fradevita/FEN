program main

  use mpi
  use precision      , only : dp
  use class_Grid     , only : base_grid, bc_type
  use constants      , only : gravity
  use multiphase     , only : rho_0, rho_1, mu_0, mu_1
  use volume_of_fluid, only : distance, vof
  use solver         , only : init_solver, advance_solution, print_solver_status, destroy_solver
  use navier_stokes  , only : set_timestep, g, v
  use io

  implicit none

  integer           :: ierror, Nx, Ny, Nz, step, ne
  real(dp)          :: Lx, Ly, Lz, time, dt, Tmax, time_event, origin(3)
  type(bc_type)     :: bc(4)
  character(len=7)  :: ss
  character(len=16) :: filename

  ! Initialize the MPI library
  call mpi_init(ierror)

  ! Create the grid
  Nx = 128
  Ny = 512
  Nz = 1
  Lx = 1.0_dp
  Ly = 4.0_dp
  Lz = Lx*float(Nz)/float(Nx)
  origin = [0.0_dp, 0.0_dp, 0.0_dp]
  bc(1)%s = 'Periodic'
  bc(2)%s = 'Periodic'
  bc(3)%s = 'Wall'
  bc(4)%s = 'Wall'
  call base_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 1, 1, bc)

  ! Setup multiphase parameter
  rho_0 = 0.1694_dp
  rho_1 = 1.2250_dp
  mu_0 = 3.13d-3
  mu_1 = mu_0

  ! Set gravity
  g(2) = -gravity

  ! Set the initial vof from a distance function
  distance => wave
  
  ! Initialize the solver
  call init_solver
  step = 0
  time = 0.0d0
  call set_timestep(dt, 1.0_dp)

  ! Set maximum time
  Tmax = 1.0_dp

  ! Save initial configuration
  write(ss,'(I0.1)') step
  filename = 'data/vof_'//ss
  call vof%write(filename)

  ! Save the vof from t = 0.7
  time_event = 0.7_dp
  ne = 1
  !==== Start Time loop ===================================================================
  do while(time < Tmax)

     step = step + 1
     time = time + dt

     ! Advance in time the solution
     call advance_solution(step, dt)

     ! Print solver status
     call print_solver_status(stdout, step, time, dt)

     ! Write output
     if (time >= time_event) then
       write(ss,'(I0.1)') ne
       filename = 'data/vof_'//ss
       call vof%write(filename)
       time_event = time_event + 0.1d0
       ne = ne + 1
     endif

  end do
  
  ! Free the memory allocated by the solver
  call destroy_solver()
  
  ! Finalize the simulation
  call MPI_FINALIZE(ierror)

contains

  !==============================================================================
  function wave(x, y) result (d)

    use constants, only : pi

    implicit none

    ! In/Out variables
    real(dp), intent(in) :: x, y

    ! Local variables
    real(dp) :: d

    d = (y - 0.05_dp*cos(2.0_dp*pi*x/Lx) - Ly/2.0_dp)

  end function wave
  !==============================================================================

end program main
