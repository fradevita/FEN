program main

  use mpi
  use precision    , only : dp
  use class_Grid   , only : base_grid, bc_type
  use class_Scalar                    
  use class_Vector
  use fields       , only : curl
  use solver       , only : init_solver, advance_solution, destroy_solver, print_solver_status
  use navier_stokes, only : viscosity, v, set_timestep, constant_CFL, CFL
  use json
  use io
  
  implicit none

  integer       :: ierror, Nx, Ny, Nz, step
  real(dp)      :: Lx, Ly, Lz, time, dt, e, origin(3)
  type(bc_type) :: bc(4)
  type(scalar)  :: uo
  type(vector)  :: omega

  ! Initialize MPI
  call mpi_init(ierror)

  ! The domain is a unit squared box
  Lx = 1.0_dp
  Ly = 1.0_dp
  
  ! Set the resolution
  Nx = 64
  Ny = 64
  Nz = 1
  Lz = Lx*float(Nz)/float(Nx)
  origin = [0.0_dp, 0.0_dp, 0.0_dp]
  bc(1)%s = 'Wall'
  bc(2)%s = 'Wall'
  bc(3)%s = 'Wall'
  bc(4)%s = 'Wall'

  ! Create the grid
  call base_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 4, 1, bc)

  ! Set the viscosity equal to 1/Re
  viscosity = 1.0e-3_dp

  ! Initialize the solver
  call init_solver
  step = 0
  time = 0.0_dp
  call set_timestep(dt, 1.0_dp)
  constant_CFL = .true.
  CFL = 0.8_dp

  ! If you want to use constant CFL uncomment the following two lines
  !constant_CFL = .true.
  !CFL = 0.9_dp
  
  ! Create a scalar field to store the old velocity for steady state check
  call uo%allocate(1)

  ! Set top bc for u equal to 1
  v%x%bc%t = 1.0_dp

  call print_setup_json(dt)

  !==== Start Time loop ===================================================================
  time_loop: do

     step = step + 1
     time = time + dt

     ! Advance in time the solution
     call advance_solution(step, dt)

     ! Advance solver status to log file
     call print_solver_status(stdout, step, time, dt)

     ! Check for steady state
     e = maxval(abs(v%x%f - uo%f))
     call mpi_allreduce(mpi_in_place,e,1,mpi_real8, mpi_max,mpi_comm_world,ierror)
     if ((e < 1.0D-8) .and. step > 1) then
       exit time_loop
     endif
     uo = v%x

  end do time_loop

  call omega%allocate()
  call curl(v, omega)

  ! Print velocity field
  call v%x%write('u.dat')
  call v%y%write('v.dat')
  call omega%x%write('vort.dat')
  
  ! Free the memory allocated by the solver
  call destroy_solver()
  call uo%destroy()
  
  ! Finalize the simulation
  call MPI_FINALIZE(ierror)

end program main
