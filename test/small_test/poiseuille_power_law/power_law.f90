program poiseuille_power_law

  use mpi
  use io
  use precision    , only : dp
  use class_Grid   , only : base_grid, bc_type
  use class_Scalar
  use solver       , only : init_solver, advance_solution, print_solver_status, destroy_solver
  use navier_stokes, only : g, v, set_timestep
  use non_newtonian, only : compute_viscosity, mu_max, mu_min, power_law, n, kc
  
  implicit none

  integer          :: ierror, Nx, Ny, Nz, step
  real(dp)         :: Lx, Ly, Lz , dt, time, origin(3)
  logical          :: periodic_bc(3)
  type(scalar)     :: uo
  character(len=3) :: sn 
  character(len=5) :: filename
  type(bc_type)    :: bc(4) 

  ! Set the tolerance
  real(dp) :: tol = 1.0d-8

  ! Initialize MPI
  call mpi_init(ierror)

  ! The domain has 4 point in x and the resolution is refined in y
  Nx = 4
  Ny = 32
  Nz = 1

  ! The physical size in y is 1
  Ly = 1.0_dp
  Lx = Ly*float(Nx)/float(Ny)
  Lz = Ly*float(Nz)/float(Ny)

  ! Body force in x
  g(1) = 1.0_dp

  origin = [0.0_dp, 0.0_dp, 0.0_dp]

  bc(1)%s = 'Periodic'
  bc(2)%s = 'Periodic'
  bc(3)%s = 'Wall'
  bc(4)%s = 'Wall'

  ! Create the grid
  call base_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 1, 1, bc)
  
  ! Old u to check steady state
  call uo%allocate(1)

  ! Set the viscosity model
  mu_max = 10.0_dp
  mu_min = 0.01_dp
  n = 2.0_dp
  kc = 1.0_dp
  compute_viscosity => power_law
  
  ! Initialize the solver
  call init_solver
  step = 0
  time = 0.0_dp
  call set_timestep(dt, 1.0_dp)
 
  !==== Start Time loop ===================================================================
  time_loop: do

     step = step + 1
     time = time + dt

     ! Advance in time the solution
     call advance_solution(step, dt)

     ! Print solver status
     call print_solver_status(stdout, step, time, dt)

     ! Check for steady state
     print *, maxval(v%x%f - uo%f)
     if ((maxval(v%x%f - uo%f) < tol .and. step > 2) .or. step > 100000) then
        exit time_loop
     endif
     uo = v%x
  end do time_loop

  ! Write u
  write(sn,'(I0.3)') ny
  filename = 'u_'//sn
  call v%x%write(filename)

  ! Free the memory allocated by the solver
  call destroy_solver()
  call uo%destroy()
  
  ! Finalize the simulation
  call MPI_FINALIZE(ierror)

end program poiseuille_power_law
