program ellipse

  use mpi
  use precision_mod
  use global_mod        , only : pi, myrank
  use grid_mod
  use solver_mod
  use navier_stokes_mod , only : viscosity, set_timestep, g
  use ibm_mod
  use IO_mod

  implicit none

  ! Parameters
  ! Scales: a, Vt, rhof
  real(dp), parameter :: a = 1.0_dp     ! Ellipse major axis
  real(dp), parameter :: Vt = 1.0_dp    ! Terminal settling velocity
  real(dp), parameter :: rhof = 1.0_dp  ! Fluid density
  real(dp), parameter :: alpha = 2.0_dp ! Particle aspect ratio: a/b
  real(dp), parameter :: eta = 4.0_dp   ! Blockage ratio: L/a
  real(dp), parameter :: gamma = 1.1_dp ! Density ratio: rhos/rhof
  real(dp), parameter :: Fr = 0.126_dp  ! Froud number: Vt/sqrt(g*a)
  real(dp), parameter :: Re = 12.5_dp   ! Reynolds number: Vt*a/nu

  ! Derived parameters
  real(dp), parameter :: b = a/alpha            ! Particle minor axis
  real(dp), parameter :: L = a*eta              ! Channel width
  real(dp), parameter :: gravity = (Vt/Fr)**2/a ! Gravity 
  real(dp), parameter :: rhos = gamma*rhof      ! Particle density
  real(dp), parameter :: nu = Vt*a/Re           ! Fluid viscosity

  ! Variables
  integer  :: ierror, Nx, Ny, Nz, step
  real(dp) :: Lx, Ly, Lz, time, dt, origin(3)
  type(grid)    :: comp_grid
  type(bc_type) :: bc(4)
  type(lagrangian_solid), target :: S
  
  ! Initialize MPI
  call mpi_init(ierror)
  call mpi_comm_rank(mpi_comm_world, myrank, ierror)

  ! The domain is a square unit box
  Lx = L
  Ly = 7.0_dp*L
  
  ! Set the resolution
  Nx = 128
  Ny = Nx*7

  ! Since 2D
  Nz = 1
  Lz = Lx*float(Nz)/float(Nx)

  origin = [0.0_dp, 0.0_dp, 0.0_dp]

  ! Set boundary conditions
  bc(1)%s = 'Wall'
  bc(2)%s = 'Wall'
  bc(3)%s = 'Wall'
  bc(4)%s = 'Wall'
  
  ! Create the grid
  call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 4, 1, bc)
  call comp_grid%print_json

  ! Create the cylinder
  call S%create('mesh.txt', name = 'S')

  ! Allocate the array of solid
  allocate(solid_list(1))
  solid_list(1)%pS => S
  
  ! Set the mass of the solid body
  S%M(1:2) = rhos*pi*a*b*0.25_dp
  S%M(3) = 0.2_dp*S%M(1)*((a/2.0_dp)**2 + (b/2.0)**2)

  ! Set the viscoisty
  viscosity = nu/rhof

  ! Set body force
  g(2) = -gravity

  ! Initialize the solver
  call init_solver(comp_grid)
  step = 0
  time = 0.0_dp
  call set_timestep(comp_grid, dt, 1.0_dp)

  ! Open output file
  if (myrank == 0) then
     call solid_list(1)%pS%write_csv(time)
  endif

  !==== Start Time loop ===================================================================
  time_loop: do while(time < 22.5_dp)

     step = step + 1
     time = time + dt

     ! Advance in time the solution
     call advance_solution(comp_grid, step, dt)

     ! Advance solver status to log file
     call print_solver_status(stdout, step, time, dt)

     ! Print position
     if (mod(step,100) == 0) then
        if (myrank == 0) then
           call solid_list(1)%pS%write_csv(time)
        endif
     endif
  end do time_loop

  ! free memory
  call destroy_solver

  ! Finalize the simulation
  call MPI_FINALIZE(ierror)

end program ellipse
