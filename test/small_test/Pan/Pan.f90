program Pan

   ! Particle migration in a channel flow at Re 12

   use mpi
   use precision            , only : dp
   use constants            , only : pi
   use io                   , only : stdout
   use class_Grid           , only : base_grid, bc_type
   use class_eulerian_circle, only : Circle
   use ibm
   use eulerian_ibm
   use navier_stokes        , only : set_timestep, g, viscosity
   use solver
   use json

   implicit none

   ! Parameters
   real(dp), parameter :: radius = 0.125_dp        !< cylinder radius
   real(dp), parameter :: Deltap = 1.763e-3_dp     !< external body force
   real(dp), parameter :: mu = 3.2498036e-3_dp     !< fluid viscosity

   ! Variables
   integer              :: ierror, Nx, Ny, Nz, step
   real(dp)             :: Lx, Ly, Lz, time, dt, origin(3)
   type(bc_type)        :: bc(4)
   type(Circle), target :: C

   ! Initialize MPI
   call mpi_init(ierror)

   ! The domain is a unit square
   Lx = 1.0_dp
   Ly = 1.0_dp
   
   ! Set the resolution
   Nx = 96
   Ny = 96

   ! Since 2D
   Nz = 1
   Lz = Lx*float(Nz)/float(Nx)
   origin = [0.0_dp, 0.0_dp, 0.0_dp]

   ! Set bc
   bc(1)%s = 'Periodic'
   bc(2)%s = 'Periodic'
   bc(3)%s = 'Wall'
   bc(4)%s = 'Wall'

   ! Create the grid
   call base_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 4, 1, bc)

   ! Set body force
   g(1) = Deltap

   ! Create the solid body
   C = circle(X = [0.5_dp, 0.4_dp, 0.0_dp], R = radius, name = 'Circle')
   allocate(Eulerian_Solid_list(1))
   Eulerian_Solid_list(1)%pS => C
   call C%load_surface_points('mesh.txt')

   ! Set the viscoisty
   viscosity = mu

   ! Initialize the solver
   call init_solver
   step = 0
   time = 0.0_dp
   call set_timestep(dt, 1.0_dp)

   call print_setup_json(dt)

   !==== Start Time loop ===================================================================
   time_loop: do while(time < 1.0e+3_dp)

      step = step + 1
      time = time + dt

      ! Advance in time the solution
      call advance_solution(step, dt)

      ! Advance solver status to log file
      call print_solver_status(stdout, step, time, dt)

      ! Print position
      if (base_grid%rank == 0) then
         call C%print_csv(time)
      endif

   end do time_loop

   ! free memory
   call destroy_solver

   ! Finalize the simulation
   call MPI_FINALIZE(ierror)

end program Pan
