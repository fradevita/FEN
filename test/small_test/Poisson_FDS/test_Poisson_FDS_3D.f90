program test_Poisson_FDS_2D

   ! Thest the Fast Direct Solver for the solution of the Poisson equation.

   use mpi
   use constants , only : pi
   use precision , only : dp
   use class_Grid, only : base_grid, bc_type
   use class_Scalar
   use Poisson   , only : init_Poisson_Solver, solve_Poisson, destroy_Poisson_solver

   implicit none

   integer      :: ierror, prow, pcol, n, Nx, Ny, Nz, i, j, k
   real(dp)     :: Lx, Ly, Lz, origin(3), e, emax
   type(scalar) :: phi
   type(bc_type) :: bc(6)
   character(len=1) :: sN
   character(len=5) :: filename

   ! Initialize MPI
   call mpi_init(ierror)

   ! Use only 1 rank
   prow = 1
   pcol = 1

   ! The test is performed with several different conditions.

   ! First test is with Neumann boundary condition on all boundaries.
   ! The domain has four walls, in order to select neumann BC on the pressure
   bc(1)%s = 'Periodic'
   bc(2)%s = 'Periodic'
   bc(3)%s = 'Periodic'
   bc(4)%s = 'Periodic'
   bc(5)%s = 'Periodic'
   bc(6)%s = 'Periodic'

   ! Open the output for the first test
   open(1, file = 'error_PPP')

   ! Solve Poisson equation for 5 levels of refinement
   do n = 1,5

      ! Create the grid
      Nx = 16*2**(n-1)
      Ny = 16*2**(n-1)
      Nz = 16*2**(n-1)
      Lx = 1.0_dp
      Ly = 1.0_dp
      Lz = 1.0_dp
      origin = [0.0_dp, 0.0_dp, 0.0_dp]
      call base_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, prow, pcol, bc(1:4))

      ! Create one scalar field for the RHS of the Poisson equation
      call phi%allocate()

      ! Set Poisson RHS
      do k = base_grid%lo(3), base_grid%hi(3)
         do j = base_grid%lo(2),base_grid%hi(2)
            do i = base_grid%lo(1),base_grid%hi(1)
               phi%f(i,j,k) = -12.0_dp*pi*pi*sin(2.0_dp*pi*base_grid%x(i))* &
                    cos(2.0_dp*pi*base_grid%y(j))*sin(2.0_dp*pi*base_grid%z(k))
            end do
         end do
      end do

      ! Initialize the Poisson Solver
      call init_Poisson_Solver(phi)

      ! Solve the Poisson equation
      call solve_Poisson(phi)

      ! Compute the error
      emax = 0.0_dp
      do k = base_grid%lo(3), base_grid%hi(3)
         do j = base_grid%lo(2),base_grid%hi(2)
            do i = base_grid%lo(1),base_grid%hi(1)
               e = abs(phi%f(i,j,k) - sin(2.0_dp*pi*base_grid%x(i))* &
                    cos(2.0_dp*pi*base_grid%y(j))*sin(2.0_dp*pi*base_grid%z(k)))
               if (e > emax) emax = e
            end do
         end do
      end do

      ! Write the error
      write(1,*) Nx, emax

      ! Free the memory
      call base_grid%destroy
      call phi%destroy()
      call destroy_Poisson_Solver

   end do
   close(1)

   ! Finalize the simulation
   call MPI_FINALIZE(ierror)

end program test_Poisson_FDS_2D
