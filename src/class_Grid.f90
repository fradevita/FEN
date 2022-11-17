!> This module contains definitions and procedures for the class Grid.
module class_Grid

   ! The class grid represent a computational grid over which solving governing equations.
   ! The grid class must contains information about: nodes location (logical and physical),
   ! boundary conditions (logical and physical) and parallel distribution.

   use precision, only : dp
   use constants, only : Ndim

   implicit none

   private
   public :: bc_type, Grid, base_grid

   ! This type is necessary to have an array of strings of variable size.
   ! It is used to specify physical boundary conditions on the domain.
   type bc_type
      character(:), allocatable :: s
   end type bc_type

   type Grid
      ! Total number of points
      integer :: Nx, Ny, Nz

      ! Local size of the grid based on stencil orientation, default is x.
      integer, dimension(3) :: lo, hi, lo_y, hi_y, lo_z, hi_z

      ! Physical size
      real(dp) :: Lx, Ly, Lz

      ! Origin of the grid
      real(dp) :: origin(3)

      ! Grid spacing
      real(dp) :: delta

      ! Cell-center position array
      real(dp), dimension(:), allocatable :: x, y, z

      ! Physical boundary conditions on the grid. It can be:
      ! - Periodic
      ! - Wall
      ! - Inflow
      ! - Outflow
      type(bc_type) :: boundary_conditions(Ndim*2)

      ! Logical array for selecting periodic boundary conditions
      logical :: periodic_bc(3)

      ! Local rank (rank), total number of ranks (nranks), parallel grid 
      ! (prow x pcol)
      integer :: rank, nranks, prow, pcol

   contains
      procedure, pass(self) :: setup
      procedure, pass(self) :: closest_grid_node
      procedure, pass(self) :: destroy
   end type Grid

   ! Define a default base grid for the simulation.
   type(Grid) :: base_grid

contains

   !========================================================================================
   subroutine setup(self, Nx, Ny, Nz, Lx, Ly, Lz, x0, prow, pcol, bc)

      ! This subroutine setup all grid variables for the simulation

      use decomp_2d, only : decomp_2d_init, nrank, nproc, xstart, xend, ystart, yend, zstart, zend
      use io       , only : print_error_message

      ! In/Out variables
      class(grid)  , intent(inout) :: self
      integer      , intent(in   ) :: Nx, Ny, Nz, prow, pcol
      real(dp)     , intent(in   ) :: Lx, Ly, Lz, x0(3)
      type(bc_type), intent(in   ) :: bc(:)

      ! Local variables
      integer :: i
      logical :: pbc(3)

      ! Set physical boundary condition on the grid
      self%boundary_conditions = bc
      pbc = [.false., .false., .false.]
      if (bc(1)%s == 'Periodic' .and. bc(2)%s == 'Periodic') then
         pbc(1) = .true.
      endif
      if (bc(3)%s == 'Periodic' .and. bc(4)%s == 'Periodic') then
         pbc(2) = .true.
      endif
      pbc(3) = .true.

      ! Create the domain decomposition with 2decomp
      call decomp_2d_init(Nx, Ny, Nz, prow, pcol, pbc)

      ! Set the grid size
      self%Nx = Nx
      self%Ny = Ny
      self%Nz = Nz
      self%Lx = Lx
      self%Ly = Ly
      self%Lz = Lz

      ! Set the origin of the grid
      self%origin = x0

      ! Set the grid spacing
      self%delta = Lx/float(Nx)

      ! Check that the spacing is equal in all directions
      if (Lx/float(Nx) .ne. Ly/float(Ny)) then
         if (nrank == 0) then
            call print_error_message('The grid spacing must be equal in all directions')
            call print_error_message('Aborting simulation ...')
         endif

         ! TODO **** call abort subroutine ****

      endif

      ! Set the cell-center grid value
      allocate(self%x(0:Nx+1))
      do i = 0,Nx+1
         self%x(i) = x0(1) + (i - 0.5_dp)*self%delta
      end do
      allocate(self%y(0:Ny+1))
      do i = 0,Ny+1
         self%y(i) = x0(2) + (i - 0.5_dp)*self%delta
      end do
      allocate(self%z(0:Nz+1))
      do i = 0,Nz+1
         self%z(i) = x0(3) + (i - 0.5_dp)*self%delta
      end do

      ! Set the local grid size
      self%lo = xstart
      self%hi = xend
      self%lo_y = ystart
      self%hi_y = yend
      self%lo_z = zstart
      self%hi_z = zend

      ! Set the periodicity
      self%periodic_bc = pbc

      ! Set the parallel variables
      self%rank = nrank
      self%nranks = nproc
      self%prow = prow
      self%pcol = pcol

   end subroutine setup
   !========================================================================================

   !========================================================================================
   function closest_grid_node(self, xl, ind) result(ie)

      ! This function returns indexes ie of the closest Eulerian point to the point xl.
      ! Ind is an integer index used to select the node location:
      ! 0 = cell center
      ! 1 = x face
      ! 2 = y face
      ! 3 = z face

      use io, only : print_error_message

      ! In/Out variables
      class(Grid), intent(in) :: self
      integer    , intent(in) :: ind
      real(dp)   , intent(in) :: xl(:)

      ! Local variables
      integer  :: ie(Ndim)
      real(dp) :: stagger(Ndim)

      select case (ind)
       case(0)
         stagger = 0.0_dp
       case(1)
#if DIM==3
         stagger = [self%delta*0.5_dp, 0.0_dp, 0.0_dp]
#else
         stagger = [self%delta*0.5_dp, 0.0_dp]
#endif
       case(2)
#if DIM==3
         stagger = [0.0_dp, self%delta*0.5_dp, 0.0_dp]
#else
         stagger = [0.0_dp, self%delta*0.5_dp]
#endif

#if DIM==3
       case(3)
         stagger = [0.0_dp, 0.0_dp, self%delta*0.5_dp]
#endif
       case default
         call print_error_message('ERROR IN CLOSEST GRID NODE')
      end select

      ie(1) = minloc(abs(self%x + stagger(1) - xl(1)),1)
      ie(2) = minloc(abs(self%y + stagger(2) - xl(2)),1)
#if DIM==3
      ie(3) = minloc(abs(self%z + stagger(3) - xl(3)),1)
#endif

   end function closest_grid_node
   !========================================================================================

   !========================================================================================
   subroutine destroy(self)

      ! Free the memory allocated by the grid obj.

      use decomp_2d, only : decomp_2d_finalize

      ! In/Out variables
      class(Grid), intent(inout) :: self

      ! Free the memory allocate by the grid
      deallocate(self%x, self%y, self%z)

      ! Free the memory of the lib 2decomp
      call decomp_2d_finalize

   end subroutine destroy
   !========================================================================================

end module class_Grid