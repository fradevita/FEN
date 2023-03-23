module class_Scalar

    ! This module contains the scalar type definition and procedures.
    ! The scalar field type has a 3D array, an integer for the number of ghost
    ! nodes and a field_bc for the boundary condition values.

    use precision, only : dp

    implicit none

    private
    public :: scalar

    ! Boundary Condition structure
    ! BC type can be: Periodic, Dirichlet, Neumann
    ! Legend:
    ! 0: Periodic
    ! 1: Dirichlet
    ! 2: Neumann
    type field_bc
        integer                               :: type_l, type_r, type_t, type_b
        real(dp), dimension(:,:), allocatable :: l, r, t, b
#if DIM==3
        integer                               :: type_f, type_e
        real(dp), dimension(:,:), allocatable :: f, e
#endif
    end type field_bc

    type scalar
        ! Array for the field values
        real(dp), dimension(:,:,:), allocatable :: f

        ! Number of ghost nodes per side (default is 0)
        integer :: gl = 0

        ! Boundary condition data structure
        type(field_bc) :: bc
    contains
        procedure, pass(self) :: allocate
        procedure, pass(self) :: set_from_function
        procedure, pass(self) :: apply_bc
        procedure, pass(self) :: max_value
        procedure, pass(self) :: integral
        procedure, pass(self) :: write
        procedure, pass(self) :: destroy
    end type scalar

contains

   !========================================================================================
   subroutine allocate(self, l)

        ! Allocate all variables for the scalar s.
        ! l: level of ghost nodes (optional).

        use class_Grid, only : base_grid

        ! In/Out variables
        class(scalar), intent(inout)           :: self
        integer      , intent(in   ), optional :: l

        ! Local variables
        integer :: lo(3), hi(3)

        lo = base_grid%lo
        hi = base_grid%hi

        ! By default the field has zero ghost nodes, if l is present set the number of
        ! ghost nodes per side equal to l
        if (present(l)) self%gl = l
        allocate(self%f(lo(1)-self%gl:hi(1)+self%gl,lo(2)-self%gl:hi(2)+self%gl,lo(3)-self%gl:hi(3)+self%gl))

        ! Set the field to zero
        self%f = 0.0_dp

        ! If using ghost node
        if (self%gl > 0) then
            ! Allocate the boundary conditions array
            allocate(self%bc%l(lo(2)-self%gl:hi(2)+self%gl,lo(3)-self%gl:hi(3)+self%gl))
            allocate(self%bc%r(lo(2)-self%gl:hi(2)+self%gl,lo(3)-self%gl:hi(3)+self%gl))
            allocate(self%bc%t(lo(1)-self%gl:hi(1)+self%gl,lo(3)-self%gl:hi(3)+self%gl))
            allocate(self%bc%b(lo(1)-self%gl:hi(1)+self%gl,lo(3)-self%gl:hi(3)+self%gl))
#if DIM==3
            allocate(self%bc%f(lo(1)-self%gl:hi(1)+self%gl,lo(2)-self%gl:hi(2)+self%gl))
            allocate(self%bc%e(lo(1)-self%gl:hi(1)+self%gl,lo(2)-self%gl:hi(2)+self%gl))
#endif

            ! Set bc to zero by default
            self%bc%l = 0.0_dp
            self%bc%r = 0.0_dp
            self%bc%t = 0.0_dp
            self%bc%b = 0.0_dp
#if DIM==3
            self%bc%f = 0.0_dp
            self%bc%e = 0.0_dp
#endif

            ! Set bc type periodic by default
            self%bc%type_l = 0
            self%bc%type_r = 0
            self%bc%type_t = 0
            self%bc%type_b = 0
#if DIM==3
            self%bc%type_f = 0
            self%bc%type_e = 0
#endif
        endif

   end subroutine allocate
   !========================================================================================

   !========================================================================================
   subroutine set_from_function(self, function)

        use class_Grid, only : base_grid
        use functions, only : function_type

        ! In/Out variables
        class(scalar)      , intent(inout) :: self
        type(function_type), intent(in   ) :: function

        ! Local variables
        integer :: i, j, k

        do k = base_grid%lo(3),base_grid%hi(3)
            do j = base_grid%lo(2),base_grid%hi(2)
                do i = base_grid%lo(1),base_grid%hi(1)
#if DIM==3
                    self%f(i,j,k) = function%f([base_grid%x(i), base_grid%y(j), base_grid%z(k)], function%args)
#else
                    self%f(i,j,k) = function%f([base_grid%x(i), base_grid%y(j)], function%args)
#endif
                end do
            end do
        end do

   end subroutine set_from_function
   !========================================================================================

   !========================================================================================
   subroutine destroy(self)

      class(scalar), intent(inout) :: self

      ! Free all memory allocated by the scalar s
      deallocate(self%f)
      if (self%gl > 0) deallocate(self%bc%l, self%bc%r, self%bc%t, self%bc%b)
#if DIM==3
      if (self%gl > 0) deallocate(self%bc%f, self%bc%e)
#endif

   end subroutine destroy
   !========================================================================================

   !========================================================================================
   real(dp) function max_value(self)

      use mpi
      use class_Grid, only : base_grid

      ! In/Out variables
      class(scalar), intent(in) :: self

      ! Local variables
      integer :: lo(3), hi(3), error

      lo = base_grid%lo
      hi = base_grid%hi

      ! Compute the maximum value of a scalar field
      max_value = maxval(self%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

      if (base_grid%nranks > 1) then
         call mpi_allreduce(mpi_in_place,max_value,1,mpi_real8,mpi_max,mpi_comm_world,error)
      endif

   end function max_value
   !========================================================================================

   !========================================================================================
   real(dp) function integral(self)

      use mpi
      use class_Grid, only : base_grid

      ! In/Out variables
      class(scalar), intent(in) :: self

      ! Local variables
      integer :: lo(3), hi(3), error
      real(dp) :: dV

      lo = base_grid%lo
      hi = base_grid%hi

      ! Compute the maximum value of a scalar field
      integral = sum(self%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

      if (base_grid%nranks > 1) then
         call mpi_allreduce(mpi_in_place,integral,1,mpi_real8,mpi_sum,mpi_comm_world,error)
      endif

      dV = base_grid%delta**3
      integral = integral*dV

   end function integral
   !========================================================================================

   !========================================================================================
   subroutine apply_bc(self)

      ! This subroutine update ghost nodes for the scalar s.

      use mpi
      use io        , only : print_error_message
      use class_grid, only : base_grid
      use halo      , only : field_halo_update

      ! In/Out variables
      class(scalar), intent(inout) :: self

      ! Local variables
      integer  :: ierror, lo(3), hi(3), l

      lo = base_grid%lo
      hi = base_grid%hi

      ! First update halos
      call field_halo_update(self%f, self%gl)

      ! Then apply phiscal bc

      ! Left Boundary
      if (self%bc%type_l == 0) then ! Periodic
         do l = 1,self%gl
            self%f(lo(1)-l,:,:) = self%f(hi(1)-l+1,:,:)
         end do
      elseif (self%bc%type_l == 1) then ! Dirichlet
         self%f(base_grid%lo(1)-1,:,:) = 2.0_dp*self%bc%l(:,:) - self%f(base_grid%lo(1),:,:)
      elseif (self%bc%type_l == 2) then ! Neumann
         self%f(base_grid%lo(1)-1,:,:) = self%f(base_grid%lo(1),:,:)
      else
         call print_error_message('ERROR: wrong left boundary condition type for scalar s')
      endif

      ! Right Boundary
      if (self%bc%type_r == 0) then ! Periodic
         do l = 1,self%gl
            self%f(hi(1)+l,:,:) = self%f(lo(1)+l-1,:,:)
         end do
      elseif (self%bc%type_r == 1) then ! Dirichlet
         self%f(base_grid%hi(1)+1,:,:) = 2.0_dp*self%bc%r(:,:) - self%f(base_grid%hi(1),:,:)
      elseif (self%bc%type_r == 2) then ! Neumann
         self%f(base_grid%hi(1)+1,:,:) = self%f(base_grid%hi(1),:,:)
      else
         call print_error_message('ERROR: wrong rigth boundary condition type for scalar s')
      endif

      ! Bottom Boundary
      !if (base_grid%rank == 0) then
      if (base_grid%lo(2) == 1) then
         if (self%bc%type_b == 0) then ! Periodic
            ! Need to update this bc only if using 1 rank in y
            if (base_grid%prow == 1) then
               do l = 1,self%gl
                  self%f(:,lo(2)-1,:) = self%f(:,hi(2),:)
               end do
            endif
         elseif (self%bc%type_b == 1) then ! Dirichlet
            self%f(:,base_grid%lo(2)-1,:) = 2.0_dp*self%bc%b(:,:) - self%f(:,base_grid%lo(2),:)
         elseif (self%bc%type_b == 2) then ! Neumann
            self%f(:,base_grid%lo(2)-1,:) = self%f(:,base_grid%lo(2),:)
         else
            call print_error_message('ERROR: wrong bottom boundary condition type for scalar s')
         endif
      endif

      ! Top Boundary
      !if (base_grid%rank == base_grid%nranks - 1) then
      if (base_grid%hi(2) == base_grid%Ny) then
         if (self%bc%type_t == 0) then ! Periodic
            ! Need to update this bc only if using 1 rank in y
            if (base_grid%prow == 1) then
               do l = 1,self%gl
                  self%f(:,hi(2)+1,:) = self%f(:,lo(2),:)
               end do
            endif
         elseif (self%bc%type_t == 1) then ! Dirichlet
            self%f(:,base_grid%hi(2)+1,:) = 2.0_dp*self%bc%t(:,:) - self%f(:,base_grid%hi(2),:)
         elseif (self%bc%type_t == 2) then ! Neumann
            self%f(:,base_grid%hi(2)+1,:) = self%f(:,base_grid%hi(2),:)
         else
            call print_error_message('ERROR: wrong top boundary condition type for scalar s')
         endif
      endif

#if DIM==3
      if (base_grid%pcol == 1) then
         self%f(:,:,lo(3)-1) = self%f(:,:,hi(3))
         self%f(:,:,hi(3)+1) = self%f(:,:,lo(3))
      end if
#endif

      ! This is necessary to avoid deadlock
      call mpi_barrier(mpi_comm_world, ierror)

   end subroutine apply_bc
   !========================================================================================

   !========================================================================================
   subroutine write(self, filename)

      use class_Grid  , only : base_grid
      use decomp_2d_io, only : decomp_2d_write_one

      class(scalar)   , intent(in) :: self
      character(len=*), intent(in) :: filename

      call decomp_2d_write_one(1, self%f(base_grid%lo(1):base_grid%hi(1), &
         base_grid%lo(2):base_grid%hi(2), &
         base_grid%lo(3):base_grid%hi(3)), filename)

   end subroutine write
   !========================================================================================

end module class_Scalar
