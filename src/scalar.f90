!< Define the scalar class
module scalar_mod

    ! This module contains the scalar type definition and procedures.
    ! The scalar field type has a 3D array, an integer for the number of ghost
    ! nodes and a field_bc for the boundary condition values.

    use precision_mod, only : dp
    use grid_mod     , only : grid

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
        !< scalar class
        real(dp), allocatable :: f(:,:,:)    !< Array for the field values
        integer               :: gl = 0      !< Number of ghost nodes per side
        type(grid), pointer   :: G => Null() !< Pointer to the grid where the scalar is defined
        type(field_bc)        :: bc          !< Boundary condition data structure
    contains
        procedure, pass(self) :: allocate           !< allocate memory
        procedure, pass(self) :: set_from_function  !< set scalar from given function
        procedure, pass(self) :: update_ghost_nodes !< fill ghost nodes
        procedure, pass(self) :: max_value          !< compute maximum scalar value
        procedure, pass(self) :: integral           !< compute scalar integral
        procedure, pass(self) :: write              !< write scalar to file
        procedure, pass(self) :: destroy
    end type scalar

contains

    !==============================================================================================
    subroutine allocate(self, G, l)

        ! Allocate all variables for the scalar s.
        ! l: level of ghost nodes (optional).

        ! In/Out variables
        class(scalar), intent(inout)           :: self !< scalar object
        type(grid)   , intent(in   ), target   :: G    !< grid where the scalar is defined
        integer      , intent(in   ), optional :: l    !< ghost node level

        ! Set the pointer to the grid
        self%G => G

        ! By default the field has zero ghost nodes, if l is present set the number of
        ! ghost nodes per side equal to l
        if (present(l)) self%gl = l
        allocate(self%f(G%lo(1) - self%gl:G%hi(1) + self%gl, &
                        G%lo(2) - self%gl:G%hi(2) + self%gl, &
                        G%lo(3) - self%gl:G%hi(3) + self%gl))

        ! Set the field to zero
        self%f = 0.0_dp

        ! If using ghost node
        if (self%gl > 0) then
            ! Allocate the boundary conditions array
            allocate(self%bc%l(G%lo(2)-self%gl:G%hi(2)+self%gl,G%lo(3)-self%gl:G%hi(3)+self%gl))
            allocate(self%bc%r(G%lo(2)-self%gl:G%hi(2)+self%gl,G%lo(3)-self%gl:G%hi(3)+self%gl))
            allocate(self%bc%t(G%lo(1)-self%gl:G%hi(1)+self%gl,G%lo(3)-self%gl:G%hi(3)+self%gl))
            allocate(self%bc%b(G%lo(1)-self%gl:G%hi(1)+self%gl,G%lo(3)-self%gl:G%hi(3)+self%gl))
#if DIM==3
            allocate(self%bc%f(G%lo(1)-self%gl:G%hi(1)+self%gl,G%lo(2)-self%gl:G%hi(2)+self%gl))
            allocate(self%bc%e(G%lo(1)-self%gl:G%hi(1)+self%gl,G%lo(2)-self%gl:G%hi(2)+self%gl))
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
    !==============================================================================================

    !==============================================================================================
    subroutine set_from_function(self, f)

        use function_mod, only : function_type

        ! In/Out variables
        class(scalar)      , intent(inout) :: self !< input scalar field 
        type(function_type), intent(in   ) :: f    !< input function 

        ! Local variables
        integer :: i, j, k

        do k = self%G%lo(3),self%G%hi(3)
            do j = self%G%lo(2),self%G%hi(2)
                do i = self%G%lo(1),self%G%hi(1)
#if DIM==3
                    self%f(i,j,k) = f%fp([self%G%x(i), self%G%y(j), self%G%z(k)], f%args)
#else
                    self%f(i,j,k) = f%fp([self%G%x(i), self%G%y(j)], f%args)
#endif

                end do
            end do
        end do

        ! if the scalar has ghost nodes, update them
        if (self%gl > 0) call self%update_ghost_nodes()

    end subroutine set_from_function
    !==============================================================================================

    !==============================================================================================
    real(dp) function max_value(self)
#ifdef MPI
        use mpi
        use global_mod, only : ierror
#endif
        
        ! In/Out variables
        class(scalar), intent(in) :: self

        ! Compute the maximum value of a scalar field
        max_value = maxval(self%f(self%G%lo(1):self%G%hi(1), &
                                  self%G%lo(2):self%G%hi(2), &
                                  self%G%lo(3):self%G%hi(3)))

#ifdef MPI
        call mpi_allreduce(mpi_in_place,max_value,1,mpi_real8,mpi_max,mpi_comm_world,ierror)
#endif

    end function max_value
    !==============================================================================================

    !==============================================================================================
    real(dp) function integral(self)
#ifdef MPI
        use mpi
        use global_mod, only : ierror
#endif
        ! In/Out variables
        class(scalar), intent(in) :: self !< input scalar

        ! Compute the maximum value of a scalar field
        integral = sum(self%f(self%G%lo(1):self%G%hi(1), &
                              self%G%lo(2):self%G%hi(2), &
                              self%G%lo(3):self%G%hi(3)))
#ifdef MPI
        call mpi_allreduce(mpi_in_place,integral,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
#endif

        integral = integral*self%G%delta**3

    end function integral
    !==============================================================================================

    !==============================================================================================
    subroutine update_ghost_nodes(self)

        ! This subroutine update ghost nodes for the scalar s.
#ifdef MPI
        use mpi
        use global_mod, only : ierror
        use halo_mod  , only : update_halos

#endif
        use IO_mod  , only : print_error_message

        ! In/Out variables
        class(scalar), intent(inout) :: self !< input scalar

        ! Local variables
        integer  :: lo(3), hi(3), l

        ! To shorten notation
        lo = self%G%lo
        hi = self%G%hi

#ifdef MPI
        ! If using MPI, first update halos
        call update_halos(self%f, self%G, self%gl)
#endif

        ! Left Boundary
        if (self%bc%type_l == 0) then ! Periodic
            do l = 1,self%gl
                self%f(lo(1)-l,:,:) = self%f(hi(1)-l+1,:,:)
            end do
        elseif (self%bc%type_l == 1) then ! Dirichlet
            self%f(lo(1)-1,:,:) = 2.0_dp*self%bc%l(:,:) - self%f(lo(1),:,:)
        elseif (self%bc%type_l == 2) then ! Neumann
            self%f(lo(1)-1,:,:) = self%f(lo(1),:,:)
        else
            call print_error_message('ERROR: wrong left boundary condition type for scalar s')
        endif

        ! Right Boundary
        if (self%bc%type_r == 0) then ! Periodic
            do l = 1,self%gl
                self%f(hi(1)+l,:,:) = self%f(lo(1)+l-1,:,:)
            end do
        elseif (self%bc%type_r == 1) then ! Dirichlet
            self%f(hi(1)+1,:,:) = 2.0_dp*self%bc%r(:,:) - self%f(hi(1),:,:)
        elseif (self%bc%type_r == 2) then ! Neumann
            self%f(hi(1)+1,:,:) = self%f(hi(1),:,:)
        else
            call print_error_message('ERROR: wrong rigth boundary condition type for scalar s')
        endif

        ! Bottom Boundary
        if (lo(2) == 1) then
            if (self%bc%type_b == 0) then ! Periodic
                ! When using 2decomp decomposition this ghost nodes is updated automatically
                ! Must be overwritten only when using 1 proc
#ifdef MPI
                if (self%G%prow == 1) then
                    do l = 1,self%gl
                        self%f(:,lo(2)-1,:) = self%f(:,hi(2),:)
                    end do
                endif
#else
                do l = 1,self%gl
                    self%f(:,lo(2)-1,:) = self%f(:,hi(2),:)
                end do
#endif
            elseif (self%bc%type_b == 1) then ! Dirichlet
                self%f(:,lo(2)-1,:) = 2.0_dp*self%bc%b(:,:) - self%f(:,lo(2),:)
            elseif (self%bc%type_b == 2) then ! Neumann
                self%f(:,lo(2)-1,:) = self%f(:,lo(2),:)
            else
                call print_error_message('ERROR: wrong bottom boundary condition type for scalar s')
            endif
        endif

        ! Top Boundary
        if (hi(2) == self%G%Ny) then
            if (self%bc%type_t == 0) then ! Periodic
                ! When using 2decomp decomposition this ghost nodes is updated automatically
                ! Must be overwritten only when using 1 proc
#ifdef MPI
                if (self%G%prow == 1) then
                    do l = 1,self%gl
                        self%f(:,hi(2)+1,:) = self%f(:,lo(2),:)
                    end do
                endif
#else
                do l = 1,self%gl
                    self%f(:,hi(2)+1,:) = self%f(:,lo(2),:)
                end do
#endif
            elseif (self%bc%type_t == 1) then ! Dirichlet
                self%f(:,hi(2)+1,:) = 2.0_dp*self%bc%t(:,:) - self%f(:,hi(2),:)
            elseif (self%bc%type_t == 2) then ! Neumann
                self%f(:,hi(2)+1,:) = self%f(:,hi(2),:)
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

#ifdef MPI
        ! This is necessary to avoid deadlock
        call mpi_barrier(mpi_comm_world, ierror)
#endif

    end subroutine
    !==============================================================================================

    !==============================================================================================
    subroutine write(self, filename)

        use decomp_2d_io, only : decomp_2d_write_one

        class(scalar)   , intent(in) :: self     !< scalar field to write
        character(len=*), intent(in) :: filename !< output file

        call decomp_2d_write_one(1, self%f(self%G%lo(1):self%G%hi(1), &
                                           self%G%lo(2):self%G%hi(2), &
                                           self%G%lo(3):self%G%hi(3)), filename)

    end subroutine write
    !==============================================================================================

    !==============================================================================================
    subroutine destroy(self)

        class(scalar), intent(inout) :: self

        ! Free all memory allocated by the scalar s
        deallocate(self%f)

        if (self%gl > 0) deallocate(self%bc%l, self%bc%r, self%bc%t, self%bc%b)
#if DIM==3
        if (self%gl > 0) deallocate(self%bc%f, self%bc%e)
#endif

    end subroutine destroy
    !==============================================================================================

end module
