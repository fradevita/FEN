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
    ! BC type can be: Periodic, Dirichlet, Neumann, Halo
    ! Legend:
    ! -1: Halo
    !  0: Periodic
    !  1: Dirichlet
    !  2: Neumann
    type field_bc
        integer               :: type_left   !< Left boundary bc type (x = 0)
        integer               :: type_right  !< Right boundary bc type (x = Lx)
        integer               :: type_bottom !< Bottom boundary bc type (y = 0)
        integer               :: type_top    !< Top boundary bc type (y = Ly)
        real(dp), allocatable :: left(:,:)   !< Left boundary bc value
        real(dp), allocatable :: right(:,:)  !< Rright boundary bc value
        real(dp), allocatable :: bottom(:,:) !< Bottom boundary bc value
        real(dp), allocatable :: top(:,:)    !< Top boundary bc value
#if DIM==3
        integer               :: type_front  !< Front boundary bc type (z = 0)
        integer               :: type_back   !< Back boundary bc type (z = Lz)
        real(dp), allocatable :: front(:,:)  !< Front boundary bc value
        real(dp), allocatable :: back(:,:)   !< Back boundary bc value
#endif
    end type field_bc

    type scalar
        !< scalar class
        real(dp), allocatable :: f(:,:,:)       !< Array for the field values
        integer               :: gl = 0         !< Number of ghost nodes per side
        type(grid), pointer   :: G => Null()    !< Pointer to the grid where the scalar is defined
        type(field_bc)        :: bc             !< Boundary condition data structure
        character(1)          :: c = 'c'        !< scalar location in the grid cell (by default cell center)
        character(len=99)     :: name = 'unset' !< scalar name
    contains
        procedure, pass(self) :: allocate           !< allocate memory
        procedure, pass(self) :: set_from_function  !< set scalar from given function
        procedure, pass(self) :: setToValue         !< set f to a given constant value
        procedure, pass(self) :: update_ghost_nodes !< fill ghost nodes
        procedure, pass(self) :: max_value          !< compute maximum scalar value
        procedure, pass(self) :: integral           !< compute scalar integral
        procedure, pass(self) :: read               !< read scalar from file
        procedure, pass(self) :: write              !< write scalar to file
        procedure, pass(self) :: destroy            !< free the memory allocated by the scalar
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

        ! Set the field value to zero
        self%f = 0.0_dp

        ! If using ghost node
        if (self%gl > 0) then
            ! Allocate the boundary conditions array
            allocate(self%bc%left(G%lo(2)-self%gl:G%hi(2)+self%gl,G%lo(3)-self%gl:G%hi(3)+self%gl))
            allocate(self%bc%right(G%lo(2)-self%gl:G%hi(2)+self%gl,G%lo(3)-self%gl:G%hi(3)+self%gl))
            allocate(self%bc%top(G%lo(1)-self%gl:G%hi(1)+self%gl,G%lo(3)-self%gl:G%hi(3)+self%gl))
            allocate(self%bc%bottom(G%lo(1)-self%gl:G%hi(1)+self%gl,G%lo(3)-self%gl:G%hi(3)+self%gl))
#if DIM==3
            allocate(self%bc%front(G%lo(1)-self%gl:G%hi(1)+self%gl,G%lo(2)-self%gl:G%hi(2)+self%gl))
            allocate(self%bc%back(G%lo(1)-self%gl:G%hi(1)+self%gl,G%lo(2)-self%gl:G%hi(2)+self%gl))
#endif

            ! Set bc to zero by default
            self%bc%left = 0.0_dp
            self%bc%right = 0.0_dp
            self%bc%top = 0.0_dp
            self%bc%bottom = 0.0_dp
#if DIM==3
            self%bc%front = 0.0_dp
            self%bc%back = 0.0_dp
#endif

            ! Set bc type periodic by default
            self%bc%type_left = 0
            self%bc%type_right = 0
            self%bc%type_top = 0
            self%bc%type_bottom = 0
#if DIM==3
            self%bc%type_front = 0
            self%bc%type_back = 0
#endif

#ifdef MPI
            ! Search for internal boundaries
            if (self%G%prow > 1) then
                if (self%G%lo(2) /= 1) self%bc%type_bottom = -1
                if (self%G%hi(2) /= self%G%Ny) self%bc%type_top = -1 
            endif
#if DIM==3
            if (self%G%pcol > 1) then
                if (self%G%lo(3) > 1) self%bc%type_front = -1
                if (self%G%hi(3) < self%G%Nz) self%bc%type_back = -1 
            endif
#endif
#endif
        endif

    end subroutine allocate
    !===============================================================================================

    !===============================================================================================
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

    !===============================================================================================
    subroutine setToValue(self, val)

        class(scalar), intent(inout) :: self
        real(dp)     , intent(in   ) :: val

        self%f = val

    end subroutine
    !=============================================================================================== 

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

        ! TODO: For now it is assumed that parallel blocks are aligned with x direction and can be 
        !       splitted in y and z (pencil decomposition with 2decomp). This should be extended
        !       to apply other kinds of domain decomposition.

#ifdef MPI
        use mpi
        use global_mod, only : ierror
        use halo_mod  , only : update_halos

#endif
        use IO_mod    , only : print_error_message

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
        if (self%bc%type_left == 0) then ! Periodic
            do l = 1,self%gl
                self%f(lo(1)-l,:,:) = self%f(hi(1)-l+1,:,:)
            end do
        elseif (self%bc%type_left == 1) then ! Dirichlet
            if (self%c == 'c' .or. self%c == 'y' .or. self%c == 'z') then
                self%f(lo(1)-1,:,:) = 2.0_dp*self%bc%left(:,:) - self%f(lo(1),:,:)
            elseif (self%c == 'x') then 
                self%f(lo(1)-1,:,:) = self%bc%left
            else
                call print_error_message('ERROR: wrong grid location type for scalar '//trim(self%name))
            endif
        elseif (self%bc%type_left == 2) then ! Neumann
            self%f(lo(1)-1,:,:) = self%f(lo(1),:,:)
        else
            call print_error_message('ERROR: wrong left boundary condition type for scalar '//trim(self%name))
        endif

        ! Right Boundary
        if (self%bc%type_right == 0) then ! Periodic
            do l = 1,self%gl
                self%f(hi(1)+l,:,:) = self%f(lo(1)+l-1,:,:)
            end do
        elseif (self%bc%type_right == 1) then ! Dirichlet
            if (self%c == 'c' .or. self%c == 'y' .or. self%c == 'z') then
                self%f(hi(1)+1,:,:) = 2.0_dp*self%bc%right(:,:) - self%f(hi(1),:,:)
            elseif (self%c == 'x') then
                self%f(hi(1)  ,:,:) = self%bc%right(:,:)
                self%f(hi(1)+1,:,:) = self%bc%right(:,:)
            else
                call print_error_message('ERROR: wrong grid location type for scalar '//trim(self%name))
            endif
        elseif (self%bc%type_right == 2) then ! Neumann
            self%f(hi(1)+1,:,:) = self%f(hi(1),:,:)
        else
            call print_error_message('ERROR: wrong rigth boundary condition type for scalar s')
        endif

        ! Bottom Boundary
        if (self%bc%type_bottom == 0) then ! Periodic
            ! When using 2decomp decomposition this ghost nodes are updated automatically
            ! Must be overwritten only when using 1 proc
            if (self%G%prow == 1) then
                do l = 1,self%gl
                    self%f(:,lo(2)-1,:) = self%f(:,hi(2),:)
                end do
            endif
        elseif (self%bc%type_bottom == 1) then ! Dirichlet
            if (self%c == 'c' .or. self%c == 'x' .or. self%c == 'z') then
                self%f(:,lo(2)-1,:) = 2.0_dp*self%bc%bottom(:,:) - self%f(:,lo(2),:)
            elseif (self%c == 'y') then
                self%f(:,lo(2)-1,:) = self%bc%bottom
            else
                call print_error_message('ERROR: wrong grid location type for scalar '//trim(self%name))
            endif
        elseif (self%bc%type_bottom == 2) then ! Neumann
                self%f(:,lo(2)-1,:) = self%f(:,lo(2),:)
        elseif (self%bc%type_bottom == -1) then ! internal
            ! do nothing
        else    
            call print_error_message('ERROR: wrong bottom boundary condition type for scalar s')
        endif
        
        ! Top Boundary
        if (self%bc%type_top == 0) then ! Periodic
            ! When using 2decomp decomposition this ghost nodes is updated automatically
            ! Must be overwritten only when using 1 proc
            if (self%G%prow == 1) then
                do l = 1,self%gl
                    self%f(:,hi(2)+1,:) = self%f(:,lo(2),:)
                end do
            endif
        elseif (self%bc%type_top == 1) then ! Dirichlet
            if (self%c == 'c' .or. self%c == 'x' .or. self%c == 'z') then
                self%f(:,hi(2)+1,:) = 2.0_dp*self%bc%top(:,:) - self%f(:,hi(2),:)
            elseif (self%c == 'y') then
               self%f(:,hi(2)  ,:) = self%bc%top
               self%f(:,hi(2)+1,:) = self%bc%top
            else
                call print_error_message('ERROR: wrong grid location type for scalar '//trim(self%name)) 
            endif
        elseif (self%bc%type_top == 2) then ! Neumann
            self%f(:,hi(2)+1,:) = self%f(:,hi(2),:)
        elseif (self%bc%type_top == -1) then ! internal
            ! do nothing
        else
            call print_error_message('ERROR: wrong top boundary condition type for scalar s')
        endif

#if DIM==3
        if (self%bc%type_front == 0) then ! Periodic 
            ! When using 2decomp decomposition this ghost nodes is updated automatically
            ! Must be overwritten only when using 1 proc
            if (self%G%pcol == 1) then
                self%f(:,:,lo(3)-1) = self%f(:,:,hi(3))
            endif
        elseif (self%bc%type_front == 1) then ! Dirichlet
            if (self%c == 'c' .or. self%c == 'x' .or. self%c == 'y') then
                self%f(:,:,lo(3)-1) = 2.0_dp*self%bc%front(:,:) - self%f(:,:,lo(3))
            elseif (self%c == 'z') then
                   self%f(:,:,lo(3)-1) = self%bc%front
            else
                call print_error_message('ERROR: wrong grid location type for scalar '//trim(self%name)) 
            endif
        elseif (self%bc%type_front == 2) then ! Neumann
            self%f(:,:,lo(3)-1) = self%f(:,:,lo(3))
        elseif(self%bc%type_front == -1) then ! internal
            ! do nothing
        else 
            call print_error_message('ERROR: wrong front boundary condition type for scalar s')
        endif

        if (self%bc%type_back == 0) then ! Periodic 
            ! When using 2decomp decomposition this ghost nodes is updated automatically
           ! Must be overwritten only when using 1 proc
           if (self%G%pcol == 1) then
               self%f(:,:,hi(3)+1) = self%f(:,:,lo(3))
           endif
       elseif (self%bc%type_back == 1) then ! Dirichlet
           if (self%c == 'c' .or. self%c == 'x' .or. self%c == 'y') then
               self%f(:,:,hi(3)+1) = 2.0_dp*self%bc%back(:,:) - self%f(:,:,hi(3))
           elseif (self%c == 'z') then
                self%f(:,:,hi(3)  ) = self%bc%back
                self%f(:,:,hi(3)+1) = self%bc%back
           else
               call print_error_message('ERROR: wrong grid location type for scalar '//trim(self%name)) 
           endif
       elseif (self%bc%type_back == 2) then ! Neumann
           self%f(:,:,hi(3)+1) = self%f(:,:,hi(3))
       elseif(self%bc%type_back == -1) then ! internal
           ! do nothing
       else
           call print_error_message('ERROR: wrong back boundary condition type for scalar s')
       endif
#endif

#ifdef MPI
        ! This is necessary to avoid deadlock
        call mpi_barrier(mpi_comm_world, ierror)
#endif

    end subroutine
    !==============================================================================================

    !==============================================================================================
    subroutine read(self, filename)

        use decomp_2d_io, only : decomp_2d_read_one

        class(scalar)   , intent(inout) :: self     !< scalar field to read
        character(len=*), intent(in)    :: filename !< input file

        integer :: unit_id, reclen
#ifdef MPI
        call decomp_2d_read_one(1, self%f(self%G%lo(1):self%G%hi(1), &
                                          self%G%lo(2):self%G%hi(2), &
                                          self%G%lo(3):self%G%hi(3)), filename)
#else
        inquire(iolength=reclen) self%f(self%G%lo(1):self%G%hi(1), &
                                        self%G%lo(2):self%G%hi(2), &
                                        self%G%lo(2):self%G%hi(3))
        open(newunit = unit_id, file = filename, form='unformatted', status='unknown', &
                access='direct', action='read', recl=reclen)
        read (unit_id, rec=1) self%f(self%G%lo(1):self%G%hi(1), &
                                     self%G%lo(2):self%G%hi(2), &
                                     self%G%lo(2):self%G%hi(3))
        close(unit_id)
#endif

    end subroutine read
    !==============================================================================================

    !==============================================================================================
    subroutine write(self, filename)

        use decomp_2d_io, only : decomp_2d_write_one

        class(scalar)   , intent(in) :: self     !< scalar field to write
        character(len=*), intent(in) :: filename !< output file

        integer :: unit_id, reclen
#ifdef MPI
        call decomp_2d_write_one(1, self%f(self%G%lo(1):self%G%hi(1), &
                                           self%G%lo(2):self%G%hi(2), &
                                           self%G%lo(3):self%G%hi(3)), filename)
#else
        inquire(iolength=reclen) self%f(self%G%lo(1):self%G%hi(1), &
                                        self%G%lo(2):self%G%hi(2), &
                                        self%G%lo(2):self%G%hi(3))
        open(newunit = unit_id, file = filename, form='unformatted', status='unknown', &
                access='direct', action='write', recl=reclen)
        write (unit_id, rec=1) self%f(self%G%lo(1):self%G%hi(1), &
                                        self%G%lo(2):self%G%hi(2), &
                                        self%G%lo(2):self%G%hi(3))
        close(unit_id)
#endif

    end subroutine write
    !==============================================================================================

    !==============================================================================================
    subroutine destroy(self)

        class(scalar), intent(inout) :: self

        ! Free all memory allocated by the scalar s
        deallocate(self%f)

        if (self%gl > 0) deallocate(self%bc%left, self%bc%right, self%bc%top, self%bc%bottom)
#if DIM==3
        if (self%gl > 0) deallocate(self%bc%front, self%bc%back)
#endif

        self%G => Null()

    end subroutine destroy
    !==============================================================================================

end module
