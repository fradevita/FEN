!> This module contains definitions and procedures for the class Grid.
module grid_mod

    ! The class grid represent a computational grid over which solving governing equations.
    ! The grid class must contains information about: nodes location (logical and physical),
    ! boundary conditions (logical and physical) and parallel distribution.

    use precision_mod, only : dp
    use global_mod   , only : Ndim

    implicit none

    private
    public :: bc_type, grid

    ! This type is necessary to have an array of strings of variable size.
    ! It is used to specify physical boundary conditions on the domain.
    type bc_type
        character(:), allocatable :: s
    end type bc_type

    type grid
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
        integer :: rank = 0, nranks = 1, prow = 1, pcol = 1

        character(len=99) :: name = 'grid'

    contains
        procedure, pass(self) :: setup
        procedure, pass(self) :: closest_grid_node
        procedure, pass(self) :: print_json
        procedure, pass(self) :: destroy
    end type grid

contains

    !==============================================================================================
    subroutine setup(self, Nx, Ny, Nz, Lx, Ly, Lz, x0, prow, pcol, bc)

        ! This subroutine setup all grid variables.
        use global_mod, only : myrank
#ifdef MPI
        use decomp_2d , only : decomp_2d_init, nrank, nproc, xstart, xend, ystart, yend, zstart, zend
#endif
        use IO_mod    , only : print_error_message

        ! In/Out variables
        class(grid)  , intent(inout)           :: self
        integer      , intent(in   )           :: Nx
        integer      , intent(in   )           :: Ny
        integer      , intent(in   )           :: Nz
        integer      , intent(in   )           :: prow
        integer      , intent(in   )           :: pcol
        real(dp)     , intent(in   )           :: Lx
        real(dp)     , intent(in   )           :: Ly
        real(dp)     , intent(in   )           :: Lz
        real(dp)     , intent(in   )           :: x0(3)
        type(bc_type), intent(in   ), optional :: bc(:)

        ! Local variables
        integer :: i
        logical :: pbc(3)

        ! Set physical boundary condition on the grid
        ! By default use periodic BC on all boundarys
        self%boundary_conditions(1)%s = 'Periodic'
        self%boundary_conditions(2)%s = 'Periodic'
        self%boundary_conditions(3)%s = 'Periodic'
        self%boundary_conditions(4)%s = 'Periodic'
#if DIM==3
        self%boundary_conditions(5)%s = 'Periodic'
        self%boundary_conditions(6)%s = 'Periodic'
#endif
        if (present(bc)) self%boundary_conditions = bc
        pbc = [.false., .false., .true.]
        if (self%boundary_conditions(1)%s == 'Periodic' .and. &
            self%boundary_conditions(2)%s == 'Periodic') then
            pbc(1) = .true.
        endif
        if (self%boundary_conditions(3)%s == 'Periodic' .and. &
            self%boundary_conditions(4)%s == 'Periodic') then
            pbc(2) = .true.
        endif
#if DIM==3
        if (self%boundary_conditions(5)%s == 'Periodic' .and. &
            self%boundary_conditions(6)%s == 'Periodic') then
            pbc(3) = .true.
        else        
            pbc(3) = .false.
        endif
#endif

#ifdef MPI
        ! Create the domain decomposition with 2decomp
        call decomp_2d_init(Nx, Ny, Nz, prow, pcol, pbc)
#endif

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
            if (myrank == 0) then
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
#ifdef MPI
        self%lo = xstart
        self%hi = xend
        self%lo_y = ystart
        self%hi_y = yend
        self%lo_z = zstart
        self%hi_z = zend
#else
        self%lo = [1,1,1]
        self%hi = [Nx,Ny,Nz]
        self%lo_y = self%lo
        self%hi_y = self%hi
        self%lo_z = self%lo
        self%hi_z = self%hi
#endif
        ! Set the periodicity
        self%periodic_bc = pbc

        ! Set the parallel variables
#ifdef MPI
        self%rank = nrank
        self%nranks = nproc
        self%prow = prow
        self%pcol = pcol
#else
        self%rank = 0
        self%nranks = 1
        self%prow = 1
        self%pcol = 1
#endif

        call self%print_json

    end subroutine setup
    !==============================================================================================

    !==============================================================================================
    function closest_grid_node(self, xl, ind) result(ie)

        ! This function returns indexes ie of the closest grid point to the point xl.
        ! Ind is an integer index used to select the node location in the cell grid:
        ! 0 = cell center
        ! 1 = x face
        ! 2 = y face
        ! 3 = z face

        use global_mod, only : stagger

        ! In/Out variables
        class(grid), intent(in) :: self  !< Grid object
        integer    , intent(in) :: ind   !< location on the cell
        real(dp)   , intent(in) :: xl(3) !< given points in the grid.
        integer                 :: ie(3) !< Output indexes

        ie(1) = minloc(abs(self%x(1:self%Nx) + stagger(1, ind)*self%delta - xl(1)),1)
        ie(2) = minloc(abs(self%y(1:self%Ny) + stagger(2, ind)*self%delta - xl(2)),1)
        ie(3) = 1
#if DIM==3
        ie(3) = minloc(abs(self%z(1:self%Nz) + stagger(3, ind)*self%delta - xl(3)),1)
#endif

    end function closest_grid_node
    !==============================================================================================

    !==============================================================================================
    subroutine print_json(self)
#ifdef MPI
        use global_mod, only : myrank
#endif
        class(grid), intent(in) :: self

        integer           :: out_id
        character(len=99) :: filename

        filename = trim(self%name)//'.json'
#ifdef MPI
        if (myrank == 0) then
#endif
        open(newunit = out_id, file = filename)
        write(out_id,'(A1)') '{'
        write(out_id,'(4x,A9)') '"Grid": {'
        write(out_id,'(8x,A6,1x,I7,A1)') '"Nx": ', self%Nx, ','
        write(out_id,'(8x,A6,1x,I7,A1)') '"Ny": ', self%Ny, ','
        write(out_id,'(8x,A6,1x,I7,A1)') '"Nz": ', self%Nz, ','
        write(out_id,'(8x,A11,1x,E16.8,A1,E16.8,A1,E16.8,A2)') '"origin": [', &
           self%origin(1), ',', self%origin(2), ',', self%origin(3), '],'
        write(out_id,'(8x,A6,1x,E16.8,A1)') '"Lx": ', self%Lx, ','
        write(out_id,'(8x,A6,1x,E16.8,A1)') '"Ly": ', self%Ly, ','
        write(out_id,'(8x,A6,1x,E16.8)'   ) '"Lz": ', self%Lz
        write(out_id,'(4x,A3)') '}'
        write(out_id,'(A1)') '}'

        call flush(out_id)
#ifdef MPI
        end if
#endif

    end subroutine
    !==============================================================================================

    !==============================================================================================
    subroutine destroy(self)

        ! Free the memory allocated by the grid obj.
#ifdef MPI
        use decomp_2d, only : decomp_2d_finalize
#endif
        ! In/Out variables
        class(grid), intent(inout) :: self

        ! Free the memory allocate by the grid
        deallocate(self%x, self%y, self%z)
        deallocate(self%boundary_conditions(1)%s)
        deallocate(self%boundary_conditions(2)%s)
        deallocate(self%boundary_conditions(3)%s)
        deallocate(self%boundary_conditions(4)%s)
#if DIM==3
        deallocate(self%boundary_conditions(5)%s)
        deallocate(self%boundary_conditions(6)%s)
#endif
#ifdef MPI
        ! Free the memory of the lib 2decomp
        call decomp_2d_finalize
#endif
    end subroutine destroy
    !==============================================================================================

end module