module vector_mod

    ! This module contains the vector type definition and procedures.
    ! The vector field type is a type with three scalar fields.

    use precision_mod, only : dp
    use grid_mod     , only : grid
    use scalar_mod   , only : scalar

    implicit none

    private
    public :: vector

    type vector
        type(grid  ), pointer :: G              !< grid where the vector is defined
        type(scalar)          :: x              !< x scalar component
        type(scalar)          :: y              !< y scalar component
        type(scalar)          :: z              !< z scalar component
        character(len=99)     :: name = 'unset' !< vector name
    contains
        procedure, pass(self) :: allocate         
        procedure, pass(self) :: update_ghost_nodes !<
        procedure, pass(self) :: destroy
    end type vector

contains

    !========================================================================================
    subroutine allocate(self, G, l)

        ! Procedure to allocate memory for the vector v,
        ! l is the number of ghost node per side.


        ! In/Out variables
        class(vector), intent(inout)           :: self
        type(grid)   , intent(in   ), target   :: G
        integer      , intent(in   ), optional :: l

        ! By default the field has zero ghost nodes, if l is present set the number of
        ! ghost nodes per side equal to l
        if (present(l)) then
            self%x%gl = l
            self%y%gl = l
            self%z%gl = l
        endif

        self%G => G

        ! Allocate memory
        call self%x%allocate(G, self%x%gl)
        call self%y%allocate(G, self%y%gl)
#if DIM == 3
        call self%z%allocate(G, self%z%gl)
#endif

        ! Set vector compontent location on the grid
        self%x%c = 'x'
        self%y%c = 'y'
#if DIM==3
        self%z%c = 'z'
#endif
    end subroutine allocate
    !==============================================================================================

    !==============================================================================================
     subroutine update_ghost_nodes(self)

        ! Apply boundary condition on scalar components of the vector.
       
        ! In/Out variables
        class(vector), intent(inout) :: self !< input vector field

        call self%x%update_ghost_nodes()
        call self%y%update_ghost_nodes()
#if DIM==3
        call self%z%update_ghost_nodes()
#endif

    end subroutine
    !==============================================================================================
    
    !==============================================================================================
    subroutine destroy(self)

        ! Free all memory allocated by vector field v.

        ! In/Out variables
        class(vector), intent(inout) :: self !< input vector field

        call self%x%destroy()
        call self%y%destroy()
#if DIM == 3
        call self%z%destroy()
#endif

        self%G => Null()

    end subroutine destroy
    !==============================================================================================

end module
