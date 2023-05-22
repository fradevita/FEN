module tensor_mod

    ! This module contains the tensor type definition and procedures.
    ! Thetensor type is a type with three vectors.

    use precision_mod, only : dp
    use vector_mod
    use grid_mod

    implicit none

    type tensor
        type(grid  ), pointer :: G              !< grid where the tensor is defined
        type(vector)          :: x              !< x vector component
        type(vector)          :: y              !< y vector component
        type(vector)          :: z              !< z vector component
        integer               :: gl             !< ghost node levels
        character(len=99)     :: name = 'unset' !< tensor name
    contains
        procedure, pass(self) :: allocate
        procedure, pass(self) :: update_ghost_nodes
        procedure, pass(self) :: symmetric
        procedure, pass(self) :: destroy     
    end type tensor

contains

    !========================================================================================
    subroutine allocate(self, G, l)

        ! In/Out variables
        class(tensor), intent(inout)           :: self !< input tensor
        type(grid)   , intent(in   ), target   :: G    !< computational grid
        integer      , intent(in   ), optional :: l    !< number of ghost nodes

        ! By default the field has zero ghost nodes, if l is present set the number of
        ! ghost nodes per side equal to l
        if (present(l)) then
            self%gl = l
            self%x%x%gl = l
            self%x%y%gl = l
            self%x%z%gl = l
            self%y%x%gl = l
            self%y%y%gl = l
            self%y%z%gl = l
            self%z%x%gl = l
            self%z%y%gl = l
            self%z%z%gl = l
        endif

        self%G => G

        call self%x%allocate(G, self%gl)
        call self%y%allocate(G, self%gl)
#if DIM == 3
        call self%z%allocate(G, self%gl)
#endif

    end subroutine allocate
    !========================================================================================

    !========================================================================================
    subroutine update_ghost_nodes(self)

        ! In/Out variables
        class(tensor), intent(inout) :: self

        call self%x%update_ghost_nodes
        call self%x%update_ghost_nodes
#if DIM==3
        call self%x%update_ghost_nodes
#endif

    end subroutine 
    !========================================================================================

    !========================================================================================
    subroutine symmetric(self)

        ! Replace T with its symmetric tensor

        ! In/Out variables
        class(tensor), intent(inout) :: self

        ! Local variables
        type(tensor) :: temp

        ! Allocate memory
        call temp%allocate(self%G, self%gl)

        ! Save T
        temp = self

        ! Replace with its symmetric tensors
        self%x%y%f = 0.5_dp*(temp%x%y%f + temp%y%x%f)
        self%y%x%f = self%x%y%f
#if DIM==3
        self%x%z%f = 0.5_dp*(temp%x%z%f + temp%z%x%f)
        self%z%x%f = T%x%z%f
        self%y%z%f = 0.5_dp*(temp%y%z%f + temp%z%y%f)
        self%z%y%f = self%y%z%f
#endif

        call temp%destroy
        if (self%gl > 0) call self%update_ghost_nodes()
    
    end subroutine symmetric
    !========================================================================================
    
    !========================================================================================
    subroutine destroy(self)

        class(tensor), intent(inout) :: self

        ! Free all memory allocated by vector field v
        call self%x%destroy()
        call self%y%destroy()
#if DIM == 3
        call self%z%destroy()
#endif

    end subroutine destroy
    !========================================================================================

end module