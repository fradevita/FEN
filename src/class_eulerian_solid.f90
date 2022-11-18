!> ! Define the eulerian solid objects for the Immersed Boundary Method
module class_eulerian_solid

    use precision, only : dp

    implicit none
    private
    public :: eulerian_solid
    public :: eulerian_solid_pointer

    !> Base object for eulerian IBM. Must define only a distance function,
    !> a norm function and a velocity function for the Dirichlet BC.
    type, abstract :: eulerian_solid
        contains
            procedure(distance_interface), pass(self), deferred :: distance
            procedure(norm_interface    ), pass(self), deferred :: norm
            procedure(velocity_interface), pass(self), deferred :: velocity
    end type eulerian_solid

    abstract interface
        !< Abastract interface for eulerian solid object
        function distance_interface(self, x) result(f)
            !< Function to compute distance from solid surface at given location x.
            use precision, only : dp
            import eulerian_solid
            class(eulerian_solid), intent(in) :: self    !< eulerian solid object
            real(dp)             , intent(in) :: x(3)    !< location where to compute distance
            real(dp)                          :: f       !< distance result
        end function distance_interface

        function velocity_interface(self, x, dir) result(v)
            !< Function to compute velocity from rigid body motion at given location x.
            use precision, only : dp
            import :: eulerian_solid
            class(eulerian_solid), intent(in) :: self !< eulerian solid object
            real(dp)             , intent(in) :: x(3) !< location where to compute distance
            integer              , intent(in) :: dir  !< velocity component
            real(dp)                          :: v    !< velocity value
        end function velocity_interface

        function norm_interface(self, x) result(n)
            !< Function to compute local norm with respect to solid surface at given location x.
            use precision, only : dp
            import :: eulerian_solid
            class(eulerian_solid), intent(in) :: self !< eulerian solid object
            real(dp)             , intent(in) :: x(3) !< location where to compute normal vector
            real(dp)                          :: n(3) !< local norm vector
        end function norm_interface
    end interface

    type eulerian_solid_pointer
        class(eulerian_solid), pointer :: pS
    end type eulerian_solid_pointer

end module class_eulerian_solid

module class_eulerian_circle

    use precision           , only : dp
    use class_eulerian_solid, only : eulerian_solid

    ! Below, some usable solids are defined.
    type, extends(eulerian_solid) :: circle
        real(dp) :: X(3) = 0.0_dp
        real(dp) :: R = 1.0_dp
        contains
            procedure, pass(self) :: distance
            procedure, pass(self) :: norm
            procedure, pass(self) :: velocity
    end type circle

contains

!=====================================================================================
    function distance(self, X) result(d)

        use precision, only : dp

        class(circle), intent(in) :: self
        real(dp)     , intent(in) :: X(3)
        real(dp)                  :: d

        d = sqrt((self%X(1) - X(1))**2 + (self%X(2) - X(2))**2 + (self%X(3) - X(3))**2) - &
            self%R

    end function distance
    !=====================================================================================

    !=====================================================================================
    function norm(self, X) result(n)

        use precision, only : dp

        class(circle), intent(in) :: self
        real(dp)     , intent(in) :: X(3)
        real(dp)                  :: n(3)

        n = 0.0_dp

    end function norm
    !=====================================================================================

    !=====================================================================================
    function velocity(self, X, dir) result(v)

        use precision, only : dp

        class(circle), intent(in) :: self
        real(dp)     , intent(in) :: X(3)
        integer      , intent(in) :: dir
        real(dp)                  :: v

#ifdef DEBUG
        ! Check that the location is on the solid boundary
        if (abs(self%distance(X)) > 1.0e-12_dp) then
            print *, 'ERROR, the point does not lies on the solid boundary.'
        endif
#endif

        v = 0.0_dp

    end function velocity
    !=====================================================================================

end module class_eulerian_circle