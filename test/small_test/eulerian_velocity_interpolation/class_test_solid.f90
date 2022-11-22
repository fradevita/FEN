!> ! Define the eulerian solid objects for the Immersed Boundary Method
module class_test_solid

    use precision           , only : dp
    use class_eulerian_solid, only : eulerian_solid

    implicit none

    ! Below, some usable solids are defined.
    type, extends(eulerian_solid) :: test_solid
        real(dp) :: X(3) = 0.0_dp
        real(dp) :: R = 1.0_dp
        contains
            procedure, pass(self) :: distance
            procedure, pass(self) :: norm
            procedure, pass(self) :: velocity
    end type test_solid

contains

    !=====================================================================================
    function distance(self, X) result(d)

        use precision , only : dp
        use class_Grid, only : base_grid

        class(test_solid), intent(in) :: self
        real(dp)     , intent(in) :: X(3)
        real(dp)                  :: d

        ! Local variables
        real(dp) :: dxp, dxm, dyp, dym

#if DIM==3
        d = sqrt((self%X(1) - X(1))**2 + (self%X(2) - X(2))**2 + (self%X(3) - X(3))**2) - self%R
#else
        d = sqrt((self%X(1) - X(1))**2 + (self%X(2) - X(2))**2 ) - self%R
#endif

    end function distance
    !=====================================================================================

    !=====================================================================================
    function norm(self, X) result(n)

        use precision, only : dp
        use constants, only : small

        class(test_solid), intent(in) :: self
        real(dp)     , intent(in) :: X(3)
        real(dp)                  :: n(3)

        ! Local variables
        integer :: i
        real(dp) :: modn

        do i = 1,2
            n(i) = X(i) - self%X(i)
        end do
#if DIM==3
        n(3) = X(3) - self%X(3)
#else
        n(3) = 0.0_dp
#endif     
        modn = sqrt(n(1)**2 + n(2)**2 + n(3)**2) + small
        n = n/modn

    end function norm
    !=====================================================================================

    !=====================================================================================
    function velocity(self, X, dir) result(v)

        use precision, only : dp
        use constants, only : pi

        class(test_solid), intent(in) :: self
        real(dp)     , intent(in) :: X(3)
        integer      , intent(in) :: dir
        real(dp)                  :: v

#ifdef DEBUG
        ! Check that the location is on the solid boundary
        if (abs(self%distance(X)) > 1.0e-12_dp) then
            print *, 'ERROR, the point does not lies on the solid boundary.'
        endif
#endif
        if (dir == 1) then
            v = -cos(2*pi*X(1))*sin(2*pi*X(2))
        elseif (dir == 2) then
            v =  sin(2*pi*X(1))*cos(2*pi*X(2))
        elseif (dir == 3) then
            v = cos(2*pi*X(1))*sin(2*pi*X(3))
        endif

    end function velocity
    !=====================================================================================

end module class_test_solid
