!> ! Define the eulerian solid objects for the Immersed Boundary Method
module class_test_solid

    use precision           , only : dp
    use class_eulerian_solid, only : eulerian_solid

    implicit none

    type, extends(eulerian_solid) :: test_solid
        real(dp) :: R = 1.0_dp
        contains
            procedure, pass(self) :: distance
            procedure, pass(self) :: norm
            procedure, pass(self) :: velocity
            procedure, pass(self) :: volume 
            procedure, pass(self) :: rotational_inertia
            end type test_solid

contains

    !=====================================================================================
    function distance(self, X) result(d)

        use class_Grid
        use precision , only : dp

        class(test_solid), intent(in) :: self
        real(dp)         , intent(in) :: X(3)
        real(dp)                      :: d

        ! Local variables
        real(dp) :: temp(3)
#if DIM==3
        d = sqrt((X(1) - self%X(1))**2 + (X(2) - self%X(2))**2 + (self%X(3) - X(3))**2) - self%R
#else
        d = sqrt((X(1) - self%X(1))**2 + (X(2) - self%X(2))**2 ) - self%R
        if (base_grid%periodic_bc(1)) then
            temp(1) = d
            temp(2) = sqrt((X(1) - self%X(1) + base_grid%Lx)**2 + (X(2) - self%X(2))**2) - self%R
            temp(3) = sqrt((X(1) - self%X(1) - base_grid%Lx)**2 + (X(2) - self%X(2))**2) - self%R
            d = temp(minloc(abs(temp), 1))
            d = minval(temp)
        endif
#endif

    end function distance
    !=====================================================================================

    !=====================================================================================
    function norm(self, X) result(n)

        use class_Grid
        use precision, only : dp
        use constants, only : small

        class(test_solid), intent(in) :: self
        real(dp)     , intent(in) :: X(3)
        real(dp)                  :: n(3)

        ! Local variables
        real(dp) :: d(3), temp(3), nx, ny, nz
        real(dp) :: modnorm

        if (base_grid%periodic_bc(1)) then
            d(1) = sqrt((X(1) - self%X(1))**2                + (X(2) - self%X(2))**2) - self%R
            d(2) = sqrt((X(1) - self%X(1) + base_grid%Lx)**2 + (X(2) - self%X(2))**2) - self%R
            d(3) = sqrt((X(1) - self%X(1) - base_grid%Lx)**2 + (X(2) - self%X(2))**2) - self%R
            temp(1) = X(1) - self%X(1)
            temp(2) = X(1) - self%X(1) + base_grid%Lx
            temp(3) = X(1) - self%X(1) - base_grid%Lx
            nx = temp(minloc(d,1))
        else
            nx = X(1) - self%X(1)
        endif
        if (base_grid%periodic_bc(2)) then
            d(1) = sqrt((X(1) - self%X(1))**2 + (X(2) - self%X(2))**2               ) - self%R
            d(2) = sqrt((X(1) - self%X(1))**2 + (X(2) - self%X(2))**2 + base_grid%Ly) - self%R
            d(3) = sqrt((X(1) - self%X(1))**2 + (X(2) - self%X(2))**2 - base_grid%Ly) - self%R
            temp(1) = X(2) - self%X(2)
            temp(2) = X(2) - self%X(2) + base_grid%Ly
            temp(3) = X(2) - self%X(2) - base_grid%Ly
            ny = temp(minloc(d,1))
        else
            ny = X(2) - self%X(2)
        endif
#if DIM==3
        nz = X(3) - self%X(3)
#else
        nz = 0.0_dp
#endif
        modnorm = sqrt(nx**2 + ny**2 + nz**2) + small
        n(1) = nx/modnorm
        n(2) = ny/modnorm
        n(3) = nz/modnorm
  
    end function norm
    !=====================================================================================

    !=====================================================================================
    function velocity(self, X, dir) result(v)

        ! Override velocity function to force velocity on the solid boundary.

        use precision, only : dp
        use constants, only : pi

        class(test_solid), intent(in) :: self
        real(dp)     , intent(in)     :: X(3)
        integer      , intent(in)     :: dir
        real(dp)                      :: v

        if (dir == 1) then
            v = -cos(2*pi*X(1))*sin(2*pi*X(2))
        elseif (dir == 2) then
            v =  sin(2*pi*X(1))*cos(2*pi*X(2))
        elseif (dir == 3) then
            v = cos(2*pi*X(1))*sin(2*pi*X(3))
        endif

    end function velocity
    !=====================================================================================

    !=====================================================================================
    function volume(self) result(v)

        class(test_solid), intent(in) :: self
        real(dp)                      :: v

    end function volume
    !=====================================================================================

    !=====================================================================================
    function rotational_inertia(self) result(v)

        class(test_solid), intent(in) :: self
        real(dp)                      :: v(3)

    end function rotational_inertia
    !=====================================================================================

end module class_test_solid
