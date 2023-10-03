!> ! Define the eulerian solid objects for the Immersed Boundary Method
module test_solid_mod

    use precision_mod     , only : dp
    use eulerian_solid_mod, only : eulerian_solid

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

        class(test_solid), intent(in) :: self
        real(dp)         , intent(in) :: X(3)
        real(dp)                      :: d

        ! Local variables
        real(dp) :: temp(3)
#if DIM==3
        d = sqrt((X(1) - self%center_of_mass%X(1))**2 + (X(2) - self%center_of_mass%X(2))**2 + &
            (self%center_of_mass%X(3) - X(3))**2) - self%R
#else
        d = sqrt((X(1) - self%center_of_mass%X(1))**2 + (X(2) - self%center_of_mass%X(2))**2 ) - self%R
        if (self%G%periodic_bc(1)) then
            temp(1) = d
            temp(2) = sqrt((X(1) - self%center_of_mass%X(1) + self%G%Lx)**2 + (X(2) - self%center_of_mass%X(2))**2) - self%R
            temp(3) = sqrt((X(1) - self%center_of_mass%X(1) - self%G%Lx)**2 + (X(2) - self%center_of_mass%X(2))**2) - self%R
            d = temp(minloc(abs(temp), 1))
            d = minval(temp)
        endif
#endif

    end function distance
    !=====================================================================================

    !=====================================================================================
    function norm(self, X) result(n)

        use global_mod, only : small

        class(test_solid), intent(in) :: self
        real(dp)     , intent(in) :: X(3)
        real(dp)                  :: n(3)

        ! Local variables
        real(dp) :: d(3), temp(3), nx, ny, nz
        real(dp) :: modnorm

        if (self%G%periodic_bc(1)) then
            d(1) = sqrt((X(1) - self%center_of_mass%X(1))**2             + (X(2) - self%center_of_mass%X(2))**2) - self%R
            d(2) = sqrt((X(1) - self%center_of_mass%X(1) + self%G%Lx)**2 + (X(2) - self%center_of_mass%X(2))**2) - self%R
            d(3) = sqrt((X(1) - self%center_of_mass%X(1) - self%G%Lx)**2 + (X(2) - self%center_of_mass%X(2))**2) - self%R
            temp(1) = X(1) - self%center_of_mass%X(1)
            temp(2) = X(1) - self%center_of_mass%X(1) + self%G%Lx
            temp(3) = X(1) - self%center_of_mass%X(1) - self%G%Lx
            nx = temp(minloc(d,1))
        else
            nx = X(1) - self%center_of_mass%X(1)
        endif
        if (self%G%periodic_bc(2)) then
            d(1) = sqrt((X(1) - self%center_of_mass%X(1))**2 + (X(2) - self%center_of_mass%X(2))**2               ) - self%R
            d(2) = sqrt((X(1) - self%center_of_mass%X(1))**2 + (X(2) - self%center_of_mass%X(2))**2 + self%G%Ly) - self%R
            d(3) = sqrt((X(1) - self%center_of_mass%X(1))**2 + (X(2) - self%center_of_mass%X(2))**2 - self%G%Ly) - self%R
            temp(1) = X(2) - self%center_of_mass%X(2)
            temp(2) = X(2) - self%center_of_mass%X(2) + self%G%Ly
            temp(3) = X(2) - self%center_of_mass%X(2) - self%G%Ly
            ny = temp(minloc(d,1))
        else
            ny = X(2) - self%center_of_mass%X(2)
        endif
#if DIM==3
        nz = X(3) - self%center_of_mass%X(3)
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

        use global_mod , only : pi

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

        v = 0.0_dp

    end function volume
    !=====================================================================================

    !=====================================================================================
    function rotational_inertia(self) result(v)

        class(test_solid), intent(in) :: self
        real(dp)                      :: v(3)

        v = 0._dp

    end function rotational_inertia
    !=====================================================================================

end module
