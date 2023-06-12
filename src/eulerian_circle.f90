module eulerian_circle_mod

    ! Implementation for the cylinder type of solid.

    use precision_mod      , only : dp
    use eulerian_solid_mod , only : eulerian_solid

    implicit none

    type, extends(eulerian_solid) :: circle
        real(dp) :: R = 1.0_dp !< radius
        contains
            procedure, pass(self) :: distance
            procedure, pass(self) :: norm
            procedure, pass(self) :: volume
            procedure, pass(self) :: rotational_inertia
    end type circle

contains

    !=====================================================================================
    function distance(self, X) result(d)

        use precision_mod, only : dp
        use marker_mod

        ! In/Out variables
        class(circle), intent(in) :: self
        real(dp)     , intent(in) :: X(3)
        real(dp)                  :: d

        ! Local variables
        real(dp) :: dd(3)

        d = sqrt((X(1) - self%center_of_mass%X(1))**2 + (X(2) - self%center_of_mass%X(2))**2) - self%R
        if (self%G%periodic_bc(1)) then
            dd(1) = d
            dd(2) = sqrt((X(1) - self%center_of_mass%X(1) + self%G%Lx)**2 + (X(2) - self%center_of_mass%X(2))**2) - self%R
            dd(3) = sqrt((X(1) - self%center_of_mass%X(1) - self%G%Lx)**2 + (X(2) - self%center_of_mass%X(2))**2) - self%R
            d = dd(minloc(abs(dd), 1))
        endif
        if (self%G%boundary_conditions(3)%s == 'Periodic') then
            dd(1) = sqrt((X(1) - self%center_of_mass%X(1))**2 + (X(2) - self%center_of_mass%X(2) + self%G%Ly)**2) - self%R
            dd(2) = sqrt((X(1) - self%center_of_mass%X(1))**2 + (X(2) - self%center_of_mass%X(2) - self%G%Ly)**2) - self%R
            d = min(d, dd(1), dd(2))
        endif
        
    end function distance
    !=====================================================================================

    !=====================================================================================
    function norm(self, X) result(n)

        use precision_mod, only : dp
        use global_mod   , only : small
        
        class(circle), intent(in) :: self
        real(dp)     , intent(in) :: X(3)
        real(dp)                  :: n(3)
        
        ! Local variables
        real(dp) :: d(3), temp(3), nx, ny
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
            d(1) = sqrt((X(1) - self%center_of_mass%X(1))**2 + (X(2) - self%center_of_mass%X(2))**2            ) - self%R
            d(2) = sqrt((X(1) - self%center_of_mass%X(1))**2 + (X(2) - self%center_of_mass%X(2))**2 + self%G%Ly) - self%R
            d(3) = sqrt((X(1) - self%center_of_mass%X(1))**2 + (X(2) - self%center_of_mass%X(2))**2 - self%G%Ly) - self%R
            temp(1) = X(2) - self%center_of_mass%X(2)
            temp(2) = X(2) - self%center_of_mass%X(2) + self%G%Ly
            temp(3) = X(2) - self%center_of_mass%X(2) - self%G%Ly
            ny = temp(minloc(d,1))
         else
            ny = X(2) - self%center_of_mass%X(2)
         endif
         modnorm = sqrt(nx**2 + ny**2) + small
         n(1) = nx/modnorm
         n(2) = ny/modnorm
  
    end function norm
    !=====================================================================================

    !=====================================================================================
    function volume(self) result(V)

        use precision_mod, only : dp
        use global_mod   , only : pi

        class(circle), intent(in) :: self
        real(dp)                  :: V

        V = pi*self%R**2
        
    end function volume
    !=====================================================================================

    !=====================================================================================
    function rotational_inertia(self) result(I)

        use precision_mod, only : dp
        use global_mod   , only : pi

        class(circle), intent(in) :: self
        real(dp)                  :: I(3)

        I = [1.0_dp, 1.0_dp, 0.5_dp*self%mass*self%R**2]
        
    end function rotational_inertia
    !=====================================================================================

end module