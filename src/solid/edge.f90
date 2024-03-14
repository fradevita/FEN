module edge_mod

    ! This module contains the definition of edge class.

    use precision_mod, only : dp
    use marker_mod   , only : marker
    use euclidean_mod, only : distance

    implicit none
    private

    public :: edge

    type edge
        integer               :: index           !< edge index
        type(marker), pointer :: m1, m2          !< markers connected to the edge
        type(marker)          :: C               !< edge centroid
        real(dp)              :: l0 = 0.0_dp     !< Stress-free length
        real(dp)              :: l  = 0.0_dp     !< Actual length
        real(dp)              :: ks = 0.0_dp     !< Spring constant
        real(dp)              :: theta0 = 0.0_dp !< stress-free angle
        real(dp)              :: n(2)            !< normal vector to edge
        character(len=99)     :: name = 'unset'  !< edge name
    contains
        procedure, pass(self) :: update_length
        procedure, pass(self) :: update_norm
        procedure, pass(self) :: update_centroid
        procedure, pass(self) :: get_interpolated_mass_values
        procedure, pass(self) :: print
    end type edge

    interface edge
        procedure constructor
    end interface

contains

    !===============================================================================================
    function constructor(m1, m2, name)

        type(edge)                             :: constructor
        type(marker)    , intent(in), target   :: m1, m2
        character(len=*), intent(in), optional :: name

        constructor%m1 => m1
        constructor%m2 => m2
        call constructor%update_length()
        if (present(name)) constructor%name = name
        constructor%l0 = distance(m1%X, m2%X)

    end function constructor
    !===============================================================================================
    
    !===============================================================================================
    subroutine update_length(self)

        ! In/Out variables
        class(edge), intent(inout) :: self

#if DIM==3
        self%l = distance(self%m1%X, self%m2%X)
#else
        self%l = distance([self%m1%X, 0.0_dp], [self%m2%X, 0.0_dp])
#endif

    end subroutine
    !===============================================================================================
       
    !===============================================================================================
    subroutine update_norm(self)

        ! In/Out variables
        class(edge)     , intent(inout) :: self
#if DIM==3
#else       
        ! Local variables
        real(dp) :: t(2), mod_t

        ! Tangent vector
        t = self%m2%X - self%m1%X
        mod_t = sqrt(t(1)**2 + t(2)**2 + 1.0e-14_dp)
        t = t / mod_t

        ! Get norm vector
        self%n = [t(2), -t(1)]
#endif

    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine update_centroid(self)

        ! In/Out variables
        class(edge), intent(inout) :: self
      
        ! Centroid location
        self%C%X = (self%m1%X*self%m1%m + self%m2%X*self%m2%m)/(self%m1%m + self%m2%m)

    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine get_interpolated_mass_values(self)

        ! In/Out variables
        class(edge), intent(inout) :: self

        ! Interpolate values on the lagrangian marker location
        self%C%V = (self%m1%V*self%m1%m + self%m2%V*self%m2%m)/(self%m1%m + self%m2%m)
        self%C%A = (self%m1%A*self%m1%m + self%m2%A*self%m2%m)/(self%m1%m + self%m2%m)

    end subroutine get_interpolated_mass_values
    !===============================================================================================

    !===============================================================================================
    subroutine print(self)
    
        ! In/Out variables
        class(edge), intent(in) :: self

        print *, 'edge ', trim(self%name), ' points from marker', self%m1%index, &
                                                    ' to marker', self%m2%index

    end subroutine print
    !===============================================================================================

end module