module triangle_mod

    ! This module contains the definition of triangle class.

    use precision_mod, only : dp
    use marker_mod   , only : marker
    use edge_mod     , only : edge

    implicit none
    private

    public :: triangle

    type triangle
        type(edge)  , pointer  :: e1, e2, e3     !< Three edges of the triangle
        type(marker), pointer  :: v1, v2, v3     !< triangle vertex
        type(marker)           :: c              !< centroid of the triangle
        real(dp), dimension(3) :: n              !< Unit norm vector
        real(dp)               :: A = 0.0_dp     !< Area of the triangle
        real(dp)               :: A0 = 0.0_dp    !< Equilibrium area
        real(dp)               :: ka = 0.0_dp    !< Area constrain constant
        character(len=99)      :: name = 'unset' !< Triangle name
    contains
        procedure, pass(self) :: updateArea 
        procedure, pass(self) :: updateNorm 
        procedure, pass(self) :: updateCentroid
        procedure, pass(self) :: print
        procedure, pass(self) :: SignedVolume
    end type triangle

    interface triangle
        procedure constructor
    end interface

contains

    !===============================================================================================
    function constructor(e1, e2, e3, dof, name)

        use marker_mod, only : markersAreEqual

        type(triangle)                         :: constructor
        type(edge)      , intent(in), target   :: e1, e2, e3
        integer         , intent(in), optional :: dof
        character(len=*), intent(in), optional :: name

        ! Local variables
        integer :: csize

        ! Setup edge pointers
        constructor%e1 => e1
        constructor%e2 => e2
        constructor%e3 => e3

        ! Setup vertex pointers
        constructor%v1 => e1%m1
        constructor%v2 => e1%m2
        constructor%v3 => e2%m2
        if (markersAreEqual(constructor%v1, constructor%v3) .or. &
            markersAreEqual(constructor%v2, constructor%v3)) then
            print *, 'ERROR, two vertex equal in triangle consturctor'
            stop
        endif

        ! Initialize the triangle centroid
        if (present(dof)) then
            csize = dof
        else
            csize = 3
        endif
        constructor%c = marker(csize, -1)
        constructor%c%X = (constructor%v1%X + &
                           constructor%v2%X + &
                           constructor%v3%X)/3.0_dp

        if (present(name)) constructor%name = name
        
        ! Setup area
        call constructor%updateArea()
        constructor%A0 = constructor%A
        
        ! Setup norm
        call constructor%updateNorm()
        
    end function
    !===============================================================================================

    !===============================================================================================
    subroutine updateArea(self)

        ! In/Out variables
        class(triangle), intent(inout) :: self

        ! Local variables
        real(dp) :: p
        
        p = 0.5_dp*(self%e1%l + self%e2%l + self%e3%l)
        self%A = sqrt(p*(p - self%e1%l)*(p - self%e2%l)*(p - self%e3%l))

    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine updateNorm(self)

        use euclidean_mod, only : distanceVector, crossProduct, distance, ZEROV

        ! In/Out variables
        class(triangle), intent(inout) :: self

        ! Local variables
        real(dp), dimension(3) :: r21, r31

        r21 = distanceVector(self%e1%m1%X, self%e2%m1%X)
        r31 = distanceVector(self%e1%m1%X, self%e3%m1%X)

        self%n = crossProduct(r21, r31)
        self%n = self%n/distance(self%n, ZEROV)

    end subroutine
    !===============================================================================================

    !===============================================================================================
    subroutine updateCentroid(self)

        ! In/Out variables
        class(triangle), intent(inout) :: self

        self%c%X = (self%v1%X + self%v2%X + self%v3%X)/3.0_dp

    end subroutine
    !===============================================================================================

    !===============================================================================================
    real(dp) function SignedVolume(self)

        ! In/Out variables
        class(triangle), intent(in), target :: self

        ! Local variables
        real(dp) :: v1(3), v2(3), v3(3)
        real(dp) :: v321, v231, v312, v132, v213, v123
        
        ! Select the vertex
        v1 = self%v1%X
        v2 = self%v2%X
        v3 = self%v3%X

        v321 = v3(1)*v2(2)*v1(3)
        v231 = v2(1)*v3(2)*v1(3)
        v312 = v3(1)*v1(2)*v2(3)
        v132 = v1(1)*v3(2)*v2(3)
        v213 = v2(1)*v1(2)*v3(3)
        v123 = v1(1)*v2(2)*v3(3)
        
        SignedVolume = (-v321 + v231 + v312 - v132 - v213 + v123)/6.0_dp

    end function
    !===============================================================================================

    !===============================================================================================
    subroutine print(self)

        ! In/Out variables
        class(triangle), intent(in) :: self

        print *, 'triangle ', trim(self%name), ' has the following edges '
        call self%e1%print()
        call self%e2%print()
        call self%e3%print()

    end subroutine
    !===============================================================================================

end module