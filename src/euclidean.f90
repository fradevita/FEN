module euclidean_mod

    ! This module contains usefull functions to work with cartesian arrays. 

    use precision_mod, only : dp

    implicit none

    ! Null vector
    real(dp), dimension(3), parameter :: ZEROV = [0.0_dp, 0.0_dp, 0.0_dp]

    ! Cartesian directions
    real(dp), parameter :: ex(3) = [1.0_dp, 0.0_dp, 0.0_dp]
    real(dp), parameter :: ey(3) = [0.0_dp, 1.0_dp, 0.0_dp]
    real(dp), parameter :: ez(3) = [0.0_dp, 0.0_dp, 1.0_dp]

contains

    !===============================================================================================
    function Rx(alpha)

        ! Rotation matrix around x axis

        ! In/Out variables
        real(dp), intent(in) :: alpha   !< input angle
        real(dp)             :: Rx(3,3) !< output matrix

        Rx(1,:) = [1.0_dp,     0.0_dp,      0.0_dp]
        Rx(2,:) = [0.0_dp, cos(alpha), -sin(alpha)]
        Rx(3,:) = [0.0_dp, sin(alpha),  cos(alpha)]
        
    end function Rx
    !===============================================================================================

    !===============================================================================================
    function Ry(beta)

        ! Rotation matrix around y axis

        ! In/Out variables
        real(dp), intent(in) :: beta    !< input angle
        real(dp)             :: Ry(3,3) !< output matrix

        Ry(1,:) = [ cos(beta), 0.0_dp, sin(beta)]
        Ry(2,:) = [    0.0_dp, 1.0_dp,    0.0_dp]
        Ry(3,:) = [-sin(beta), 0.0_dp, cos(beta)]
        
    end function Ry
    !===============================================================================================

    !===============================================================================================
    function Rz(gamma)

        ! Rotation matrix around z axis

        ! In/Out variables
        real(dp), intent(in) :: gamma   !< input angle
        real(dp)             :: Rz(3,3) !< output matrix

        Rz(1,:) = [cos(gamma), -sin(gamma), 0.0_dp]
        Rz(2,:) = [sin(gamma),  cos(gamma), 0.0_dp]
        Rz(3,:) = [    0.0_dp,      0.0_dp, 1.0_dp]
        
    end function Rz
    !===============================================================================================
    
    !===============================================================================================
    real(dp) function dotProduct(a, b)

        ! compute dot product between array a and b.

        ! In/Out variables
        real(dp), dimension(3), intent(in) :: a, b

        dotProduct = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

    end function 
    !===============================================================================================

    !===============================================================================================
    function crossProduct(a, b) result(c)

        ! compute cross product between array a and b.

        ! In/Out variables
        real(dp), dimension(3), intent(in) :: a, b
        real(dp), dimension(3)             :: c

        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)

    end function 
    !===============================================================================================

    !===============================================================================================
    real(dp) function distance(a, b)

        ! compute distance between position vector a and b.

        ! In/Out variables
        real(dp), dimension(3), intent(in) :: a, b

        distance = sqrt((a(1) - b(1))**2 + (a(2) - b(2))**2 + (a(3) - b(3))**2)

    end function distance
    !===============================================================================================

    !===============================================================================================
    function distanceVector(a, b) result(c)

        ! Compute vector pointing from position vector a to b

        ! In/Out variables
        real(dp), dimension(3), intent(in) :: a, b
        real(dp), dimension(3)             :: c

        c(1) = b(1) - a(1)
        c(2) = b(2) - a(2)
        c(3) = b(3) - a(3)

    end function distanceVector
    !===============================================================================================

    !===============================================================================================
    function tilda_matrix(a)

        ! In/Out variables
        real(dp), intent(in)     :: a(3)
        real(dp), dimension(3,3) :: tilda_matrix

        tilda_matrix = 0.0_dp
        tilda_matrix(1,2) =  a(3)
        tilda_matrix(1,3) = -a(2)
        tilda_matrix(2,1) = -a(3)
        tilda_matrix(2,3) =  a(1)
        tilda_matrix(3,1) =  a(2)
        tilda_matrix(3,2) = -a(1)

    end function tilda_matrix
    !===============================================================================================
    
end module euclidean_mod
