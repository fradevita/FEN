!< A set of utilites.
module utils_mod

   use precision_mod, only : dp

   implicit none

contains

    !==============================================================================================
    pure real(dp) function clamp(x, xmin, xmax)
        !< Clamp variable x betwenn the minimum xmin and maximum xmax values.

        ! In/Out variables
        real(dp), intent(in) :: x    !< input value to be clamped
        real(dp), intent(in) :: xmin !< smaller acceptable value
        real(dp), intent(in) :: xmax !< larger acceptable value

        ! Local variables
        real(dp) :: temp

        temp = merge(xmin, x, x < xmin)
        clamp = merge(xmax, temp, temp > xmax)

    end function clamp
    !==============================================================================================

    !==============================================================================================
    pure function cross_product(a, b) result(c)
        !< Perform cross product between vector a and b
    
        ! In/Out variables
        real(dp), dimension(3), intent(in) :: a !< input vector a
        real(dp), dimension(3), intent(in) :: b !< input vector b
        real(dp), dimension(3)             :: c !< output vector

        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)] 
    
    end function cross_product
    !==============================================================================================

    !==============================================================================================
    function Rx(alpha)

        ! Rotation matrix around x axis

        real(dp), intent(in) :: alpha   !< input angle
        real(dp)             :: Rx(3,3) !< output matrix

        Rx(1,:) = [1.0_dp,     0.0_dp,      0.0_dp]
        Rx(2,:) = [0.0_dp, cos(alpha), -sin(alpha)]
        Rx(3,:) = [0.0_dp, sin(alpha),  cos(alpha)]
        
    end function Rx
    !==============================================================================================
    
    !==============================================================================================
    function Ry(beta)

        ! Rotation matrix around y axis

        real(dp), intent(in) :: beta    !< input angle
        real(dp)             :: Ry(3,3) !< output matrix

        Ry(1,:) = [ cos(beta), 0.0_dp, sin(beta)]
        Ry(2,:) = [    0.0_dp, 1.0_dp,    0.0_dp]
        Ry(3,:) = [-sin(beta), 0.0_dp, cos(beta)]
        
    end function Ry
    !==============================================================================================
    
    !==============================================================================================
    function Rz(gamma)

        ! Rotation matrix around z axis

        real(dp), intent(in) :: gamma   !< input angle
        real(dp)             :: Rz(3,3) !< output matrix

        Rz(1,:) = [cos(gamma), -sin(gamma), 0.0_dp]
        Rz(2,:) = [sin(gamma),  cos(gamma), 0.0_dp]
        Rz(3,:) = [    0.0_dp,      0.0_dp, 1.0_dp]
        
    end function Rz
    !==============================================================================================

    !==============================================================================================
    pure function linear_interpolation(x, x1, x2, y1, y2) result(y)
        !< Perform linear interpolation in (x,y) between (x1,y1) and (x2,y2)
    
        ! In/out variables
        real(dp), intent(in) :: x, x1, x2, y1, y2
        real(dp)             :: y

        y = y1 + (y2 - y1)*(x - x1)/(x2 - x1)

    end function
    !==============================================================================================

    !==============================================================================================
    pure function bilinear_interpolation(x, y, f11, f21, f12, f22) result(f)
        !< Perform bilinear interpolation in (x,y) in a unit square with node value fij

        ! In/out variables
        real(dp), intent(in) :: x, y, f11, f21, f12, f22
        real(dp)             :: f

        f = f11*(1.0_dp - x)*(1.0_dp - y) + f21*x*(1.0_dp - y) + f12*(1.0_dp - x)*y + f22*x*y

    end function
    !==============================================================================================

    !==============================================================================================
    function trilinear_interpolation(x,y,z,f000,f100,f010,f001,f110,f101,f011,f111) result(f)
        ! Perform trilinear interpolation in (x,y,z) in a unit cube with node value fijk

        ! In/out variables
        real(dp), intent(in) :: x, y, z, f000, f100, f010, f001, f110, f101, f011, f111
        real(dp)             :: f

        f = ((f000*(1.0_dp - x) + f100*x)*(1.0_dp - y) + &
             (f010*(1.0_dp - x) + f110*x)*y)*(1.0_dp - z) + &
            ((f001*(1.0_dp - x) + f101*x)*(1.0_dp - y) + &
             (f011*(1.0_dp - x) + f111*x)*y)*z

    end function
    !==============================================================================================

    !==============================================================================================
    function solve_quadratic(a, b, c) result(x)
        !< Solve second order equation in the form a*x**2 + b*x + c = 0

        ! In/Out variables
        real(dp), intent(in) :: a, b, c
        real(dp) :: x(2)

        ! Local variables
        real(dp) :: delta

        delta = b**2 - 4.0_dp*a*c

        x(1) = 0.5_dp*(-b - sqrt(delta))/a
        x(2) = 0.5_dp*(-b + sqrt(delta))/a

    end function solve_quadratic
    !==============================================================================================

    !==============================================================================================
    function integrate_1D(x, f, periodic) result(g)
        !< Integrate the function g(x) = \int_x f(x)\, dx using the trapeziodal rule.

        use IO_mod

        ! In/Out variables
        real(dp), intent(in) :: x(:), f(:)
        logical , intent(in) :: periodic
        real(dp)             :: g

        ! Local variables
        integer  :: i, N

        ! Check that size of x and f are the same
        if (size(x) /= size(f)) then
            call print_error_message('ERROR in integral_1D: x and f arrays must have same size')
        endif

        N = size(x)
        g = 0.0_dp
        do i = 2,N
            g = g + 0.5_dp*(f(i) + f(i-1))*(x(i) - x(i-1))
        end do

        if (periodic) g = g + 0.5_dp*(f(N) + f(1))*(x(N) - x(N-1))

    end function integrate_1D
    !==============================================================================================

    !==============================================================================================
    function urandom_seed(n, stat) result(seed)
        
        ! Routine taken from https://cyber.dabamos.de/programming/modernfortran/random-numbers.html
        ! Returns a seed array filled with random values from `/dev/urandom`.
        
        integer, intent(in)            :: n
        integer, intent(out), optional :: stat
        integer                        :: seed(n)
        integer                        :: fu, rc

        open (access='stream', action='read', file='/dev/urandom', &
              form='unformatted', iostat=rc, newunit=fu)
        if (present(stat)) stat = rc
        if (rc == 0) read (fu) seed
        close (fu)

    end function urandom_seed
    !========================================================================================

end module