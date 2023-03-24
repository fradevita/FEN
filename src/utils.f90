module utils

   use precision, only : dp

   implicit none

contains

    !========================================================================================
    pure real(dp) function clamp(x, xmin, xmax)

        ! In/Out variables
        real(dp), intent(in) :: x, xmin, xmax

        ! Local variables
        real(dp) :: temp

        temp = merge(xmin, x, x < xmin)
        clamp = merge(xmax, temp, temp > xmax)

    end function clamp
    !========================================================================================

    !========================================================================================
    pure function cross_product(a, b) result(c)
    
        ! In/Out variables
        real(dp), dimension(3), intent(in) :: a, b
        real(dp), dimension(3)             :: c

        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)] 
    
    end function cross_product
    !========================================================================================

    !========================================================================================
    pure function linear(x, x1, x2, y1, y2) result(y)

        ! Perform linear interpolation in (x,y) between (x1,y1) and (x2,y2)
    
        ! In/out variables
        real(dp), intent(in) :: x, x1, x2, y1, y2
        real(dp)             :: y

        y = y1 + (y2 - y1)*(x - x1)/(x2 - x1)

    end function linear
    !========================================================================================

    !========================================================================================
    pure function bilinear(x, y, f11, f21, f12, f22) result(q)

        ! Perform bilinear interpolation in (x,y) in a unit square with node value fij

        ! In/out variables
        real(dp), intent(in) :: x, y, f11, f21, f12, f22
        real(dp) :: q

        q = f11*(1.0_dp - x)*(1.0_dp - y) + f21*x*(1.0_dp - y) + f12*(1.0_dp - x)*y + f22*x*y

    end function bilinear
    !========================================================================================

    !========================================================================================
    function trilinear(x,y,z,f000,f100,f010,f001,f110,f101,f011,f111) result(q)

        ! Perform trilinear interpolation in (x,y,z) in a unit cube with node value fijk

        implicit none

        real(dp), intent(in) :: x, y, z, f000, f100, f010, f001, f110, f101, f011, f111
        real(dp) :: q

        q = ((f000*(1.0_dp - x) + f100*x)*(1.0_dp - y) + &
             (f010*(1.0_dp - x) + f110*x)*y)*(1.0_dp - z) + &
            ((f001*(1.0_dp - x) + f101*x)*(1.0_dp - y) + &
             (f011*(1.0_dp - x) + f111*x)*y)*z

    end function trilinear
    !========================================================================================

    !========================================================================================
    function solve_quadratic(a, b, c) result(x)

        ! Solve second order equation in the form a*x**2 + b*x + c = 0

        ! In/Out variables
        real(dp), intent(in) :: a, b, c
        real(dp) :: x(2)

        ! Local variables
        real(dp) :: delta

        delta = b**2 - 4.0_dp*a*c

        x(1) = 0.5_dp*(-b - sqrt(delta))/a
        x(2) = 0.5_dp*(-b + sqrt(delta))/a

    end function solve_quadratic
    !========================================================================================

    !========================================================================================
    function integrate_1D(x, f, periodic) result(g)

        ! Integrate the function g(x) = \int_x f(x)\, dx using the trapeziodal rule.

        use io

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
    !========================================================================================

    !========================================================================================
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

end module utils
