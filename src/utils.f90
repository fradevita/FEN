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
  
end module utils
