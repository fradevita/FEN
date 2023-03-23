!> This module contains definitions for functions.
module functions
    
    use precision, only : dp

    implicit none

    abstract interface
        function function_procedure(x, args) result(f)
            use precision, only : dp
            real(dp), intent(in) :: x(:)    !< location where to evaluate function f
            real(dp), intent(in) :: args(:) !< function arguments
            real(dp)             :: f       !< value of the function
        end function function_procedure
    end interface

    type function_type
        real(dp), dimension(:), allocatable :: args         !< list of function arguments
        procedure(function_procedure), pointer, nopass :: f !< pointer to function definition
    end type function_type

contains

    !==============================================================================================
    real(dp) function circle(x, args)

        ! In/Out variables
        real(dp), intent(in) :: x(2)
        real(dp), intent(in) :: args(3)

        ! Local variables
        real(dp) :: xc, yc, radius

        xc = args(1)
        yc = args(2)
        radius = args(3)

        circle = radius - sqrt( (x(1) - xc)**2 + (x(1) - yc)**2 )

    end function circle
    !==============================================================================================

end module functions
