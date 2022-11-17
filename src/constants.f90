!> This module contains some predefined physical constants.
module constants

  use precision, only : dp
#ifdef DIM==3
  integer , parameter :: Ndim = 3
#else
  integer , parameter :: Ndim = 2 
#endif
  real(dp), parameter :: pi = acos(-1.0_dp)   !< Pi grek
  real(dp), parameter :: gravity = 9.80665_dp !< Acceleration of gravity
  real(dp), parameter :: small = 1.0e-16_dp   !< Small positive number used to avoid divison by zero
  
end module constants