!> This module contains some predefined physical constants.
module constants

  use precision, only : dp
#if DIM==3
  integer , parameter :: Ndim = 3             !< For 3D simulations
#else
  integer , parameter :: Ndim = 2             !< For 2D simulations
#endif
  real(dp), parameter :: pi = acos(-1.0_dp)   !< Pi grek
  real(dp), parameter :: gravity = 9.80665_dp !< Acceleration of gravity
  real(dp), parameter :: small = 1.0e-16_dp   !< Small positive number used to avoid divison by zero
  
  !< Staggering of variables: indexes are: (index, location)
  ! index = i,j,k, the value of the staggering in each direction
  ! location = 0 -> cell center, 1 -> x-face, 2 -> y-face, 3 -> z-face 
  real(dp), dimension(3,0:3), parameter :: stagger = reshape([0.5_dp, 0.5_dp, 0.5_dp, &
                                                              0.0_dp, 0.5_dp, 0.5_dp, &
                                                              0.5_dp, 0.0_dp, 0.5_dp, &
                                                              0.5_dp, 0.5_dp, 0.0_dp], shape(stagger))

end module constants