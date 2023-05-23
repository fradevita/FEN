!> This module contains some predefined physical constants.
module global_mod

    use precision_mod, only : dp

#if DIM==3
    integer , parameter :: Ndim = 3                    !< For 3D simulations
#else
    integer , parameter :: Ndim = 2                    !< For 2D simulations
#endif
    integer , parameter :: tdof = merge(3, 2, Ndim==3) !< traslational degree of freedom
    integer , parameter :: rdof = merge(3, 1, Ndim==3) !< rotational degree of freedom
    integer , parameter :: dofs = tdof + rdof          !< total degree of freedom
    real(dp), parameter :: pi = acos(-1.0_dp)          !< Pi grek
    real(dp), parameter :: gravity = 9.80665_dp        !< Acceleration of gravity
    real(dp), parameter :: small = 1.0e-14_dp          !< Small positive number

    !< Staggering of variables: indexes are: (index, location)
    ! index = i,j,k, the value of the staggering in each direction
    ! location = 0 -> cell center, 1 -> x-face, 2 -> y-face, 3 -> z-face
    ! starting from cell center adding stagger(ind,location) gives the 
    ! physical coordinates of the point.
    real(dp), dimension(3,0:3), parameter :: stagger = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                                                                0.5_dp, 0.0_dp, 0.0_dp, &
                                                                0.0_dp, 0.5_dp, 0.0_dp, &
                                                                0.0_dp, 0.0_dp, 0.5_dp], shape(stagger))

    integer :: myrank = 0
    integer :: ierror

end module