program main

    use precision_mod , only : dp
    use grid_mod
    use euclidean_mod , only : ZEROV, distance

    implicit none

    ! Parameters
    integer , parameter :: N = 16
    real(dp), parameter :: L = 1.0_dp

    ! Variables
    integer    :: ie(3)
    real(dp)   :: Xp(3)
    type(grid) :: G

    call G%setup(N, N, N, L, L, L, ZEROV, 1, 1)

    Xp = [G%delta*1.1, G%delta, G%delta*1.1]
    ie = G%closest_grid_node(Xp, 0)
    
    call G%destroy()

end program main