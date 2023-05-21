program memory

    use precision_mod
    use grid_mod
    use scalar_mod

    implicit none

    type(grid)   :: G
    type(scalar) :: s1, s2

    ! Before define a scalar we must define a grid
    call G%setup(10, 10, 10, 1.0_dp, 1.0_dp, 1.0_dp, [0.0_dp, 0.0_dp, 0.0_dp], 1, 1)

    ! Create scalar s1 without ghost node
    call s1%allocate(G)
    
    ! Create scalar s2 with 1 level of ghost node
    call s2%allocate(G, 1)
    
    ! Free memory
    call s1%destroy()
    call s2%destroy()
    call G%destroy()

end program