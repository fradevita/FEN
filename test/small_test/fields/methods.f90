program methods
#ifdef MPI
    use mpi
    use global_mod    , only : ierror
#endif
    use precision_mod , only : dp
    use grid_mod      , only : grid
    use scalar_mod    , only : scalar
    use function_mod  , only : function_type, circle

    implicit none

    type(grid  )        :: G
    type(scalar)        :: s
    type(function_type) :: my_func 

#ifdef MPI
    call mpi_init(ierror)
#endif

    ! Setup the grid
    call G%setup(50, 50, 1, 2.0_dp, 2.0_dp, 2.0_dp/50.0_dp, [-0.5_dp, -0.5_dp, 0.0_dp], 1, 1)

    ! Setup the function
    my_func%args = [0.0_dp, 0.0_dp, 1.0_dp]
    my_func%fp => circle

    ! Create the scalar
    call s%allocate(G)

    ! Setup from the function
    call s%set_from_function(my_func)

    ! Write the scalar on file
    call s%write('s.raw')

#ifdef MPI
    call mpi_finalize(ierror)
#endif

end program