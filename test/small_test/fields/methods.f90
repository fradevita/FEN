program methods

#ifdef MPI
    use mpi
    use global_mod    , only : ierror
#endif
    use global_mod    , only : small
    use precision_mod , only : dp
    use grid_mod      , only : grid
    use scalar_mod    , only : scalar
    use IO_mod

    implicit none

    integer      :: i, j, k
    real(dp)     :: x, y, z, tol = small
    type(grid  ) :: G
    type(scalar) :: s

#ifdef MPI
    call mpi_init(ierror)
#endif

    ! Setup the grid
    call G%setup(64, 64, 64, 1.0_dp, 1.0_dp, 1.0_dp, [-0.5_dp, -0.5_dp, 0.0_dp], 1, 8)

    ! Create the scalar
    call s%allocate(G, 1)

     ! Initialize the scalar field s
    do k = G%lo(3),G%hi(3)
        z = (float(k) - 0.5_dp)*G%delta
        do j = G%lo(2),G%hi(2)
           y = (float(j) - 0.5_dp)*G%delta
           do i = G%lo(1),G%hi(1)
              x = (float(i) - 0.5_dp)*G%delta
              s%f(i,j,k) = sample_function(x, y, z)
           end do
        end do
    end do
    call s%update_ghost_nodes()
   
    ! Check that the ghost nodes have been correctly filled
    call check_ghost_nodes(s)

#ifdef MPI
    call mpi_finalize(ierror)
#endif

contains

    !==============================================================================================
    pure real(dp) function sample_function(x, y, z)

        use global_mod , only : pi

        real(dp), intent(in) :: x, y, z

        sample_function = sin(2.0_dp*pi*x)*cos(2.0_dp*pi*y)*sin(2.0_dp*pi*z)

    end function sample_function
    !==============================================================================================

    !==============================================================================================
    subroutine check_ghost_nodes(s)

        real(dp)                 :: e
        type(scalar), intent(in) :: s

        ! Check that the ghost nodes have been correctly filled
        ! ** Left boundary **
        i = G%hi(1) - 1
        x = (float(i) - 0.5_dp)*G%delta
        do k = G%lo(3),G%hi(3)
        z = (float(k) - 0.5_dp)*G%delta
        do j = G%lo(2),G%hi(2)
            y = (float(j) - 0.5_dp)*G%delta
            e = abs(s%f(i,j,k) - sample_function(x,y,z))
            if (e > tol) then
                call print_error_message('Error greather than tollerance on left boundary')
            endif
        end do
        end do

        ! ** Right boundary **
        i = G%hi(1) + 1
        x = (float(i) - 0.5_dp)*G%delta
        do k = G%lo(3),G%hi(3)
        z = (float(k) - 0.5_dp)*G%delta
        do j = G%lo(2),G%hi(2)
            y = (float(j) - 0.5_dp)*G%delta
            e = abs(s%f(i,j,k) - sample_function(x,y,z))
            if (e > tol) then
                print *, i, j, k, e
                call print_error_message('Error greather than tollerance on right boundary')
            endif
        end do
        end do

        ! ** Top boundary **
        j = G%hi(2) + 1
        y = (float(j) - 0.5_dp)*G%delta
        do k = G%lo(3),G%hi(3)
        z = (float(k) - 0.5_dp)*G%delta
        do i = G%lo(1),G%hi(1)
            x = (float(i) - 0.5_dp)*G%delta
            e = abs(s%f(i,j,k) - sample_function(x,y,z))
            if (e > tol) then
                call print_error_message('Error greather than tollerance on top boundary')
            endif
        end do
        end do

        ! ** Bottom boundary **
        j = G%lo(2) - 1
        y = (float(j) - 0.5_dp)*G%delta
        do k = G%lo(3),G%hi(3)
        z = (float(k) - 0.5_dp)*G%delta 
        do i = G%lo(1),G%hi(1)
            x = (float(i) - 0.5_dp)*G%delta
            e = abs(s%f(i,j,k) - sample_function(x,y,z))
            if (e > tol) then
                call print_error_message('Error greather than tollerance on bottom boundary')
            endif
        end do
        end do

        ! ** Front boundary **
        k = G%hi(3) + 1
        z = (float(k) - 0.5_dp)*G%delta
        do j = G%lo(2),G%hi(2)
        y = (float(j) - 0.5_dp)*G%delta
        do i = G%lo(1),G%hi(1)
            x = (float(i) - 0.5_dp)*G%delta
            e = abs(s%f(i,j,k) - sample_function(x,y,z))
            if (e > tol) then
                call print_error_message('Error greather than tollerance on front boundary')
            endif
        end do
        end do

        ! ** Back boundary **
        k = G%lo(3) - 1
        z = (float(k) - 0.5_dp)*G%delta
        do j = G%lo(2),G%hi(2)
        y = (float(j) - 0.5_dp)*G%delta
        do i = G%lo(1),G%hi(1)
            x = (float(i) - 0.5_dp)*G%delta
            e = abs(s%f(i,j,k) - sample_function(x,y,z))
            if (e > tol) then
                call print_error_message('Error greather than tollerance on back boundary')
            endif
        end do
        end do
        
    end subroutine check_ghost_nodes
    !==============================================================================================

end program