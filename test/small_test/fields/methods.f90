program methods

#ifdef MPI
    use mpi
    use global_mod    , only : ierror, myrank
#endif
    use global_mod    , only : small, pi
    use precision_mod , only : dp
    use grid_mod      , only : grid
    use scalar_mod    , only : scalar
    use vector_mod    , only : vector
    use fields_mod    , only : gradient, laplacian
    use IO_mod

    implicit none

    ! Parameters
    real(dp), parameter :: L = 1.0_dp

    ! Variables
    integer         :: r, N, i, j, k, prow = 1, pcol = 8, out_id
    real(dp)        :: L_inf_x, L_inf_y, L_inf_z, L_inf_2
    real(dp)        :: L_2_x, L_2_y, L_2_z, L_2_2
    type(grid  )    :: G
    type(scalar)    :: s
    type(vector)    :: error

#ifdef MPI
    call mpi_init(ierror)
    call mpi_comm_rank(mpi_comm_world, myrank, ierror)
#endif

    open(newunit = out_id, file = 'out.csv')
    if (myrank == 0) write(out_id,'(A)') 'N,Linfx,Linfy,Linfz,L2x,L2y,L2z,Linf2,L22'

    resolution_cycle: do r = 3,6
        N = 2**r

        ! Setup the grid
        call G%setup(N, N, N, L, L, L, [-L/2.0_dp, -L/2.0_dp, -L/2.0_dp], prow, pcol)

        ! Create the scalar
        call s%allocate(G, 1)

        ! Initialize the scalar field s
        do k = G%lo(3),G%hi(3)
            do j = G%lo(2),G%hi(2)
                do i = G%lo(1),G%hi(1)
                    s%f(i,j,k) = f(G%x(i), G%y(j), G%z(k))
                end do
            end do
        end do
        call s%update_ghost_nodes()

        ! Check that the ghost nodes have been correctly filled
        call check_ghost_nodes()

        ! Create the vector for error
        call error%allocate(G)
   
        ! Evaluate gradien of scalar s
        call check_gradient(error)
        L_inf_x = error%x%max_value()
        L_inf_y = error%y%max_value()
        L_inf_z = error%z%max_value()
        L_2_x = 0.0_dp
        L_2_y = 0.0_dp
        L_2_z = 0.0_dp
        do k = G%lo(3),G%hi(3)
            do j = G%lo(2),G%hi(2)
                do i = G%lo(1),G%hi(1)
                    L_2_x = L_2_x + error%x%f(i,j,k)**2
                    L_2_y = L_2_y + error%y%f(i,j,k)**2
                    L_2_z = L_2_z + error%z%f(i,j,k)**2
                end do
            end do
        end do
        L_2_x = sqrt(L_2_x/reaL(N, dp))
        L_2_y = sqrt(L_2_y/reaL(N, dp))
        L_2_z = sqrt(L_2_z/reaL(N, dp))

        ! Evaluate laplacian of s
        call check_laplacian(error%x)
        L_inf_2 = error%x%max_value()
        L_2_2 = 0.0_dp
        do k = G%lo(3),G%hi(3)
            do j = G%lo(2),G%hi(2)
                do i = G%lo(1),G%hi(1)
                    L_2_2 = L_2_2 + error%x%f(i,j,k)**2
                end do
            end do
        end do
        L_2_2 = sqrt(L_2_2/reaL(N, dp))
        
        if (myrank == 0) then
            write(out_id,'(I03,",",*(E16.8,:,","))') N, L_inf_x, L_inf_y, L_inf_z, L_2_x, L_2_y, &
                                                     L_2_z, L_inf_2, L_2_2
        endif

        ! Free memory
        call s%destroy()
        call error%destroy()
        call G%destroy()

    end do resolution_cycle
#ifdef MPI
    call mpi_finalize(ierror)
#endif

contains

    !===============================================================================================
    pure real(dp) function f(x, y, z)

        ! In/Out variables
        real(dp), intent(in) :: x, y, z

        f = sin(2.0_dp*pi*x)*cos(2.0_dp*pi*y)*sin(2.0_dp*pi*z)

    end function f
    !===============================================================================================

    !===============================================================================================
    pure real(dp) function dfdx(x, y, z)

        ! In/Out variables
        real(dp), intent(in) :: x, y, z

        dfdx = 2.0_dp*pi*cos(2.0_dp*pi*x)*cos(2.0_dp*pi*y)*sin(2.0_dp*pi*z)

    end function dfdx
    !===============================================================================================

    !===============================================================================================
    pure real(dp) function dfdy(x, y, z)

        ! In/Out variables
        real(dp), intent(in) :: x, y, z

        dfdy = -2.0_dp*pi*sin(2.0_dp*pi*x)*sin(2.0_dp*pi*y)*sin(2.0_dp*pi*z)

    end function dfdy
    !===============================================================================================

    !===============================================================================================
    pure real(dp) function dfdz(x, y, z)

        ! In/Out variables
        real(dp), intent(in) :: x, y, z

        dfdz = 2.0_dp*pi*sin(2.0_dp*pi*x)*cos(2.0_dp*pi*y)*cos(2.0_dp*pi*z)

    end function dfdz
    !===============================================================================================

    !===============================================================================================
    pure real(dp) function d2fdx2(x, y, z)

        ! In/Out variables
        real(dp), intent(in) :: x, y, z

        d2fdx2 = -4.0_dp*pi*pi*sin(2.0_dp*pi*x)*cos(2.0_dp*pi*y)*sin(2.0_dp*pi*z)

    end function d2fdx2
    !===============================================================================================

    !===============================================================================================
    pure real(dp) function d2fdy2(x, y, z)

        ! In/Out variables
        real(dp), intent(in) :: x, y, z

        d2fdy2 = -4.0_dp*pi*pi*sin(2.0_dp*pi*x)*cos(2.0_dp*pi*y)*sin(2.0_dp*pi*z)

    end function d2fdy2
    !===============================================================================================

    !===============================================================================================
    pure real(dp) function d2fdz2(x, y, z)

        ! In/Out variables
        real(dp), intent(in) :: x, y, z

        d2fdz2 = -4.0_dp*pi*pi*sin(2.0_dp*pi*x)*cos(2.0_dp*pi*y)*sin(2.0_dp*pi*z)

    end function d2fdz2
    !===============================================================================================

    !===============================================================================================
    subroutine check_ghost_nodes()

        ! Local variables
        real(dp) :: e

        ! Check that the ghost nodes have been correctly filled
        ! **** Left boundary ***********************************************************************
        i = G%lo(1) - 1
        do k = G%lo(3),G%hi(3)
            do j = G%lo(2),G%hi(2)
                e = abs(s%f(i,j,k) - f(G%x(i),G%y(j),G%z(k)))
                if (e > small) then
                    call print_error_message('Error greather than tollerance on left boundary')
                endif
            end do
        end do

        ! **** Right boundary **********************************************************************
        i = G%hi(1) + 1
        do k = G%lo(3),G%hi(3)
            do j = G%lo(2),G%hi(2)
                e = abs(s%f(i,j,k) - f(G%x(i),G%y(j),G%z(k)))
                if (e > small) then
                    call print_error_message('Error greather than tollerance on right boundary')
                endif
            end do
        end do

        ! **** Top boundary ************************************************************************
        j = G%hi(2) + 1
        do k = G%lo(3),G%hi(3)
            do i = G%lo(1),G%hi(1)
                e = abs(s%f(i,j,k) - f(G%x(i),G%y(j),G%z(k)))
                if (e > small) then
                    call print_error_message('Error greather than tollerance on top boundary')
                endif
            end do
        end do

        ! **** Bottom boundary *********************************************************************
        j = G%lo(2) - 1
        do k = G%lo(3),G%hi(3) 
            do i = G%lo(1),G%hi(1)
                e = abs(s%f(i,j,k) - f(G%x(i),G%y(j),G%z(k)))
                if (e > small) then
                    call print_error_message('Error greather than tollerance on bottom boundary')
                endif
            end do
        end do

        ! **** Front boundary **********************************************************************
        k = G%hi(3) + 1
        do j = G%lo(2),G%hi(2)
            do i = G%lo(1),G%hi(1)
                e = abs(s%f(i,j,k) - f(G%x(i),G%y(j),G%z(k)))
                if (e > small) then
                    call print_error_message('Error greather than tollerance on front boundary')
                endif
            end do
        end do

        ! **** Back boundary ***********************************************************************
        k = G%lo(3) - 1
        do j = G%lo(2),G%hi(2)
            do i = G%lo(1),G%hi(1)
                e = abs(s%f(i,j,k) - f(G%x(i),G%y(j),G%z(k)))
                if (e > small) then
                    call print_error_message('Error greather than tollerance on back boundary')
                endif
            end do
        end do

    end subroutine check_ghost_nodes
    !===============================================================================================

    !===============================================================================================
    subroutine check_gradient(error)

        ! In/Out variables
        type(vector), intent(inout) :: error

        ! Local variables
        real(dp)     :: d2
        type(vector) :: nabla_s

        d2 = G%delta/2.0_dp

        ! Compute gradient
        call nabla_s%allocate(G)
        call gradient(s, nabla_s)

        ! Compute error
        do k = G%lo(3),G%hi(3)
            do j = G%lo(2),G%hi(2)
                do i = G%lo(1),G%hi(1)
                    error%x%f(i,j,k) = abs(nabla_s%x%f(i,j,k) - dfdx(G%x(i) + d2, G%y(j), G%z(k)))
                    error%y%f(i,j,k) = abs(nabla_s%y%f(i,j,k) - dfdy(G%x(i), G%y(j) + d2, G%z(k)))
                    error%z%f(i,j,k) = abs(nabla_s%z%f(i,j,k) - dfdz(G%x(i), G%y(j), G%z(k) + d2))                    
                end do
            end do
        end do

        call nabla_s%destroy()

    end subroutine check_gradient
    !===============================================================================================

    !===============================================================================================
    subroutine check_laplacian(error)

        ! In/Out variables
        type(scalar), intent(inout) :: error

        ! Local variables
        type(scalar) :: nabla2_s

        ! Compute laplacian
        call nabla2_s%allocate(G)
        call laplacian(s, nabla2_s)

        ! Compute error
        do k = G%lo(3),G%hi(3)
            do j = G%lo(2),G%hi(2)
                do i = G%lo(1),G%hi(1)
                    error%f(i,j,k) = abs(nabla2_s%f(i,j,k) - d2fdx2(G%x(i), G%y(j), G%z(k))  &
                                                           - d2fdy2(G%x(i), G%y(j), G%z(k))  &
                                                           - d2fdz2(G%x(i), G%y(j), G%z(k)))                    
                end do
            end do
        end do

        call nabla2_s%destroy()

    end subroutine check_laplacian
    !===============================================================================================

end program
