program main

    ! Test the polynomial reconstruction of a given fuciont

    use mpi
    use precision                , only : dp
    use constants                , only : pi
    use class_Grid               , only : base_grid, bc_type
    use polynomial_reconstruction
    
    implicit none

    integer  :: ierror, Nx, Ny, Nz, i, ie(3), j, l, res, out_id
    real(dp) :: Lx, Ly, Lz, x, y, f(-1:1,-1:1), f_reconstructed, e, e_Linf, e_L2
    type(bc_type) :: bc(4)

    call mpi_init(ierror)

    open(newunit = out_id, file = 'error')

    resolution_loop: do res = 3,7

        ! Create a grid
        Nx = 2**res
        Ny = 2**res
        Nz = 1
        Lx = 1.0_dp
        Ly = 1.0_dp
        Lz = Lx/reaL(Nx, dp)*real(Nz, dp)
        do i = 1,4
            bc(i)%s = 'Periodic'
        end do 
        call base_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, [0.0_dp, 0.0_dp, 0.0_dp], 1, 1, bc)

        ! Init the polynomial reconstruction routine
        call init(2)

        ! Evaluate reconstruction error for several location inside the unit domain
        e_Linf = 0.0_dp
        e_L2 = 0.0_dp
        do l = 1,100
            x = rand()
            y = rand()

            ! Find the closest grid cell center
            ie = base_grid%closest_grid_node([x, y, 0.0_dp], 0)

            ! Build the support kernel for the polynomail interpolation   
            do j = -1,1
                do i = -1,1
                    f(i,j) = test_function([base_grid%x(ie(1)) + i*base_grid%delta, &
                                            base_grid%y(ie(2)) + j*base_grid%delta])
                end do
            end do
        
            ! Evaluate the polynomial on (x,y)
            f_reconstructed = P2nd([x - base_grid%x(ie(1)), y - base_grid%y(ie(2))], &
                                    get_coefficients(f))

            ! Evaluate error
            e = abs(f_reconstructed - test_function([x,y]))
            if (e > e_Linf) e_Linf = e
            e_L2 = e_L2 + e**2
        end do

        write(out_id,*) Nx, e_Linf, sqrt(e_L2/l)

        call base_grid%destroy()
        call destroy()

    end do resolution_loop

   ! Finalize the simulation
   call MPI_FINALIZE(ierror)

contains

    !=====================================================================================
    function test_function(X) result(f)

        ! In/Out variables
        real(dp), intent(in) :: X(:)
        real(dp)             :: f

        f = cos(2.0_dp*pi*x(1))*sin(2.0_dp*x(2))

    end function test_function
    !=====================================================================================

end program
