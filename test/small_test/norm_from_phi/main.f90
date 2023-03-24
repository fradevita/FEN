program main

    ! Test the norm computation from a signed distance function using 2nd order central differencing 
    ! and polynomial reconstruction with different order 

    use mpi
    use precision
    use functions
    use class_Grid
    use class_Scalar
    use polynomial_reconstruction, only : ap, get_coefficients, init, destroy, P, ks

    implicit none

    ! Parameters
    real(dp), parameter :: xc = 0.525_dp, yc = 0.464_dp, r = 0.368_dp
    real(dp), parameter :: eps = 1.0e-8_dp, sd = 1.0e-14

    ! Variables
    integer :: ierror, Nx, Ny, Nz, i, j, k, res, resolutions(7), ii, jj, c, order
    real(dp) :: Lx, Ly, Lz, origin(3), epsilon, max_ex, max_ey, ex, ey, modnorm, xx, yy
    type(bc_type) :: bc(4)
    type(scalar) :: phi, norm_x, norm_y
    type(function_type) :: distance
    
    ! Initialize MPI
    call mpi_init(ierror)

    ! The domain is a unit squared box
    Lx = 1.0_dp
    Ly = 1.0_dp
    origin = [0.0_dp, 0.0_dp, 0.0_dp]
    bc(1)%s = 'Periodic'
    bc(2)%s = 'Periodic'
    bc(3)%s = 'Periodic'
    bc(4)%s = 'Periodic'

    open(10, file = 'ecd.dat')
    open(12, file = 'epr2.dat')
    open(14, file = 'epr4.dat')
   
    resolutions = [8, 16, 32, 64, 128, 256, 512]

    do res = 1,7

        ! Set resolution
        Nx = resolutions(res)
        Ny = resolutions(res)
        Nz = 1
        Lz = Lx*float(Nz)/float(Nx)
    
        ! Create the grid
        call base_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 1, 1, bc)

        epsilon = sqrt(2.0_dp)*base_grid%delta/2.0_dp

        call phi%allocate(1)
   
        distance%f => circle
        distance%args = [xc, yc, r]
        call phi%set_from_function(distance)
        call phi%apply_bc

        ! Norm from phi via 2nd order central differencing
        call norm_x%allocate()
        call norm_y%allocate()
        max_ex = 0.0_dp
        max_ey = 0.0_dp
        do k = base_grid%lo(3),base_grid%hi(3)
            do j = base_grid%lo(2),base_grid%hi(2)
                do i = base_grid%lo(1),base_grid%hi(1)
                    if ( abs(phi%f(i,j,1)) < epsilon ) then
                        norm_x%f(i,j,k) = (phi%f(i+1,j-1,k) + 2.0_dp*phi%f(i+1,j,k) + phi%f(i+1,j+1,k) - &
                                           phi%f(i-1,j-1,k) - 2.0_dp*phi%f(i-1,j,k) - phi%f(i-1,j+1,k))
                        norm_y%f(i,j,k) = (phi%f(i-1,j+1,k) + 2.0_dp*phi%f(i,j+1,k) + phi%f(i+1,j+1,k) - &
                                           phi%f(i-1,j-1,k) - 2.0_dp*phi%f(i,j-1,k) - phi%f(i+1,j-1,k))
                        modnorm = sqrt( norm_x%f(i,j,k)**2 + norm_y%f(i,j,k)**2 + sd)
                        norm_x%f(i,j,k) = norm_x%f(i,j,k)/modnorm
                        norm_y%f(i,j,k) = norm_y%f(i,j,k)/modnorm
                        ex = abs(norm_x%f(i,j,k) - an_norm_x(base_grid%x(i),base_grid%y(j)))
                        if (ex > max_ex) max_ex = ex
                        ey = abs(norm_y%f(i,j,k) - an_norm_y(base_grid%x(i),base_grid%y(j)))
                        if (ey > max_ey) max_ey = ey
                    endif
                end do
            end do
        end do
        write(10,*) Nx, max_ex, max_ey

        ! Norm from phi via polynomial reconstruction with order 2 and 4
        order_loop: do order = 2,4,2

            ! For 4th order we must increase ghost nodes of phi
            if (order == 4) then
                call phi%destroy()
                call phi%allocate(2)
                call phi%set_from_function(distance)
                call phi%apply_bc
            end if

            call init(order)
            max_ex = 0.0_dp
            max_ey = 0.0_dp
            norm_x%f = 0.0_dp
            norm_y%f = 0.0_dp
            c = 0
            do k = base_grid%lo(3),base_grid%hi(3)
                do j = base_grid%lo(2),base_grid%hi(2)
                    do i = base_grid%lo(1),base_grid%hi(1)
                        if ( abs(phi%f(i,j,1)) < epsilon ) then
                            ap = get_coefficients(phi%f(i-ks:i+ks,j-ks:j+ks,k))
                            if (order == 1) then
                                norm_x%f(i,j,k) = ap(2)
                                norm_y%f(i,j,k) = ap(3)
                            elseif (order == 2) then
                                norm_x%f(i,j,k) = ap(2)
                                norm_y%f(i,j,k) = ap(4)
                            elseif (order == 4) then
                                norm_x%f(i,j,k) = ap(2)
                                norm_y%f(i,j,k) = ap(6)
                            endif
                            modnorm = sqrt( norm_x%f(i,j,k)**2 + norm_y%f(i,j,k)**2 + sd)
                            norm_x%f(i,j,k) = norm_x%f(i,j,k)/modnorm
                            norm_y%f(i,j,k) = norm_y%f(i,j,k)/modnorm
                            ex = abs(norm_x%f(i,j,k) - an_norm_x(base_grid%x(i),base_grid%y(j)))
                            if (ex > max_ex) max_ex = ex
                            ey = abs(norm_y%f(i,j,k) - an_norm_y(base_grid%x(i),base_grid%y(j)))
                            if (ey > max_ey) max_ey = ey

                            ! Save configuration for one case
                            if (res == 2) then
                                do jj = 1,5
                                    yy = base_grid%y(j) - base_grid%delta/2.0_dp + (jj - 1.0_dp)*base_grid%delta/4.0_dp
                                    do ii = 1,5
                                        xx = base_grid%x(i) - base_grid%delta/2.0_dp + (ii - 1.0_dp)*base_grid%delta/4.0_dp
                                        write(100*order + c,*) xx, yy, P([xx - base_grid%x(i), yy - base_grid%y(j)], ap)
                                    end do
                                    write(100*order + c,*) ''
                                end do
                                c = c + 1
                            end if

                        endif
                    end do
                end do
            end do
            
            write(10 + order,*) Nx, max_ex, max_ey
            
            call destroy()

        end do order_loop

        call phi%destroy()
        call norm_x%destroy
        call norm_y%destroy
        call base_grid%destroy

    end do

contains

    !==============================================================================================
    real(dp) function an_norm_x(x,y)

        real(dp), intent(in) :: x, y
        real(dp) :: nx, ny

        nx = x - xc
        ny = y - yc

        an_norm_x = -nx/sqrt(nx**2 + ny**2 + sd)

    end function
    !==============================================================================================

    !==============================================================================================
    real(dp) function an_norm_y(x,y)

        real(dp), intent(in) :: x, y
        real(dp) :: nx, ny

        nx = x - xc
        ny = y - yc

        an_norm_y = -ny/sqrt(nx**2 + ny**2 + sd)

    end function
    !==============================================================================================

end program main