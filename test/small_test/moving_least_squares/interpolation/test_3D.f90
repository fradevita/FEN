program test_mls_interpolation_3D

    ! Program to test the interpolation of a scalar function with the 
    ! Moving Least Square method. The scalar function is defined in 
    ! test_function subroutine.

    use mpi
    use precision_mod , only : dp
    use global_mod    , only : ierror, pi, stagger, Ndim, small
    use grid_mod 
    use scalar_mod
    use mls_mod       , only : interpolate

    implicit none

    integer       :: r, Nx, Ny, Nz, np, i, j, k, ie(3)
    real(dp)      :: Lx, Ly, Lz, delta, deltal, x, y, z, solution, omega
    real(dp)      :: xl(Ndim), fl(Ndim+1), origin(3), e(7)
    type(grid)    :: comp_grid
    type(scalar)  :: f, dfdx, dfdy, dfdz

    ! Initialize MPI
    call mpi_init(ierror)

    ! Domain size is [0, 10][0, 10][0, 10]
    Lx = 10.0_dp
    Ly = 10.0_dp
    Lz = 10.0_dp

    ! Wavenumber
    omega = 2.0_dp*pi/10.0_dp

    open(1, file = 'error_3D.csv')
    write(1,'(A28)') 'N,f1,f2,f3,f4,dfdx,dfdy,dfdz'
    ! Cycle over the resolution
    resolution_cycle: do r = 4,8
     
        Nx = 2**r
        Ny = 2**r
        Nz = 2**r
        origin = [0.0_dp, 0.0_dp, 0.0_dp]
        call comp_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 1, 1)

        ! Eulerian grid spacing
        delta = comp_grid%delta

        ! Allocate memory for the scalar field f and its derivatives
        call f%allocate(comp_grid, 1)
        call dfdx%allocate(comp_grid, 1)
        call dfdy%allocate(comp_grid, 1)
        call dfdz%allocate(comp_grid, 1)

        ! Set the test field
        do k = comp_grid%lo(3),comp_grid%hi(3)
            z = comp_grid%z(k)
            do j = comp_grid%lo(2),comp_grid%hi(2)
                y = comp_grid%y(j)
                do i = comp_grid%lo(1),comp_grid%hi(1)
                    x = comp_grid%x(i)
                       f%f(i,j,k) =    test_function(x, y, z)
                    dfdx%f(i,j,k) = test_function_dx(x, y, z)
                    dfdy%f(i,j,k) = test_function_dy(x, y, z)
                    dfdz%f(i,j,k) = test_function_dz(x, y, z)
                end do
            end do
        end do
        call f%update_ghost_nodes()
        call dfdx%update_ghost_nodes()
        call dfdy%update_ghost_nodes()
        call dfdz%update_ghost_nodes()

        ! Test the value of the interpolation in a cubic box of size 
        ! [0.4, 9.4][0.4, 9.4][0.4, 9.4] with np +1 points per side
        np = 10

        ! Spacing between interpolation points
        deltal = 9.0_dp/real(np, dp)
        
        e = 0.0_dp
        do k = 1,np + 1
            xl(3) = 0.4_dp + (k - 1.0_dp)*deltal
            do j = 1,np + 1
                xl(2) = 0.4_dp + (j - 1.0_dp)*deltal
                do i = 1,np + 1
                    xl(1) = 0.4_dp + (i - 1.0_dp)*deltal
                    ! Find closest grid point
                    ie = comp_grid%closest_grid_node([xl(1), xl(2), xl(3)], 0)
            
                    ! Interpolate in xl
                    fl = interpolate(f, xl, ie, 0)

                    ! Compute the errors
                    solution = test_function(xl(1), xl(2), xl(3))
                    e(1) = e(1) + abs(fl(1) - solution)/abs(solution + small)
                    
                    solution = test_function_dx(xl(1), xl(2), xl(3))
                    e(2) = e(2) + abs(fl(2) - solution)/abs(solution + small)
                    
                    solution = test_function_dy(xl(1), xl(2), xl(3))
                    e(3) = e(3) + abs(fl(3) - solution)/abs(solution + small)
                    
                    solution = test_function_dz(xl(1), xl(2), xl(3))
                    e(4) = e(4) + abs(fl(4) - solution)/abs(solution + small)
                    
                    ! Interpolate dfdx in xl
                    fl = 0.0_dp
                    fl = interpolate(dfdx, xl, ie, 0)
                    solution = test_function_dx(xl(1), xl(2), xl(3))
                    e(5) = e(5) + abs(fl(1) - solution)/abs(solution + small)
                    
                    ! Interpolate dfdy in xl
                    fl = 0.0_dp
                    fl = interpolate(dfdy, xl, ie, 0)
                    solution = test_function_dy(xl(1), xl(2), xl(3))
                    e(6) = e(6) + abs(fl(1) - solution)/abs(solution + small)

                    ! Interpolate dfdz in xl
                    fl = 0.0_dp
                    fl = interpolate(dfdz, xl, ie, 0)
                    solution = test_function_dz(xl(1), xl(2), xl(3))
                    e(7) = e(7) + abs(fl(1) - solution)/abs(solution + small)
                end do
            end do
        end do

        ! Free the memory
        call f%destroy()
        call dfdx%destroy()
        call dfdy%destroy()
        call dfdz%destroy()
        call comp_grid%destroy()

        ! Print errors
        if (comp_grid%rank == 0) write(1,'(*(G0.7,:,","))') Nx, &
                            e/float((np+1)*(np+1)*(np+1))
     
    end do resolution_cycle

    ! Finalize the simulation
    call MPI_FINALIZE(ierror)

contains

    !===========================================================================
    function test_function(x, y, z) result(f)

        real(dp), intent(in) :: x, y, z
        real(dp) :: f

        f = 1.5_dp + sin(omega*x)*cos(omega*z) + cos(omega*y)*cos(omega*z) 

    end function test_function
    !===========================================================================

    !===========================================================================
    function test_function_dx(x, y, z) result(f)

        real(dp), intent(in) :: x, y, z
        real(dp) :: f

        f = omega*cos(omega*x)*cos(omega*z)

    end function test_function_dx
    !===========================================================================

    !===========================================================================
    function test_function_dy(x, y, z) result(f)
        
        real(dp), intent(in) :: x, y, z
        real(dp) :: f
        
        f = -omega*sin(omega*y)*cos(omega*z)
        
    end function test_function_dy
    !===========================================================================

    !===========================================================================
    function test_function_dz(x, y, z) result(f)
        
        real(dp), intent(in) :: x, y, z        
        real(dp) :: f
        
        f = -omega*sin(omega*x)*sin(omega*z) - omega*cos(omega*y)*sin(omega*z)
        
    end function test_function_dz
    !===========================================================================

end program test_mls_interpolation_3D
