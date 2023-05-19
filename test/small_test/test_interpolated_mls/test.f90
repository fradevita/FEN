program test_mls_interpolation

  ! Program to test the interpolation of a scalar function with the Moving Least Square method.
  ! The scalar function is defined in test_function

  use mpi
  use precision , only : dp
  use constants , only : pi
  use class_Grid, only : base_grid, bc_type
  use class_Scalar
  use mls       , only : interpolate

  implicit none

  integer  :: ierror, r, Nx, Ny, Nz, np, i, j, k, ie(3)
  real(dp) :: Lx, Ly, Lz, delta, deltal, xl(2), fl(3), origin(3), e(5)
  type(scalar) :: f, dfdx, dfdy
  type(bc_type) :: bc(4)

  ! Initialize MPI
  call mpi_init(ierror)

  ! Domain size is [0, 10][0, 10]
  Lx = 10.0_dp
  Ly = 10.0_dp
  bc(1)%s = 'Periodic'
  bc(2)%s = 'Periodic'
  bc(3)%s = 'Periodic'
  bc(4)%s = 'Periodic'

  open(1, file = 'error.dat')  
  
  ! Cycle over the resolution
  resolution_cycle: do r = 5,10
     
     Nx = 2**r
     Ny = 2**r
     Nz = 1
     Lz = Lx*float(Nz)/float(Nx)
     origin = [0.0_dp, 0.0_dp, 0.0_dp]
     call base_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 4, 1, bc)

     ! Eulerian grid spacing
     delta = base_grid%delta

     ! Allocate memory for the scalar field f and its derivatives
     call f%allocate(1)
     call dfdx%allocate(1)
     call dfdy%allocate(1)

     ! Set the test field
     do k = base_grid%lo(3),base_grid%hi(3)
        do j = base_grid%lo(2),base_grid%hi(2)
           do i = base_grid%lo(1),base_grid%hi(1)
              f%f(i,j,k) = test_function(base_grid%x(i),base_grid%y(j))
              dfdx%f(i,j,k) = test_function_dx(base_grid%x(i),base_grid%y(j))
              dfdy%f(i,j,k) = test_function_dy(base_grid%x(i),base_grid%y(j))
           end do
        end do
     end do
     call f%apply_bc()
     call dfdx%apply_bc()
     call dfdy%apply_bc()

     ! Test the value of the interpolation in a box [0.4, 9.4][0.4, 9.4]
     ! with Np +1 points per side
     np = 10
     ! Spacing between interpolation points
     deltal = 9.0_dp/real(np, dp)
     e = 0.0_dp
     do j = 1,np + 1
        xl(2) = 0.4_dp + (j - 1.0_dp)*deltal
        do i = 1,np + 1
           xl(1) = 0.4_dp + (i - 1.0_dp)*deltal

           ! Find closest grid point
           ie = base_grid%closest_grid_node([xl(1), xl(2), 0.0_dp], 0)
           
           ! Interpolate in xl
           fl = interpolate(f, xl, ie, 0)

           ! Compute the errors
           e(1) = e(1) + abs(fl(1) - test_function(xl(1), xl(2)))/ &
                abs(test_function(xl(1),xl(2)) + 1.0e-14_dp)
           e(2) = e(2) + abs(fl(2) - test_function_dx(xl(1), xl(2)))/ &
                abs(test_function_dx(xl(1),xl(2)) + 1.0e-14_dp)
           e(3) = e(3) + abs(fl(3) - test_function_dy(xl(1), xl(2)))/ &
                abs(test_function_dx(xl(1),xl(2)) + 1.0e-14_dp)

           ! Interpolate dfdx in xl
           fl = 0.0_dp
           fl = interpolate(dfdx, xl, ie, 0)
           e(4) = e(4) + abs(fl(1) - test_function_dx(xl(1), xl(2)))/ &
                abs(test_function_dx(xl(1),xl(2)) + 1.0e-14_dp)
           
           ! Interpolate dfdy in xl
           fl = 0.0_dp
           fl = interpolate(dfdy, xl, ie, 0)
           e(5) = e(5) + abs(fl(1) - test_function_dy(xl(1), xl(2)))/ &
                abs(test_function_dx(xl(1),xl(2)) + 1.0e-14_dp)
        end do
     end do

     call f%destroy()
     call dfdx%destroy()
     call dfdy%destroy()
     call base_grid%destroy()

     ! Print errors
     if (base_grid%rank == 0) write(1,*) Nx, e/float((np+1)*(np+1))
     
  end do resolution_cycle

  ! Finalize the simulation
  call MPI_FINALIZE(ierror)

contains

  !========================================================================================
  function test_function(x,y) result(f)

    real(dp), intent(in) :: x, y

    real(dp) :: f

    f = sin(2.0_dp*pi*x/10.0_dp)*cos(2.0*pi*y/10.0_dp) + 1.5_dp

  end function test_function
  !========================================================================================

  !========================================================================================
  function test_function_dx(x, y) result(f)

    real(dp), intent(in) :: x, y

    real(dp) :: f

    f = 2.0_dp*pi/10.0_dp*cos(2.0_dp*pi*x/10.0_dp)*cos(2.0*pi*y/10.0_dp)

  end function test_function_dx
  !========================================================================================

  !========================================================================================
  function test_function_dy(x, y) result(f)
    
    real(dp), intent(in) :: x, y
    
    real(dp) :: f
    
    f = -2.0_dp*pi/10.0_dp*sin(2.0_dp*pi*x/10.0_dp)*sin(2.0*pi*y/10.0_dp)
    
  end function test_function_dy
  !========================================================================================

end program test_mls_interpolation
