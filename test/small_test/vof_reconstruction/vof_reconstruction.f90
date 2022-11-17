program main

  use mpi
  use precision      , only : dp
  use class_Grid     , only : base_grid, bc_type
  use volume_of_fluid, only : allocate_vof_fields, get_h_from_vof, h, distance
  use volume_of_fluid, only : get_vof_from_distance
  use class_Scalar
  use io
  use decomp_2d
  
  implicit none

  integer  :: ierror, Nx, Ny, Nz
  real(dp) :: Lx, Ly, Lz, origin(3)
  type(bc_type) :: bc(4)
  type(scalar)  :: phi 
  
  ! Initialize MPI
  call mpi_init(ierror)

  ! The domain is a unit squared box
  Lx = 1.0_dp
  Ly = 1.0_dp
  
  ! Set the resolution
  Nx = 16
  Ny = 16
  Nz = 1
  Lz = Lx*float(Nz)/float(Nx)
  origin = [0.0_dp, 0.0_dp, 0.0_dp]
  bc(1)%s = 'Periodic'
  bc(2)%s = 'Periodic'
  bc(3)%s = 'Periodic'
  bc(4)%s = 'Periodic'

  ! Create the grid
  call base_grid%setup(Nx, Ny, Nz, Lx, Ly, Lz, origin, 1, 1, bc)

  ! Allocate all vof fields
  call allocate_vof_fields
  
  ! Set the initial VoF from a distance function
  distance => circle
  call get_vof_from_distance

  ! Reconstruct the interface
  call get_h_from_vof

  ! Print the hyperbolic tangent H on a file
  call h%write('h.dat')

  ! Get distance function from H
  call phi%allocate(1)
  call get_distance(phi, h)

  ! Check error on distance function
  call check(phi)

  ! Finalize the simulation
  call decomp_2d_finalize
  call MPI_FINALIZE(ierror)

contains

  !========================================================================================
  function circle(x,y) result(d)

    ! In/Out variables
    real(dp), intent(in) :: x, y

    ! Local variables
    real(dp) :: x0, y0, r, d

    x0 = 0.5_dp
    y0 = 0.5_dp
    r  = 0.433_dp
    d = sqrt((x - x0)**2 + (y - y0)**2) - r
    
  end function circle
  !========================================================================================

  !========================================================================================
  subroutine get_distance(phi, h)

    use volume_of_fluid, only : beta
    
    type(scalar), intent(in   ) :: h
    type(scalar), intent(inout) :: phi

    integer :: i, j, k
    
    do k = base_grid%lo(3),base_grid%hi(3)
       do j = base_grid%lo(2),base_grid%hi(2)
          do i = base_grid%lo(1),base_grid%hi(1)
             phi%f(i,j,k) = -1.0_dp/beta*atanh(2.0_dp*h%f(i,j,k) - 1.0_dp)*base_grid%delta
          end do
       end do
    end do

  end subroutine get_distance
  !========================================================================================

  !========================================================================================
  subroutine check(phi)

    type(scalar), intent(inout) :: phi

    integer :: i, j, k
    real(dp) :: sol, e, e_max

    e_max = 0.0_dp
    do k = base_grid%lo(3),base_grid%hi(3)
       do j = base_grid%lo(2),base_grid%hi(2)
          do i = base_grid%lo(1),base_grid%hi(1)
             sol = circle(base_grid%x(i), base_grid%y(j))
             e = abs(sol - phi%f(i,j,k))
             if (e > e_max) e_max = e
             write(99,*) i, j, sol, phi%f(i,j,k)
          end do
          write(99,*) ''
       end do
    end do

    print *, e_max
  end subroutine check
  !========================================================================================
  
end program main
